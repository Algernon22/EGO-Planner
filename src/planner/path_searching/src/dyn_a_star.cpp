#include "path_searching/dyn_a_star.h"

using namespace std;
using namespace Eigen;

// 析构函数在程序结束后自动释放内存
AStar::~AStar()
{
    for (int i = 0; i < POOL_SIZE_(0); i++)
        for (int j = 0; j < POOL_SIZE_(1); j++)
            for (int k = 0; k < POOL_SIZE_(2); k++)
                delete GridNodeMap_[i][j][k];
}

void AStar::initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size)
{
    POOL_SIZE_ = pool_size; //初始化大小为100，对应栅格地图最大索引值
    CENTER_IDX_ = pool_size / 2;

    GridNodeMap_ = new GridNodePtr **[POOL_SIZE_(0)];
    for (int i = 0; i < POOL_SIZE_(0); i++)
    {
        GridNodeMap_[i] = new GridNodePtr *[POOL_SIZE_(1)];
        for (int j = 0; j < POOL_SIZE_(1); j++)
        {
            GridNodeMap_[i][j] = new GridNodePtr[POOL_SIZE_(2)];
            for (int k = 0; k < POOL_SIZE_(2); k++)
            {
                GridNodeMap_[i][j][k] = new GridNode;
            }
        }
    }

    grid_map_ = occ_map;
}

/*  对角启发式函数  计算任意两个节点之间距离的closed-form solution */
double AStar::getDiagHeu(GridNodePtr node1, GridNodePtr node2)
{
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    double h = 0.0;
    int diag = min(min(dx, dy), dz);
    dx -= diag;
    dy -= diag;
    dz -= diag;

    if (dx == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + 1.0 * abs(dy - dz);
    }
    if (dy == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + 1.0 * abs(dx - dz);
    }
    if (dz == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + 1.0 * abs(dx - dy);
    }
    return h;
}

/*  曼哈顿启发式函数*/
double AStar::getManhHeu(GridNodePtr node1, GridNodePtr node2)
{
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    return dx + dy + dz;
}

/*  欧式距离启发式函数*/
double AStar::getEuclHeu(GridNodePtr node1, GridNodePtr node2)
{
    return (node2->index - node1->index).norm();
}

/*
* @function      AStar::retrievePath
* @brief        回溯获得最优路径
* @param        GridNodePtr current
* @return       vector<Vector3d> path     最优路径
*/
vector<GridNodePtr> AStar::retrievePath(GridNodePtr current)
{
    vector<GridNodePtr> path;
    path.push_back(current);    //提取最后openSet的节点路径，也就是终点

    // 提取闭集closeSet所存的父节点，直到最开始存放的父节点
    while (current->cameFrom != NULL)
    // while (current !=startPtr) 这种写法路径中不包含头结点
    {
        current = current->cameFrom;
        path.push_back(current);
    }

    return path;    //返回path_
}

/* 
    将将欧式空间坐标转换成grid_map索引值，并修正起点和终点位置
*/
bool AStar::ConvertToIndexAndAdjustStartEndPoints(Vector3d start_pt, Vector3d end_pt, Vector3i &start_idx, Vector3i &end_idx)
{
    // 将欧式空间坐标转换成grid_map索引值
    if (!Coord2Index(start_pt, start_idx) || !Coord2Index(end_pt, end_idx))
        return false;
    // 检查起点位置是否在障碍物里面，如果在，通过归一化始末位置差与步长的乘积修改起点位置
    if (checkOccupancy(Index2Coord(start_idx)))
    {
        ROS_WARN("Start point is insdide an obstacle.");
        do
        {
            start_pt = (start_pt - end_pt).normalized() * step_size_ + start_pt;
            if (!Coord2Index(start_pt, start_idx))
                return false;
        } while (checkOccupancy(Index2Coord(start_idx)));
    }
    // 检查终点位置是否在障碍物里面，如果在，通过归一化始末位置差与步长的乘积修改终点位置
    if (checkOccupancy(Index2Coord(end_idx)))
    {
        ROS_WARN("End point is insdide an obstacle.");
        do
        {
            end_pt = (end_pt - start_pt).normalized() * step_size_ + end_pt;
            if (!Coord2Index(end_pt, end_idx))
                return false;
        } while (checkOccupancy(Index2Coord(end_idx)));
    }

    return true;
}

bool AStar::AstarSearch(const double step_size, Vector3d start_pt, Vector3d end_pt)
{
    ros::Time time_1 = ros::Time::now();
    ++rounds_;  // 一个标志位，记录某个节点（neighborPtr）最近一次被探索的轮次

    step_size_ = step_size; // 步长：0.1
    inv_step_size_ = 1 / step_size;
    center_ = (start_pt + end_pt) / 2;

    // 起点和终点的索引向量
    Vector3i start_idx, end_idx;
    // 将将欧式空间坐标转换成grid_map索引值，并修正起点和终点位置
    if (!ConvertToIndexAndAdjustStartEndPoints(start_pt, end_pt, start_idx, end_idx))
    {
        ROS_ERROR("Unable to handle the initial or end point, force return!");
        return false;
    }

    // if ( start_pt(0) > -1 && start_pt(0) < 0 )
    //     cout << "start_pt=" << start_pt.transpose() << " end_pt=" << end_pt.transpose() << endl;
    
    // 初始化起点和终点的节点(指针) 
    GridNodePtr startPtr = GridNodeMap_[start_idx(0)][start_idx(1)][start_idx(2)];
    GridNodePtr endPtr = GridNodeMap_[end_idx(0)][end_idx(1)][end_idx(2)];

    /*
    * @brief 定义优先级队列
    * @param   GridNodePtr：队列元素的类型
    * @param   std::vector<GridNodePtr>：底层容器（默认用 vector 实现堆）
    * @param   NodeComparator：自定义比较器，决定优先级规则
    */
    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> empty;
    // 清空 openSet_
    openSet_.swap(empty); 

    GridNodePtr neighborPtr = NULL;
    // current代表openSet中f(n)最小的节点
    GridNodePtr current = NULL; 

    startPtr->index = start_idx;
    startPtr->rounds = rounds_;
    startPtr->gScore = 0;
    startPtr->fScore = getHeu(startPtr, endPtr);    //获取启发式函数值，采用对角函数
    startPtr->state = GridNode::OPENSET; 
    startPtr->cameFrom = NULL;
    openSet_.push(startPtr); 
    // openSet.insert(make_pair(startPtr->fScore, startPtr));   //将起点加入openset中

    endPtr->index = end_idx;
    double tentative_gScore;

    int num_iter = 0;
    // 开始进入循环
    while (!openSet_.empty())
    {
        num_iter++;
        current = openSet_.top();   // 访问最小的f(n)cost节点
        openSet_.pop(); // 堆顶元素出队也就是优先级最高的

        // if ( num_iter < 10000 )
        //     cout << "current=" << current->index.transpose() << endl;

        // 判断是否到达终点
        if (current->index(0) == endPtr->index(0) && current->index(1) == endPtr->index(1) && current->index(2) == endPtr->index(2))
        {
            // ros::Time time_2 = ros::Time::now();
            // printf("\033[34mA star iter:%d, time:%.3f\033[0m\n",num_iter, (time_2 - time_1).toSec()*1000);
            // if((time_2 - time_1).toSec() > 0.1)
            // ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );

            gridPath_ = retrievePath(current);  // 返回path
            return true;
        }
        // 将当前节点从开集移到闭集中去
        current->state = GridNode::CLOSEDSET; 
        // 八连通搜索 
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    // 搜索自己时，终止当前循环，开启下一次循环
                    if (dx == 0 && dy == 0 && dz == 0)
                        continue;   

                    Vector3i neighborIdx;
                    neighborIdx(0) = (current->index)(0) + dx;
                    neighborIdx(1) = (current->index)(1) + dy;
                    neighborIdx(2) = (current->index)(2) + dz;
                    
                    // 避免访问网格边缘, 超出时跳出循环 
                    if (
                        neighborIdx(0) < 1 || neighborIdx(0) >= POOL_SIZE_(0) - 1 || 
                        neighborIdx(1) < 1 || neighborIdx(1) >= POOL_SIZE_(1) - 1 || 
                        neighborIdx(2) < 1 || neighborIdx(2) >= POOL_SIZE_(2) - 1)
                    {
                        continue;
                    }

                    // 初始化 拓展邻近节点
                    neighborPtr = GridNodeMap_[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];
                    neighborPtr->index = neighborIdx;

                    // 若相等 → 节点已在当前轮次被探索过（flag_explored = true）; 若不等 → 节点未被当前轮次探索
                    bool flag_explored = neighborPtr->rounds == rounds_;

                    // 判断拓展的节点是否在闭集中，如果已经在闭集中，则跳出循环
                    if (flag_explored && neighborPtr->state == GridNode::CLOSEDSET)
                    {
                        continue; 
                    }

                    neighborPtr->rounds = rounds_;

                    // 检查拓展的邻近节点索引是否在障碍物里面，如果在，则跳出循环
                    if (checkOccupancy(Index2Coord(neighborPtr->index)))
                    {
                        continue;
                    }

                    // 计算从当前节点到邻居节点的欧几里得距离
                    double static_cost = sqrt(dx * dx + dy * dy + dz * dz); 
                    // 从起点到当前邻居节点的预估代价
                    tentative_gScore = current->gScore + static_cost;   

                    // 如果第一次探索，发现新的节点（没有在闭集中），加入开集openSet
                    if (!flag_explored)
                    {
                        neighborPtr->state = GridNode::OPENSET;
                        neighborPtr->cameFrom = current;    
                        neighborPtr->gScore = tentative_gScore;
                        neighborPtr->fScore = tentative_gScore + getHeu(neighborPtr, endPtr);
                        openSet_.push(neighborPtr); 
                    }
                    // 找到了一条更短路径到该邻居节点
                    else if (tentative_gScore < neighborPtr->gScore)
                    { 
                        neighborPtr->cameFrom = current;
                        neighborPtr->gScore = tentative_gScore;
                        neighborPtr->fScore = tentative_gScore + getHeu(neighborPtr, endPtr);   //更新完成
                    }
                }
        ros::Time time_2 = ros::Time::now();
        if ((time_2 - time_1).toSec() > 0.2)
        {
            ROS_WARN("Failed in A star path searching !!! 0.2 seconds time limit exceeded.");
            return false;
        }
    }

    ros::Time time_2 = ros::Time::now();

    if ((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in A star path finding is %.3fs, iter=%d", (time_2 - time_1).toSec(), num_iter);

    return false;
}

vector<Vector3d> AStar::getPath()
{
    vector<Vector3d> path;

    for (auto ptr : gridPath_)
        path.push_back(Index2Coord(ptr->index));

    reverse(path.begin(), path.end());  // 由于先进去为终点，需要反一下顺序
    return path;
}
