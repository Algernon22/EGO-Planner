#ifndef _DYN_A_STAR_H_
#define _DYN_A_STAR_H_

#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include <plan_env/grid_map.h>
#include <queue>

constexpr double inf = 1 >> 20;
struct GridNode;
typedef GridNode *GridNodePtr;

struct GridNode
{
	enum enum_state
	{
		OPENSET = 1,
		CLOSEDSET = 2,
		UNDEFINED = 3
	};

	int rounds{0}; // Distinguish every call
	enum enum_state state
	{
		UNDEFINED
	};
	Eigen::Vector3i index;

	double gScore{inf}, fScore{inf};
	GridNodePtr cameFrom{NULL};
};

class NodeComparator
{
	public:
		bool operator()(GridNodePtr node1, GridNodePtr node2)
		{
			//从堆底到堆顶 降序排序 ,也就是小的优先级高，优先级队列
			return node1->fScore > node2->fScore;
		}
};

class AStar
{
	private:
		GridMap::Ptr grid_map_;

		inline void coord2gridIndexFast(const double x, const double y, const double z, int &id_x, int &id_y, int &id_z);

		double getDiagHeu(GridNodePtr node1, GridNodePtr node2);
		double getManhHeu(GridNodePtr node1, GridNodePtr node2);
		double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
		inline double getHeu(GridNodePtr node1, GridNodePtr node2);

		bool ConvertToIndexAndAdjustStartEndPoints(const Eigen::Vector3d start_pt, const Eigen::Vector3d end_pt, Eigen::Vector3i &start_idx, Eigen::Vector3i &end_idx);

		inline Eigen::Vector3d Index2Coord(const Eigen::Vector3i &index) const;
		inline bool Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const;

		//bool (*checkOccupancyPtr)( const Eigen::Vector3d &pos );
		/*
		* @function   checkOccupancy
		* @brief        判断某一栅格是否占用
		* @param       Eigen::Vector3d &pos   要查询的节点位置
		* @return       bool isOccupied         是否占用
		*/
		inline bool checkOccupancy(const Eigen::Vector3d &pos) { return (bool)grid_map_->getInflateOccupancy(pos); }

		std::vector<GridNodePtr> retrievePath(GridNodePtr current);

		double step_size_, inv_step_size_;
		Eigen::Vector3d center_;
		Eigen::Vector3i CENTER_IDX_, POOL_SIZE_;
		const double tie_breaker_ = 1.0 + 1.0 / 10000;

		std::vector<GridNodePtr> gridPath_;

		GridNodePtr ***GridNodeMap_;	//三维指针，指向一个指向指针的指针
		std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> openSet_;

		int rounds_{0};

	public:
		typedef std::shared_ptr<AStar> Ptr;

		AStar(){};
		~AStar();

		void initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size);

		bool AstarSearch(const double step_size, Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);

		std::vector<Eigen::Vector3d> getPath();
};

inline double AStar::getHeu(GridNodePtr node1, GridNodePtr node2)
{
	return tie_breaker_ * getDiagHeu(node1, node2);
}

/*
* @function   AStar::Index2Coord
* @brief        将某个栅格的索引转化为欧式空间中的坐标(取栅格的中心位置)
* @param        Vector3i & index        栅格的索引向量
* @return       Vector3d pt             对应的欧式空间坐标向量
*/
inline Eigen::Vector3d AStar::Index2Coord(const Eigen::Vector3i &index) const
{
	return ((index - CENTER_IDX_).cast<double>() * step_size_) + center_;
};

/*
* @定义内联函数的形式，不用初始化运行速度快 
* @function     AStar::Coord2Index
* @brief        将某个欧式空间中的坐标转化为栅格的索引(向无穷小取整)
* @param        Vector3d & pt           对应的欧式空间坐标向量
* @     				 Vector3i idx            栅格的索引向量
*/
inline bool AStar::Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const
{
	idx = ((pt - center_) * inv_step_size_ + Eigen::Vector3d(0.5, 0.5, 0.5)).cast<int>() + CENTER_IDX_;

	if (idx(0) < 0 || idx(0) >= POOL_SIZE_(0) || idx(1) < 0 || idx(1) >= POOL_SIZE_(1) || idx(2) < 0 || idx(2) >= POOL_SIZE_(2))
	{
		ROS_ERROR("Ran out of pool, index=%d %d %d", idx(0), idx(1), idx(2));
		return false;
	}

	return true;
};

#endif
