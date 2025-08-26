// #include <fstream>
#include <plan_manage/planner_manager.h>
#include <thread>
#include "visualization_msgs/Marker.h" // zx-todo

namespace ego_planner
{

  // SECTION interfaces for setup and query

  EGOPlannerManager::EGOPlannerManager() {}

  EGOPlannerManager::~EGOPlannerManager() {}

  void EGOPlannerManager::initPlanModules(ros::NodeHandle &nh, PlanningVisualization::Ptr vis)
  {
    /* 读取运动约束参数 */
    nh.param("manager/max_vel", pp_.max_vel_, -1.0);
    nh.param("manager/max_acc", pp_.max_acc_, -1.0);
    nh.param("manager/max_jerk", pp_.max_jerk_, -1.0);

    nh.param("manager/feasibility_tolerance", pp_.feasibility_tolerance_, 0.0);   // 允许轨迹违反约束的阈值
    nh.param("manager/control_points_distance", pp_.ctrl_pt_dist, -1.0);    // B样条控制点间距
    nh.param("manager/planning_horizon", pp_.planning_horizen_, 5.0);     // 规划视野范围(米)
    nh.param("manager/use_distinctive_trajs", pp_.use_distinctive_trajs, false);
    nh.param("manager/drone_id", pp_.drone_id, -1);

    local_data_.traj_id_ = 0;   // 初始化轨迹ID计数器
    grid_map_.reset(new GridMap);
    grid_map_->initMap(nh);

    // obj_predictor_.reset(new fast_planner::ObjPredictor(nh));
    // obj_predictor_->init();
    // obj_pub_ = nh.advertise<visualization_msgs::Marker>("/dynamic/obj_prdi", 10); // zx-todo

    /* B样条优化器初始化 */
    bspline_optimizer_.reset(new BsplineOptimizer);
    bspline_optimizer_->setParam(nh);
    bspline_optimizer_->setEnvironment(grid_map_, obj_predictor_);

    /* A*路径搜索初始化 */
    bspline_optimizer_->a_star_.reset(new AStar);
    bspline_optimizer_->a_star_->initGridMap(grid_map_, Eigen::Vector3i(100, 100, 100));

    visualization_ = vis;
  }

  /**
   * @brief 进行无人机路径的重规划。
   *
   * 该函数使用B样条优化技术生成无碰撞路径，并在必要时进行时间重分配。
   * 它根据当前状态和目标状态尝试生成一条可行的路径。
   * 
   * 后端轨迹优化核心
   * 
   * @param start_pt 起始位置向量。
   * @param start_vel 起始速度向量。
   * @param start_acc 起始加速度向量。
   * @param local_target_pt 局部目标位置向量。
   * @param local_target_vel 局部目标速度向量。
   * @param flag_polyInit 是否使用多项式轨迹初始化。
   * @param flag_randomPolyTraj 是否使用随机多项式轨迹。
   *
   * @return 布尔值，表示重规划是否成功。
   *         - `true` 表示规划成功。
   *         - `false` 表示规划失败或已接近目标。
   */
  bool EGOPlannerManager::reboundReplan(
      Eigen::Vector3d start_pt, 
      Eigen::Vector3d start_vel,
      Eigen::Vector3d start_acc, 
      Eigen::Vector3d local_target_pt,
      Eigen::Vector3d local_target_vel, 
      bool flag_polyInit, 
      bool flag_randomPolyTraj
    )
  {
    // 打印重规划的次数
    static int count = 0;
    printf("\033[47;30m\n[drone %d replan %d]==============================================\033[0m\n", pp_.drone_id, count++);
    
    // cout.precision(3);
    // cout << "start: " << start_pt.transpose() << ", " << start_vel.transpose() << "\ngoal:" << local_target_pt.transpose() << ", " << local_target_vel.transpose()
    //      << endl;

    // 如果当前位置和局部目标位置欧式距离相差小于0.2，判定接近目标，返回false，不进行重规划
    if ((start_pt - local_target_pt).norm() < 0.2)
    {
      cout << "Close to goal" << endl;
      continous_failures_count_++;
      return false;
    }

    bspline_optimizer_->setLocalTargetPt(local_target_pt);

    ros::Time t_start = ros::Time::now();
    ros::Duration t_init, t_opt, t_refine;

    /*** STEP 1: INIT ***/
    // 根据距离设计步长
    double ts = (start_pt - local_target_pt).norm() > 0.1 ? 
      pp_.ctrl_pt_dist / pp_.max_vel_ * 1.5 : 
      pp_.ctrl_pt_dist / pp_.max_vel_ * 5; 
    
    vector<Eigen::Vector3d> point_set, start_end_derivatives;
    static bool flag_first_call = true, flag_force_polynomial = false;
    bool flag_regenerate = false;
    do
    {
      point_set.clear();  // 清空路径点集
      start_end_derivatives.clear();  // 清空边界条件（速度/加速度）
      flag_regenerate = false;

      // 多项式初始化路径--无人机从静止启动飞向目标
      if (flag_first_call || flag_polyInit || flag_force_polynomial /*|| ( start_pt - local_target_pt ).norm() < 1.0*/)   
      {
        flag_first_call = false;
        flag_force_polynomial = false;
        PolynomialTraj gl_traj;

        // 轨迹时间计算
        double dist = (start_pt - local_target_pt).norm();
        double time = pow(pp_.max_vel_, 2) / pp_.max_acc_ > dist ?  // pow(pp_.max_vel_, 2) / pp_.max_acc_: 仅加速和减速即可达到最大速度的最小距离
          sqrt(dist / pp_.max_acc_) :   // 短距离,只有加速和减速  依据 dist = 2 * 0.5at^2 求解t
          (dist - pow(pp_.max_vel_, 2) / pp_.max_acc_) / pp_.max_vel_ + 2 * pp_.max_vel_ / pp_.max_acc_; // 加速 → 匀速 → 减速 

        // 非随机模式
        if (!flag_randomPolyTraj)
        {
          // 生成单段多项式轨迹：满足起点/终点的位置、速度、加速度约束（末端加速度设为0
          gl_traj = PolynomialTraj::one_segment_traj_gen(
            start_pt, start_vel, start_acc, 
            local_target_pt, local_target_vel, Eigen::Vector3d::Zero(), time);
        }
        // 随机模式
        else
        {
          // 通过叉积构造第三个正交方向，形成局部坐标系
          Eigen::Vector3d horizen_dir = ((start_pt - local_target_pt).cross(Eigen::Vector3d(0, 0, 1))).normalized();
          Eigen::Vector3d vertical_dir = ((start_pt - local_target_pt).cross(horizen_dir)).normalized();

          // 生成随机中间点
          Eigen::Vector3d random_inserted_pt = (start_pt + local_target_pt) / 2 +
                                               (((double)rand()) / RAND_MAX - 0.5) * (start_pt - local_target_pt).norm() * horizen_dir * 0.8 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989) +
                                               (((double)rand()) / RAND_MAX - 0.5) * (start_pt - local_target_pt).norm() * vertical_dir * 0.4 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989);
          Eigen::MatrixXd pos(3, 3);
          pos.col(0) = start_pt;
          pos.col(1) = random_inserted_pt;
          pos.col(2) = local_target_pt;

          // 构建最小加加速度（Min-Snap）轨迹
          Eigen::VectorXd t(2); // 时间分配（两段轨迹）
          t(0) = t(1) = time / 2;
          gl_traj = PolynomialTraj::minSnapTraj(
            pos, 
            start_vel, 
            local_target_vel, 
            start_acc, 
            Eigen::Vector3d::Zero(), 
            t
          );
        }

        /*
            得到多项式轨迹后,通过动态调整 ts, 从连续轨迹中提取离散点，确保点距合理、数量足够，并保存关键运动信息,使得这些离散点尽可能的描述曲线的形状
        */
        double t;
        bool flag_too_far;
        ts *= 1.5; 

        /* 
            通过循环动态调整时间步长 ts，直到满足以下两个条件：
              相邻采样点距离不超过 pp_.ctrl_pt_dist * 1.5
              采样点数量至少为7（point_set.size() < 7）
          */
        do
        {
          ts /= 1.5;
          point_set.clear();
          flag_too_far = false;
          Eigen::Vector3d last_pt = gl_traj.evaluate(0);  // 记录起点
          for (t = 0; t < time; t += ts)
          {
            // 采样当前点
            Eigen::Vector3d pt = gl_traj.evaluate(t);
            // 若相邻点距离超过阈值（ctrl_pt_dist * 1.5），立即终止当前采样，缩小 ts 后重试
            if ((last_pt - pt).norm() > pp_.ctrl_pt_dist * 1.5)
            {
              flag_too_far = true;
              break;
            }
            last_pt = pt;
            point_set.push_back(pt);
          }
        } while (flag_too_far || point_set.size() < 7);

        t -= ts;
        start_end_derivatives.push_back(gl_traj.evaluateVel(0));
        start_end_derivatives.push_back(local_target_vel);
        start_end_derivatives.push_back(gl_traj.evaluateAcc(0));
        start_end_derivatives.push_back(gl_traj.evaluateAcc(t));
      }
      else  // 从历史轨迹初始化路径--无人机正在执行已有轨迹时收到新的目标点，需将历史轨迹与新目标平滑衔接
      {
        double t;
        double t_cur = (ros::Time::now() - local_data_.start_time_).toSec();  // 计算当前时刻相对于轨迹起始时间的时间差

        vector<double> pseudo_arc_length; // 存储累积路径长度
        vector<Eigen::Vector3d> segment_point;  // 存储采样点
        pseudo_arc_length.push_back(0.0); // 初始化累积长度为0

        // 将连续轨迹离散化后，用折线段长度之和近似真实弧长, 从t_cur开始，按时间步长ts采样历史轨迹
        for (t = t_cur; t < local_data_.duration_ + 1e-3 /*local_data_.duration表示当前轨迹的总时间长度*/; t += ts)
        {
          segment_point.push_back(local_data_.position_traj_.evaluateDeBoorT(t));  // 获取t时刻的位置
          // 计算相邻点距离并累积到pseudo_arc_length
          if (t > t_cur)
          {
            pseudo_arc_length.push_back(
              (segment_point.back() - segment_point[segment_point.size() - 2]).norm() + 
              pseudo_arc_length.back()
            );
          }
        }
        // 回退到最后一个有效采样点的时间
        t -= ts;

        // 计算从当前点到目标点的多项式轨迹时间
        double poly_time = (local_data_.position_traj_.evaluateDeBoorT(t) - local_target_pt).norm() / pp_.max_vel_ * 2;
        // 如果当前点到目标点的距离需要更长时间（poly_time > ts），则生成一条多项式过渡轨迹（gl_traj）连接历史轨迹和目标点
        if (poly_time > ts)
        {
          // 生成从当前状态到目标点的多项式轨迹
          PolynomialTraj gl_traj = PolynomialTraj::one_segment_traj_gen(local_data_.position_traj_.evaluateDeBoorT(t),
                                                                        local_data_.velocity_traj_.evaluateDeBoorT(t),
                                                                        local_data_.acceleration_traj_.evaluateDeBoorT(t),
                                                                        local_target_pt, local_target_vel, Eigen::Vector3d::Zero(), poly_time);
          // 采样多项式轨迹的点并添加到segment_point
          for (t = ts; t < poly_time; t += ts)
          {
            if (!pseudo_arc_length.empty())
            {
              segment_point.push_back(gl_traj.evaluate(t));
              pseudo_arc_length.push_back(
                (segment_point.back() - segment_point[segment_point.size() - 2]).norm() + 
                pseudo_arc_length.back()
              );
            }
            else
            {
              ROS_ERROR("pseudo_arc_length is empty, return!"); // 错误处理
              continous_failures_count_++;
              return false;
            }
          }
        }

        double sample_length = 0;
        double cps_dist = pp_.ctrl_pt_dist * 1.5; 
        size_t id = 0;
        // 调整控制点间距，直到点数≥7
        do
        {
          cps_dist /= 1.5;
          point_set.clear();
          sample_length = 0;
          id = 0;
          while ((id <= pseudo_arc_length.size() - 2) && sample_length <= pseudo_arc_length.back())
          {
            if (sample_length >= pseudo_arc_length[id] && sample_length < pseudo_arc_length[id + 1])
            {
              point_set.push_back((sample_length - pseudo_arc_length[id]) / (pseudo_arc_length[id + 1] - pseudo_arc_length[id]) * segment_point[id + 1] +
                                  (pseudo_arc_length[id + 1] - sample_length) / (pseudo_arc_length[id + 1] - pseudo_arc_length[id]) * segment_point[id]);
              sample_length += cps_dist;
            }
            else
              id++;
          }
          point_set.push_back(local_target_pt);
        } while (point_set.size() < 7);

        start_end_derivatives.push_back(local_data_.velocity_traj_.evaluateDeBoorT(t_cur));
        start_end_derivatives.push_back(local_target_vel);
        start_end_derivatives.push_back(local_data_.acceleration_traj_.evaluateDeBoorT(t_cur));
        start_end_derivatives.push_back(Eigen::Vector3d::Zero());

        if (point_set.size() > pp_.planning_horizen_ / pp_.ctrl_pt_dist * 3) 
        {
          flag_force_polynomial = true;
          flag_regenerate = true;
        }
      }
    } while (flag_regenerate);  // 若路径过长则重新生成

    // 定义控制点矩阵
    Eigen::MatrixXd ctrl_pts, ctrl_pts_temp;
    // 将离散路径点（point_set）转化为B样条曲线的控制点（ctrl_pts），同时考虑时间分配和边界约束
    UniformBspline::parameterizeToBspline(ts, point_set, start_end_derivatives, ctrl_pts/*生成的B样条控制点矩阵*/); 
    // 存储分段信息 pair<int, int>表示一个片段的起始和结束 控制点 索引
    vector<std::pair<int, int>> segments;

    /* 贼长，通过A star 搜索来使轨迹无碰撞   对应Fast_Planner的ESDF提供的障碍物距离作用，和check_collision_and_rebound()很像*/
    /*--------------贼重要---------------*/
    segments = bspline_optimizer_->initControlPoints(ctrl_pts, true);

    t_init = ros::Time::now() - t_start;
    t_start = ros::Time::now();

    /*** STEP 2: OPTIMIZE ***/
    bool flag_step_1_success = false;
    vector<vector<Eigen::Vector3d>> vis_trajs;

    // use_distinctive_trajs=true 
    if (pp_.use_distinctive_trajs)
    {
      // cout << "enter" << endl;
      std::vector<ControlPoints> trajs = bspline_optimizer_->distinctiveTrajs(segments);
      cout << "\033[1;33m"
           << "multi-trajs=" << trajs.size() << "\033[1;0m" << endl;

      double final_cost, min_cost = 999999.0;
      // 多条轨迹,优化多条轨迹，使用最小代价的一条
      for (int i = trajs.size() - 1; i >= 0; i--)
      {
        // B样条优化 
        if (bspline_optimizer_->BsplineOptimizeTrajRebound(ctrl_pts_temp, final_cost, trajs[i], ts))
        {
          /*优化成功*/
          cout << "traj " << trajs.size() - i << " success." << endl;

          flag_step_1_success = true;
          /*得到最小cost的轨迹*/
          if (final_cost < min_cost)
          {
            min_cost = final_cost;
            ctrl_pts = ctrl_pts_temp;
          }

          // visualization
          point_set.clear();
          for (int j = 0; j < ctrl_pts_temp.cols(); j++)
          {
            point_set.push_back(ctrl_pts_temp.col(j));
          }
          vis_trajs.push_back(point_set);
        }
        else
        {
          cout << "traj " << trajs.size() - i << " failed." << endl;
        }
      }

      t_opt = ros::Time::now() - t_start;

      visualization_->displayMultiInitPathList(vis_trajs, 0.2); // This visuallization will take up several milliseconds.
    }else
    {
      flag_step_1_success = bspline_optimizer_->BsplineOptimizeTrajRebound(ctrl_pts, ts);
      t_opt = ros::Time::now() - t_start;
      //static int vis_id = 0;
      visualization_->displayInitPathList(point_set, 0.2, 0);
    }
    cout << "plan_success=" << flag_step_1_success << endl;

    // 规划失败
    if (!flag_step_1_success)
    {
      visualization_->displayOptimalList(ctrl_pts, 0);
      continous_failures_count_++;
      return false;
    }

    t_start = ros::Time::now();

    // 用优化后的控制点生成实际执行的轨迹，并设置物理约束,reboundReplan()生成的轨迹不一定是最优路径，所以只是临时保存
    UniformBspline pos = UniformBspline(ctrl_pts, 3, ts);
    pos.setPhysicalLimits(pp_.max_vel_, pp_.max_acc_, pp_.feasibility_tolerance_);

    /*** STEP 3: 在单机模式下时间重分配，多机不需要 ***/
    // 单机模式
    if (pp_.drone_id <= 0)
    {
      double ratio;                     // 时间重分配后的比例系数
      bool flag_step_2_success = true;
      
      // 首先类似 fast_planner 进行动力学合理性检查求并解得到 ratio
      if (!pos.checkFeasibility(ratio, false))
      {
        // 有超限现象，则需要重新分配时间
        cout << "Need to reallocate time." << endl;
        Eigen::MatrixXd optimal_control_points;

        /* 时间重分配 */
        flag_step_2_success = refineTrajAlgo(pos, start_end_derivatives, ratio, ts, optimal_control_points);
        if (flag_step_2_success)
          pos = UniformBspline(optimal_control_points, 3, ts);
      }

      if (!flag_step_2_success)
      {
        printf("\033[34mThis refined trajectory hits obstacles. It doesn't matter if appeares occasionally. But if continously appearing, Increase parameter \"lambda_fitness\".\n\033[0m");
        continous_failures_count_++;
        return false;
      }
    }
    else
    {
      static bool print_once = true;
      if (print_once)
      {
        print_once = false;
        ROS_ERROR("IN SWARM MODE, REFINE DISABLED!");
      }
    }

    t_refine = ros::Time::now() - t_start;

    // save planned results
    /* 保存最终更新的轨迹*/
    updateTrajInfo(pos, ros::Time::now());

    static double sum_time = 0;
    static int count_success = 0;
    sum_time += (t_init + t_opt + t_refine).toSec();
    count_success++;
    cout << "total time:\033[42m" << (t_init + t_opt + t_refine).toSec() << "\033[0m,optimize:" << (t_init + t_opt).toSec() << ",refine:" << t_refine.toSec() << ",avg_time=" << sum_time / count_success << endl;

    // success. YoY
    continous_failures_count_ = 0;
    return true;
  }

  bool EGOPlannerManager::EmergencyStop(Eigen::Vector3d stop_pos)
  {
    Eigen::MatrixXd control_points(3, 6);
    for (int i = 0; i < 6; i++)
    {
      control_points.col(i) = stop_pos;
    }

    updateTrajInfo(UniformBspline(control_points, 3, 1.0), ros::Time::now());

    return true;
  }

  bool EGOPlannerManager::checkCollision(int drone_id)
  {
    if (local_data_.start_time_.toSec() < 1e9) // It means my first planning has not started
      return false;

    double my_traj_start_time = local_data_.start_time_.toSec();
    double other_traj_start_time = swarm_trajs_buf_[drone_id].start_time_.toSec();

    double t_start = max(my_traj_start_time, other_traj_start_time);
    double t_end = min(my_traj_start_time + local_data_.duration_ * 2 / 3, other_traj_start_time + swarm_trajs_buf_[drone_id].duration_);

    for (double t = t_start; t < t_end; t += 0.03)
    {
      if ((local_data_.position_traj_.evaluateDeBoorT(t - my_traj_start_time) - swarm_trajs_buf_[drone_id].position_traj_.evaluateDeBoorT(t - other_traj_start_time)).norm() < bspline_optimizer_->getSwarmClearance())
      {
        return true;
      }
    }

    return false;
  }

  bool EGOPlannerManager::planGlobalTrajWaypoints(const Eigen::Vector3d &start_pos, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                                  const std::vector<Eigen::Vector3d> &waypoints, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc)
  {

    // generate global reference trajectory

    vector<Eigen::Vector3d> points;
    points.push_back(start_pos);

    for (size_t wp_i = 0; wp_i < waypoints.size(); wp_i++)
    {
      points.push_back(waypoints[wp_i]);
    }

    double total_len = 0;
    total_len += (start_pos - waypoints[0]).norm();
    for (size_t i = 0; i < waypoints.size() - 1; i++)
    {
      total_len += (waypoints[i + 1] - waypoints[i]).norm();
    }

    // insert intermediate points if too far
    vector<Eigen::Vector3d> inter_points;
    double dist_thresh = max(total_len / 8, 4.0);

    for (size_t i = 0; i < points.size() - 1; ++i)
    {
      inter_points.push_back(points.at(i));
      double dist = (points.at(i + 1) - points.at(i)).norm();

      if (dist > dist_thresh)
      {
        int id_num = floor(dist / dist_thresh) + 1;

        for (int j = 1; j < id_num; ++j)
        {
          Eigen::Vector3d inter_pt =
              points.at(i) * (1.0 - double(j) / id_num) + points.at(i + 1) * double(j) / id_num;
          inter_points.push_back(inter_pt);
        }
      }
    }

    inter_points.push_back(points.back());

    // for ( int i=0; i<inter_points.size(); i++ )
    // {
    //   cout << inter_points[i].transpose() << endl;
    // }

    // write position matrix
    int pt_num = inter_points.size();
    Eigen::MatrixXd pos(3, pt_num);
    for (int i = 0; i < pt_num; ++i)
      pos.col(i) = inter_points[i];

    Eigen::Vector3d zero(0, 0, 0);
    Eigen::VectorXd time(pt_num - 1);
    for (int i = 0; i < pt_num - 1; ++i)
    {
      time(i) = (pos.col(i + 1) - pos.col(i)).norm() / (pp_.max_vel_);
    }

    time(0) *= 2.0;
    time(time.rows() - 1) *= 2.0;

    PolynomialTraj gl_traj;
    if (pos.cols() >= 3)
      gl_traj = PolynomialTraj::minSnapTraj(pos, start_vel, end_vel, start_acc, end_acc, time);
    else if (pos.cols() == 2)
      gl_traj = PolynomialTraj::one_segment_traj_gen(start_pos, start_vel, start_acc, pos.col(1), end_vel, end_acc, time(0));
    else
      return false;

    auto time_now = ros::Time::now();
    global_data_.setGlobalTraj(gl_traj, time_now);

    return true;
  }

  /*
      生成从起点到终点的全局参考轨迹
  */
  bool EGOPlannerManager::planGlobalTraj(const Eigen::Vector3d &start_pos, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                         const Eigen::Vector3d &end_pos, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc)
  {
    vector<Eigen::Vector3d> points;   // 定义容器用来存放始末位置点
    points.push_back(start_pos);
    points.push_back(end_pos);

    vector<Eigen::Vector3d> inter_points;   // 存储插值后的路径点
    const double dist_thresh = 4.0;

    for (size_t i = 0; i < points.size() - 1; ++i)
    {
      inter_points.push_back(points.at(i));   // 加入当前点
      double dist = (points.at(i + 1) - points.at(i)).norm();   // 计算两点间距离: 始末xyz位置差平方和开根号

      // 如果距离超过阈值，需要插入中间点 
      if (dist > dist_thresh) 
      {
        // 计算需要插入的点数, floor(x)为取整函数，返回小于等于x的整数
        int id_num = floor(dist / dist_thresh) + 1; 

        // 线性插值插入中间点
        for (int j = 1; j < id_num; ++j)
        {
          /* 
              这里面参考了线性贝塞尔曲线 Bezier Curve  B(t)=P0+(P1-P0)*t=(1-t)*P0+t*P1 ，t 取[0,1] 
              P_0 是起点（points.at(i)）; P_1 是终点（points.at(i+1)）; t 是比例参数，t=0 时为 P_0，t=1 时为 P_1; B(t) 是插值得到的中间点
              t 被离散化为 j / id_num

              假设 P_0 = (0, 0, 0)，P_1 = (9, 0, 0)，dist_thresh = 4.0：
              计算 id_num = floor(9/4) + 1 = 3
              j=1 → t=1/3 → 插值点位于 1/3 处
              j=2 → t=2/3 → 插值点位于 2/3 处
              插入点：
              t=1/3 → (0,0,0)*2/3 + (9,0,0)*1/3 = (3,0,0)
              t=2/3 → (0,0,0)*1/3 + (9,0,0)*2/3 = (6,0,0)
              最终路径点：[ (0,0,0), (3,0,0), (6,0,0) ]
          */
          Eigen::Vector3d inter_pt = 
              points.at(i) * (1.0 - double(j) / id_num) + 
              points.at(i + 1) * double(j) / id_num;
          inter_points.push_back(inter_pt); // 插点
        }
      }
    }

    // 将终点位置传入新的waypoint点中
    inter_points.push_back(points.back());  

    // 构建位置矩阵
    int pt_num = inter_points.size();     // 共有多少个waypoints
    Eigen::MatrixXd pos(3, pt_num);     // 定义一个3 x pt_num的态大小的双精度浮点数矩阵用来存放位置信息
    
    // 将路径点存入矩阵
    for (int i = 0; i < pt_num; ++i)
      pos.col(i) = inter_points[i]; 

    Eigen::Vector3d zero(0, 0, 0);
    Eigen::VectorXd time(pt_num - 1);   // 每段路径的时间分配
    
    // 每一段时间点等于两点间欧式距离除以之前设置好的最大速度
    for (int i = 0; i < pt_num - 1; ++i)
    {
      time(i) = (pos.col(i + 1) - pos.col(i)).norm() / (pp_.max_vel_);  
    }

    // 首尾段时间加倍（可能是为了加减速预留时间）
    time(0) *= 2.0;
    time(time.rows() - 1) *= 2.0;

    PolynomialTraj gl_traj;  // 多项式轨迹对象

    // 根据路径点数选择不同的全局路线生成方式
    if (pos.cols() >= 3)
    {
      // 多点路径：使用minimum snap方法生成平滑轨迹
      gl_traj = PolynomialTraj::minSnapTraj(pos, start_vel, end_vel, start_acc, end_acc, time);
    }
    else if (pos.cols() == 2)
      // 两点路径：生成单段轨迹
      gl_traj = PolynomialTraj::one_segment_traj_gen(start_pos, start_vel, start_acc, end_pos, end_vel, end_acc, time(0));
    else
      return false; // 路径点不足，返回失败

    auto time_now = ros::Time::now();

    // 传入全局轨迹和总时间
    global_data_.setGlobalTraj(gl_traj, time_now);

    return true;
  }

  /**
   * @brief 优化B样条轨迹以满足动力学约束（时间重分配+控制点优化）
   * @param[in/out] traj                   待优化的B样条轨迹（输入不可行轨迹，输出可行轨迹）
   * @param[in]     start_end_derivative   起点/终点的速度、加速度约束（vector<Eigen::Vector3d>）
   * @param[in]     ratio                  时间调整比例（由checkFeasibility计算得出）
   * @param[in/out] ts                     时间步长（输入原始值，输出调整后值）
   * @param[out]    optimal_control_points 优化后的控制点（N×3矩阵）
   * @return bool                          优化是否成功（true：成功，false：失败）
   */
  bool EGOPlannerManager::refineTrajAlgo(UniformBspline &traj, vector<Eigen::Vector3d> &start_end_derivative, double ratio, double &ts, Eigen::MatrixXd &optimal_control_points)
  {
    double t_inc;             // 时间增量
    Eigen::MatrixXd ctrl_pts; // 存储调整后的控制点

    // std::cout << "ratio: " << ratio << std::endl;
    
    // 更新B样条参数
    reparamBspline(traj, start_end_derivative, ratio, ctrl_pts, ts, t_inc);

    traj = UniformBspline(ctrl_pts, 3, ts);

    // 生成参考路径点
    double t_step = traj.getTimeSum() / (ctrl_pts.cols() - 3);
    bspline_optimizer_->ref_pts_.clear();
    for (double t = 0; t < traj.getTimeSum() + 1e-4; t += t_step)
      bspline_optimizer_->ref_pts_.push_back(traj.evaluateDeBoorT(t));
    
      // 调用优化器优化控制点
    bool success = bspline_optimizer_->BsplineOptimizeTrajRefine(ctrl_pts, ts, optimal_control_points);

    return success;
  }

  void EGOPlannerManager::updateTrajInfo(const UniformBspline &position_traj, const ros::Time time_now)
  {
    local_data_.start_time_ = time_now;
    local_data_.position_traj_ = position_traj;
    local_data_.velocity_traj_ = local_data_.position_traj_.getDerivative();
    local_data_.acceleration_traj_ = local_data_.velocity_traj_.getDerivative();
    local_data_.start_pos_ = local_data_.position_traj_.evaluateDeBoorT(0.0);
    local_data_.duration_ = local_data_.position_traj_.getTimeSum();
    local_data_.traj_id_ += 1;
  }

  void EGOPlannerManager::reparamBspline(UniformBspline &bspline, vector<Eigen::Vector3d> &start_end_derivative, double ratio,
                                         Eigen::MatrixXd &ctrl_pts, double &dt, double &time_inc)
  {
    double time_origin = bspline.getTimeSum();// 轨迹总时间
    int seg_num = bspline.getControlPoint().cols() - 3;// 样条轨迹段数
    // double length = bspline.getLength(0.1);
    // int seg_num = ceil(length / pp_.ctrl_pt_dist);

    bspline.lengthenTime(ratio);  //对节点向量重分配
    double duration = bspline.getTimeSum();
    dt = duration / double(seg_num);
    time_inc = duration - time_origin;

    vector<Eigen::Vector3d> point_set;
    for (double time = 0.0; time <= duration + 1e-4; time += dt)
    {
      point_set.push_back(bspline.evaluateDeBoorT(time));//更新控制点所对应的位置
    }
    UniformBspline::parameterizeToBspline(dt, point_set, start_end_derivative, ctrl_pts);//重新得到B样条控制点
  }

} // namespace ego_planner
