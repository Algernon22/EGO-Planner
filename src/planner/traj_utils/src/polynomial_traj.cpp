#include <iostream>
#include <traj_utils/polynomial_traj.h>

/**
 * 基于 minimum snap 的多段多项式轨迹生成
 * 
 * Pos: 3×N矩阵，表示N个路径点的坐标(x,y,z)
 * start_vel, end_vel: 起点和终点的速度约束
 * start_acc, end_acc: 起点和终点的加速度约束
 * Time: 每段轨迹的时间分配
 */
PolynomialTraj PolynomialTraj::minSnapTraj(const Eigen::MatrixXd &Pos, const Eigen::Vector3d &start_vel,
                                           const Eigen::Vector3d &end_vel, const Eigen::Vector3d &start_acc,
                                           const Eigen::Vector3d &end_acc, const Eigen::VectorXd &Time)
{
  int seg_num = Time.size();  // 轨迹段数
  Eigen::MatrixXd poly_coeff(seg_num, 3 * 6); // 存储每段轨迹的五次多项式系数(6个x系数+6个y系数+6个z系数)
  Eigen::VectorXd Px(6 * seg_num), Py(6 * seg_num), Pz(6 * seg_num);  // Px:存储所有段在该轴上的多项式系数 ......

  int num_f, num_p; // number of fixed and free variables
  int num_d;        // number of all segments' derivatives

  // 阶乘计算函数
  const static auto Factorial = [](int x) 
  {
    int fac = 1;
    for (int i = x; i > 0; i--)
      fac = fac * i;
    return fac;
  };

  // 端点及端点微分约束, 每个轨迹段有6个约束(起始位置、起始速度、起始加速度)。
  Eigen::VectorXd Dx = Eigen::VectorXd::Zero(seg_num * 6);
  Eigen::VectorXd Dy = Eigen::VectorXd::Zero(seg_num * 6);
  Eigen::VectorXd Dz = Eigen::VectorXd::Zero(seg_num * 6);
  // Dx存储了所有轨迹段的在x轴上的约束，对任意一段轨迹，Dx上的元素依次代表：起始位置、终止位置、起始速度、终止速度、起始加速度、终止加速度 ......
  // 如Dx0 - Dx5 依次代表第一段轨迹在x轴上的约束：起始位置、终止位置、起始速度、终止速度、起始加速度、终止加速度
  for (int k = 0; k < seg_num; k++)
  {
    // 设置每段轨迹起点和终点的位置约束
    Dx(k * 6) = Pos(0, k);
    Dx(k * 6 + 1) = Pos(0, k + 1);
    Dy(k * 6) = Pos(1, k);
    Dy(k * 6 + 1) = Pos(1, k + 1);
    Dz(k * 6) = Pos(2, k);
    Dz(k * 6 + 1) = Pos(2, k + 1);

    /* 对初始速度、加速度进行约束*/
    if (k == 0)
    {
      Dx(k * 6 + 2) = start_vel(0);
      Dy(k * 6 + 2) = start_vel(1);
      Dz(k * 6 + 2) = start_vel(2);

      Dx(k * 6 + 4) = start_acc(0);
      Dy(k * 6 + 4) = start_acc(1);
      Dz(k * 6 + 4) = start_acc(2);
    }
    /* 对终点速度、加速度进行约束*/
    else if (k == seg_num - 1)
    {
      Dx(k * 6 + 3) = end_vel(0);
      Dy(k * 6 + 3) = end_vel(1);
      Dz(k * 6 + 3) = end_vel(2);

      Dx(k * 6 + 5) = end_acc(0);
      Dy(k * 6 + 5) = end_acc(1);
      Dz(k * 6 + 5) = end_acc(2);
    }
  }

  // 映射矩阵A的构建: 建立多项式系数 P 和导数约束 D 之间的关系 D = AP   A的详细由来见飞书minimum snap笔记中的引用部分
  Eigen::MatrixXd Ab;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(seg_num * 6, seg_num * 6);

  /* K代表段数*/
  for (int k = 0; k < seg_num; k++)
  {
    // Ab 是一个 6×6 的矩阵, 每一行对应一个约束条件
    Ab = Eigen::MatrixXd::Zero(6, 6);
    for (int i = 0; i < 3; i++)
    {
      Ab(2 * i, i) = Factorial(i);
      for (int j = i; j < 6; j++)
        Ab(2 * i + 1, j) = Factorial(j) / Factorial(j - i) * pow(Time(k), j - i);
    }
    // 将局部约束矩阵 Ab 填充到全局约束矩阵 A 的对应位置
    A.block(k * 6, k * 6, 6, 6) = Ab;
  }

  // 选择矩阵C的构建: 分离 固定约束（fixed） 和 自由变量（free） 这里的设计我没看懂，暂且认为选择矩阵这样设计
  Eigen::MatrixXd Ct, C;
  // Fixed Constraints: 起点：位置、速度、加速度; 终点：位置、速度、加速度; 中间点位置连续性：(seg_num - 1) * 2; 3 + 3 + (seg_num - 1) * 2 = 2m + 4
  num_f = 2 * seg_num + 4; 
  // Free Variables: 每个中间点有速度和加速度 2 个自由度，共 (seg_num - 1) * 2 个; (seg_num - 1) * 2 = 2m - 2
  num_p = 2 * seg_num - 2; 
  num_d = 6 * seg_num;  
  Ct = Eigen::MatrixXd::Zero(num_d, num_f + num_p);
  // 设置起点约束
  Ct(0, 0) = 1;
  Ct(2, 1) = 1;
  Ct(4, 2) = 1; 
  // 设置终点约束
  Ct(1, 3) = 1;
  Ct(3, 2 * seg_num + 4) = 1;
  Ct(5, 2 * seg_num + 5) = 1;

  Ct(6 * (seg_num - 1) + 0, 2 * seg_num + 0) = 1;
  Ct(6 * (seg_num - 1) + 1, 2 * seg_num + 1) = 1; // Stack the end point
  Ct(6 * (seg_num - 1) + 2, 4 * seg_num + 0) = 1;
  Ct(6 * (seg_num - 1) + 3, 2 * seg_num + 2) = 1; // Stack the end point
  Ct(6 * (seg_num - 1) + 4, 4 * seg_num + 1) = 1;
  Ct(6 * (seg_num - 1) + 5, 2 * seg_num + 3) = 1; // Stack the end point
  // 设置中间点约束
  for (int j = 2; j < seg_num; j++)
  {
    Ct(6 * (j - 1) + 0, 2 + 2 * (j - 1) + 0) = 1;
    Ct(6 * (j - 1) + 1, 2 + 2 * (j - 1) + 1) = 1;
    Ct(6 * (j - 1) + 2, 2 * seg_num + 4 + 2 * (j - 2) + 0) = 1;
    Ct(6 * (j - 1) + 3, 2 * seg_num + 4 + 2 * (j - 1) + 0) = 1;
    Ct(6 * (j - 1) + 4, 2 * seg_num + 4 + 2 * (j - 2) + 1) = 1;
    Ct(6 * (j - 1) + 5, 2 * seg_num + 4 + 2 * (j - 1) + 1) = 1;
  }

  C = Ct.transpose();

  Eigen::VectorXd Dx1 = C * Dx;
  Eigen::VectorXd Dy1 = C * Dy;
  Eigen::VectorXd Dz1 = C * Dz;

  // 构建Minimum Snap 的 Hessian 矩阵 Q   学习一下写法
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(seg_num * 6, seg_num * 6);

  for (int k = 0; k < seg_num; k++)
  {
    for (int i = 3; i < 6; i++)
    {
      for (int j = 3; j < 6; j++)
      {
        Q(k * 6 + i, k * 6 + j) =
            i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) / (i + j - 5) * pow(Time(k), (i + j - 5));
      }
    }
  }

  // 将原始优化问题转化为带约束的二次规划（QP）形式
  Eigen::MatrixXd R = C * A.transpose().inverse() * Q * A.inverse() * Ct;

  // 为每个坐标轴（x, y, z）初始化固定约束向量 Dxf, Dyf, Dzf，维度为 2 * seg_num + 4: 起点约束（位置、速度、加速度）：3 个; 终点约束（位置）：1 个;中间点连续性约束（位置、速度）：每段 2 个（共 2 * seg_num 个）。
  Eigen::VectorXd Dxf(2 * seg_num + 4), Dyf(2 * seg_num + 4), Dzf(2 * seg_num + 4);
  // 提取固定约束部分,从完整约束向量 Dx1, Dy1, Dz1 中提取前 2 * seg_num + 4 个元素，赋值给固定约束向量 Dxf, Dyf, Dzf
  Dxf = Dx1.segment(0, 2 * seg_num + 4);
  Dyf = Dy1.segment(0, 2 * seg_num + 4);
  Dzf = Dz1.segment(0, 2 * seg_num + 4);

  // 将大型矩阵 R 分解为固定约束（fixed）和自由变量（free）对应的子矩阵，以便后续高效求解二次规划（QP）问题
  Eigen::MatrixXd Rff(2 * seg_num + 4, 2 * seg_num + 4);  // 固定约束部分的 Hessian 子矩阵
  Eigen::MatrixXd Rfp(2 * seg_num + 4, 2 * seg_num - 2);  // 固定约束对自由变量的耦合子矩阵
  Eigen::MatrixXd Rpf(2 * seg_num - 2, 2 * seg_num + 4);  // 自由变量对固定约束的耦合子矩阵
  Eigen::MatrixXd Rpp(2 * seg_num - 2, 2 * seg_num - 2);  // 自由变量部分的 Hessian 子矩阵

  // block() 方法参数：block(start_row, start_col, rows, cols) 从 (start_row, start_col) 开始提取 rows × cols 的子矩阵。
  Rff = R.block(0, 0, 2 * seg_num + 4, 2 * seg_num + 4);
  Rfp = R.block(0, 2 * seg_num + 4, 2 * seg_num + 4, 2 * seg_num - 2);
  Rpf = R.block(2 * seg_num + 4, 0, 2 * seg_num - 2, 2 * seg_num + 4);
  Rpp = R.block(2 * seg_num + 4, 2 * seg_num + 4, 2 * seg_num - 2, 2 * seg_num - 2);

  /* 
      闭式求
      通过解析方法直接求解自由变量（Dxp, Dyp, Dzp）的最优值
      具体数学流程见飞书笔记引用部分
  */
  Eigen::VectorXd Dxp(2 * seg_num - 2), Dyp(2 * seg_num - 2), Dzp(2 * seg_num - 2);
  Dxp = -(Rpp.inverse() * Rfp.transpose()) * Dxf;
  Dyp = -(Rpp.inverse() * Rfp.transpose()) * Dyf;
  Dzp = -(Rpp.inverse() * Rfp.transpose()) * Dzf;

  // 将优化后的自由变量（Dxp, Dyp, Dzp）回填到完整的轨迹参数向量（Dx1, Dy1, Dz1）中，从而构建最终的轨迹解。
  Dx1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dxp;
  Dy1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dyp;
  Dz1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dzp;

  /* 
      将优化后的轨迹参数（Dx1, Dy1, Dz1）转换为多项式轨迹的系数矩阵（Px, Py, Pz），从而生成可执行的连续轨迹。
      P = A^-1 · D
      实际优化中，D 可能包含冗余或约束（如连续性），需通过 Ct 矩阵将 D 投影到最小参数空间：
      P = A^-1 · Ct · D
  */
  Px = (A.inverse() * Ct) * Dx1;
  Py = (A.inverse() * Ct) * Dy1;
  Pz = (A.inverse() * Ct) * Dz1;

  // 将分段多项式系数（Px, Py, Pz）按轨迹段顺序重组到统一的系数矩阵 poly_coeff 中
  for (int i = 0; i < seg_num; i++)
  {
    poly_coeff.block(i, 0, 1, 6) = Px.segment(i * 6, 6).transpose();
    poly_coeff.block(i, 6, 1, 6) = Py.segment(i * 6, 6).transpose();
    poly_coeff.block(i, 12, 1, 6) = Pz.segment(i * 6, 6).transpose();
  }

  // 从系数矩阵构建多项式轨迹
  PolynomialTraj poly_traj;
  for (int i = 0; i < poly_coeff.rows(); ++i) // poly_coeff.rows()代表段数
  {
    vector<double> cx(6), cy(6), cz(6);
    for (int j = 0; j < 6; ++j)
    {
      cx[j] = poly_coeff(i, j), cy[j] = poly_coeff(i, j + 6), cz[j] = poly_coeff(i, j + 12);
    }
    reverse(cx.begin(), cx.end());
    reverse(cy.begin(), cy.end());
    reverse(cz.begin(), cz.end());
    double ts = Time(i);
    poly_traj.addSegment(cx, cy, cz, ts); // 将所求得的每一段轨迹系数加入到多项式轨迹中
  }

  return poly_traj;
}

/* 一段minSnap轨迹生成 */
PolynomialTraj PolynomialTraj::one_segment_traj_gen(const Eigen::Vector3d &start_pt, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                                    const Eigen::Vector3d &end_pt, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc,
                                                    double t)
{
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 6), Crow(1, 6);
  Eigen::VectorXd Bx(6), By(6), Bz(6);

  /*
  C = [ 1,  0,      0,       0,        0,        0;
      1, dt,    dt2,     dt3,      dt4,      dt5;
      0,  1,      0,       0,        0,        0;
      0,  1, 2 * dt, 3 * dt2,  4 * dt3,  5 * dt4;
      0,  0,      2,       0,        0,        0;
      0,  0,      2,  6 * dt, 12 * dt2, 20 * dt3];
      注：dtn=t*t*...*t,t的个数为n
  */
  C(0, 5) = 1;
  C(1, 4) = 1;
  C(2, 3) = 2;
  Crow << pow(t, 5), pow(t, 4), pow(t, 3), pow(t, 2), t, 1;
  C.row(3) = Crow;
  Crow << 5 * pow(t, 4), 4 * pow(t, 3), 3 * pow(t, 2), 2 * t, 1, 0;
  C.row(4) = Crow;
  Crow << 20 * pow(t, 3), 12 * pow(t, 2), 6 * t, 2, 0, 0;
  C.row(5) = Crow;

  // 始末位置、速度、加速度
  Bx << start_pt(0), start_vel(0), start_acc(0), end_pt(0), end_vel(0), end_acc(0);
  By << start_pt(1), start_vel(1), start_acc(1), end_pt(1), end_vel(1), end_acc(1);
  Bz << start_pt(2), start_vel(2), start_acc(2), end_pt(2), end_vel(2), end_acc(2);
  
  // 使用QR分解法求解 CX = B，得到三轴的系数向量
  Eigen::VectorXd Cofx = C.colPivHouseholderQr().solve(Bx);
  Eigen::VectorXd Cofy = C.colPivHouseholderQr().solve(By);
  Eigen::VectorXd Cofz = C.colPivHouseholderQr().solve(Bz);

  vector<double> cx(6), cy(6), cz(6);
  for (int i = 0; i < 6; i++)
  {
    cx[i] = Cofx(i);
    cy[i] = Cofy(i);
    cz[i] = Cofz(i);
  }

  PolynomialTraj poly_traj;
  poly_traj.addSegment(cx, cy, cz, t);  //将所求得的每一段轨迹系数加入到多项式轨迹中

  return poly_traj;
}
