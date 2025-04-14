#include "diffusion.hpp"

Diffusion::Diffusion(const Parameter input_para, const double input_equ_ave_freq, 
  const int input_sep, const std::vector<double>& input_ave_freq, 
  const std::vector<std::vector<double>>& input_freq_each){

  para = input_para;
  equ_ave_freq = input_equ_ave_freq;
  sep = input_sep;

  ave_freq = input_ave_freq;
  freq_each = input_freq_each;

  sum_sojourn = 0.0;
  convergence = 1;
}

void Diffusion::calculate_diffusion(){
  std::vector<double> xs;
  std::vector<double> mx_s;
  std::vector<double> vx_s;
  std::vector<double> l1_length;
  std::vector<double> list_eigen;

  std::vector<double> ret_ps_zero, extinct_prob;

  double mm_v_zero = eigen_at_zero();
  double mm_v_one = eigen_at_one(); 
  eigen_at_x(xs, mx_s, vx_s, l1_length, list_eigen);

  for(int i = 1; i < static_cast<int>(xs.size()) - 1; i++){
    if(l1_length.at(i) < 1e-8 && mx_s.at(i - 1) * mx_s.at(i + 1) < 0.0){
      mx_s.at(i) = (mx_s.at(i - 1) + mx_s.at(i + 1)) / 2.0;
      vx_s.at(i) = (vx_s.at(i - 1) + vx_s.at(i + 1)) / 2.0;
      list_eigen.at(i) = (list_eigen.at(i - 1) + list_eigen.at(i + 1)) / 2.0;
    }

    if(std::isnan(vx_s.at(i)) && (!std::isnan(vx_s.at(i - 1)) && !std::isnan(vx_s.at(i + 1)))){
      mx_s.at(i) = (mx_s.at(i - 1) + mx_s.at(i + 1)) / 2.0;
      vx_s.at(i) = (vx_s.at(i - 1) + vx_s.at(i + 1)) / 2.0;
      list_eigen.at(i) = (list_eigen.at(i - 1) + list_eigen.at(i + 1)) / 2.0;
    }
  }

  max_eigen = list_eigen.at(0);
  for(int i = 0; i < static_cast<int>(xs.size()); i++){
    if(list_eigen.at(i) > 0.0){
      convergence = 0;
    }

    if(max_eigen < list_eigen.at(i)){
      max_eigen = list_eigen.at(i);
    }
  }

  double x0 = equ_ave_freq;

  std::vector<double> ys;
  std::vector<double> vys;
  std::vector<double> log_gys;
  std::vector<double> log_int_gys_0_y;
  std::vector<double> log_int_gys_y_1;
  std::vector<double> tmps;

  double int_mm_v = 0.0;
  double log_int_g = std::log(0.0);

  double h1 = 1.0 / 2.0 / sep;
  double h2 = 1.0 / sep;

  for(int i = 0; i < sep; i++){
    double y4 = 1.0 * (i + 1) / sep;

    double mm_v0, mm_v1, mm_v2, mm_v3, mm_v4;
    if(i == 0){
      mm_v0 = mm_v_zero;
      mm_v1 = 2.0 * mx_s.at(4 * i) / vx_s.at(4 * i);
      mm_v2 = 2.0 * mx_s.at(4 * i + 1) / vx_s.at(4 * i + 1);
      mm_v3 = 2.0 * mx_s.at(4 * i + 2) / vx_s.at(4 * i + 2);
      mm_v4 = 2.0 * mx_s.at(4 * i + 3) / vx_s.at(4 * i + 3);
    }else if(i == sep - 1){
      mm_v0 = 2.0 * mx_s.at(4 * i - 1) / vx_s.at(4 * i - 1);
      mm_v1 = 2.0 * mx_s.at(4 * i) / vx_s.at(4 * i);
      mm_v2 = 2.0 * mx_s.at(4 * i + 1) / vx_s.at(4 * i + 1);
      mm_v3 = 2.0 * mx_s.at(4 * i + 2) / vx_s.at(4 * i + 2);
      mm_v4 = mm_v_one;
    }else{
      mm_v0 = 2.0 * mx_s.at(4 * i - 1) / vx_s.at(4 * i - 1);
      mm_v1 = 2.0 * mx_s.at(4 * i) / vx_s.at(4 * i);
      mm_v2 = 2.0 * mx_s.at(4 * i + 1) / vx_s.at(4 * i + 1);
      mm_v3 = 2.0 * mx_s.at(4 * i + 2) / vx_s.at(4 * i + 2);
      mm_v4 = 2.0 * mx_s.at(4 * i + 3) / vx_s.at(4 * i + 3);
    }

    double log_g0 = -int_mm_v;
    int_mm_v += h1 / 6.0 * (mm_v0 + 4.0 * mm_v1 + mm_v2);
    double log_g2 = -int_mm_v;
    int_mm_v += h1 / 6.0 * (mm_v2 + 4.0 * mm_v3 + mm_v4);
    double log_g4 = -int_mm_v;
    
    {
      double max_log = std::max({log_int_g, log_g0, log_g2, log_g4});
      log_int_g = max_log + std::log(
        std::exp(log_int_g - max_log) + h2 / 6.0 * (
          std::exp(log_g0 - max_log) + 4.0 * std::exp(log_g2 - max_log) + std::exp(log_g4 - max_log)
        )
      );
    }
    
    for(auto& j: log_int_gys_y_1){
      double max_log = std::max({j, log_g0, log_g2, log_g4});

      j = max_log + std::log(
        std::exp(j - max_log) + h2 / 6.0 * (
          std::exp(log_g0 - max_log) + 4.0 * std::exp(log_g2 - max_log) + std::exp(log_g4 - max_log)
        )
      );
    }

    if(i < sep - 1){
      ys.push_back(y4);
      vys.push_back(vx_s.at(4 * i + 3));
      log_gys.push_back(log_g4);
      log_int_gys_0_y.push_back(log_int_g);
      log_int_gys_y_1.push_back(std::log(0.0));
      tmps.push_back(mx_s.at(4 * i + 3));
    }
  }

  double ux, vx;

  if(x0 < ys.at(0)){
    double x_lower = 0.0;
    double x_upper = ys.at(0);

    double u_lower = 0.0;
    double u_upper = std::exp(log_int_gys_0_y.at(0) - log_int_g);

    double v_lower = 1.0;
    double v_upper = std::exp(log_int_gys_y_1.at(0) - log_int_g);

    ux = u_lower + (x0 - x_lower) / (x_upper - x_lower) * (u_upper - u_lower);
    vx = v_lower + (x0 - x_lower) / (x_upper - x_lower) * (v_upper - v_lower);
  }else if(x0 > ys.back()){
    double x_lower = ys.back();
    double x_upper = 1.0;

    double u_lower = std::exp(log_int_gys_0_y.back() - log_int_g);
    double u_upper = 1.0;

    double v_lower = std::exp(log_int_gys_y_1.back() - log_int_g);
    double v_upper = 0.0;

    ux = u_lower + (x0 - x_lower) / (x_upper - x_lower) * (u_upper - u_lower);
    vx = v_lower + (x0 - x_lower) / (x_upper - x_lower) * (v_upper - v_lower);
  }else{
    int index = 0;
    while(ys.at(index) < x0){
      index++;
    }

    double x_lower = ys.at(index - 1);
    double x_upper = ys.at(index);

    double u_lower = std::exp(log_int_gys_0_y.at(index - 1) - log_int_g);
    double u_upper = std::exp(log_int_gys_0_y.at(index) - log_int_g);

    double v_lower = std::exp(log_int_gys_y_1.at(index - 1) - log_int_g);
    double v_upper = std::exp(log_int_gys_y_1.at(index) - log_int_g);

    ux = u_lower + (x0 - x_lower) / (x_upper - x_lower) * (u_upper - u_lower);
    vx = v_lower + (x0 - x_lower) / (x_upper - x_lower) * (v_upper - v_lower);
  }


  for(int i = 0; i < static_cast<int>(ys.size()); i++){
    if(x0 > ys.at(i)){
      double ret = 2.0 * vx / vys.at(i) * std::exp(log_int_gys_0_y.at(i) - log_gys.at(i));
      sum_sojourn += ret / sep;
    }else{
      double ret = 2.0 * ux / vys.at(i) * std::exp(log_int_gys_y_1.at(i) - log_gys.at(i));
      sum_sojourn += ret / sep;
    }
  }
}

double Diffusion::eigen_at_zero(){
  int pop_num = static_cast<int>(para.ss.size());

  // calculate leading eigenvector
  Eigen::MatrixXd MM(pop_num, pop_num);

  if(pop_num == 1){
    MM(0, 0) = para.ss.at(0);
  }else{
    for(int i = 0; i < pop_num; i++){
      for(int j = 0; j < pop_num; j++){
        if(i == j){
          MM(i, i) = para.ss.at(i) + para.ms.at(i).at(i);
        }else{
          MM(i, j) = para.ms.at(i).at(j);
        }
      }
    }
  }

  Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(MM);
  double max = eigensolver.eigenvalues()(0).real();
  int max_index = 0;
  std::vector<double> eigen_values;

  for(int j = 1; j < pop_num; j++){
    if(max < eigensolver.eigenvalues()(j).real()){
      max = eigensolver.eigenvalues()(j).real();
      max_index = j;
    }
  }

  Eigen::VectorXd L = eigensolver.eigenvectors().col(max_index).real();
  L = L.normalized();

  // constitute normalized orthogonal basis
  Eigen::MatrixXd U(pop_num, pop_num);
  U.col(0) = L;
  for(int i = 1; i < pop_num; i++){
    Eigen::VectorXd Vi(pop_num);
    Vi.setZero();
    Vi(i) = 1;
    for(int j = 0; j < i; j++){
      double inner_product_UV = U.col(j).dot(Vi);
      Vi -= inner_product_UV * U.col(j);
    }
    Vi = Vi.normalized();
    U.col(i) = Vi;
  }

  // constitute approximate space in the standard coordinates
  Eigen::MatrixXd MM2(pop_num, pop_num);

  if(pop_num == 1){
    MM2(0, 0) = (1.0 - L(0) * L(0)) * para.ss.at(0);
  }else{
    for(int i = 0; i < pop_num; i++){
      for(int j = 0; j < pop_num; j++){
        if(i == j){
          double m_sum = 0.0;
          double llm_sum = 0.0;
          for(int k = 0; k < pop_num; k++){
            if(k != i){
              m_sum += para.ms.at(i).at(k);
              llm_sum += L(i) * L(k) * para.ms.at(k).at(i);
            }
          }
          MM2(i, i) = (1.0 - L(i) * L(i)) * 
            (para.ss.at(i) - m_sum) - llm_sum;
        }else{
          double m_sum = 0.0;
          double llm_sum = 0.0;

          for(int k = 0; k < pop_num; k++){
            if(k != j){
              m_sum += para.ms.at(j).at(k);
            }
            if(k != i && k != j){
              llm_sum += L(i) * L(k) * para.ms.at(k).at(j);
            }
          }
          MM2(i, j) = (1.0 - L(i) * L(i)) * para.ms.at(i).at(j) - 
            L(i) * L(j) * para.ss.at(j) + L(i) * L(j) * m_sum - llm_sum;
        }
      }
    }
  }

  // change basis
  Eigen::MatrixXd MM3(pop_num, pop_num);
  MM3 = U.transpose() * MM2 * U;

  Eigen::MatrixXd MM4(pop_num - 1, pop_num - 1);
  for(int i = 1; i < pop_num; i++){
    for(int j = 1; j < pop_num; j++){
      MM4(i - 1, j - 1) = MM3(i, j);
    }
  }

  // calculate eigenvalues
  Eigen::EigenSolver<Eigen::MatrixXd> eigensolver2(MM4);
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> P_tmp(pop_num - 1, pop_num - 1);
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> P(pop_num, pop_num);
  P.setZero();
  P_tmp = eigensolver2.eigenvectors();

  P(0, 0) = 1.0;
  for(int i = 0; i < pop_num - 1; i++){
    for(int j = 0; j < pop_num - 1; j++){
      P(i + 1, j + 1) = P_tmp(i, j);
    }
  }

  for(int i = 0; i < pop_num - 1; i++){
    if(eigensolver2.eigenvalues()(i).real() > 0.0){
      convergence = 0;
    }
  }

  // change basis in the direction of eigen vector
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> U2(pop_num, pop_num);
  U2 = U * P;

  std::vector<std::complex<double>> alphas = {U2.col(0).dot(L)};
  for(int i = 1; i < pop_num; i++){
    std::complex<double> sum = 0.0;
    for(int j = 0; j < pop_num; j++){
      std::complex<double> mig_sum = 0.0;
      for(int k = 0; k < pop_num; k++){
        if(j != k){
          mig_sum += para.ms.at(j).at(k) * (U2(k, i) - U2(j, i));
        }
      }
      sum += L(j) * (para.ss.at(j) * U2(j, i) + mig_sum);
    }
    alphas.push_back(sum);
  }

  std::vector<std::complex<double>> weights_complex(pop_num);
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> U2_inv(pop_num, pop_num);
  U2_inv = U2.inverse();

  for(int i = 0; i < pop_num; i++){
    Eigen::Vector<std::complex<double>, Eigen::Dynamic> Tmp(pop_num);
    Tmp.setZero();
    Tmp(i) = 1.0;
    
    Eigen::Vector<std::complex<double>, Eigen::Dynamic> A(pop_num);
    A = U2_inv * Tmp;

    for(int j = 0; j < pop_num; j++){
      if(j != 0){
        weights_complex.at(i) -= A(j) * alphas.at(j) / eigensolver2.eigenvalues()(j - 1);
      }else{
        weights_complex.at(i) += A(j) * alphas.at(0);
      }
    }
  }

  std::vector<double> weights(pop_num);
  for(int i = 0; i < pop_num; i++){
    weights.at(i) = weights_complex.at(i).real();
  }

  int total_pop_size = 0;
  for(int i = 0; i < pop_num; i++){
    total_pop_size += para.pop_sizes.at(i);
  }

  double sum_l = 0.0;
  for(int i = 0; i < pop_num; i++){
    sum_l += L(i) * para.pop_sizes.at(i);
  }

  double mx = max;
  double vx = 0.0;

  for(int i = 0; i < pop_num; i++){
    vx += total_pop_size * L(i) / sum_l / para.pop_sizes.at(i) * weights.at(i) * weights.at(i) *
      sum_l * sum_l / total_pop_size / total_pop_size;
  }

  return(2.0 * mx / vx);
}

double Diffusion::eigen_at_one(){
  int pop_num = static_cast<int>(para.ss.size());

  // calculate leading eigenvector
  Eigen::MatrixXd MM(pop_num, pop_num);

  if(pop_num == 1){
    MM(0, 0) = -para.ss.at(0);
  }else{
    for(int i = 0; i < pop_num; i++){
      for(int j = 0; j < pop_num; j++){
        if(i == j){
          MM(i, i) = -para.ss.at(i) + para.ms.at(i).at(i);
        }else{
          MM(i, j) = para.ms.at(i).at(j);
        }
      }
    }
  }

  Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(MM);
  double max = eigensolver.eigenvalues()(0).real();
  int max_index = 0;
  std::vector<double> eigen_values;

  for(int j = 1; j < pop_num; j++){
    if(max < eigensolver.eigenvalues()(j).real()){
      max = eigensolver.eigenvalues()(j).real();
      max_index = j;
    }
  }

  Eigen::VectorXd L = eigensolver.eigenvectors().col(max_index).real();
  L = L.normalized();

  // constitute normalized orthogonal basis
  Eigen::MatrixXd U(pop_num, pop_num);
  U.col(0) = L;
  for(int i = 1; i < pop_num; i++){
    Eigen::VectorXd Vi(pop_num);
    Vi.setZero();
    Vi(i) = 1;
    for(int j = 0; j < i; j++){
      double inner_product_UV = U.col(j).dot(Vi);
      Vi -= inner_product_UV * U.col(j);
    }
    Vi = Vi.normalized();
    U.col(i) = Vi;
  }

  // constitute approximate space in the standard coordinates
  Eigen::MatrixXd MM2(pop_num, pop_num);

  if(pop_num == 1){
    MM2(0, 0) = (1.0 - L(0) * L(0)) * (-para.ss.at(0));
  }else{
    for(int i = 0; i < pop_num; i++){
      for(int j = 0; j < pop_num; j++){
        if(i == j){
          double m_sum = 0.0;
          double llm_sum = 0.0;
          for(int k = 0; k < pop_num; k++){
            if(k != i){
              m_sum += para.ms.at(i).at(k);
              llm_sum += L(i) * L(k) * para.ms.at(k).at(i);
            }
          }
          MM2(i, i) = (1.0 - L(i) * L(i)) * 
            (-para.ss.at(i) - m_sum) - llm_sum;
        }else{
          double m_sum = 0.0;
          double llm_sum = 0.0;

          for(int k = 0; k < pop_num; k++){
            if(k != j){
              m_sum += para.ms.at(j).at(k);
            }
            if(k != i && k != j){
              llm_sum += L(i) * L(k) * para.ms.at(k).at(j);
            }
          }
          MM2(i, j) = (1.0 - L(i) * L(i)) * para.ms.at(i).at(j) - 
            L(i) * L(j) * (-para.ss.at(j)) + L(i) * L(j) * m_sum - llm_sum;
        }
      }
    }
  }

  // change basis
  Eigen::MatrixXd MM3(pop_num, pop_num);
  MM3 = U.transpose() * MM2 * U;

  Eigen::MatrixXd MM4(pop_num - 1, pop_num - 1);
  for(int i = 1; i < pop_num; i++){
    for(int j = 1; j < pop_num; j++){
      MM4(i - 1, j - 1) = MM3(i, j);
    }
  }

  // calculate eigenvalues
  Eigen::EigenSolver<Eigen::MatrixXd> eigensolver2(MM4);
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> P_tmp(pop_num - 1, pop_num - 1);
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> P(pop_num, pop_num);
  P.setZero();
  P_tmp = eigensolver2.eigenvectors();

  P(0, 0) = 1.0;
  for(int i = 0; i < pop_num - 1; i++){
    for(int j = 0; j < pop_num - 1; j++){
      P(i + 1, j + 1) = P_tmp(i, j);
    }
  }

  for(int i = 0; i < pop_num - 1; i++){
    if(eigensolver2.eigenvalues()(i).real() > 0.0){
      convergence = 0;
    }
  }

  // change basis in the direction of eigen vector
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> U2(pop_num, pop_num);
  U2 = U * P;

  std::vector<std::complex<double>> alphas = {U2.col(0).dot(L)};
  for(int i = 1; i < pop_num; i++){
    std::complex<double> sum = 0.0;
    for(int j = 0; j < pop_num; j++){
      std::complex<double> mig_sum = 0.0;
      for(int k = 0; k < pop_num; k++){
        if(j != k){
          mig_sum += para.ms.at(j).at(k) * (U2(k, i) - U2(j, i));
        }
      }
      sum += L(j) * (-para.ss.at(j) * U2(j, i) + mig_sum);
    }
    alphas.push_back(sum);
  }

  std::vector<std::complex<double>> weights_complex(pop_num);
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> U2_inv(pop_num, pop_num);
  U2_inv = U2.inverse();

  for(int i = 0; i < pop_num; i++){
    Eigen::Vector<std::complex<double>, Eigen::Dynamic> Tmp(pop_num);
    Tmp.setZero();
    Tmp(i) = 1.0;
    
    Eigen::Vector<std::complex<double>, Eigen::Dynamic> A(pop_num);
    A = U2_inv * Tmp;

    for(int j = 0; j < pop_num; j++){
      if(j != 0){
        weights_complex.at(i) -= A(j) * alphas.at(j) / eigensolver2.eigenvalues()(j - 1);
      }else{
        weights_complex.at(i) += A(j) * alphas.at(0);
      }
    }
  }

  std::vector<double> weights(pop_num);
  for(int i = 0; i < pop_num; i++){
    weights.at(i) = weights_complex.at(i).real();
  }

  int total_pop_size = 0;
  for(int i = 0; i < pop_num; i++){
    total_pop_size += para.pop_sizes.at(i);
  }

  double sum_l = 0.0;
  for(int i = 0; i < pop_num; i++){
    sum_l += L(i) * para.pop_sizes.at(i);
  }

  double mx = max;
  double vx = 0.0;

  for(int i = 0; i < pop_num; i++){
    vx += total_pop_size * L(i) / sum_l / para.pop_sizes.at(i) * weights.at(i) * weights.at(i) *
      sum_l * sum_l / total_pop_size / total_pop_size;
  }

  return(-2.0 * mx / vx);
}

void Diffusion::eigen_at_x(std::vector<double>& xs,
  std::vector<double>& mx_s, std::vector<double>& vx_s, 
  std::vector<double>& l1_length, std::vector<double>& list_eigen){

  xs.clear();
  mx_s.clear();
  vx_s.clear();
  l1_length.clear();
  list_eigen.clear();

  int pop_num = static_cast<int>(para.ss.size());

  for(int ite = 1; ite < static_cast<int>(ave_freq.size()) - 1; ite++){
    // calculate evolutionary direction
    Eigen::VectorXd L(pop_num);

    for(int i = 0; i < pop_num; i++){
      L(i) = para.ss.at(i) * freq_each.at(ite).at(i) * (1.0 - freq_each.at(ite).at(i));

      for(int j = 0; j < pop_num; j++){
        L(i) += para.ms.at(i).at(j) * freq_each.at(ite).at(j);
      }
    }

    l1_length.push_back(L.norm());

    L = L.normalized();

    // constitute normalized orthogonal basis
    Eigen::MatrixXd U(pop_num, pop_num);
    U.col(0) = L;
    for(int i = 1; i < pop_num; i++){
      Eigen::VectorXd Vi(pop_num);
      Vi.setZero();
      Vi(i) = 1;
      for(int j = 0; j < i; j++){
        double inner_product_UV = U.col(j).dot(Vi);
        Vi -= inner_product_UV * U.col(j);
      }
      Vi = Vi.normalized();
      U.col(i) = Vi;
    }

    // constitute approximate space in the standard coordinates
    Eigen::MatrixXd MM2(pop_num, pop_num);

    if(pop_num == 1){
      MM2(0, 0) = (1.0 - L(0) * L(0)) * para.ss.at(0) * 
        (1.0 - 2.0 * freq_each.at(ite).at(0));
    }else{
      for(int i = 0; i < pop_num; i++){
        for(int j = 0; j < pop_num; j++){
          if(i == j){
            double m_sum = 0.0;
            double llm_sum = 0.0;
            for(int k = 0; k < pop_num; k++){
              if(k != i){
                m_sum += para.ms.at(i).at(k);
                llm_sum += L(i) * L(k) * para.ms.at(k).at(i);
              }
            }
            MM2(i, i) = (1.0 - L(i) * L(i)) * 
              (para.ss.at(i) * (1.0 - 2.0 * freq_each.at(ite).at(i)) - m_sum) - llm_sum;
          }else{
            double m_sum = 0.0;
            double llm_sum = 0.0;

            for(int k = 0; k < pop_num; k++){
              if(k != j){
                m_sum += para.ms.at(j).at(k);
              }
              if(k != i && k != j){
                llm_sum += L(i) * L(k) * para.ms.at(k).at(j);
              }
            }
            MM2(i, j) = (1.0 - L(i) * L(i)) * para.ms.at(i).at(j) - 
              L(i) * L(j) * para.ss.at(j) * (1.0 - 2.0 * freq_each.at(ite).at(j)) + L(i) * L(j) * m_sum - llm_sum;
          }
        }
      }
    }

    // change basis
    Eigen::MatrixXd MM3(pop_num, pop_num);
    MM3 = U.transpose() * MM2 * U;

    Eigen::MatrixXd MM4(pop_num - 1, pop_num - 1);
    for(int i = 1; i < pop_num; i++){
      for(int j = 1; j < pop_num; j++){
        MM4(i - 1, j - 1) = MM3(i, j);
      }
    }

    // calculate eigenvalues
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver2(MM4);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> P_tmp(pop_num - 1, pop_num - 1);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> P(pop_num, pop_num);
    P.setZero();
    P_tmp = eigensolver2.eigenvectors();

    P(0, 0) = 1.0;
    for(int i = 0; i < pop_num - 1; i++){
      for(int j = 0; j < pop_num - 1; j++){
        P(i + 1, j + 1) = P_tmp(i, j);
      }
    }

    list_eigen.push_back(eigensolver2.eigenvalues()(0).real());
    for(int i = 1; i < pop_num - 1; i++){
      if(list_eigen.back() < eigensolver2.eigenvalues()(i).real()){
        list_eigen.back() = eigensolver2.eigenvalues()(i).real();
      }
    }

    // change basis in the direction of eigen vector
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> U2(pop_num, pop_num);
    U2 = U * P;

    std::vector<std::complex<double>> alphas = {U2.col(0).dot(L)};
    for(int i = 1; i < pop_num; i++){
      std::complex<double> sum = 0.0;
      for(int j = 0; j < pop_num; j++){
        std::complex<double> mig_sum = 0.0;
        for(int k = 0; k < pop_num; k++){
          if(j != k){
            mig_sum += para.ms.at(j).at(k) * (U2(k, i) - U2(j, i));
          }
        }
        sum += L(j) * (para.ss.at(j) * (1.0 - 2.0 * freq_each.at(ite).at(j)) * U2(j, i) + mig_sum);
      }
      alphas.push_back(sum);
    }

    std::vector<std::complex<double>> weights_complex(pop_num);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> U2_inv(pop_num, pop_num);
    U2_inv = U2.inverse();

    for(int i = 0; i < pop_num; i++){
      Eigen::Vector<std::complex<double>, Eigen::Dynamic> Tmp(pop_num);
      Tmp.setZero();
      Tmp(i) = 1.0;

      Eigen::Vector<std::complex<double>, Eigen::Dynamic> A(pop_num);
      A = U2_inv * Tmp;

      for(int j = 0; j < pop_num; j++){
        if(j != 0){
          weights_complex.at(i) -= A(j) * alphas.at(j) / eigensolver2.eigenvalues()(j - 1);
        }else{
          weights_complex.at(i) += A(j) * alphas.at(0);
        }
      }
    }

    std::vector<double> weights(pop_num);
    for(int i = 0; i < pop_num; i++){
      weights.at(i) = weights_complex.at(i).real();

      if(weights_complex.at(i).imag() > 1e-6){
        std::cout << "Large imaginary" << std::endl;
      }
    }

    double sum_l = 0.0;
    for(int i = 0; i < pop_num; i++){
      sum_l += L(i) * para.pop_sizes.at(i);
    }

    int total_pop_size = 0;
    for(int i = 0; i < pop_num; i++){
      total_pop_size += para.pop_sizes.at(i);
    }

    double mx = 0.0;
    for(int i = 0; i < pop_num; i++){
      mx += para.ss.at(i) * freq_each.at(ite).at(i) * (1.0 - freq_each.at(ite).at(i)) * para.pop_sizes.at(i);

      for(int j = 0; j < pop_num; j++){
        mx += para.ms.at(i).at(j) * freq_each.at(ite).at(j) * para.pop_sizes.at(i);
      }
    }
    mx /= total_pop_size;

    double vx = 0.0;

    for(int i = 0; i < pop_num; i++){
      vx += freq_each.at(ite).at(i) * (1.0 - freq_each.at(ite).at(i)) / para.pop_sizes.at(i) * 
        weights.at(i) * weights.at(i) * sum_l * sum_l / total_pop_size / total_pop_size;
    }

    xs.push_back(ave_freq.at(ite));
    mx_s.push_back(mx);
    vx_s.push_back(vx);
  }
}