#ifndef DIFFUSION
#define DIFFUSION

#include <fstream>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <numeric>
#include <vector>
#include <complex>
#include "parameter.hpp"

class Diffusion{
private:
  Parameter para;
  int sep;
  int ini_pop;
  bool invasion;
  double sum_sojourn;
  bool convergence;
  double est_prob;
  double max_eigen;
  double x0;

  double abs_time;
  double abs_time_sq;
  double cond_time1;
  double cond_time2;
  double fix_prob;

  std::vector<double> ave_freq;
  std::vector<std::vector<double>> freq_each;

public:
  Diffusion(const Parameter input_para, const int input_ini_pop, 
    const int input_sep, const std::vector<double>& input_ave_freq, 
    const std::vector<std::vector<double>>& input_freq_each);
  void calculate_diffusion();
  double eigen_at_zero();
  double eigen_at_one();
  void eigen_at_x(std::vector<double>& xs,
    std::vector<double>& mx_s, std::vector<double>& vx_s, 
    std::vector<double>& l1_length, std::vector<double>& list_eigen);
  void fixation_prob_branching(std::vector<double>& extinct_prob);

  void calculate_var_abs_time(const double x0, const std::vector<double>& xs,
    const std::vector<double>& mx_s, const std::vector<double>& vx_s);

  bool return_invasion() const{ return(invasion); };
  double return_sum_sojourn() const{ return(sum_sojourn); };
  int return_convergence() const{ return(convergence); };
  double return_est_prob() const{ return(est_prob); };
  double return_max_eigen() const{ return(max_eigen); };
  double return_x0() const{ return(x0); };

  double return_fix_prob() const{ return(fix_prob); };
  double return_cond_time1() const{ return(cond_time1); };
  double return_cond_time2() const{ return(cond_time2); };
  double return_abs_time() const{ return(abs_time); };
  double return_var_abs_time() const{ return(std::pow(abs_time_sq - abs_time * abs_time, 0.5)); };
};


#endif