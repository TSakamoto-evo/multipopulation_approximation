#ifndef DETSIMU
#define DETSIMU

#include "parameter.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <cmath>
#include <chrono>

class Detsimu{
private:
  Parameter para;
  double inv_fit1;
  double inv_fit2;
  double fin_ave_freq;
  std::vector<double> fin_freq;

public:
  Detsimu(const Parameter input_para);
  double calculate_det_trajectory(const int sep, std::vector<double>& ave_freq, 
    std::vector<std::vector<double>>& freq_each, const double conv_crit);
  double ret_invasion1() const{ return(inv_fit1); };
  double ret_invasion2() const{ return(inv_fit2); };
  double return_fin_ave_freq() const{ return(fin_ave_freq); };
  std::vector<double> return_fin_freq() const{ return(fin_freq); };
};

#endif