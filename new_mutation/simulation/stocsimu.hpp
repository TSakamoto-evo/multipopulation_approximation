#ifndef STOCSIMU
#define STOCSIMU

#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include "parameter.hpp"

class Stocsimu{
private:
  Parameter para;
  int ini_pop;
  int rep;

  double total_seg_time;
  double total_seg_time_sq;
  bool too_long;
  long long int gen;
  std::vector<int> allele_num;

public:
  Stocsimu(const Parameter input_para, const int input_ini_pop, const int input_rep);
  void run_simulation();
};

#endif