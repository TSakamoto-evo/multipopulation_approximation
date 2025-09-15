#ifndef STOCSIMU
#define STOCSIMU

#include <fstream>
#include <vector>
#include <random>
#include <boost/random.hpp>
#include <string>
#include <filesystem>
#include <sstream>
#include <chrono>
#include "parameter.hpp"

class Stocsimu{
private:
  Parameter para;
  int ini_pop;
  int rep;

  int fix1;
  int fix2;
  double total_seg_time1;
  double total_seg_time2;
  double total_seg_time1_sq;
  double total_seg_time2_sq;
  double total_seg_time;
  double total_seg_time_sq;
  int run_num;
  bool too_long;
  long long int gen;
  std::vector<int> allele_num;


public:
  Stocsimu(const Parameter input_para, const int input_ini_pop, const int input_rep);
  void run_simulation();
};

#endif