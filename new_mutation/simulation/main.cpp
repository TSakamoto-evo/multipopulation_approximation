#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include "parameter.hpp"
#include "stocsimu.hpp"

int main(){
  // population into which a new mutation is introdiced
  int ini_pop = 0;

  // number of replicates
  int rep = 100000;

  // specify parameters as a vector (could be arbitrary dimension)
  double s1 = 0.02;
  double s2 = -0.05;
  int n1 = 2 * 1000;
  int n2 = 2 * 1500;
  double m1 = 0.02;
  double m2 = 0.03;

  std::vector<int> pop_sizes = {n1, n2};
  std::vector<double> ss = {s1, s2};
  std::vector<std::vector<double>> ms(2, std::vector<double>(2, 0.0));
  ms.at(0).at(0) = -m1;
  ms.at(0).at(1) = m1;
  ms.at(1).at(0) = m2;
  ms.at(1).at(1) = -m2;

  Parameter para;
  para.ss = ss;
  para.ms = ms;
  para.pop_sizes = pop_sizes;

  Stocsimu stoc(para, ini_pop, rep);
  stoc.run_simulation();
}
