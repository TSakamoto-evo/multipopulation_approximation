#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include "parameter.hpp"
#include "stocsimu.hpp"

int main(int argc, char* argv[]){
  double s1, s2, m;
  int rep;

  if(argc == 5){
    sscanf(argv[1], "%lf", &s1);
    sscanf(argv[2], "%lf", &s2);
    sscanf(argv[3], "%lf", &m);
    sscanf(argv[4], "%d", &rep);
  }else{
    std::cerr << "error" << std::endl;
    std::exit(1);
  }

  int pop_size = 200;
  int ini_pop = 0;

  std::vector<int> pop_sizes = {pop_size, pop_size};
  std::vector<double> ss = {s1, s2};
  std::vector<std::vector<double>> ms(2, std::vector<double>(2, 0.0));
  ms.at(0).at(0) = -m;
  ms.at(0).at(1) = m;
  ms.at(1).at(0) = m;
  ms.at(1).at(1) = -m;

  Parameter para;
  para.ss = ss;
  para.ms = ms;
  para.pop_sizes = pop_sizes;

  Stocsimu stoc(para, ini_pop, rep);
  stoc.run_simulation();
}
