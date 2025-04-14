#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include "parameter.hpp"
#include "detsimu.hpp"
#include "stocsimu.hpp"

int main(){
  // number of replicates
  int rep = 100000;

  // number of the points in discretization
  int sep = 1000;

  // specify parameters as a vector (could be arbitrary dimension)
  // DO NOT assume neutral cases of si = 0.0 for all i
  double s1 = 0.02;
  double s2 = -0.05;
  int n1 = 2 * 200;
  int n2 = 2 * 500;
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


  // calculate deterministic trajectory
  std::vector<double> ave_freq;
  std::vector<std::vector<double>> freq_each;

  Detsimu det(para);

  double conv = 1.0;
  double conv_crit = 1.0;

  // repeat calculation until the two trajectories reach the same state
  while(conv > 1.5 / 4.0 / sep){
    conv = det.calculate_det_trajectory(4 * sep, ave_freq, freq_each, conv_crit);
    conv_crit /= 10.0;
  }

  if(conv < -0.5){
    // if the trajectories do not reach the same state
    std::cout << "error: fail to construct deterministic trajectory" << std::endl;
  }else{
    // invasion fitness of alleles
    double inv1 = det.ret_invasion1();
    double inv2 = det.ret_invasion2();

    // allele frequency at the internal equilibrium
    std::vector<double> equ_freq = det.return_fin_freq();

    if(inv1 < 1e-8 || inv2 < 1e-8){
      // if either of alleles A or a is not invasive
      std::cout << "either of the alleles is not invasive" << std::endl;
    }else{
      Stocsimu stoc(para, equ_freq, rep);
      stoc.run_simulation();
    }
  }
}
