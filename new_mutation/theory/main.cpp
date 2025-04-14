#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include <random>
#include "parameter.hpp"
#include "detsimu.hpp"
#include "diffusion.hpp"

int main(){
  // population into which a new mutation is introdiced
  int ini_pop = 0;

  // number of the points in discretization
  int sep = 1000;

  // specify parameters as a vector (could be arbitrary dimension)
  // DO NOT assume neutral cases of si = 0.0 for all i
  double s1 = 0.1;
  double s2 = -0.085;
  int n1 = 2 * 100;
  int n2 = 2 * 100;
  double m1 = 0.002;
  double m2 = 0.002;

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
    // invasion fitness of allele A
    double inv1 = det.ret_invasion1();

    if(inv1 < 1e-8){
      // if allele A is not invasive
      std::cout << "allele A is not invasive" << std::endl;
    }else{
      // calculate diffusion coefficients
      Diffusion diff(para, ini_pop, sep, ave_freq, freq_each);
      diff.calculate_diffusion();

      bool invasion = diff.return_invasion();
      double sum_sojourn = diff.return_sum_sojourn();
      bool convergence = diff.return_convergence();
      double max_eigen = diff.return_max_eigen();

      std::cout << "invasive or not   : " << invasion << std::endl;
      std::cout << "applicable or not : " << convergence << std::endl;
      std::cout << "max lambda        : " << max_eigen << std::endl;
      std::cout << "absorption time   : " << sum_sojourn << std::endl;
    }
  }

  return(0);
}
