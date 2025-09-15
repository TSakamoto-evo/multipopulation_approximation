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
  int pop_size = 200;
  double m = 0.01;

  int ini_pop = 0;
  int sep = 1000;

  std::ofstream ofs1("segregation_time_theory.txt");
  ofs1 << "s1\ts2\tm\tinvasion\tconvergence\teig\tx0\tT\tu1\tu2\tT1\tT2\tT3\tS" << std::endl;

  std::ofstream ofs2("error_theory.txt");
  ofs2 << "s1\ts2\tm" << std::endl;

  for(int rep1 = -20; rep1 <= 20; rep1++){
    for(int rep2 = -20; rep2 <= 20; rep2++){
      double s1 = 0.05 * rep1 / 20;
      double s2 = 0.05 * rep2 / 20;

      if(rep1 == 0 && rep2 == 0){
        ofs1 << s1 << "\t" << s2 << "\t" << m << "\t" << 0 << 
          "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" << std::endl;
        continue;
      }

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

      std::vector<double> ave_freq;
      std::vector<std::vector<double>> freq_each;

      Detsimu det(para);

      double conv = 1.0;
      double conv_crit = 1.0;
      while(conv > 1.5 / 4.0 / sep){
        conv = det.calculate_det_trajectory(4 * sep, ave_freq, freq_each, conv_crit);
        conv_crit /= 10.0;
      }

      if(conv < -0.5){
        ofs2 << s1 << "\t" << s2 << "\t" << m << std::endl;
      }else{
        double inv1 = det.ret_invasion1();
        double inv2 = det.ret_invasion2();
        double equ_ave_freq = det.return_fin_ave_freq();
        std::vector<double> equ_freq = det.return_fin_freq();

        if(inv1 < 1e-8){
          ofs1 << s1 << "\t" << s2 << "\t" << m << "\t" << 0 << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" << std::endl;
        }else{
          Diffusion diff(para, ini_pop, sep, ave_freq, freq_each);
          diff.calculate_diffusion();

          bool invasion = diff.return_invasion();
          double sum_sojourn = diff.return_sum_sojourn();
          bool convergence = diff.return_convergence();
          double max_eigen = diff.return_max_eigen();
          double x0 = diff.return_x0();

          double est_prob = diff.return_est_prob();
          double fix_prob = diff.return_fix_prob();

          double cond_time1 = diff.return_cond_time1();
          double cond_time2 = diff.return_cond_time2();
          double abs_time = diff.return_abs_time();
          double var_abs_time = diff.return_var_abs_time();

          ofs1 << s1 << "\t" << s2 << "\t" << m << "\t"
            << invasion << "\t" << convergence << "\t"
            << max_eigen << "\t" << x0 << "\t" << sum_sojourn << "\t"
            << est_prob << "\t" << fix_prob << "\t"
            << cond_time1 << "\t" << cond_time2 << "\t"
            << abs_time << "\t" << var_abs_time << std::endl;
        }
      }
    }
  }

  return(0);
}
