#include "stocsimu.hpp"

Stocsimu::Stocsimu(const Parameter input_para, const int input_ini_pop, const int input_rep){
  para = input_para;
  ini_pop = input_ini_pop;
  rep = input_rep;

  total_seg_time = 0.0;
  total_seg_time_sq = 0.0;
  too_long = 0;
  gen = 0;

  int pop_num = static_cast<int>(para.ss.size());
  std::vector<int> tmp(pop_num, 0);
  allele_num = tmp;
  allele_num.at(ini_pop) = 1;
}

void Stocsimu::run_simulation(){
  std::random_device seed;
  std::mt19937 mt(seed());

  int pop_num = static_cast<int>(para.ss.size());
  int total_pop_size = 0;
  for(int i = 0; i < pop_num; i++){
    total_pop_size += para.pop_sizes.at(i);
  }

  if(too_long == 0){
    for(long long int i = 1; i <= rep; i++){
      int total_num = 0;
      for(const auto& j: allele_num){
        total_num += j;
      }

      while(total_num > 0 && total_num < total_pop_size){
        // if still segregating
        std::vector<int> parent(allele_num);
        total_num = 0;

        for(int j = 0; j < pop_num; j++){
          double p = 1.0 * parent.at(j) / para.pop_sizes.at(j);
          double expect = p + para.ss.at(j) * p * (1.0 - p);

          for(int k = 0; k < pop_num; k++){
            expect += para.ms.at(j).at(k) * parent.at(k) / para.pop_sizes.at(k);
          }

          std::binomial_distribution<> det(para.pop_sizes.at(j), expect);
          allele_num.at(j) = det(mt);
          total_num += allele_num.at(j);
        }

        gen++;

        if(gen > 1e+9){
          too_long = 1;
          break;
        }
      }

      total_seg_time += gen;
      total_seg_time_sq += gen * gen;

      // initialize
      std::vector<int> tmp(pop_num, 0);
      allele_num = tmp;
      allele_num.at(ini_pop) = 1;

      gen = 0;

      if(too_long){
        break;
      }
    }
  }

  std::cout << "# of replicates          : " << rep << std::endl;
  std::cout << "mean absorption time     : " << total_seg_time / rep << std::endl;
  std::cout << "sd of absorption time    : " << std::sqrt(total_seg_time_sq / rep - (total_seg_time / rep) * (total_seg_time / rep)) << std::endl;
  std::cout << "existence of long run    : " << too_long << std::endl;
}
