#include "stocsimu.hpp"

Stocsimu::Stocsimu(const Parameter input_para, const int input_ini_pop, const int input_rep){
  para = input_para;
  ini_pop = input_ini_pop;
  rep = input_rep;

  fix1 = 0;
  fix2 = 0;
  total_seg_time = 0.0;
  total_seg_time_sq = 0.0;

  total_seg_time1 = 0.0;
  total_seg_time2 = 0.0;
  total_seg_time1_sq = 0.0;
  total_seg_time2_sq = 0.0;

  run_num = 0;
  too_long = 0;
  gen = 0;

  int pop_num = static_cast<int>(para.ss.size());
  std::vector<int> tmp(pop_num, 0);
  allele_num = tmp;
  allele_num.at(ini_pop) = 1;

  if(std::filesystem::is_regular_file("count.txt")){
    std::ifstream ifs("count.txt");
    if(!ifs){
      std::cerr << "Fail to open the file!" << std::endl;
      std::exit(1);
    }

    std::string line;
    while (getline(ifs, line)){
      std::istringstream iss(line);
      std::string tmp_list;
      std::vector<std::string> list;

      while(getline(iss, tmp_list, '\t')){
        list.push_back(tmp_list);
      }

      run_num = std::stoi(list.at(0));

      fix1 = std::stoi(list.at(1));
      fix2 = std::stoi(list.at(2));
      total_seg_time = std::stod(list.at(3));
      total_seg_time_sq = std::stod(list.at(4));

      total_seg_time1 = std::stod(list.at(5));
      total_seg_time1_sq = std::stod(list.at(6));
      total_seg_time2 = std::stod(list.at(7));
      total_seg_time2_sq = std::stod(list.at(8));

      too_long = std::stoi(list.at(9));

      gen = std::stoll(list.at(10));

      for(int i = 0; i < pop_num; i++){
        allele_num.at(i) = std::stoi(list.at(11 + i));
      }
    }
  }
}

void Stocsimu::run_simulation(){
  std::random_device seed;
  boost::random::mt19937 mt(seed());

  int pop_num = static_cast<int>(para.ss.size());
  int total_pop_size = 0;
  for(int i = 0; i < pop_num; i++){
    total_pop_size += para.pop_sizes.at(i);
  }

  if(too_long == 0){
    int start_num = run_num;
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    for(long long int i = start_num + 1; i <= rep; i++){
      int total_num = 0;
      for(const auto& j: allele_num){
        total_num += j;
      }

      while(total_num > 0 && total_num < total_pop_size){
        std::vector<int> parent(allele_num);
        total_num = 0;

        for(int j = 0; j < pop_num; j++){
          double p = 1.0 * parent[j] / para.pop_sizes.at(j);
          double expect = p + para.ss[j] * p * (1.0 - p);

          for(int k = 0; k < pop_num; k++){
            expect += para.ms[j][k] * parent[k] / para.pop_sizes.at(k);
          }

          boost::random::binomial_distribution<> det(para.pop_sizes.at(j), expect);
          allele_num[j] = det(mt);
          total_num += allele_num[j];
        }

        gen++;

        if(gen > 1e+9){
          too_long = 1;
          break;
        }

        if(gen % 100000 == 1){
          std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
          auto time = now - start;
          double minute = std::chrono::duration_cast<std::chrono::minutes>(time).count();

          if(minute > 5){
            std::ofstream ofs("count.txt");
            ofs << run_num << "\t" << fix1 << "\t" << fix2 << "\t" << total_seg_time << "\t" << total_seg_time_sq << "\t" <<
              total_seg_time1 << "\t" << total_seg_time1_sq << "\t" <<
              total_seg_time2 << "\t" << total_seg_time2_sq << "\t" << too_long << "\t" << gen;

            for(const auto& j: allele_num){
              ofs << "\t" << j;
            }

            ofs << std::endl;
            ofs.close();

            start = std::chrono::system_clock::now();
          }
        }
      }

      if(total_num == 0){
        fix1++;
        total_seg_time1 += gen;
        total_seg_time1_sq += gen * gen;
      }else if(total_num == total_pop_size){
        fix2++;
        total_seg_time2 += gen;
        total_seg_time2_sq += gen * gen;
      }

      total_seg_time += gen;
      total_seg_time_sq += gen * gen;
      run_num++;

      // initialize
      std::vector<int> tmp(pop_num, 0);
      allele_num = tmp;
      allele_num.at(ini_pop) = 1;

      gen = 0;

      if(i == rep || too_long == 1){
        std::ofstream ofs("count.txt");
        ofs << run_num << "\t" << fix1 << "\t" << fix2 << "\t" << total_seg_time << "\t" << total_seg_time_sq << "\t" <<
              total_seg_time1 << "\t" << total_seg_time1_sq << "\t" <<
              total_seg_time2 << "\t" << total_seg_time2_sq << "\t" << too_long << "\t" << gen;

        for(const auto& j: allele_num){
          ofs << "\t" << j;
        }

        ofs << std::endl;
        ofs.close();

        if(too_long == 1){
          break;
        }
      }
    }
  }
}
