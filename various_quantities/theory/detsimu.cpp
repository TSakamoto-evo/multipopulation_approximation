#include "detsimu.hpp"

Detsimu::Detsimu(const Parameter input_para){
  para = input_para;
}

double Detsimu::calculate_det_trajectory(const int sep, 
  std::vector<double>& ave_freq, 
  std::vector<std::vector<double>>& freq_each, 
  const double conv_crit){

  int pop_num = static_cast<int>(para.ss.size());

  std::vector<double> list_ave_freq1 = {0.0};
  std::vector<std::vector<double>> list_freqs1(1, std::vector<double>(pop_num, 0.0));

  std::vector<double> list_ave_freq2 = {1.0};
  std::vector<std::vector<double>> list_freqs2(1, std::vector<double>(pop_num, 1.0));

  double convergence = 1e-6 * conv_crit;
  std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

  double fin1 = 0.0;
  double fin2 = 1.0;

  std::vector<double> fin1_freq(pop_num, 0.0);
  std::vector<double> fin2_freq(pop_num, 1.0);

  {
    Eigen::MatrixXd MM(pop_num, pop_num);

    if(pop_num == 1){
      MM(0, 0) = para.ss.at(0);
    }else{
      for(int i = 0; i < pop_num; i++){
        for(int j = 0; j < pop_num; j++){
          if(i == j){
            MM(i, i) = para.ss.at(i) + para.ms.at(i).at(i);
          }else{
            MM(i, j) = para.ms.at(i).at(j);
          }
        }
      }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(MM);
    double max = eigensolver.eigenvalues()(0).real();
    int max_index = 0;
    std::vector<double> eigen_values;

    for(int j = 1; j < pop_num; j++){
      if(max < eigensolver.eigenvalues()(j).real()){
        max = eigensolver.eigenvalues()(j).real();
        max_index = j;
      }
    }

    inv_fit1 = eigensolver.eigenvalues()(max_index).real();

    if(eigensolver.eigenvalues()(max_index).real() > 0.0){
      Eigen::VectorXd L0 = eigensolver.eigenvectors().col(max_index).real();

      int now_sep = 1;
      
      double sum = 0.0;
      int total_pop_size = 0;
      for(int i = 0; i < pop_num; i++){
        sum += L0(i) * para.pop_sizes.at(i);
        total_pop_size += para.pop_sizes.at(i);
      }

      std::vector<double> ps(pop_num);

      for(int i = 0; i < pop_num; i++){
        ps.at(i) = L0(i) / sum * total_pop_size * 1e-6;
      }

      {
        std::vector<double> p_parent(pop_num);
        double ave_p_parent = 0.0;
        double ave_p = 0.0;

        for(int i = 0; i < pop_num; i++){
          ave_p_parent += p_parent.at(i) * para.pop_sizes.at(i);
          ave_p += ps.at(i) * para.pop_sizes.at(i);
        }

        ave_p_parent /= total_pop_size;
        ave_p /= total_pop_size;

        while(now_sep < ave_p * sep){
          double x = 1.0 * now_sep / sep;

          std::vector<double> tmp;
          for(int i = 0; i < pop_num; i++){
            double p_lower = p_parent.at(i);
            double p_higher = ps.at(i);
            double p = p_lower + (p_higher - p_lower) * ((x - ave_p_parent) / (ave_p - ave_p_parent));

            tmp.push_back(p);
          }

          list_ave_freq1.push_back(x);
          list_freqs1.push_back(tmp);

          now_sep++;
        }
      }

      double change = 1.0;
      int rep_step = 0;

      while(change > convergence){
        std::vector<double> p_parent(ps);
        for(int i = 0; i < pop_num; i++){
          ps.at(i) += para.ss.at(i) * p_parent.at(i) * (1.0 - p_parent.at(i));
          for(int j = 0; j < pop_num; j++){
            ps.at(i) += para.ms.at(i).at(j) * p_parent.at(j);
          }
        }

        change = 0.0;
        for(int i = 0; i < pop_num; i++){
          change += std::abs(ps.at(i) - p_parent.at(i));
        }

        {
          double ave_p_parent = 0.0;
          double ave_p = 0.0;

          for(int i = 0; i < pop_num; i++){
            ave_p_parent += p_parent.at(i) * para.pop_sizes.at(i);
            ave_p += ps.at(i) * para.pop_sizes.at(i);
          }

          ave_p_parent /= total_pop_size;
          ave_p /= total_pop_size;

          while(now_sep < ave_p * sep){
            double x = 1.0 * now_sep / sep;

            std::vector<double> tmp;
            for(int i = 0; i < pop_num; i++){
              double p_lower = p_parent.at(i);
              double p_higher = ps.at(i);
              double p = p_lower + (p_higher - p_lower) * ((x - ave_p_parent) / (ave_p - ave_p_parent));

              tmp.push_back(p);
            }

            list_ave_freq1.push_back(x);
            list_freqs1.push_back(tmp);

            now_sep++;
          }

          fin1 = ave_p;
          fin1_freq = ps;
        }

        rep_step++;
        if(rep_step % 10000 == 0){
          std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
          auto time = now - start;
          double minute = std::chrono::duration_cast<std::chrono::minutes>(time).count();

          if(minute > 10.0){
            return(-1.0);
          }
        }
      }
    }
  }

  {
    Eigen::MatrixXd MM(pop_num, pop_num);

    if(pop_num == 1){
      MM(0, 0) = -para.ss.at(0);
    }else{
      for(int i = 0; i < pop_num; i++){
        for(int j = 0; j < pop_num; j++){
          if(i == j){
            MM(i, i) = -para.ss.at(i) + para.ms.at(i).at(i);
          }else{
            MM(i, j) = para.ms.at(i).at(j);
          }
        }
      }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(MM);
    double max = eigensolver.eigenvalues()(0).real();
    int max_index = 0;
    std::vector<double> eigen_values;

    for(int j = 1; j < pop_num; j++){
      if(max < eigensolver.eigenvalues()(j).real()){
        max = eigensolver.eigenvalues()(j).real();
        max_index = j;
      }
    }

    inv_fit2 = eigensolver.eigenvalues()(max_index).real();

    if(eigensolver.eigenvalues()(max_index).real() > 0.0){
      Eigen::VectorXd L1 = eigensolver.eigenvectors().col(max_index).real();

      int now_sep = 1;
      
      double sum = 0.0;
      int total_pop_size = 0;
      for(int i = 0; i < pop_num; i++){
        sum += L1(i) * para.pop_sizes.at(i);
        total_pop_size += para.pop_sizes.at(i);
      }

      std::vector<double> qs(pop_num);

      for(int i = 0; i < pop_num; i++){
        qs.at(i) = L1(i) / sum * total_pop_size * 1e-6;
      }
      
      {
        std::vector<double> q_parent(pop_num);
        double ave_q_parent = 0.0;
        double ave_q = 0.0;

        for(int i = 0; i < pop_num; i++){
          ave_q_parent += q_parent.at(i) * para.pop_sizes.at(i);
          ave_q += qs.at(i) * para.pop_sizes.at(i);
        }

        ave_q_parent /= total_pop_size;
        ave_q /= total_pop_size;

        while(now_sep < ave_q * sep){
          double x = 1.0 * now_sep / sep;

          std::vector<double> tmp;
          for(int i = 0; i < pop_num; i++){
            double q_lower = q_parent.at(i);
            double q_higher = qs.at(i);
            double q = q_lower + (q_higher - q_lower) * ((x - ave_q_parent) / (ave_q - ave_q_parent));

            tmp.push_back(1.0 - q);
          }

          list_ave_freq2.push_back(1.0 - x);
          list_freqs2.push_back(tmp);

          now_sep++;
        }
      }

      double change = 1.0;
      int rep_step = 0;

      while(change > convergence){
        std::vector<double> q_parent(qs);
        for(int i = 0; i < pop_num; i++){
          qs.at(i) += -para.ss.at(i) * q_parent.at(i) * (1.0 - q_parent.at(i));
          for(int j = 0; j < pop_num; j++){
            qs.at(i) += para.ms.at(i).at(j) * q_parent.at(j);
          }
        }

        change = 0.0;
        for(int i = 0; i < pop_num; i++){
          change += std::abs(qs.at(i) - q_parent.at(i));
        }

        {
          double ave_q_parent = 0.0;
          double ave_q = 0.0;

          for(int i = 0; i < pop_num; i++){
            ave_q_parent += q_parent.at(i) * para.pop_sizes.at(i);
            ave_q += qs.at(i) * para.pop_sizes.at(i);
          }

          ave_q_parent /= total_pop_size;
          ave_q /= total_pop_size;

          while(now_sep < ave_q * sep){
            double x = 1.0 * now_sep / sep;

            std::vector<double> tmp;
            for(int i = 0; i < pop_num; i++){
              double q_lower = q_parent.at(i);
              double q_higher = qs.at(i);
              double q = q_lower + (q_higher - q_lower) * ((x - ave_q_parent) / (ave_q - ave_q_parent));
              tmp.push_back(1.0 - q);
            }

            list_ave_freq2.push_back(1.0 - x);
            list_freqs2.push_back(tmp);

            now_sep++;
          }

          fin2 = 1.0 - ave_q;
          fin2_freq = qs;
          for(auto& i: fin2_freq){
            i = 1.0 - i;
          }
        }

        rep_step++;
        if(rep_step % 10000 == 0){
          std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
          auto time = now - start;
          double minute = std::chrono::duration_cast<std::chrono::minutes>(time).count();

          if(minute > 10.0){
            return(-1.0);
          }
        }
      }
    }
  }

  ave_freq.clear();
  freq_each.clear();

  double equ1 = list_ave_freq1.back();
  double equ2 = list_ave_freq2.back();

  fin_ave_freq = (fin1 + fin2) / 2.0;
  fin_freq.clear();
  for(int i = 0; i < pop_num; i++){
    double tmp = (fin1_freq.at(i) + fin2_freq.at(i)) / 2.0;
    fin_freq.push_back(tmp);
  }

  for(int i = 0; i < static_cast<int>(list_ave_freq1.size()); i++){
    ave_freq.push_back(list_ave_freq1.at(i));
    freq_each.push_back(list_freqs1.at(i));
  }

  if(std::floor(fin2 * sep) - std::floor(fin1 * sep) && fin2 - fin1 < 1e-6 &&
    (std::floor(fin1 * sep) > 0) && (std::floor(fin2 * sep) < sep)){
    double x = 1.0 * std::floor(fin2 * sep) / sep;

    std::vector<double> tmp;
    for(int i = 0; i < pop_num; i++){
      double p_lower = fin1_freq.at(i);
      double p_higher = fin2_freq.at(i);
      double p = p_lower + (p_higher - p_lower) * ((x - fin1) / (fin2 - fin1));
      tmp.push_back(p);
    }

    ave_freq.push_back(x);
    freq_each.push_back(tmp);
    equ1 = x;
  }

  std::reverse(list_ave_freq2.begin(), list_ave_freq2.end());
  std::reverse(list_freqs2.begin(), list_freqs2.end());

  for(int i = 0; i < static_cast<int>(list_ave_freq2.size()); i++){
    ave_freq.push_back(list_ave_freq2.at(i));
    freq_each.push_back(list_freqs2.at(i));
  }

  return(equ2 - equ1);
}