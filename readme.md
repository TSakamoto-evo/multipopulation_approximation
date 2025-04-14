Code for Sakamoto and Yeaman (2025).

Boost and Eigen libraries are needed. 

## Folder "new_mutation"
This folder contains the codes to calculate the absorption time starting from a single mutation. 

### theory folder
#### How to use
Please do not assume purely neutral cases.

(Step 1) Specify parameters in the main file:  
ini_pop [int]: the id of subpopulation where the new mutation arises  
sep [int]: number of discretization points in the numerical calculation  
pop_sizes [1D vector]: haploid population size of each subpopulation  
ss [1D vector]: selection coefficient on the mutant allele in each subpopulation  
ms [2D vector]: backward migration rate. ms.at(i).at(j) is $m_{ij}$ (i.e., the rate from population $j$ to $i$). ms.at(i).at(i) should equal to $-\sum_{j\ne i} m_{ij}$.

(Step 2) Compile the c++ codes by the command:
```
  g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out
```

This code first calculate the deterministic trajectory from the state where each of alleles is rare.
Although the two trajectories mostly reach the same state, in some exceptional cases, they do not reach the same state within a reasonable time mainly due to too weak selection. If so, the program outputs the error message. Otherwise, the program next check if the new mutation is invasive or not. If invasive, the evolutionary directions around the trajectories are investigated. Then, the program outputs if the method is applicable (1) or not (0), in addition to the absorption time. Although the absorption time is output also for non-applicable cases, it is out of the theoretical applicability and may be untrustable. 

#### Brief explanation for each file
```parameter.hpp``` defines the class containing all relevant parameters.  

```detsimu.hpp/.cpp``` calculates the deterministic trajectory. This function is implemented in the function calculate_det_trajectory. Starting from a initial frequency, this file calculate the frequency following:
$$p_i^* = p_i + s_i p_i (1-p_i) + \sum_{j\ne i} m_{ij}(p_j - p_i)$$
The frequency trajectory is returned after the discretization. 

```diffusion.hpp/.cpp``` calculates the diffusion coefficients. Function eigen_at_x calculates the convergence tendency $\lambda(\bar{p})$ and diffusion coefficient $M(\bar{p}), V(\bar{p})$ for each point. These values are used in function calculate_diffusion to calculate the absorption time, following diffusion theory on one population. Function fixation_prob_branching calculates the extinction probability of the mutation by using multitype branching process, which is used to determine the effective initial frequency in the function calculate_diffusion.

### simulation folder
#### How to use
(Step 1) Specify parameters in the main file:  
ini_pop [int]: the id of subpopulation where the new mutation arises  
rep [int]: number of replicates 
pop_sizes [1D vector]: haploid population size of each subpopulation  
ss [1D vector]: selection coefficient on the mutant allele in each subpopulation  
ms [2D vector]: backward migration rate. ms.at(i).at(j) is $m_{ij}$ (i.e., the rate from population $j$ to $i$). ms.at(i).at(i) should equal to $-\sum_{j\ne i} m_{ij}$.

(Step 2) Compile the c++ codes by the command:
```
  g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out
```

Each replicate ends when either of the alleles is lost from the entire population. 
The program outputs mean absoption time and standard deviation of absorption time. 
If some replicates continue very long (in the code, $T=10^9$ generations is assumed), the program stops and returns ```long_run = 1```.


## Folder "new_mutation"
This folder contains the codes to calculate the absorption time starting from polymorphic equilibrium. The codes and usage are almost similar to the new mutation case.

#### How to use
Please do not assume purely neutral cases.

(Step 1) Specify parameters in the main file:  
sep [int]: number of discretization points in the numerical calculation  
pop_sizes [1D vector]: haploid population size of each subpopulation  
ss [1D vector]: selection coefficient on the mutant allele in each subpopulation  
ms [2D vector]: backward migration rate. ms.at(i).at(j) is $m_{ij}$ (i.e., the rate from population $j$ to $i$). ms.at(i).at(i) should equal to $-\sum_{j\ne i} m_{ij}$.

(Step 2) Compile the c++ codes by the command:
```
  g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out
```

### simulation folder
#### How to use
(Step 1) Specify parameters in the main file:  
sep [int]: number of discretization points in the numerical calculation  
rep [int]: number of replicates 
pop_sizes [1D vector]: haploid population size of each subpopulation  
ss [1D vector]: selection coefficient on the mutant allele in each subpopulation  
ms [2D vector]: backward migration rate. ms.at(i).at(j) is $m_{ij}$ (i.e., the rate from population $j$ to $i$). ms.at(i).at(i) should equal to $-\sum_{j\ne i} m_{ij}$.

(Step 2) Compile the c++ codes by the command:
```
  g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out
```