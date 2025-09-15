import os
import sys
import subprocess
import shutil
import glob
import time

args = sys.argv

time.sleep(10)

if os.path.isfile("../regi_mean_var.txt"):
    with open("../regi_mean_var.txt", mode="r") as f:
        for line in f:
            line = line.rstrip("\n")
            list = line.split(" ")

            if list[0] == args[1]:
                time.sleep(10)
                exit()


para_list = open("../list_para.txt", mode="r")

line_no = 1
for lines in para_list:
    if line_no == int(args[1]):
        line = lines.rstrip("\n")
        list_all = line
    line_no += 1

para_list.close()

subprocess.run("g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o test.out", shell=True)
subprocess.run("./test.out "+list_all, shell=True)

with open("count.txt", mode="r") as f:
    with open("../regi_mean_var.txt", mode="a") as g:
        for line in f:
            line = line.rstrip("\n")
            list = line.split("\t")

            run_num = int(list[0])
            
            fix1 = int(list[1])
            fix2 = int(list[2])
            total = float(list[3])
            total_sq = float(list[4])

            total1 = float(list[5])
            total1_sq = float(list[6])
            total2 = float(list[7])
            total2_sq = float(list[8])

            too_long = int(list[9])

            mean_time = total / run_num
            var_time = total_sq / run_num - mean_time ** 2
            if var_time < 0.0:
                var_time = 0.0
            else:
                var_time = var_time ** 0.5

            if fix1 > 0:
                mean_time1 = total1 / fix1
                var_time1 = total1_sq / fix1 - mean_time1 ** 2
                if var_time1 < 0.0:
                    var_time1 = 0.0
                else:
                    var_time1 = var_time1 ** 0.5
            else:
                mean_time1 = 0.0
                var_time1 = 0.0
        
            if fix2 > 0:
                mean_time2 = total2 / fix2
                var_time2 = total2_sq / fix2 - mean_time2 ** 2
                if var_time2 < 0.0:
                    var_time2 = 0.0
                else:
                    var_time2 = var_time2 ** 0.5
            else:
                mean_time2 = 0.0
                var_time2 = 0.0

            print(args[1], list_all, run_num, fix1, fix2, mean_time, var_time, 
                  mean_time1, var_time1, mean_time2, var_time2, too_long, file=g)

time.sleep(10)

exit()
