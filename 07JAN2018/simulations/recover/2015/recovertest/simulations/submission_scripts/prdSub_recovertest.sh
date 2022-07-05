#!/bin/bash

trial=(0 1 2 3)
namestring=(B0_30_S0_3)
Temp=("500000 490000 480000 470000 460000 450000 440000 430000" "500000 490000 480000 470000 460000 450000 440000 430000" "500000 490000 480000 470000 460000 450000 440000 430000" "500000 490000 480000 470000 460000 450000 440000 430000")

cd /Users/SJH/Dropbox/UA_Research/Codes/cppWork/cppWork/simulations/recover/2015/recovertest/simulations

for i in ${trial[@]}; do
  for ii in ${namestring[@]}; do
    for iii in ${Temp[$i]}; do
      qsub_File=./production/submission_files/prd_recovertest_00${i}_${ii}_T${iii}.qsub
      test -e ${qsub_File}
      if [ "$?" -eq "0" ]; then
        qsub ${qsub_File}
      else
        continue
      fi
    done
  done
done
