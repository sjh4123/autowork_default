#!/bin/bash

trial=(0 1 2 3)
namestring=(B0_30_S0_3)

cd /Users/SJH/Dropbox/UA_Research/Codes/cppWork/cppWork/simulations/recover/2015/recovertest/simulations

for i in ${trial[@]}; do
  for ii in ${namestring[@]}; do
      qsub_File=./quench/submission_files/qch_recovertest_00${i}_${ii}_Regime0.qsub
      test -e ${qsub_File}
      if [ "$?" -eq "0" ]; then
        qsub ${qsub_File}
      else
        continue
      fi
  done
done
