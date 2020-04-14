#!/bin/bash

cp /work3/artirain/Cloud_seeding/*.ncl .
sleep 3
ncl cal_cloud.ncl > info1.txt &
ncl cal_CTH.ncl > info2.txt &
for job in `jobs -p`
do
    echo "Wait job: ${job}"
    #wait $job
done
wait
ncl assess.ncl > info3.txt &
PID=$! ; wait $PID
echo "assess has done!!"
ncl plt_score_series.ncl ; wait

FileName=ser_fig.m
cp /work3/artirain/Cloud_seeding/$FileName .
sleep 1 ; chmod +x $FileName 
matlab -nodesktop -nosplash -nojvm -r "run ./$FileName;quit;"&

