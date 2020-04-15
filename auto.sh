#!/bin/bash

cp /work3/artirain/Cloud_seeding/*.ncl .
cp /work3/artirain/Cloud_seeding/*.m .
sleep 3 ; chmod +x *.m
ncl plt_rain.ncl &
PID=$! 
ncl cal_cloud.ncl > info1.txt &
ncl cal_CTH.ncl > info2.txt &
for job in `jobs -p`
do
    echo "Wait job: ${job}"
    #wait $job
done
wait $PID
matlab -nodesktop -nosplash -nojvm -r "run ./rain_fig.m;quit;"&

wait #cal_cloud.ncl & cal_CTH.ncl
ncl assess.ncl > info3.txt & ; wait
echo "assess has done!!"
ncl plt_score_series.ncl ; wait

FileName=ser_fig.m
#cp /work3/artirain/Cloud_seeding/$FileName .
matlab -nodesktop -nosplash -nojvm -r "run ./$FileName;quit;"&

