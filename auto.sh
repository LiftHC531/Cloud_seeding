#!/bin/bash
if [ 1 ]
then
   cp /work3/artirain/Cloud_seeding/*.ncl .
   cp /work3/artirain/Cloud_seeding/*.m .
   sleep 3 ; chmod +x *.m
   ncl plt_rain.ncl &
   PID=$! ; sleep 0.5
   wait $PID
   matlab -nodesktop -nosplash -nojvm -r "run ./rain_fig.m;quit;"&
   ncl cal_cloud.ncl >& info1.txt ; sleep 5 
   ncl cal_CTH.ncl >& info2.txt ; sleep 30
   for job in `jobs -p`
   do
       echo "Wait job: ${job}"
       #wait $job
   done
fi

wait #cal_cloud.ncl & cal_CTH.ncl
ncl assess.ncl >& info3.txt ; wait
echo "assess has done!!"
ncl plt_score_series.ncl ; wait

FileName=ser_fig.m
#cp /work3/artirain/Cloud_seeding/$FileName .
matlab -nodesktop -nosplash -nojvm -r "run ./$FileName;quit;"&
wait
echo "All Done!!!!! get score_series_panel.png"
