#!/bin/bash

start_time=$(date +"%s") 
source activate met_env
if (( 1 )); then
   cp /work3/artirain/Cloud_seeding/*.ncl .
   cp /work3/artirain/Cloud_seeding/*.m .
   cp /work3/artirain/Cloud_seeding/PyNGL/wrf_dim_info.py .
   cp /work3/artirain/Cloud_seeding/PyNGL/cal_cloud.py .
   cp /work3/artirain/Cloud_seeding/PyNGL/cal_CTH.py .
   sleep 3 ; chmod +x *.m
   ncl plt_rain.ncl &
   PID=$! ; sleep 0.5
   wait $PID
   matlab -nodesktop -nosplash -nojvm -r "run ./rain_fig.m;quit;"&
   python cal_CTH.py &  
   python cal_cloud.py &  
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
end_time=$(date +"%s"); end_time=$[end_time-start_time] 
end_time=$(echo "scale=2;$end_time/60" | bc)
echo -e "auto.sh has done!\nTime elapsed: ${end_time} mins."
echo "All Done!!!!! get score_series_panel.png"
