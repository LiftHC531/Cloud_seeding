#!/bin/bash
NC='\033[0m' # No Color
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'

start_time=$(date +"%s") 
source activate met_env
if (( 1 )); then
   \cp /work3/artirain/Cloud_seeding/plt_rain.ncl   .
   \cp /work3/artirain/Cloud_seeding/assess.ncl     .
   \cp /work3/artirain/Cloud_seeding/plt_score_series.ncl       .
   \cp /work3/artirain/Cloud_seeding/get_forecast_sheet.ncl     .
   \cp /work3/artirain/Cloud_seeding/*.m .
   \cp /work3/artirain/Cloud_seeding/PyNGL/cal_cloud.py .
   \cp /work3/artirain/Cloud_seeding/PyNGL/cal_CTH.py   .
   sleep 1; chmod +x *.m
   ncl plt_rain.ncl &
   PID=$! ; sleep 0.5
   wait $PID
   matlab -nodesktop -nosplash -nojvm -r "run ./rain_fig.m;quit;"&
   python cal_CTH.py &  
   python cal_cloud.py &  
   for job in `jobs -p`
   do
       echo -e "${BLUE}job ID: ${job}${NC}"
       #wait $job
   done
fi

wait 
ncl assess.ncl >& info3.txt ; wait
echo "assess has done!!"
ncl plt_score_series.ncl ; wait
ncl get_forecast_sheet.ncl

FileName=ser_fig.m
#cp /work3/artirain/Cloud_seeding/$FileName .
matlab -nodesktop -nosplash -nojvm -r "run ./$FileName;quit;"&
wait
end_time=$(date +"%s"); end_time=$[end_time-start_time] 
end_time=$(echo "scale=2;$end_time/60" | bc)
echo -e "${GREEN}auto.sh has done!\nTime elapsed: ${end_time} mins.${NC}"
echo "All Done!!!!! get score_series_panel.png"
