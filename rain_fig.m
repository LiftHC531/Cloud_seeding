clear
close all
%%
%[s,figname] = system('ls -t  *score_series.png')
dirOutput=dir(fullfile('./','*12hr_QPF*.png'));
%dirOutput.name
m = size(dirOutput,1);
%fileNames={dirOutput.name};
LL=400;
RR=500;
UU=200;
DD=30;
clear newfig;
for i = 1 : m
    fig = zeros(256,256);
    fig = imread(dirOutput(i).name);
    ss = size(fig);
    xx = ss(2)-RR-LL;
    yy = ss(1)-UU-DD;
    newfig(1:yy , xx*(i-1)+1:i*xx , :) ...
          =fig(DD+1:end-UU,LL+1:end-RR,:);
   clear fig;
   clear ss ;
end;
imwrite(newfig,['QPF_panel','.png'],'png')
