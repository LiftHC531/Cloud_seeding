clear
close all
%%
%[s,figname] = system('ls -t  *score_series.png')
dirOutput=dir(fullfile('./','*score_series.png'));
m = size(dirOutput,1);
%fileNames={dirOutput.name};
LL=10;
RR=10;
UU=500;
DD=500;
clear newfig;
for i = 1 : m
    fig = zeros(256,256);
    fig = imread(dirOutput(i).name);
    ss = size(fig);
    xx = ss(2)-RR-LL;
    yy = ss(1)-UU-DD;
    newfig(yy*(i-1)+1:i*yy , 1:xx , :) ...
          =fig(DD+1:end-UU,LL+1:end-RR,:);
   clear fig;
   clear ss ;
end;
imwrite(newfig,['score_series_panel','.png'],'png')
