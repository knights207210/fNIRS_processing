function HRF=BlockAverage(start_plot, end_plot,hb_data,s)
% start_plot=Plot Start Time (in seconds) relative to stimulus onset - most
% likely, this value should be negative
% end_plot=Plot End Time (in seconds)relative to stimulus onset
% hb_data= data from ntbx converted into oxy and deoxy hemoglobin values


jobs = nirs.modules.Run_HOMER2();
jobs.fcn='hmrBlockAvg';  
job.keepoutputs = true; 

jobs.vars.trange=[start_plot end_plot];

HRF = jobs.run(hb_data);

[t,z]=size(HRF.probe.link);
[tt,zz,dd]=size(HRF.data);
time=linspace(start_plot,end_plot,tt);
time=time';
t=t/2;
num=sqrt(t);
num=round(num);
m=num^2;
r=t-m;
i=1;
if r>0
    num2=num+1;
else 
    num2=num;
end 
p=1;
q=1;
[c,k]=size(hb_data.stimulus.keys);
y=1;
optodes=table2array(HRF.probe.link(:,1:2));

mkdir '/Users/yichen/Downloads/fNIRS/untitled folder/fnirsdata/' BlockAverage_rawdata
root_dir = '/Users/yichen/Downloads/fNIRS/untitled folder/fnirsdata/'
while(y<(k+1))
figure 
hold on
suptitle(hb_data.stimulus.keys(1,y))
while(i<=t)
        subplot(num2,num,p) 
        hold on
        plot(time,HRF.data(:,q,y),'r');
        plot(time,HRF.data(:,(q+1),y),'c');
        title([ 'S',num2str(optodes(q,1)),' D',num2str(optodes(q,2))])
        if y < 2
            ta_hbo = table(time, HRF.data(:,q,y));
            ta_hbr = table(time, HRF.data(:,(q+1),y));
            writetable(ta_hbo,[root_dir,'/BlockAverage_rawdata/sub',num2str(s),'_hbo_infantCry_source',num2str(optodes(q,1)),'_detector',num2str(optodes(q,2)),'.txt'], 'Delimiter', ' ');
            writetable(ta_hbr,[root_dir,'/BlockAverage_rawdata/sub',num2str(s),'_hbr_infantCry_source',num2str(optodes(q,1)),'_detector',num2str(optodes(q,2)),'.txt'], 'Delimiter', ' ');
        else
            ta_hbo = table(time, HRF.data(:,q,y));
            ta_hbr = table(time, HRF.data(:,(q+1),y));
            writetable(ta_hbo,[root_dir,'/BlockAverage_rawdata/sub',num2str(s),'_hbo_infantNoise_source',num2str(optodes(q,1)),'_detector',num2str(optodes(q,2)),'.txt'], 'Delimiter', ' ');
            writetable(ta_hbr,[root_dir,'/BlockAverage_rawdata/sub',num2str(s),'_hbr_infantNoise_source',num2str(optodes(q,1)),'_detector',num2str(optodes(q,2)),'.txt'], 'Delimiter', ' ');
        end       
        hold off
        p=p+1;
        q=q+2;      
i=i+1;        
end
q=1;
p=1;
i=1;
hold off
y=y+1;
end 



end 