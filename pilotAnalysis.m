clear;
clc;

%% group number
group_number = 1;

%% load data
root_dir = '/Users/yichen/Downloads/2021-05-18'
raw = nirs.io.loadDirectory(root_dir, {'group','subject','session'});

%% Create demographics table

demographics = nirs.createDemographicsTable(raw);
disp(demographics)

%% Manually add events info
temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable_remember', [21], [43], [1])
temp('open_agreeable_remember')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [73], [114], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_agreeable_rest', [200], [17], [1])
temp('open_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'open_neutral_remember', [349], [36], [1]);
temp('open_neutral_remember')=temp_se;
temp_se = nirs.design.StimulusEvents( 'open_neutral', [401], [99], [1]);
temp('open_neutral')=temp_se;
temp_se = nirs.design.StimulusEvents( 'open_neutral_rest', [513], [17], [1]);
temp('open_neutral_rest')=temp_se;

temp_se = nirs.design.StimulusEvents( 'choice_agreeable_remember', [652], [33], [1]);
temp('choice_agreeable_remember')=temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [707], [84], [1]);
temp('choice_agreeable')=temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_agreeable_rest', [798], [16], [1]);
temp('choice_agreeable_rest')=temp_se;

temp_se = nirs.design.StimulusEvents( 'choice_neutral_remember', [936], [34], [1]);
temp('choice_neutral_remember')=temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_neutral', [982], [78], [1]);
temp('choice_neutral')=temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_neutral_rest', [1067], [18], [1]);
temp('choice_neutral_rest')=temp_se;

raw(1,1).stimulus = temp



temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable_remember', [17], [37], [1])
temp('open_agreeable_remember')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [70], [121], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_agreeable_rest', [197], [16], [1])
temp('open_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'open_neutral_remember', [333], [34], [1])
temp('open_neutral_remember')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_neutral', [379], [109], [1])
temp('open_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_neutral_rest', [492], [16], [1])
temp('open_neutral_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'choice_agreeable_remember', [621], [32], [1])
temp('choice_agreeable_remember')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [664], [79], [1])
temp('choice_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable_rest', [750], [18], [1])
temp('choice_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'choice_neutral_remember', [882], [34], [1])
temp('choice_neutral_remember')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_neutral', [927], [72], [1])
temp('choice_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_neutral_rest', [1006], [17], [1])
temp('choice_neutral_rest')=temp_se

raw(2,1).stimulus = temp



temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable_remember', [17], [30], [1])
temp('open_agreeable_remember')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [58], [100], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_agreeable_rest', [168], [16], [1])
temp('open_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'open_neutral_remember', [302], [35], [1])
temp('open_neutral_remember')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_neutral', [348], [93], [1])
temp('open_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_neutral_rest', [447], [17], [1])
temp('open_neutral_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'choice_agreeable_remember', [581], [39], [1])
temp('choice_agreeable_remember')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [632], [71], [1])
temp('choice_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable_rest', [712], [20], [1])
temp('choice_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'choice_neutral_remember', [846], [31], [1])
temp('choice_neutral_remember')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_neutral', [886], [74], [1])
temp('choice_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_neutral_rest', [965], [18], [1])
temp('choice_neutral_rest')=temp_se

raw(3,1).stimulus = temp


%% Space Registration
% add an anchor to attach the probe to the center of the forehead (FpZ)
% TODO: The following code cannot be used for further visualizations yet.
% Just used for a simple demostration on how to generate 10-20
% need space registration information here
% ---------------------------------------------------------------------
%Name{1}='FpZ';
%xyz(1,:)=[0 0 0];
%Type{1}='FID-anchor';  % This is an anchor point
%Units{1}='mm';

%Now add a few more
%Name{2}='Cz';
%xyz(2,:)=[0 100 0];
%Type{2}='FID-attractor';  % This is an attractor
%Units{2}='mm';

%Name{3}='T7';
%xyz(3,:)=[-200 0 0];
%Type{3}='FID-attractor';  % This is an attractor
%Units{3}='mm';

%Name{4}='T8';
%xyz(4,:)=[200 0 0];
%Type{4}='FID-attractor';  % This is an attractor
%Units{4}='mm';

% now add these points to the optodes field in the probe. 
%fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    %'VariableNames',{'Name','X','Y','Z','Type','Units'});
%raw(1).probe.optodes=[raw(1).probe.optodes; fid];

% Use the default head size
%probe1020=nirs.util.registerprobe1020(raw(1).probe);
%probe1020.defaultdrawfcn='10-20';
%probe1020.draw;
%close;

%% Quality checking 
mkdir '/Users/yichen/Downloads/2021-05-18/' rawdata
mkdir '/Users/yichen/Downloads/2021-05-18/' qualitycheck
for s = 1:size(raw)
%for s =1:1
% view each subject's graph
    data = raw(s);
    data.draw;
    saveas(gcf,[root_dir,'/rawdata/sub',num2str(s),'.png'])
    close;

% turncate for quality checking
    pre_Baseline = 0.5;
    post_Baseline = 0.5;
    fs = data.Fs;
    onset_firststi = data.stimulus.values{1,1}.onset; 
    onset_laststi = data.stimulus.values{1,12}.onset+data.stimulus.values{1,12}.dur;
    
    data = raw(s).data(round((onset_firststi-pre_Baseline)*fs):round((onset_laststi+post_Baseline)*fs),:);

% quality checking using cv, cv defined as dev/m
    m = mean(data);
    dev = std(data);
    cvThresh = 0.15;  % can be changed
    cvMask = ones(size(m));
    cvOutput = ones(size(m));
    for i = 1:40
        cv = dev(i)/m(i);
        %disp(cv)
        if cv < cvThresh
            cvMask(i) = 1;
        else
            cvMask(i) = 0;
            disp('bad link:');
            disp(raw(s).probe.link(i,:));
        end
    end
    T = table(table2array(raw(s).probe.link(:,1)), table2array(raw(s).probe.link(:,2)),table2array(raw(s).probe.link(:,3)), transpose(cvMask));
    xlswrite([root_dir,'/qualitycheck/sub',num2str(s),'.csv'], table2array(T));

end

%% with or without rest
%% with rest
mkdir '/Users/yichen/Downloads/2021-05-18/' data_cleaning
for s = 1:size(raw)
    
    data = raw(s);
    fs = data.Fs;
    onset_open_agreeable = data.stimulus.values{1,2}.onset;
    stop_open_agreeable = data.stimulus.values{1,2}.onset+data.stimulus.values{1,2}.dur;
    
    onset_open_agreeable_rest = data.stimulus.values{1,3}.onset;
    stop_open_agreeable_rest = data.stimulus.values{1,3}.onset+data.stimulus.values{1,3}.dur;
    
    onset_open_neutral = data.stimulus.values{1,5}.onset;
    stop_open_neutral = data.stimulus.values{1,5}.onset+data.stimulus.values{1,5}.dur;
    
    onset_open_neutral_rest = data.stimulus.values{1,6}.onset;
    stop_open_neutral_rest = data.stimulus.values{1,6}.onset+data.stimulus.values{1,6}.dur;
    
    onset_choice_agreeable = data.stimulus.values{1,8}.onset;
    stop_choice_agreeable = data.stimulus.values{1,8}.onset+data.stimulus.values{1,8}.dur;
    
    onset_choice_agreeable_rest = data.stimulus.values{1,9}.onset;
    stop_choice_agreeable_rest = data.stimulus.values{1,9}.onset+data.stimulus.values{1,9}.dur;
    
    onset_choice_neutral = data.stimulus.values{1,11}.onset;
    stop_choice_neutral = data.stimulus.values{1,11}.onset+data.stimulus.values{1,11}.dur;
    
    onset_choice_neutral_rest = data.stimulus.values{1,12}.onset;
    stop_choice_neutral_rest = data.stimulus.values{1,12}.onset+data.stimulus.values{1,12}.dur;
    
    data = cat(1, raw(s).data(round(onset_open_agreeable*fs):round(stop_open_agreeable*fs),:), raw(s).data(round(onset_open_agreeable_rest*fs):round(stop_open_agreeable_rest*fs),:));
    data = cat(1, data, raw(s).data(round(onset_open_neutral*fs):round(stop_open_neutral*fs),:));
    data = cat(1, data, raw(s).data(round(onset_open_neutral_rest*fs):round(stop_open_neutral_rest*fs),:));
    data = cat(1, data, raw(s).data(round(onset_choice_agreeable*fs):round(stop_choice_agreeable*fs),:));
    data = cat(1, data, raw(s).data(round(onset_choice_agreeable_rest*fs):round(stop_choice_agreeable_rest*fs),:));
    data = cat(1, data, raw(s).data(round(onset_choice_neutral*fs):round(stop_choice_neutral*fs),:));
    data = cat(1, data, raw(s).data(round(onset_choice_neutral_rest*fs):round(stop_choice_neutral_rest*fs),:));
   
    
    time = raw(s).time(1:length(data),:);
    
    
    raw(s).data = data;
    raw(s).time = time;
    
end

temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [0], [114], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_agreeable_rest', [114], [17], [1])
temp('open_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'open_neutral', [131], [99], [1]);
temp('open_neutral')=temp_se;
temp_se = nirs.design.StimulusEvents( 'open_neutral_rest', [230], [17], [1]);
temp('open_neutral_rest')=temp_se;

temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [247], [84], [1]);
temp('choice_agreeable')=temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_agreeable_rest', [331], [16], [1]);
temp('choice_agreeable_rest')=temp_se;

temp_se = nirs.design.StimulusEvents( 'choice_neutral', [347], [78], [1]);
temp('choice_neutral')=temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_neutral_rest', [425], [18], [1]);
temp('choice_neutral_rest')=temp_se;

raw(1,1).stimulus = temp



temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [0], [121], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_agreeable_rest', [121], [16], [1])
temp('open_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'open_neutral', [137], [109], [1])
temp('open_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_neutral_rest', [246], [16], [1])
temp('open_neutral_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [262], [79], [1])
temp('choice_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable_rest', [341], [18], [1])
temp('choice_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'choice_neutral', [359], [72], [1])
temp('choice_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_neutral_rest', [431], [17], [1])
temp('choice_neutral_rest')=temp_se

raw(2,1).stimulus = temp



temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [0], [100], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_agreeable_rest', [100], [16], [1])
temp('open_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'open_neutral', [116], [93], [1])
temp('open_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_neutral_rest', [209], [17], [1])
temp('open_neutral_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [226], [71], [1])
temp('choice_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable_rest', [297], [20], [1])
temp('choice_agreeable_rest')=temp_se

temp_se = nirs.design.StimulusEvents( 'choice_neutral', [317], [74], [1])
temp('choice_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_neutral_rest', [391], [18], [1])
temp('choice_neutral_rest')=temp_se

raw(3,1).stimulus = temp

for s=1:size(raw)
    raw(s).draw;
    saveas(gcf,[root_dir,'/data_cleaning/sub',num2str(s),'.png'])
    close;
end

%% without rest
mkdir '/Users/yichen/Downloads/2021-05-18/' data_cleaning
for s = 1:size(raw)
    
    data = raw(s);
    fs = data.Fs;
    onset_open_agreeable = data.stimulus.values{1,2}.onset;
    stop_open_agreeable = data.stimulus.values{1,2}.onset+data.stimulus.values{1,2}.dur;
    
    onset_open_neutral = data.stimulus.values{1,5}.onset;
    stop_open_neutral = data.stimulus.values{1,5}.onset+data.stimulus.values{1,5}.dur;
    
    onset_choice_agreeable = data.stimulus.values{1,8}.onset;
    stop_choice_agreeable = data.stimulus.values{1,8}.onset+data.stimulus.values{1,8}.dur;
    
    onset_choice_neutral = data.stimulus.values{1,11}.onset;
    stop_choice_neutral = data.stimulus.values{1,11}.onset+data.stimulus.values{1,11}.dur;
    
    data = cat(1, raw(s).data(round(onset_open_agreeable*fs):round(stop_open_agreeable*fs),:), raw(s).data(round(onset_open_neutral*fs):round(stop_open_neutral*fs),:));
    data = cat(1, data, raw(s).data(round(onset_choice_agreeable*fs):round(stop_choice_agreeable*fs),:));
    data = cat(1, data, raw(s).data(round(onset_choice_neutral*fs):round(stop_choice_neutral*fs),:));
    
    time = raw(s).time(1:length(data),:);
    
    
    raw(s).data = data;
    raw(s).time = time;
    
end

temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [0], [114], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_neutral', [114], [99], [1]);
temp('open_neutral')=temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [213], [84], [1]);
temp('choice_agreeable')=temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_neutral', [297], [78], [1]);
temp('choice_neutral')=temp_se;

raw(1,1).stimulus = temp



temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [0], [121], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_neutral', [121], [109], [1])
temp('open_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [230], [79], [1])
temp('choice_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_neutral', [309], [72], [1])
temp('choice_neutral')=temp_se

raw(2,1).stimulus = temp



temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [0], [100], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'open_neutral', [100], [93], [1])
temp('open_neutral')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [193], [71], [1])
temp('choice_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_neutral', [264], [74], [1])
temp('choice_neutral')=temp_se
raw(3,1).stimulus = temp

for s=1:size(raw)
    raw(s).draw;
    saveas(gcf,[root_dir,'/data_cleaning/sub',num2str(s),'.png'])
    close;
end

%% free to pick
mkdir '/Users/yichen/Downloads/2021-05-18/' data_cleaning
%% with rest
for s = 1:size(raw)
    
    data = raw(s);
    fs = data.Fs;
    onset_open_agreeable = data.stimulus.values{1,2}.onset;
    stop_open_agreeable = data.stimulus.values{1,2}.onset+data.stimulus.values{1,2}.dur;
    
    onset_choice_agreeable = data.stimulus.values{1,8}.onset;
    stop_choice_agreeable = data.stimulus.values{1,8}.onset+data.stimulus.values{1,8}.dur;
    
    data = cat(1, raw(s).data(round(onset_open_agreeable*fs):round(stop_open_agreeable*fs),:), raw(s).data(round(onset_choice_agreeable*fs):round(stop_choice_agreeable*fs),:));
    
    time = raw(s).time(1:length(data),:);
    
    
    raw(s).data = data;
    raw(s).time = time;
    
end

temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [0], [114], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [114], [84], [1]);
temp('choice_agreeable')=temp_se;

raw(1,1).stimulus = temp



temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [0], [121], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [121], [79], [1])
temp('choice_agreeable')=temp_se

raw(2,1).stimulus = temp



temp=Dictionary()
temp_se = nirs.design.StimulusEvents( 'open_agreeable', [0], [100], [1])
temp('open_agreeable')=temp_se
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', [100], [71], [1])
temp('choice_agreeable')=temp_se

raw(3,1).stimulus = temp

for s=1:size(raw)
    raw(s).draw;
    saveas(gcf,[root_dir,'/data_cleaning/sub',num2str(s),'.png'])
    close;
end

%% Prepare the data for preprocessing
%% Resample
mkdir '/Users/yichen/Downloads/2021-05-18/' preparation
j = nirs.modules.RemoveStimless();
% Just keep stim events for infant cry (these will change depending on the
% task being analyzed)

j = nirs.modules.KeepStims( j);
%j.listOfStims = {'open_agreeable', 'choice_agreeable', 'open_neutral', 'choice_neutral','open_agreeable_rest', 'choice_agreeable_rest', 'open_neutral_rest', 'choice_neutral_rest'}; 
j.listOfStims = {'open_agreeable', 'choice_agreeable', 'open_neutral', 'choice_neutral'};
%j.listOfStims = {'open_agreeable', 'choice_agreeable'}
% resample to 1Hz 
j = nirs.modules.Resample( j);
j.Fs = 1; 
j = nirs.modules.TrimBaseline( j);
j.preBaseline = 0.5; %can change these values
j.postBaseline = 0.5;
raw = j.run(raw);

for s = 1:size(raw)
%for s = 1
% view each subject's graph
    raw(s).draw;
    saveas(gcf,[root_dir,'/preparation/sub',num2str(s),'.png'])
    close;
end



%% Run the conversions 
mkdir '/Users/yichen/Downloads/2021-05-18/' conversion
j = nirs.modules.OpticalDensity();
%% the DFP parameters can be modified in nirs.modules.BeerLamberLaw.m
%% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);    
hb = j.run(raw);
% Results visualization
for s = 1:size(hb)
%for s = 1:1
% view each subject's graph
    hb(s).draw;
    saveas(gcf,[root_dir,'/conversion/sub',num2str(s),'.png'])
    close;
end

%% to decide on best parameters for the pipeline, currently we have following parameters to decide: 
% 1. with signal preprocessing or not, 2. PCA or not, 3. if PCA, what is best ncomp
% This can take around several minutes for each subject

% job1: GLM only
job1 = nirs.modules.GLM();
job1.type = 'AR-IRLS';
job1.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);
job1.basis = Dictionary();
job1.basis('default') = nirs.design.basis.Canonical(); 
job1.basis('open_agreeable') = nirs.design.basis.Canonical();
job1.basis('choice_agreeable') = nirs.design.basis.Canonical();
job1.AddShortSepRegressors = false;

% job2: GLM without PCA
% motion correction to remove DC-shifts
job2 = nirs.modules.BaselineCorrection(); 
job2 = nirs.modules.GLM(job2);
job2.type = 'AR-IRLS';
job2.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);
job2.basis = Dictionary();
job2.basis('default') = nirs.design.basis.Canonical(); 
job2.basis('open_agreeable') = nirs.design.basis.Canonical();
job2.basis('choice_agreeable') = nirs.design.basis.Canonical();
job2.AddShortSepRegressors = true;

% job3: GLM with PCA ncomp=1
job3 = nirs.modules.PCAFilter(); 
job3.ncomp = 1;
job3 = nirs.modules.GLM(job3);
job3.type = 'AR-IRLS';
job3.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);
job3.basis = Dictionary();
job3.basis('default') = nirs.design.basis.Canonical(); 
job3.basis('open_agreeable') = nirs.design.basis.Canonical();
job3.basis('choice_agreeable') = nirs.design.basis.Canonical();
job3.AddShortSepRegressors = true;

% job4: GLM with PCA ncomp=2
job4 = nirs.modules.PCAFilter(); 
job4.ncomp = 2;
job4 = nirs.modules.GLM(job4);
job4.type = 'AR-IRLS';
job4.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);
job4.basis = Dictionary();
job4.basis('default') = nirs.design.basis.Canonical(); 
job4.basis('open_agreeable') = nirs.design.basis.Canonical();
job4.basis('choice_agreeable') = nirs.design.basis.Canonical();
job4.AddShortSepRegressors = true;

% job5: GLM with PCA ncomp=3
job5 = nirs.modules.PCAFilter(); 
job5.ncomp = 3;
job5 = nirs.modules.GLM(job5);
job5.type = 'AR-IRLS';
job5.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);
job5.basis = Dictionary();
job5.basis('default') = nirs.design.basis.Canonical(); 
job5.basis('open_agreeable') = nirs.design.basis.Canonical();
job5.basis('infantNoise') = nirs.design.basis.Canonical();
job5.AddShortSepRegressors = true;

% job6: GLM with PCA ncomp=4
job6 = nirs.modules.PCAFilter(); 
job6.ncomp = 4;
job6 = nirs.modules.GLM(job6);
job6.type = 'AR-IRLS';
job6.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);
job6.basis = Dictionary();
job6.basis('default') = nirs.design.basis.Canonical(); 
job6.basis('open_agreeable') = nirs.design.basis.Canonical();
job6.basis('choice_agreeable') = nirs.design.basis.Canonical();
job6.AddShortSepRegressors = true;

%% draw ROC to select
mkdir '/Users/yichen/Downloads/2021-05-18/' evalwithROC
%for s=1:size(hb)
for s=1:1
    ROC=nirs.testing.ChannelStatsROC;
    ROC.simfunc=@()nirs.testing.simData(hb(s));
    ROC.pipeline={job1,job2,job3,job4,job5,job6};
    ROC = ROC.run(10);
    ROC.draw('hbo');
        for i = 1:2
            saveas(gcf,[root_dir,'/evalwithROC/sub',num2str(s),'_', num2str(i),'.png'])
            close;
        end
end

%% Generate 'pipeline.xlsx'
T_pipeline = table(demographics.subject, 'VariableNames', {'sID'});
writetable(T_pipeline, [root_dir, '/pipeline.xlsx'])

% after we have the blank 'pipeline.xlsx', we need to manually start a new
% column called 'pipeline', then fill in which job you hope to apply on
% that specific subject, don't need to worry about the order, just follow
% the order in evalwithROC folder.
%% Signal Preprocessing + GLM
mkdir '/Users/yichen/Downloads/2021-05-18/' preprocessing
%read fine-tune information from xlsx file.
table_pipeline = readtable('/Users/yichen/Downloads/2021-05-18/pipeline.xlsx');
for s = 1:size(hb)
%for s=1:1
    if strcmp(table_pipeline.sID(s),hb(s).demographics.values(9))
        disp(table_pipeline.sID(s));
        disp(table_pipeline.pipeline(s));
    else
        msg = 'Error occured.';
        error(msg);
    end
    
    if table_pipeline.pipeline(s) == 1
        SubjStats(s) = job1.run(hb(s));
    elseif table_pipeline.pipeline(s) == 2
        SubjStats(s) = job2.run(hb(s));
    elseif table_pipeline.pipeline(s) == 3
        SubjStats(s) = job3.run(hb(s));
    elseif table_pipeline.pipeline(s) == 4
        SubjStats(s) = job4.run(hb(s));
    elseif table_pipeline.pipeline(s) == 5
        SubjStats(s) = job5.run(hb(s)); 
    elseif table_pipeline.pipeline(s) == 6
        SubjStats(s) = job6.run(hb(s)); 
    end
end
       
%% Block Average for visualizations only
mkdir '/Users/yichen/Downloads/2021-05-18/' blockAverage
for s = 1:size(hb)
%for s = 1:1
    HRF=BlockAverage(-1, 19, hb(s),s);  %the parameters of start and end here is based on experimental protocols:Number of Trials: 20 trials
                              %Trials by Stimuli: 10 trials control cry, 10 trials control noise ? 10s block
                              %Pre/post durations: 10s pre-stim, 10s
                              %post-stim--import block average function and
                              %homer as well
               
    for i = 1:2
        saveas(gcf,[root_dir,'/blockAverage/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
end


%% save GLM model
mkdir '/Users/yichen/Downloads/2021-05-18/' GLM
for s = 1:size(transpose(SubjStats))
%for s = 1:1
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.txt'], 'Delimiter', ' '); 
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.xls']);
end


%% export for 2-nd level analysis
mkdir 'rootdir/' secondlevelAnalysis_exported
T_temp = table(repmat(SubjStats(1).demographics.values(9),80,1), repmat(SubjStats(1).demographics.values(10),80,1), repmat(SubjStats(1).demographics.values(8),80,1), SubjStats(1).variables.source, SubjStats(1).variables.detector, SubjStats(1).variables.type, SubjStats(1).variables.cond, SubjStats(1).beta, 'VariableNames', {'subjectID','session', 'group', 'source','detector', 'type', 'cond', 'beta'});
for sID = 2:size(transpose(SubjStats)) 
    T_temp_loop = table(repmat(SubjStats(sID).demographics.values(9),80,1), repmat(SubjStats(sID).demographics.values(10),80,1), repmat(SubjStats(sID).demographics.values(8),80,1), SubjStats(sID).variables.source, SubjStats(sID).variables.detector, SubjStats(sID).variables.type, SubjStats(sID).variables.cond, SubjStats(sID).beta, 'VariableNames', {'subjectID','session', 'group', 'source','detector', 'type', 'cond', 'beta'}); 
    T_temp = [T_temp ; T_temp_loop]; 
end
writetable(T_temp, [root_dir, '/secondlevelAnalysis_exported/exportedforSPSS.xls'])

%% data view on individual level
mkdir '/Users/yichen/Downloads/2021-05-18/' data_individual
for s = 1:size(transpose(SubjStats))
%for s = 1:1
    %SubjStats(s).probe.defaultdrawfcn='3D mesh'/'10-20';  % cannot work, error message 'No public field defaultdrawfcn exists for class nirs.core.Probe.'
    SubjStats(s).probe.defaultdrawfcn='3D-mesh';
    SubjStats(s).draw('tstat', [], 'q < 0.05');
    for i = 1:4
        saveas(gcf,[root_dir,'/data_individual/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
    writetable(SubjStats(s).table,[root_dir,'/data_individual/SubjStatssub',num2str(s),'.xls']);
end

%% contrast on individual level and ROI this is only IC and Infant Noise
%% display conditions
disp(SubjStats(1).conditions)
%% Define some contrasts
c = [-1 1]

mkdir '/Users/yichen/Downloads/2021-05-18/' contrast_individual
mkdir '/Users/yichen/Downloads/2021-05-18/' data_individual_ROI
for s = 1:size(transpose(SubjStats))
    ContrastStats = SubjStats(s).ttest(c);
    % Display the contrasts
    ContrastStats.draw('tstat', [], 'q < 0.05');
    
    for i = 1:2
        saveas(gcf,[root_dir,'/contrast_individual/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
    %writetable(ContrastStats.table,[root_dir,'/contrast_individual/ContrastStatssub',num2str(s),'.xls'], 'Delimiter', ' ');

%% ROI Analysis Individual
%Src - Det to include in the ROI, need to be manually defined by ourself,
%below is just sample Src-Det set

    MeasList=[4 3;...
            4 2;...
            2 2;...
            2 3];
      
    Region{1} = table(MeasList(:,1),MeasList(:,2),'VariableNames',{'source','detector'});
    ROItable=nirs.util.roiAverage(ContrastStats,Region,{'region1'});
    disp(ROItable);
%%
    job_ROI = nirs.modules.ApplyROI();
    job_ROI.listOfROIs = table(MeasList(1,1),MeasList(1,2),{'Edge1'},'VariableNames',{'source','detector','name'});
    job_ROI.listOfROIs(1,:) = table(MeasList(1,1),MeasList(1,2),{'Edge1'},'VariableNames',{'source','detector','name'});
    job_ROI.listOfROIs(2,:) = table(MeasList(2,1),MeasList(2,2),{'Edge2'},'VariableNames',{'source','detector','name'});
    job_ROI.listOfROIs(3,:) = table(MeasList(3,1),MeasList(3,2),{'Edge3'},'VariableNames',{'source','detector','name'});
    job_ROI.listOfROIs(4,:) = table(MeasList(4,1),MeasList(4,2),{'Edge4'},'VariableNames',{'source','detector','name'});
    job.weighted = false;
    dataROI = job_ROI.run( ContrastStats );
%% visualize ROI analysis stats

    
    dataROI.draw('tstat', [], 'p < 0.05');
    
    for i = 1:2
        saveas(gcf,[root_dir,'/data_individual_ROI/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
    %writetable(dataROI.table,[root_dir,'/data_individual_ROI/dataROIsub',num2str(s),'.xls'], 'Delimiter', ' ');
end

%% Group Analysis
%% remove outlier
job=nirs.modules.RemoveOutlierSubjects;
SubjStats = job.run(SubjStats);
%tbl=nirs.util.grouplevelleveragestats(SubjStats); % This reports the leverage per subject, channel, condition

%% group models - mixedEffects model with ANOVA contrast
% ANOVA first and then mixed effects to check detailed model information 
job = nirs.modules.Anova();
%job.formula = 'beta ~ -1 + group:cond:session + (1|subject)'; 
% if to go with both main effects and interactions, do the following
% formula 
job.formula = 'beta ~ cond + (1|subject)'
job.dummyCoding = 'effects';
% if to go with both main effects and interactions, do the following
% dummyCoding
% "effects"
GroupStats_ANOVA = job.run(SubjStats);

mkdir '/Users/yichen/Downloads/2021-05-18/' group_stats_ANOVA
save([root_dir,'/group_stats_ANOVA/group_stats_ANOVA.mat'],'GroupStats_ANOVA');
writetable(GroupStats_ANOVA.table,'/Users/yichen/Downloads/2021-05-18/ANOVATable.xls');
%%
disp(GroupStats_ANOVA.conditions);
%% Visualize the group level stats - ANOVA
mkdir '/Users/yichen/Downloads/2021-05-18/' group_ANOVA
GroupStats_ANOVA.draw([], 'p < 0.05')

for i = 1:4
    saveas(gcf,[root_dir,'/group_ANOVA/',num2str(i),'.png'])
    close;
end
writetable(GroupStats_ANOVA.table,'/Users/yichen/Downloads/2021-05-18/ANOVATable2.xls');

%% FDR corrected (q < 0.05) - ANOVA
mkdir '/Users/yichen/Downloads/2021-05-18/' FDRcorrected_group_ANOVA
GroupStats_ANOVA.draw([], 'q < 0.05')

for i = 1:8
    saveas(gcf,[root_dir,'/FDRcorrected_group_ANOVA/',num2str(i),'.png'])
    close;
 writetable(GroupStats_ANOVA.table,'/Users/yichen/Downloads/2021-05-18/fdrANOVATable.xls');
end

%% Use MixedEffects to check details 
job = nirs.modules.MixedEffects();
job.formula = 'beta ~ cond + (1|subject)';
% if to go with both main effects and interactions, do the following
% formula 
% job.formula = 'beta ~ cond + group + session + cond:group + cond:session +
% group:session + cond:group:session + (1|subject)'
job.dummyCoding = 'effects';
% if to go with both main effects and interactions, do the following
% dummyCoding
% "effects"
job.include_diagnostics=true;
GroupStats = job.run(SubjStats);

mkdir '/Users/yichen/Downloads/2021-05-18/' group_stats
save([root_dir,'/group_stats/group_stats.mat'],'GroupStats');
%%
disp(GroupStats.conditions);
%% Visualize the group level stats
%mkdir 'Y:\2_Individual_Projects\RISE_PV1_CannabisCry\06_22_2020_Prenatal_Processing\PV1\02.16.2021_finalprocessing\nocovariates' data_group
mkdir '/Users/yichen/Downloads/2021-05-18/' data_group
GroupStats.draw('tstat', [], 'p < 0.05')
for i = 1:4
    saveas(gcf,[root_dir,'/data_group/',num2str(i),'.png'])
    close;
end
StatsTable=GroupStats.table;
writetable(StatsTable,'GroupLevelResults.xls');
%% FDR corrected (q < 0.05)
%mkdir 'Y:\2_Individual_Projects\RISE_PV1_CannabisCry\06_22_2020_Prenatal_Processing\PV1\02.16.2021_finalprocessing\nocovariates' data_FDRcorrected_group
mkdir '/Users/yichen/Downloads/2021-05-18/' data_FDRcorrected_group
GroupStats.draw('tstat', [], 'q < 0.05')
for i = 1:4
    saveas(gcf,[root_dir,'/data_FDRcorrected_group/',num2str(i),'.png'])
    close;
end

%% examine detailed models
s1d1_effects = GroupStats.variables(GroupStats.variables.source == 1 & GroupStats.variables.detector == 1, :);
s1d1_effects.model{1}
%directly go to vairables to look for evaluation criteria, such as BIC,
%etc.

%% ROI Analysis Group
%Src - Det to include in the ROI MeasuresList, first column is source
%second is detector
MeasList=[4 4;...
          5 4;...
          5 6;...
          6 6;...
          6 5;...
          4 5];
      
Region{1} = table(MeasList(:,1),MeasList(:,2),'VariableNames',{'source','detector'});
ROItable=nirs.util.roiAverage(GroupStats_ANOVA,Region,{'region1'});
disp(ROItable);

job_ROI = nirs.modules.ApplyROI();
job_ROI.listOfROIs = table(MeasList(1,1),MeasList(1,2),{'Edge1'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(1,:) = table(MeasList(1,1),MeasList(1,2),{'Edge1'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(2,:) = table(MeasList(2,1),MeasList(2,2),{'Edge2'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(3,:) = table(MeasList(3,1),MeasList(3,2),{'Edge3'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(4,:) = table(MeasList(4,1),MeasList(4,2),{'Edge4'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(3,:) = table(MeasList(3,1),MeasList(5,2),{'Edge5'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(4,:) = table(MeasList(4,1),MeasList(6,2),{'Edge6'},'VariableNames',{'source','detector','name'});
job.weighted = false;
dataROI = job_ROI.run( GroupStats_ANOVA );

% visualize ROI analysis stats 
mkdir 'rootdir/' data_group_ROI
for s = 1:size(dataROI)
    dataROI(s).draw([], 'q < 0.05');
     for i = 1:7*group_number
        saveas(gcf,[root_dir,'/data_group_ROI/',num2str(s),'_', num2str(i),'.png'])
        close;
    end
    writetable(dataROI(s).table,'rootdir/ROIFPTable.xls');
end