% 1. change the pipeline order
% 2. use pre-whitening/DCT trending
% 3. use FIR as basis function
% 4. use mixed-effects group level to time-average subject-level 
% 5. ROI analysis
% 6. Haven't done ANCOVA yet, need to ask about ANCOVA first
clear;
clc;

%% group number
group_number = 1;

%% load data
root_dir = '/Users/hanxu/Downloads/fnirs_test'
raw = nirs.io.loadDirectory(root_dir, {'group', 'subject'});

%% Space Registration
% add an anchor to attach the probe to the center of the forehead (FpZ)

Name{1}='FpZ';
xyz(1,:)=[0 0 0];
Type{1}='FID-anchor';  % This is an anchor point
Units{1}='mm';

%Now add a few more
Name{2}='Cz';
xyz(2,:)=[0 100 0];
Type{2}='FID-attractor';  % This is an attractor
Units{2}='mm';

Name{3}='T7';
xyz(3,:)=[-200 0 0];
Type{3}='FID-attractor';  % This is an attractor
Units{3}='mm';

Name{4}='T8';
xyz(4,:)=[200 0 0];
Type{4}='FID-attractor';  % This is an attractor
Units{4}='mm';

% now add these points to the optodes field in the probe. 
fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});
raw.probe.optodes=[raw.probe.optodes; fid];

% Use the default head size
probe1020=nirs.util.registerprobe1020(raw.probe);
probe1020.defaultdrawfcn='10-20';
probe1020.draw;

%% Quality checking -> bad channel rejecting
mkdir /Users/hanxu/Downloads/fnirs_test rawdata
mkdir /Users/hanxu/Downloads/fnirs_test qualitycheck
for s = 1:size(raw)
% view each subject's graph
    data = raw(s);
    data.draw;
    saveas(gcf,[root_dir,'/rawdata/sub',num2str(s),'.png'])
    close;
    
% Change stimuli duration
    raw(s) = nirs.viz.StimUtil(raw(s));

% turncate for quality checking
    pre_Baseline = 10;
    post_Baseline = 10;
    fs = data.Fs;
    [m,n] = size(data.stimulus.keys);
    val6 = 6;
    val7 = 7;
    for a = 1:n
        if strcmp(data.stimulus.keys{1,a}, 'stim_channel6')
            val6 = a;
        end
        %if data.stimulus.keys{1,a} == 'stim_channel11'
        if strcmp(data.stimulus.keys{1,a}, 'stim_channel7')
            val7 = a;
            break;
        end
    end
  
   
    onset_firststi = data.stimulus.values{1,val6}.onset(1); 
    disp(onset_firststi);
    onset_laststi = data.stimulus.values{1,val7}.onset(10);
    
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

%% Prepare the data for preprocessing
%% Remove stimless files
mkdir /Users/hanxu/Downloads/fnirs_test stimremove_turncate
j = nirs.modules.RemoveStimless();

% Just keep stim events for infant cry (these will change depending on the
% task being analyzed)

j = nirs.modules.KeepStims( j);
j.listOfStims = {'stim_channel6', 'stim_channel7'};  % take a look at the raw data containing condition information, to see whether we can fix the condition information

% Trim pre and post baseline
j = nirs.modules.TrimBaseline( j);
j.preBaseline = 10; %can change these values
j.postBaseline = 10;
raw = j.run(raw);

for s = 1:size(raw)
% view each subject's graph
    raw(s).draw;
    saveas(gcf,[root_dir,'/stimremove_turncate/sub',num2str(s),'.png'])
    close;
end

%% Signal Preprocessing -> PCA
mkdir /Users/hanxu/Downloads/fnirs_test preprocessing
% motion correction to remove DC-shifts
jobs = nirs.modules.BaselineCorrection(); 
jobs.PCA=true;
filtered_raw = jobs.run( raw );
% visualize results
for s = 1:size(filtered_raw)
% view each subject's graph
    filtered_raw(s).draw;
    saveas(gcf,[root_dir,'/preprocessing/sub',num2str(s),'.png'])
    close;
end
%% Run the conversions 
mkdir /Users/hanxu/Downloads/fnirs_test conversion
j = nirs.modules.OpticalDensity();
%% the DFP parameters can be modified in nirs.modules.BeerLamberLaw.m
%% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);    
hb = j.run(filtered_raw);
% Results visualization
for s = 1:size(hb)
% view each subject's graph
    hb(s).draw;
    saveas(gcf,[root_dir,'/conversion/sub',num2str(s),'.png'])
    close;
end

%% Block Average
mkdir /Users/hanxu/Downloads/fnirs_test blockAverage
HRF=BlockAverage(-5, 15, hb)
for i = 1:2
    saveas(gcf,[root_dir,'/blockAverage/sub',num2str(s),'_', num2str(i),'.png'])
    close;
end
%% Run the GLM
j = nirs.modules.GLM();
j.verbose = true;

%choose the AR-IRLS model
j.type = 'AR-IRLS';

% DCT terms with a frequency cutoff of 0.08.
j.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);

%Canonical basis function
j.basis = Dictionary();
j.basis('default') = nirs.design.basis.Canonical(); 
j.basis('stim_channel6') = nirs.design.basis.Canonical();
j.basis('stim_channel7') = nirs.design.basis.Canonical();
SubjStats = j.run( hb );

% The FIR basis function
%j.basis = Dictionary();
%j.basis('default') = nirs.design.basis.FIR(); 
%j.basis('stim_channel6') = nirs.design.basis.FIR();
%j.basis('stim_channel7') = nirs.design.basis.FIR();
%SubjStats = j.run( hb );
%% the parameters can be modified in nirs.design.basis.Canonical
%% save GLM model
mkdir /Users/hanxu/Downloads/fnirs_test GLM
for s = 1:size(transpose(SubjStats))
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.txt'], 'Delimiter', ' '); 
end

%% data view on individual level
mkdir /Users/hanxu/Downloads/fnirs_test data_individual
for s = 1:size(transpose(SubjStats))
    %SubjStats(s).probe.defaultdrawfcn='3D mesh';
    SubjStats(s).draw('tstat', [], 'p < 0.05');
    for i = 1:4
        saveas(gcf,[root_dir,'/data_individual/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
end

%% contrast on individual level
% Define some contrasts
c = [-1 1]

mkdir /Users/hanxu/Downloads/fnirs_test contrast_individual
for s = 1:size(transpose(SubjStats))
    ContrastStats = SubjStats(s).ttest(c);
    % Display the contrasts
    ContrastStats.draw('tstat', [], 'p < 0.05');
    for i = 1:2
        saveas(gcf,[root_dir,'/contrast_individual/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
end

%% ROI Analysis Individual
%Src - Det to include in the ROI
MeasList=[1 1;...
          1 2;...
          1 6;...
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
job.weighted = true;
dataROI = job_ROI.run( ContrastStats );
%% visualize ROI analysis stats
mkdir /Users/hanxu/Downloads/fnirs_test data_individual_ROI
for s = 1:size(transpose(dataROI))
    dataROI(s).draw('tstat', [], 'p < 0.05');
    for i = 1:2
        saveas(gcf,[root_dir,'/data_individual_ROI/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
end

%% Group Analysis
%% remove outlier
job=nirs.modules.RemoveOutlierSubjects;
SubjStats = job.run(SubjStats);
%tbl=nirs.util.grouplevelleveragestats(SubjStats); % This reports the leverage per subject, channel, condition

%% group models
% Uses Wilkinson notation for specifying models
job = nirs.modules.MixedEffects();
job.formula = 'beta ~ -1 + group:cond + (1|subject)';
% this would model the two conditions (A and B) per group treating subject as a random variable.
job.dummyCoding = 'full';
job.include_diagnostics=true;
GroupStats = job.run(SubjStats);
%%%job=nirs.modules.AnovaN;
%%%job.variables={'cond','group'};
%%%GroupStats = job.run(SubjStats);

mkdir /Users/hanxu/Downloads/fnirs_test group_stats
save([root_dir,'/group_stats/group_stats.mat'],'GroupStats');

%% Visualize the group level stats
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test data_group
GroupStats.draw('tstat', [-10 10], 'p < 0.05')
for i = 1:4*group_number
    saveas(gcf,[root_dir,'/data_group/',num2str(i),'.png'])
    close;
end
%% FDR corrected (q < 0.05)
mkdir /Users/hanxu/Downloads/fnirs_test data_FDRcorrected_group
GroupStats.draw('tstat', [-5 5], 'q < 0.05')
for i = 1:4*group_number
    saveas(gcf,[root_dir,'/data_FDRcorrected_group/',num2str(i),'.png'])
    close;
end

%% group contrast
%% display condition
disp(GroupStats.conditions);
%% run contrast
c = [eye(4);  % all 5 of the original variables
     1 0 -1 0; % X - Y for group 1
     0 1 0 -1; % X - Y for group 2
     1 -1 0 0; % G1 - G2 for X
     0 0 1 -1]; % G1 - G2 for Y
%c = [1 -1 0 0]
size_c = size(c)
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat', [-5 5], 'p < 0.05');
mkdir /Users/hanxu/Downloads/fnirs_test contrast_group
for i = 1:2*size_c(1)
    saveas(gcf,[root_dir,'/contrast_group/',num2str(i),'.png'])
    close;
end
%% FDR corrected (q < 0.05)
ContrastStats.draw('tstat', [-5 5], 'q < 0.05');
mkdir /Users/hanxu/Downloads/fnirs_test contrast_FDRcorrected_group
for i = 1:2*size_c(1)
    saveas(gcf,[root_dir,'/contrast_FDRcorrected_group/',num2str(i),'.png'])
    close;
end

%% ROI Analysis Group
%Src - Det to include in the ROI
MeasList=[1 1;...
          1 2;...
          1 6;...
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
job.weighted = true;
dataROI = job_ROI.run( ContrastStats );

%% visualize ROI analysis stats
mkdir /Users/hanxu/Downloads/fnirs_test data_group_ROI
for s = 1:size(transpose(SubjStats))
    dataROI(s).draw('tstat', [], 'p < 0.05');
    for i = 1:2
        saveas(gcf,[root_dir,'/data_individual_ROI/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
end
%% Use mixed-effects group level to time-average subject-level
%HRF = GroupStats.HRF;
%nirs.viz.plot2D(HRF);

% find HRF peak via visual inspection

%SubjStats_averaged_stimchannel6 = mean(SubjStats{'stim_channel6[3:8s]'})
%SubjStats_averaged_stimchannel7 = mean(SubjStats{'stim_channel7[3:8s]'})

