% 1. change the pipeline order
% 2. use pre-whitening/DCT trending
% 3. use FIR as basis function
% 4. use mixed-effects group level to time-average subject-level 
% 5. ROI analysis
% 6. Haven't done ANCOVA yet, need to ask about ANCOVA first
% 7. raw.probe.defaultdrawfcn='10-20';
% 8. find(conStats.ttest.p<0.05)
%figure
%raw.probe.draw
%raw.probe.defaultdrawfcn='3D Mesh (top)'; 
% nirs.io.loadNIRx
%figure
%raw.probe.draw
clear;
clc;

%% group number
group_number = 2;

%% load data
%root_dir = '/Users/hanxu/Downloads/fnirs_test'
%raw = nirs.io.loadDirectory(root_dir, {'group', 'subject'});
%raw = nirs.io.loadNIRx(root_dir)

%for visualization
%root_dir = '/Users/hanxu/Downloads/prelim_data/5084/2019-05-10_002';  % The path here needs to be the path where it stores the
%source files of fnirs data, rather than the coverted .nirs
%raw = nirs.io.loadNIRx(root_dir);
%raw.probe.defaultdrawfcn='3D Mesh (top)';
%raw.probe.draw

root_dir = '/Users/hanxu/Downloads/fnirs_test'
raw = nirs.io.loadDirectory(root_dir, {'group', 'subject'});

%% Create demographics table

demographics = nirs.createDemographicsTable(raw);
disp(demographics)
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
mkdir /Users/hanxu/Downloads/fnirs_test rawdata
mkdir /Users/hanxu/Downloads/fnirs_test qualitycheck
for s = 1:size(raw)
% view each subject's graph
    data = raw(s);
    data.draw;
    saveas(gcf,[root_dir,'/rawdata/sub',num2str(s),'.png'])
    close;
    
% Change stimuli duration
    %raw(s) = nirs.viz.StimUtil(raw(s));
    raw(s)=nirs.design.change_stimulus_duration(raw(s),'channel_6',10);
    raw(s)=nirs.design.change_stimulus_duration(raw(s),'channel_7',10);
    %raw(s).gui

% turncate for quality checking
    pre_Baseline = 10;
    post_Baseline = 10;
    fs = data.Fs;
    [m,n] = size(data.stimulus.keys);
    val6 = 6;
    val7 = 7;
    for a = 1:n
        if strcmp(data.stimulus.keys{1,a}, 'channel_6')
            val6 = a;
        end
        %if data.stimulus.keys{1,a} == 'stim_channel11'
        if strcmp(data.stimulus.keys{1,a}, 'channel_7')
            val7 = a;
            break;
        end
    end
  
   
    onset_firststi_val6 = data.stimulus.values{1,val6}.onset(1); 
    onset_firststi_val7 = data.stimulus.values{1,val7}.onset(1); 
    
    if onset_firststi_val6 < onset_firststi_val7
        onset_firststi = onset_firststi_val6
    else
        onset_firststi = onset_firststi_val7
    end
    disp(onset_firststi);
    %onset_laststi = data.stimulus.values{1,val7}.onset(10);
    [z1_val6,z2_val6] = size(data.stimulus.values{1,val6}.onset);
    [z1_val7,z2_val7] = size(data.stimulus.values{1,val7}.onset);

    
    onset_laststi_val6 = data.stimulus.values{1,val6}.onset(z1_val6);
    onset_laststi_val7 = data.stimulus.values{1,val7}.onset(z1_val7);
   
    if onset_laststi_val6 > onset_laststi_val7
        display(onset_laststi_val6);
        data = raw(s).data(round((onset_firststi-pre_Baseline)*fs):round((onset_laststi_val6+post_Baseline)*fs),:);
    else
        display(onset_laststi_val7);
        data = raw(s).data(round((onset_firststi-pre_Baseline)*fs):round((onset_laststi_val7+post_Baseline)*fs),:);
    end  
    
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
j.listOfStims = {'channel_6', 'channel_7'};  % take a look at the raw data containing condition information, to see whether we can fix the condition information

% rename stims
j = nirs.modules.RenameStims( j );
j.listOfChanges = {
    'channel_6', 'infantCry'; 
    'channel_7', 'infantNoise'};

% resample to 4Hz
j = nirs.modules.Resample( j );
j.Fs = 4;  % Sets the new sample rate to 2Hz (was 10Hz).

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

%% Run the conversions 
mkdir /Users/hanxu/Downloads/fnirs_test conversion
j = nirs.modules.OpticalDensity();
%% the DFP parameters can be modified in nirs.modules.BeerLamberLaw.m
%% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);    
hb = j.run(raw);
% Results visualization
for s = 1:size(hb)
% view each subject's graph
    hb(s).draw;
    saveas(gcf,[root_dir,'/conversion/sub',num2str(s),'.png'])
    close;
end

%% Signal Preprocessing -> PCA
mkdir /Users/hanxu/Downloads/fnirs_test preprocessing
% motion correction to remove DC-shifts
jobs = nirs.modules.PCAFilter(); 
%jobs.PCA=true;
hb = jobs.run( hb );
% visualize results
for s = 1:size(hb)
% view each subject's graph
    hb(s).draw;
    saveas(gcf,[root_dir,'/preprocessing/sub',num2str(s),'.png'])
    close;
end


%% Block Average for visualizations only
mkdir /Users/hanxu/Downloads/fnirs_test blockAverage
for s = 1:size(hb)
    HRF=BlockAverage(-1, 19, hb(s),s);  %the parameters of start and end here is based on experimental protocols:Number of Trials: 20 trials
                              %Trials by Stimuli: 10 trials control cry, 10 trials control noise ? 10s block
                              %Pre/post durations: 10s pre-stim, 10s post-stim
               
    for i = 1:2
        saveas(gcf,[root_dir,'/blockAverage/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
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
j.basis('infantCry') = nirs.design.basis.Canonical();
j.basis('infantNoise') = nirs.design.basis.Canonical();
SubjStats = j.run( hb );
%% the parameters can be modified in nirs.design.basis.Canonical
%% save GLM model
mkdir /Users/hanxu/Downloads/fnirs_test GLM
for s = 1:size(transpose(SubjStats))
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.txt'], 'Delimiter', ' '); 
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.xls'], 'Delimiter', ' ');
end

%% Add demographic information
j_demo = nirs.modules.AddDemographics();
j_demo.demoTable = readtable( [root_dir '/demographics.csv'] );

% We are going to match the subject column in the above table
j_demo.varToMatch = 'subject';  
% Note "subject" (or whatever you are matching based on) needs to be an
% entry in BOTH the CSV file and the nirs.core.<> class that you are
% loading the demographics to.  

% This is a job, so we need to run it like the others.  In this case, let's
% run the job on the SubjStats variable we just created.  We could, of
% course, have ran this job earlier on the raw data
SubjStats = j_demo.run(SubjStats);

% Just like before, we can check the demographics by using
disp(nirs.createDemographicsTable( SubjStats ));
%% data view on individual level
mkdir /Users/hanxu/Downloads/fnirs_test data_individual
for s = 1:size(transpose(SubjStats))
    %SubjStats(s).probe.defaultdrawfcn='3D mesh'/'10-20';  % cannot work, error message 'No public field defaultdrawfcn exists for class nirs.core.Probe.'
    SubjStats(s).probe.defaultdrawfcn='3D mesh';
    SubjStats(s).draw('tstat', [], 'p < 0.05');
    for i = 1:4
        saveas(gcf,[root_dir,'/data_individual/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
    %writetable(SubjStats(s).table,[root_dir,'/data_individual/SubjStatssub',num2str(s),'.xls'], 'Delimiter', ' ');
end

%% contrast on individual level and ROI
% Define some contrasts
c = [-1 1]

mkdir /Users/hanxu/Downloads/fnirs_test contrast_individual
mkdir /Users/hanxu/Downloads/fnirs_test data_individual_ROI
for s = 1:size(transpose(SubjStats))
    ContrastStats = SubjStats(s).ttest(c);
    % Display the contrasts
    ContrastStats.draw('tstat', [], 'p < 0.05');
    
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

%% group models
% Uses Wilkinson notation for specifying models
job = nirs.modules.MixedEffects();
job.formula = 'beta ~ -1 + group:cond + demo1 + (1|subject)';
% this would model the two conditions (A and B) per group treating subject as a random variable.
%job.dummyCoding = 'full';
job.dummyCoding = 'full';
job.include_diagnostics=true;
GroupStats_ME = job.run(SubjStats);
%%%job=nirs.modules.AnovaN;
%%%job.variables={'cond','group'};
%%%GroupStats = job.run(SubjStats);

mkdir /Users/hanxu/Downloads/fnirs_test group_stats
save([root_dir,'/group_stats/group_stats.mat'],'GroupStats_ME');

%% Visualize the group level stats
mkdir /Users/hanxu/Downloads/fnirs_test data_group
GroupStats_ME.draw('tstat', [-10 10], 'p < 0.05')

for i = 1:5*group_number
    saveas(gcf,[root_dir,'/data_group/',num2str(i),'.png'])
    close;
end
%writetable(GroupStats.table,[root_dir,'/data_group/GroupStats','.xls'], 'Delimiter', ' ');

%% FDR corrected (q < 0.05)
mkdir /Users/hanxu/Downloads/fnirs_test data_FDRcorrected_group
GroupStats_ME.draw('tstat', [-5 5], 'q < 0.05')

for i = 1:5*group_number
    saveas(gcf,[root_dir,'/data_FDRcorrected_group/',num2str(i),'.png'])
    close;
end

%% group contrast -- compare between two conditions or groups
%% display condition
disp(GroupStats_ME.conditions);
%% run contrast
%c = [eye(4);  % all 5 of the original variables
     %1 0 -1 0; % X - Y for group 1
     %0 1 0 -1; % X - Y for group 2
     %1 -1 0 0; % G1 - G2 for X
     %0 0 1 -1]; % G1 - G2 for Y
c = [0 1 -1 0 0] % need to adjust 'c' based on how many covariates are we going to include
size_c = size(c);
mkdir /Users/hanxu/Downloads/fnirs_test contrast_group
ContrastStats = GroupStats_ME.ttest(c);
ContrastStats.draw('tstat', [-5 5], 'p < 0.05');



for i = 1:2*size_c(1)
    saveas(gcf,[root_dir,'/contrast_group/',num2str(i),'.png'])
    close;
end
%writetable(ContrastStats.table,[root_dir,'/contrast_group/ContrastStats','.xls'], 'Delimiter', ' ');
%% FDR corrected (q < 0.05)
ContrastStats.draw('tstat', [-5 5], 'q < 0.05');
mkdir /Users/hanxu/Downloads/fnirs_test contrast_FDRcorrected_group
for i = 1:2*size_c(1)
    saveas(gcf,[root_dir,'/contrast_FDRcorrected_group/',num2str(i),'.png'])
    close;
end

%% ROI Analysis Group
%Src - Det to include in the ROI
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
mkdir /Users/hanxu/Downloads/fnirs_test data_group_ROI
for s = 1:size(dataROI)
    dataROI(s).draw('tstat', [], 'p < 0.05');
    for i = 1:2*size_c(1)
        saveas(gcf,[root_dir,'/data_group_ROI/',num2str(s),'_', num2str(i),'.png'])
        close;
    end
    %writetable(dataROI(s).table,[root_dir,'/data_group_ROI/dataROI',num2str(s), '.xls'], 'Delimiter', ' ');
end


%% group contrast -- compare among multiple variables - using ANCOVA
job = nirs.modules.Anova();
job.formula = 'beta ~ group*cond + demo1 + (1|subject)';
job.dummyCoding = 'effects';
GroupStats_ANOVA = job.run(SubjStats);

mkdir /Users/hanxu/Downloads/fnirs_test group_stats_ANOVA
save([root_dir,'/group_stats_ANOVA/group_stats_ANOVA.mat'],'GroupStats_ANOVA');
%% Visualize the group level stats - ANOVA
mkdir /Users/hanxu/Downloads/fnirs_test group_ANOVA
GroupStats_ANOVA.draw(100, 'p < 0.05')

for i = 1:5*group_number
    saveas(gcf,[root_dir,'/group_ANOVA/',num2str(i),'.png'])
    close;
end
%writetable(GroupStats.table,[root_dir,'/data_group/GroupStats','.xls'], 'Delimiter', ' ');

%% FDR corrected (q < 0.05) - ANOVA
mkdir /Users/hanxu/Downloads/fnirs_test FDRcorrected_group_ANOVA
GroupStats_ANOVA.draw(100, 'q < 0.05')

for i = 1:5*group_number
    saveas(gcf,[root_dir,'/FDRcorrected_group_ANOVA/',num2str(i),'.png'])
    close;
end