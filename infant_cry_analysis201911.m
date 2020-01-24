% need explore:
% 1. ANOVA group analysis
% 2. the formula
% FDR corrected Q value


clear;
clc;

%% load data
root_dir = '/Users/hanxu/Downloads/prelim_data_nirs_test'
raw = nirs.io.loadDirectory(root_dir, {'group', 'subject'});

%% Quality checking
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test rawdata
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test qualitycheck
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
    onset_firststi = data.stimulus.values{1,6}.onset(1);
    disp(onset_firststi);
    onset_laststi = data.stimulus.values{1,7}.onset(10);
    data = raw(s).data(round((onset_firststi-pre_Baseline)*fs):round((onset_laststi+post_Baseline)*fs),:);

% quality checking using cv, cv defined as dev/m
    m = mean(data);
    dev = std(data);
    cvThresh = 0.15;  % can be changed
    cvMask = ones(size(m));
    cvOutput = ones(size(m));
    for i = 1:40
        cv = dev(i)/m(i);
        disp(cv)
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
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test stimremove_turncate
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

%% Signal Preprocessing
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test preprocessing
% motion correction to remove DC-shifts
jobs = nirs.modules.BaselineCorrection(); 
jobs.PCA=true;
BaselineCorrection_raw = jobs.run( raw ); 

% Filter to remove outliers (motion) and low freq characteristics
jobs = nirs.modules.WaveletFilter(); 
filtered_raw = jobs.run(BaselineCorrection_raw); 

% visualize results
for s = 1:size(filtered_raw)
% view each subject's graph
    filtered_raw(s).draw;
    saveas(gcf,[root_dir,'/preprocessing/sub',num2str(s),'.png'])
    close;
end

%% Run the conversions 
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test conversion
j = nirs.modules.OpticalDensity();
%% the DFP parameters can be modified in nirs.modules.BeerLamberLaw.m
%% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);    
hb = j.run( filtered_raw );

% Results visualization
for s = 1:size(hb)
% view each subject's graph
    hb(s).draw;
    saveas(gcf,[root_dir,'/conversion/sub',num2str(s),'.png'])
    close;
end

%% Run the GLM
j = nirs.modules.GLM();
j.verbose = true;
% This example specifies DCT terms with a frequency cutoff of 0.08.
j.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);
j.basis = Dictionary();
j.basis('default') = nirs.design.basis.Canonical(); 
j.basis('stim_channel6') = nirs.design.basis.Canonical();
j.basis('stim_channel7') = nirs.design.basis.Canonical();
SubjStats = j.run( hb );
%% the DFP parameters can be modified in nirs.design.basis.Canonical
%% save GLM model
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test GLM
for s = 1:size(transpose(SubjStats))
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.txt'], 'Delimiter', ' '); 
end

%% data view on individual level
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test data_individual
for s = 1:size(transpose(SubjStats))
    SubjStats(s).draw('tstat', [-10 10], 'p < 0.05');
    for i = 1:4
        saveas(gcf,[root_dir,'/data_individual/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
end

%% contrast on individual level
% Define some contrasts
c = [-1 1]

mkdir /Users/hanxu/Downloads/prelim_data_nirs_test contrast_individual
for s = 1:size(transpose(SubjStats))
    ContrastStats = SubjStats(s).ttest(c);
    % Display the contrasts
    ContrastStats(s).draw('tstat', [-10 10], 'p < 0.05');
    for i = 1:2
        saveas(gcf,[root_dir,'/contrast_individual/sub',num2str(s),'_', num2str(i),'.png'])
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
job = nirs.modules.MixedEffects;
job.formula = 'beta ~ -1 + group:cond';
% this would model the two conditions (A and B) per group treating subject as a random variable.
GroupStats = job.run(SubjStats);
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test group_stats
save([root_dir,'/group_stats/group_stats.mat'],'GroupStats');

%% Visualize the group level stats
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test data_group
GroupStats.draw('tstat', [-10 10], 'p < 0.05')
for i = 1:4
    saveas(gcf,[root_dir,'/data_group/',num2str(i),'.png'])
    close;
end
% FDR corrected (q < 0.05)
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test data_FDRcorrected_group
GroupStats.draw('tstat', [-10 10], 'q < 0.05')
for i = 1:4
    saveas(gcf,[root_dir,'/data_FDRcorrected_group/',num2str(i),'.png'])
    close;
end

%% group contrast
c = [-1 1];
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat', [-10 10], 'p < 0.05');
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test contrast_group
for i = 1:2
    saveas(gcf,[root_dir,'/contrast_group/',num2str(i),'.png'])
    close;
end
%% FDR corrected (q < 0.05)
ContrastStats.draw('tstat', [-10 10], 'q < 0.05');
mkdir /Users/hanxu/Downloads/prelim_data_nirs_test contrast_FDRcorrected_group
for i = 1:2
    saveas(gcf,[root_dir,'/contrast_FDRcorrected_group/',num2str(i),'.png'])
    close;
end

%%
% This would be the Fixed Effects version which will run much faster
%job.formula = 'beta ~ -1 + group:cond';
%GroupStatsFE = job.run(SubjStats);   
% ANOVA
%job=nirs.modules.AnovaN;
%job.depvar = 'tstat';
%job.variables={'cond','group'};
%GroupTStats = job.run(SubjStats);
%GroupTStats.draw([],'p<0.05');




