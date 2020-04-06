% 1. the order of stim6 and stim 7 are different
% 2. the number of stim6 in subj 13 is 12, in subj 18 is 9
% 3. use the sample data to test whether the script it robust
clear;
clc;

%% group number
group_number = 2;

%% load data
root_dir = '/Users/hanxu/Downloads/fnirs_test'
raw = nirs.io.loadDirectory(root_dir, {'group', 'subject'});

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

%% Signal Preprocessing -> PCA vs. Wavelet vs. Spline
mkdir /Users/hanxu/Downloads/fnirs_test preprocessing
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
mkdir /Users/hanxu/Downloads/fnirs_test conversion
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
%% the parameters can be modified in nirs.design.basis.Canonical
%% save GLM model
mkdir /Users/hanxu/Downloads/fnirs_test GLM
for s = 1:size(transpose(SubjStats))
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.txt'], 'Delimiter', ' '); 
end

%% data view on individual level
mkdir /Users/hanxu/Downloads/fnirs_test data_individual
for s = 1:size(transpose(SubjStats))
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

%% Group Analysis -> wired results, take a look at this
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




