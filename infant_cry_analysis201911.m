% need explore:
% 1. ANOVA group analysis
% 2. the formula
% 3. GroupStatsME and GroupStatsFE, and how to draw them

clear;
clc;

%% load data
root_dir = '/Users/hanxu/Downloads/prelim_data_nirs'
raw = nirs.io.loadDirectory(root_dir, {'group', 'subject'});

%% Prepare the data for preprocessing
% Remove stimless files
j = nirs.modules.RemoveStimless();

% Just keep stim events for infant cry (these will change depending on the
% task being analyzed)

j = nirs.modules.KeepStims( j);
j.listOfStims = {'stim_channel6', 'stim_channel7'};

% Trim pre and post baseline
j = nirs.modules.TrimBaseline( j);
j.preBaseline = 1; %can change these values
j.postBaseline = 1;

%% Run the conversions 
j = nirs.modules.OpticalDensity( j);

% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);

% Run the GLM

hb = j.run( raw );
j = nirs.modules.GLM();
j.verbose = true;
% This example specifies DCT terms with a frequency cutoff of 0.08.
j.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);
j.basis = Dictionary();
j.basis('default') = nirs.design.basis.Canonical();
%j.basis('X') = nirs.design.basis.Canonical();
%j.basis('Y') = nirs.design.basis.Canonical();
j.basis('stim_channel6') = nirs.design.basis.Canonical();
j.basis('stim_channel7') = nirs.design.basis.Canonical();
SubjStats = j.run( hb );


%% Create a table of the output of the first level GLM for each participant
hb_table1 = SubjStats(1).table;
hb_table2 = SubjStats(2).table;
hb_table3 = SubjStats(3).table;
hb_table4 = SubjStats(4).table;
hb_table5 = SubjStats(5).table;
hb_table6 = SubjStats(6).table;
hb_table7 = SubjStats(7).table;
hb_table8 = SubjStats(8).table;
hb_table9 = SubjStats(9).table;
hb_table10 = SubjStats(10).table;
hb_table11 = SubjStats(11).table;
hb_table12 = SubjStats(12).table;
hb_table13 = SubjStats(13).table;
hb_table14 = SubjStats(14).table;
hb_table15 = SubjStats(15).table;
hb_table16 = SubjStats(16).table;
hb_table17 = SubjStats(17).table;
hb_table18 = SubjStats(18).table;
hb_table19 = SubjStats(19).table;
hb_table20 = SubjStats(20).table;

%% Display the contrasts (individual level)
SubjStats(4).draw('tstat', [-10 10], 'p < 0.05');

%% Group Analysis
job = nirs.modules.MixedEffects;
job.formula = 'beta ~ -1 + group:cond';
% this would model the two conditions (A and B) per group treating subject
% as a random variable. Visit number is ignored.  
GroupStats = job.run(SubjStats);
GroupStats.draw('tstat', [-10 10], 'p < 0.05');

% This would be the Fixed Effects version which will run much faster
%job.formula = 'beta ~ -1 + group:cond';
%GroupStatsFE = job.run(SubjStats);   
% ANOVA
%job=nirs.modules.AnovaN;
%job.depvar = 'tstat';
%job.variables={'cond','group'};
%GroupTStats = job.run(SubjStats);
%GroupTStats.draw([],'p<0.05');




