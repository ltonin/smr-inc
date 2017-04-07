clearvars; clc;

sublist     = {'b4', 'a1', 'b9', 'c2', 'c4', 'c5', 'c6', 'c7'};
nsubjects   = length(sublist);

pattern      = 'mi_hfr';
extension    = '.mat';
analysisdir  = 'analysis/';
selectedtask = [771 773];
NumSelFeatures = 6;

% eegc3 style settings
eegc3.modules.smr.gau.somunits 	= [1 1]; % QDA-style
eegc3.modules.smr.gau.sharedcov = 'f'; % No difference anyway
eegc3.modules.smr.gau.epochs 	= 4;
eegc3.modules.smr.gau.mimean	= 0.01;
eegc3.modules.smr.gau.micov		= 0.001;
eegc3.modules.smr.gau.th        = 0.70;
eegc3.modules.smr.gau.terminate	= true;


%% Training classifier for all subjects
for sId = 1:nsubjects
    util_bdisp(['[io] - Subject ' num2str(sId) '/' num2str(nsubjects) ': ' sublist{sId}]);
    
    % Concatenating data
    [psd, events, labels, settings] = smrinc_concatenate_data(analysisdir, sublist{sId}, pattern, extension);
    F = log(psd);
    nsamples = size(psd, 1);
    
    % Generate event labels
    [~, CfbEvents] = proc_get_event2(781, nsamples, events.POS, events.TYP, events.DUR);
    [~, CueEvents] = proc_get_event2([771 773 783], nsamples, events.POS, events.TYP, events.DUR);
    ntrials = length(CfbEvents.TYP);
    
    TrialLb = zeros(nsamples, 1);
    for eId = 1:ntrials
        cstart = CfbEvents.POS(eId);
        cstop  = cstart + CfbEvents.DUR(eId) - 1;
        TrialLb(cstart:cstop) = CueEvents.TYP(eId);
    end
    
    % Select the classes
    TaskIdx = false(nsamples, 1);
    for cId = 1:length(selectedtask)
        TaskIdx = TaskIdx | TrialLb == selectedtask(cId);
    end
    
    % Select the modality (offline)
    ModalityIdx = labels.Mk == 0;
    
    % General indexes
    Indexes = ModalityIdx & TaskIdx;
    
    % Compute discriminancy over the selected tasks
    util_bdisp('[proc] - Computing fisher score');
    fs = proc_fisher2(F(Indexes, :, :), TrialLb(Indexes));
    
    % Rank fs and select features
    util_bdisp(['[proc] - Select the best ' num2str(NumSelFeatures) ' features']);
    [~, FeaturesRankIdx] =  sort(fs, 'descend');
    FeatureIdx = FeaturesRankIdx(1:NumSelFeatures);
    
    % Reshape data
    D  = proc_reshape_ts_bc(F);
    P  = D(Indexes, FeatureIdx);
    Pk = TrialLb(Indexes);
    T  = D(:, FeatureIdx);
    
    % Train classifier
    util_bdisp('[proc] - Train classifier');
    [gau, performace] = eegc3_train_gau(eegc3, P, Pk);
    
    % Evaluate classifier
    util_bdisp('[proc] - Evaluate classifier');
    postprob = zeros(nsamples, 2);
    for tId = 1:nsamples
        [~, postprob(tId, :)] = gauClassifier(gau.M, gau.C, T(tId, :));
    end
    
    % Save analysis
    Mk = labels.Mk;
    Ck = TrialLb;
    Sl = sublist{sId};
    filename = [analysisdir sublist{sId} '.probabilities.raw.mat'];
    util_bdisp(['[out] - Saving raw probabilities in: ' filename]);
    save(filename, 'postprob', 'Mk', 'Ck', 'Sl');
    
    
%     Tasks = [773 783 771]; 
%     TaskLb = {'Left', 'Rest', 'Right'};
%     h = zeros(length(Tasks), 1);
% 
%     for i = 1:length(Tasks);
%         h(i) = subplot(1, 3, i); 
%         hist(postprob(Ck == Tasks(i), 1), 100); 
%         title([Sl '-' TaskLb{i}]);
%     end; 
%     plot_set_limits(h, 'y', 'minmax');
        
end