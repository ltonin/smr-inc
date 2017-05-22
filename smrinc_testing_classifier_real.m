clearvars; clc; 

sublist     = {'b4', 'a1', 'b9', 'c2', 'c4', 'c5', 'c6', 'c7'};
subdates    = {'20160301', '20160310', '20160316', '20160321', '20160324', '20160329', '20160330', '20160331'};
nsubjects   = length(sublist);

pattern      = 'mi_hfr';
extension    = '.mat';
analysisdir  = 'analysis/';
selectedtask = [771 773];
FreqGrid     = 4:2:48;
ChanGrid     = 1:16;

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
    
    % Get classifier
    classifpath = [analysisdir sublist{sId} '_bhbf_' subdates{sId} '.mat'];
    disp(['     - Loading classifier: ' classifpath])
    cclassif = load(classifpath);
    
    % Concatenating data
    [psd, events, labels, settings] = smrinc_concatenate_data(analysisdir, sublist{sId}, pattern, extension);
    F = log(psd);
    nsamples = size(psd, 1);
   
    % Generate event labels
    [~, CfbEvents] = proc_get_event2(781, nsamples, events.POS, events.TYP, events.DUR);
    [~, CueEvents] = proc_get_event2([771 773 783], nsamples, events.POS, events.TYP, events.DUR);
    ntrials = length(CfbEvents.TYP);
    
    TrialLb = zeros(nsamples, 1);
    TrialId = zeros(nsamples, 1);
    for eId = 1:ntrials
        cstart = CfbEvents.POS(eId);
        cstop  = cstart + CfbEvents.DUR(eId) - 1;
        TrialLb(cstart:cstop) = CueEvents.TYP(eId);
        TrialId(cstart:cstop) = eId;
    end
    
    % Import classifier settings
    util_bdisp('[proc] - Importing classifier parameters');
    FeatureIdx = smrinc_eegc32index_features(cclassif.analysis.tools.features.bands, 4:2:48, 1:16);
    gau = cclassif.analysis.tools.net.gau;
    
    % Reshaping data
    D  = proc_reshape_ts_bc(F);
    T  = D(:, FeatureIdx);
    
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
    Tk = TrialId;

    filename = [analysisdir sublist{sId} '.probabilities.raw.real.mat'];
    util_bdisp(['[out] - Saving raw probabilities in: ' filename]);
    save(filename, 'postprob', 'Mk', 'Ck', 'Sl', 'Tk');
    
    
%     Tasks = [773 783 771]; 
%     TaskLb = {'Left', 'Rest', 'Right'};
%     h = zeros(length(Tasks), 1);
% 
%     for i = 1:length(Tasks);
%         h(i) = subplot(nsubjects, 3, i); 
%         hist(postprob(Ck == Tasks(i), 1), 100); 
%         title([Sl '-' TaskLb{i}]);
%     end; 
%     plot_set_limits(h, 'y', 'minmax');
        
end