clearvars; clc;

sublist     = {'b4', 'a1', 'b9', 'c2', 'c4', 'c5', 'c6', 'c7'};
nsubjects   = length(sublist);
analysisdir = 'analysis/';

for sId = 1:nsubjects
    cfilename = [sublist{sId} '.probabilities.raw'];
    util_bdisp(['[io] - Importing file ' num2str(sId) '/' num2str(nsubjects) ': ' cfilename '.mat']);
    cdata = load([analysisdir cfilename '.mat']);
    
    util_bdisp(['[io] - Exporting text file: ' cfilename '.txt']);
    edata = [cdata.postprob(:, 1), cdata.Ck, cdata.Mk];
    fId = fopen([analysisdir cfilename '.txt'],'w');
    fprintf(fId,'%7s\t\t%4s\t%8s\n','pp', 'task', 'modality');
    fprintf(fId,'%1.5f\t\t%4d\t%8d\n',edata');
    fclose(fId);
    
end