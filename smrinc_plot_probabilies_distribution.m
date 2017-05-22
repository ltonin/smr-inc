clearvars; clc;

sublist     = {'b4', 'a1', 'b9', 'c2', 'c4', 'c5', 'c6', 'c7'};
nsubjects   = length(sublist);
analysisdir  = 'analysis/';
SelTasks   = [771 783 773];
SelTaskLbs = {'BothFeet', 'Rest', 'BothHands'};
ColorTasks = {'b', 'g', 'r'};
WidthModality = {1, 2};
StyleModality = {':', '-'};
NumCols = 3;

fig1 = figure;
fig_set_position(fig1, 'All');

for sId = 1:nsubjects
    util_bdisp(['[io] - Subject ' num2str(sId) '/' num2str(nsubjects) ': ' sublist{sId}]);
    
    % Load raw probabilities for current subject
    cdata = load([analysisdir '/' sublist{sId} '.probabilities.raw.mat']);
    
    Ck = cdata.Ck;
    Mk = cdata.Mk;
    pp = cdata.postprob;
    
    Mods    = unique(Mk);
    NumMods = length(Mods);
    
    subplot(ceil(nsubjects/NumCols), NumCols, sId);
    hold on;
    for mId = 1:NumMods
        for tId = 1:length(SelTasks)
            cpp = pp(Ck == SelTasks(tId) & Mk == Mods(mId), 2);
            [b, x] = hist(cpp, 50);
            plot(x, 100*b./sum(b), ColorTasks{tId}, 'LineStyle', StyleModality{mId}, 'LineWidth', WidthModality{mId});
        end
    end
    legend(SelTaskLbs);
    hold off;
    %ylim([0 100]);
    grid on;
    title(['Subject ' cdata.Sl]);
    
    
    
end