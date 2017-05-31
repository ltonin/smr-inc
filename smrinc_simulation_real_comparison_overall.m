clearvars; clc;

sublist = {'a1', 'b4', 'b9', 'c4', 'c5', 'c6'};
nsubjects = length(sublist);
datadir  = 'analysis/';

eVAL =[]; eDUR = []; eTYP = [];
dVAL =[]; dDUR = []; dTYP = [];
aVAL =[]; aDUR = []; aTYP = [];
TYP = []; POS = []; DUR = []; ODT = []; MOD = []; SUB = [];
for sId = 1:nsubjects
    csubject = sublist{sId};
    util_bdisp(['[io] - Loading and concatenating data for subject: ' csubject]);
    
    filename = [datadir '/' csubject '.integrator.comparison.mat'];
    disp(['       Filename: ' filename]);
    cdata = load(filename);
    
    
    eVAL = cat(1, eVAL, cdata.cmd.exp.VAL);
    eDUR = cat(1, eDUR, cdata.cmd.exp.DUR);
    eTYP = cat(1, eTYP, cdata.cmd.exp.TYP);
    
    dVAL = cat(1, dVAL, cdata.cmd.dyn.VAL);
    dDUR = cat(1, dDUR, cdata.cmd.dyn.DUR);
    dTYP = cat(1, dTYP, cdata.cmd.dyn.TYP);
    
    aVAL = cat(1, aVAL, cdata.cmd.alp.VAL);
    aDUR = cat(1, aDUR, cdata.cmd.alp.DUR);
    aTYP = cat(1, aTYP, cdata.cmd.alp.TYP);
    
    TYP = cat(1, TYP, cdata.evt.TYP);
    POS = cat(1, POS, cdata.evt.POS);
    DUR = cat(1, DUR, cdata.evt.DUR);
    ODT = cat(1, ODT, cdata.evt.ODT);
    MOD = cat(1, MOD, cdata.evt.MOD);
    SUB = cat(1, SUB, sId*ones(length(MOD), 1));
end

cmd.exp.VAL = eVAL;
cmd.exp.DUR = eDUR;
cmd.exp.TYP = eTYP;

cmd.dyn.VAL = dVAL;
cmd.dyn.DUR = dDUR;
cmd.dyn.TYP = dTYP;

cmd.alp.VAL = aVAL;
cmd.alp.DUR = aDUR;
cmd.alp.TYP = aTYP;

evt.TYP = TYP;
evt.POS = POS;
evt.DUR = DUR;
evt.ODT = ODT;
evt.MOD = MOD;
evt.SUB = SUB;

icIdx   = evt.MOD == 1 & evt.TYP ~= 783;
incIdx  = evt.MOD == 1 & evt.TYP == 783;

%% Correct delivery with rejection
% Exponential
nd_e  = isnan(cmd.exp.TYP);
cor_e = 100*sum(cmd.exp.TYP(icIdx & nd_e == false) == evt.TYP(icIdx & nd_e == false))./sum(icIdx & nd_e == false);
rej_e = 100*sum(nd_e(icIdx))./length(nd_e(icIdx));

% Dynamic
nd_d  = isnan(cmd.dyn.TYP);
cor_d = 100*sum(cmd.dyn.TYP(icIdx & nd_d == false) == evt.TYP(icIdx & nd_d == false))./sum(icIdx & nd_d == false);
rej_d = 100*sum(nd_d(icIdx))./length(nd_d(icIdx));

% Alpha
nd_a  = isnan(cmd.alp.TYP);
cor_a = 100*sum(cmd.alp.TYP(icIdx & nd_a == false) == evt.TYP(icIdx & nd_a == false))./sum(icIdx & nd_a == false);
rej_a = 100*sum(nd_a(icIdx))./length(nd_a(icIdx));

%% Durations

nd_e = isnan(cmd.exp.TYP);
nd_d = isnan(cmd.dyn.TYP);
nd_a = isnan(cmd.alp.TYP);

dur_e = cmd.exp.DUR;
dur_d = cmd.dyn.DUR;
dur_a = cmd.alp.DUR;

%% Plots
fig1 = figure;
fig_set_position(fig1, 'Top');
dt = 0.0625;

subplot(1, 3, 1);
bar([cor_e rej_e; cor_d rej_d; cor_a rej_a]);
grid on;
plot_hline(100, 'k-');
ylim([0 120]);
ylabel('%');
grid on;
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Exp', 'Dyn', 'Alp'});
legend('Correct delivered', 'Not delivered');
title('IC trials');


subplot(1, 3, 2);
icdur_d = [dur_e(icIdx & nd_e == false); dur_d(icIdx & nd_d == false); dur_a(icIdx & nd_a == false)];
icdur_l = [1*ones(sum(icIdx & nd_e == false), 1); 2*ones(sum(icIdx & nd_d == false), 1); 3*ones(sum(icIdx & nd_a == false), 1)];
boxplot(icdur_d*dt, icdur_l, 'DataLim', [0 10], 'extrememode', 'compress');
ylabel('Time [s]');
grid on;
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Exp', 'Dyn', 'Alp'});
title('IC trials');

subplot(1, 3, 3);
incdur_d = [dur_e(incIdx); dur_d(incIdx); dur_a(incIdx)];
incdur_l = [1*ones(sum(incIdx), 1); 2*ones(sum(incIdx), 1); 3*ones(sum(incIdx), 1)];
boxplot(incdur_d*dt, incdur_l, 'DataLim', [0 10], 'extrememode', 'compress');
ylabel('Time [s]');
grid on;
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Exp', 'Dyn', 'Alp'});
title('INC trials');

suptitle(['Overall']);