clearvars; clc;

subject = 'c6';
datadir  = 'analysis/';
datapath = [datadir subject '.probabilities.raw.real.mat'];

%% Force function definition

sx = [0   0.1     0.2      0.4   0.5     0.6     0.8     0.9     1];
sy = [0  -0.03      0     0.01     0   -0.01       0    0.03     0];

degrees = 10;
scoeff  = polyfit(sx, sy, degrees);

%% Exponential smoothing definition
alpha = 0.97;
rejection = 0.55;


%% Loading posterior probabilities
data = load(datapath);

pp = data.postprob;
nsamples = size(pp, 1);
Tk = data.Tk;
Ck = data.Ck;
Mk = data.Mk;

trialid = setdiff(unique(Tk), 0);
ntrials = max(unique(Tk));

Yd = 0.5*ones(nsamples, 1);
Ye = 0.5*ones(nsamples, 1);

dt      = 0.0625;
for trId = 1:ntrials
    cindex    = Tk == trialid(trId);
    cnsamples = sum(cindex);

    cstart    = find(Tk == trialid(trId), 1, 'first');
    cstop     = find(Tk == trialid(trId), 1, 'last');
    
    
    for nId = cstart+1:cstop
        currx = pp(nId, 1);
        
        prevy_d = Yd(nId-1);                                               % Previous integrated values
        prevy_e = Ye(nId-1);
        
                                                         % Current random input from distribution
     
        
        curry_d = smrinc_integrator_dynamic(currx, prevy_d, scoeff, dt);    % Current integrated Y 
        curry_e = alpha*prevy_exp + (1-alpha)*currx_exp;
        Ydyn(nId) = curry_dyn;
        Yexp(nId) = curry_exp;
    end
end

thexp = 0.7;
thdyn = 0.7;
Dexp = nan(ntrials, 1);
Ddyn = nan(ntrials, 1);
TCk = zeros(ntrials, 1);
TMk = zeros(ntrials, 1);
for trId = 1:ntrials
    cindex    = Tk == trialid(trId);
    
    cyexp = Yexp(cindex);
    cydyn = Ydyn(cindex);
    
    
    othexp = find(cyexp >= thexp | cyexp <= 1-thexp, 1, 'first');
    othdyn = find(cydyn >= thdyn | cydyn <= 1-thdyn, 1, 'first');
    
    if(isempty(othexp) == false)
        Dexp(trId) = othexp;
%     else
%         [~, Dexp(trId)] = max(abs(cyexp - 0.5));
    end
    
    if(isempty(othdyn) == false)
        Ddyn(trId) = othdyn;
%     else
%         [~, Ddyn(trId)] = max(abs(cydyn - 0.5));
    end
    
    TMk(trId) = unique(Mk(cindex));
    TCk(trId) = unique(Ck(cindex));
end
    

ICIdx  = TCk == 771 | TCk == 773;
INCIdx = TCk == 783;

DTexp = Dexp*dt;
DTdyn = Ddyn*dt;

MdExp_IC = nanmedian(DTexp(ICIdx & TMk == 1));
SeExp_IC = nanstd(DTexp(ICIdx & TMk == 1))./sqrt(sum(ICIdx & TMk == 1));
MdDyn_IC = nanmedian(DTdyn(ICIdx & TMk == 1));
SeDyn_IC = nanstd(DTdyn(ICIdx & TMk == 1))./sqrt(sum(ICIdx & TMk == 1));

MdExp_INC = nanmedian(DTexp(INCIdx & TMk == 1));
SeExp_INC = nanstd(DTexp(INCIdx & TMk == 1))./sqrt(sum(INCIdx & TMk == 1));
MdDyn_INC = nanmedian(DTdyn(INCIdx & TMk == 1));
SeDyn_INC = nanstd(DTdyn(INCIdx & TMk == 1))./sqrt(sum(INCIdx & TMk == 1));

subplot(1, 2, 1);
plot_barerrors([MdExp_IC MdDyn_IC], [SeExp_IC SeDyn_IC]);
h(1) = gca;
grid on;
ylabel('Time [s]');
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Exp', 'Dyn'});
title('IC trials');

subplot(1, 2, 2);
plot_barerrors([MdExp_INC MdDyn_INC], [SeExp_INC SeDyn_INC]);
h(2) = gca;
grid on;
ylabel('Time [s]');
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Exp', 'Dyn'});
title('INC trials');

plot_set_limits(h, 'y', [0 5]);
suptitle(['Subject ' subject]);


