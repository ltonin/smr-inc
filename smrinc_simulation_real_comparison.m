clearvars; clc;

subject = 'b4';
datadir  = 'analysis/';
datapath = [datadir subject '.probabilities.raw.real.mat'];

%% Dynamic integrator parameters
% sx = [0   0.1     0.2      0.3   0.5     0.7     0.8     0.9     1];
% sy = [0  -0.03      0     0.01     0   -0.01       0    0.03     0];
% 
% degrees = 8;
% 
% th = 0.7;
% x1 = (1-th):0.01:th;
% y1 = -0.1*(x1-0.5);
% 
% sx = [0  (1-th)/2 x1 th+(1-th)/2  1];
% sy = [0 -0.05 y1 0.05 0];
% degrees = 10;
% 
% 
% % sx = [0   0.1     0.15     0.2  0.35  0.425  0.5   0.575    0.65   0.8     0.85     0.9     1];
% % sy = [0  -0.03      0      0.02  0.01   0.005  0   -0.005   -0.01  -0.02        0    0.03     0];
% % 
% % degrees = 10;
% 
% scoeff  = polyfit(sx, sy, degrees);

inclim = 0.6;
nrpt   = 0.8;
[scoeff, support] = smrinc_get_forceprofile(inclim, nrpt);

switch(subject)
    case 'a1'
        phi  = 0.8;
        rho   = 0.5;
        gamma = 0.6;
    case 'b4'
        phi = 0.8;
        rho   = 0.1;
        gamma = 0.2;
    case 'b9'
        phi = 0.8;
        rho   = 0.2;
        gamma = 0.3;
    case 'c4'
        phi = 0.75;
        rho   = 0.2;
        gamma = 0.25;
    case 'c5'
        phi = 0.8;
        rho   = 0.3;
        gamma = 0.3;
    case 'c6'
        phi = 0.5;
        rho   = 0.3;
        gamma = 0.2;
end

%% Exponential integrator parameters
alpha     = 0.97;
if strcmp(subject, 'b4')
    alpha = 0.98;
end
rejection = 0.55;

%% Alpha integrator parameters
% rho   = 0.3;
% gamma = 0.2;

%% Loading posterior probabilities
data = load(datapath);

pp = data.postprob;
trialId = setdiff(unique(data.Tk), 0);
ntrials = max(unique(data.Tk));

dt      = 0.0625;
timeout = 10;

%% Extract event info
evt.TYP = zeros(ntrials, 1);
evt.POS = zeros(ntrials, 1);
evt.DUR = zeros(ntrials, 1);
evt.ODT = zeros(ntrials, 1);
evt.MOD = zeros(ntrials, 1);
rpp = [];
for trId = 1:ntrials
    cstart = find(data.Tk == trialId(trId), 1, 'first');
    cstop  = find(data.Tk == trialId(trId), 1, 'last');
    cCk    = unique(data.Ck(cstart:cstop));
    cMk    = unique(data.Mk(cstart:cstop));
    tdur   = length(cstart:cstop); 
    tpp    = pp(cstart:cstop, 1);
    
    app    = repmat(tpp(randperm(tdur)), [ceil(floor(timeout/dt)/tdur) 1]);
   % app    = app(randperm(length(app)));    
    cpp    = cat(1, tpp, app);
    
    evt.POS(trId) = length(rpp) + 1;
    evt.DUR(trId) = length(cpp);
    evt.TYP(trId) = cCk;
    evt.ODT(trId) = tdur;
    evt.MOD(trId) = cMk;
    rpp    = cat(1, rpp, cpp);
end


%% Simulate integration frameworks
nsamples = size(rpp, 1);
Yd = 0.5*ones(nsamples, 1);
Ye = 0.5*ones(nsamples, 1);
Ya = 0.5*ones(nsamples, 1);

trialId = 1;
for sId = 1:nsamples
    
    if(trialId < ntrials)
        npos = evt.POS(trialId);
        if(sId == npos) 
            prevy_e = 0.5;
            prevy_d = 0.5;
            prevy_a = 0.5;
            trialId = trialId+1;
        end
    end
    
    currx = rpp(sId);
    curry_e = smrinc_integrator_exponential(currx, prevy_e, alpha, rejection);
    curry_d = smrinc_integrator_dynamic(currx, prevy_d, scoeff, phi, dt);  
    curry_a = smrinc_integrator_alpha(currx, prevy_a, rho, gamma, dt);
    Ye(sId) = curry_e;
    Yd(sId) = curry_d;
    Ya(sId) = curry_a;
    
    prevy_e = curry_e;
    prevy_d = curry_d;
    prevy_a = curry_a;
end




The = 0.7;
Thd = 0.7;
Tha = 0.7;

cmd.exp = [];
cmd.dyn = [];
cmd.alp = [];
classes = [771 773];
for trId = 1:ntrials
    
    cstart = evt.POS(trId);
    cstop  = cstart + evt.DUR(trId) - 1;
    
    cinte = Ye(cstart:cstop, 1);
    cintd = Yd(cstart:cstop, 1);
    cinta = Ya(cstart:cstop, 1);
    
    % Exponential
    id_e  = find(cinte <= (1-The) | cinte >= The, 1, 'first');
    if(isempty(id_e))
        cmd.exp.VAL(trId, 1) = nan;
        cmd.exp.TYP(trId, 1) = nan;
        cmd.exp.DUR(trId, 1) = length(cinte);
    else
        cmd.exp.VAL(trId, 1) = cinte(id_e);
        cmd.exp.TYP(trId, 1) = classes(double(cinte(id_e)-0.5>0)+1);
        cmd.exp.DUR(trId, 1) = id_e;
    end
    
    % Dynamic
    id_d  = find(cintd <= (1-Thd) | cinte >= Thd, 1, 'first');
    if(isempty(id_d))
        cmd.dyn.VAL(trId, 1) = nan;
        cmd.dyn.TYP(trId, 1) = nan;
        cmd.dyn.DUR(trId, 1) = length(cintd);
    else
        cmd.dyn.VAL(trId, 1) = cintd(id_d);
        cmd.dyn.TYP(trId, 1) = classes(double(cintd(id_d)-0.5>0)+1);
        cmd.dyn.DUR(trId, 1) = id_d;
    end
    
    % Alpha
    id_a  = find(cinta <= (1-Tha) | cinta >= Tha, 1, 'first');
    if(isempty(id_a))
        cmd.alp.VAL(trId, 1) = nan;
        cmd.alp.TYP(trId, 1) = nan;
        cmd.alp.DUR(trId, 1) = length(cinta);
    else
        cmd.alp.VAL(trId, 1) = cinta(id_a);
        cmd.alp.TYP(trId, 1) = classes(double(cinta(id_a)-0.5>0)+1);
        cmd.alp.DUR(trId, 1) = id_a;
    end
 
end

icIdx   = evt.MOD == 1 & evt.TYP ~= 783;
incIdx  = evt.MOD == 1 & evt.TYP == 783;

%% Check original duration w/r exponential smoothing
fprintf('Check duration ratio: %3.2f%s%3.2f\n', mean(cmd.exp.DUR(icIdx)./evt.ODT(icIdx)),char(177), std(cmd.exp.DUR(icIdx)./evt.ODT(icIdx)));

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

cmd.exp.th = The;
cmd.dyn.th = Thd;
cmd.alp.th = Tha;
cmd.exp.alpha = alpha;
cmd.exp.rejection = rejection;
cmd.dyn.phi = phi;
cmd.dyn.coeff = scoeff;
cmd.alp.rho = rho;
cmd.alp.gamma = gamma;

filename = [datadir '/' subject '.integrator.comparison.mat'];
util_bdisp(['Saving results in: ' filename]);
save(filename, 'cmd', 'evt', 'subject');


%% Plots
fig1 = figure;
fig_set_position(fig1, 'Top');

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

suptitle(['Subject ' subject]);
% subplot(1, 2, 2);
% bar(incratios);
% grid on;
% plot_hline(100, 'k-');
% ylim([0 200]);

% fig1 = figure;
% subplot(1, 2, 1);
% cindex = Mk == 1 & Ck ~= 783;
% icdata = 100*[sum(Re(cindex) == false)./sum(cindex) sum(Rd(cindex) == false)./sum(cindex) sum(Ra(cindex) == false)./sum(cindex)];
% plot_barerrors(mean(icdata, 1), std(icdata, [], 1)./sqrt(size(icdata,1)));
% grid on;
% h(1) = gca;
% ylabel('%');
% set(gca, 'XTick', [1 2 3]);
% set(gca, 'XTickLabel', {'Exp', 'Dyn', 'Alpha'});
% title('IC trials');
% 
% 
% subplot(1, 2, 2);
% cindex = Mk == 1 & Ck == 783;
% incdata = 100*[sum(Re(cindex) == false)./sum(cindex) sum(Rd(cindex) == false)./sum(cindex) sum(Ra(cindex) == false)./sum(cindex)];
% plot_barerrors(mean(incdata, 1), std(incdata, [], 1)./sqrt(size(incdata,1)));
% grid on;
% h(2) = gca;
% ylabel('%');
% set(gca, 'XTick', [1 2 3]);
% set(gca, 'XTickLabel', {'Exp', 'Dyn', 'Alpha'});
% title('INC trials');
% 
% plot_set_limits(h, 'y', [0 120]);

% fig1 = figure;
% subplot(1, 2, 1);
% cindex = evt.MOD == 1 & evt.TYP ~= 783;
% icdata = 100*[Dd(cindex)./De(cindex) Da(cindex)./De(cindex)];
% plot_barerrors(nanmean(icdata, 1), nanstd(icdata, [], 1)./sqrt(size(icdata,1)));
% plot_hline(100, 'k-');
% grid on;
% h(1) = gca;
% ylabel('% vs. Exponential');
% set(gca, 'XTick', [1 2]);
% set(gca, 'XTickLabel', {'Dynamic', 'Alpha'});
% title('IC trials');
% 
% subplot(1, 2, 2);
% cindex = evt.MOD == 1 & evt.TYP == 783;
% incdata = 100*[Dd(cindex)./De(cindex) Da(cindex)./De(cindex)];
% plot_barerrors(nanmean(incdata, 1), nanstd(incdata, [], 1)./sqrt(size(incdata,1)));
% plot_hline(100, 'k-');
% grid on;
% h(2) = gca;
% ylabel('% vs. Exponential');
% set(gca, 'XTick', [1 2]);
% set(gca, 'XTickLabel', {'Dynamic', 'Alpha'});
% title('INC trials');
% 
% plot_set_limits(h, 'y', [0 200]);

% thexp = 0.7;
% thdyn = 0.7;
% Dexp = nan(ntrials, 1);
% Ddyn = nan(ntrials, 1);
% TCk = zeros(ntrials, 1);
% TMk = zeros(ntrials, 1);
% for trId = 1:ntrials
%     cindex    = Tk == trialid(trId);
%     
%     cyexp = Yexp(cindex);
%     cydyn = Ydyn(cindex);
%     
%     
%     othexp = find(cyexp >= thexp | cyexp <= 1-thexp, 1, 'first');
%     othdyn = find(cydyn >= thdyn | cydyn <= 1-thdyn, 1, 'first');
%     
%     if(isempty(othexp) == false)
%         Dexp(trId) = othexp;
% %     else
% %         [~, Dexp(trId)] = max(abs(cyexp - 0.5));
%     end
%     
%     if(isempty(othdyn) == false)
%         Ddyn(trId) = othdyn;
% %     else
% %         [~, Ddyn(trId)] = max(abs(cydyn - 0.5));
%     end
%     
%     TMk(trId) = unique(Mk(cindex));
%     TCk(trId) = unique(Ck(cindex));
% end
%     
% 
% ICIdx  = TCk == 771 | TCk == 773;
% INCIdx = TCk == 783;
% 
% DTexp = Dexp*dt;
% DTdyn = Ddyn*dt;
% 
% MdExp_IC = nanmedian(DTexp(ICIdx & TMk == 1));
% SeExp_IC = nanstd(DTexp(ICIdx & TMk == 1))./sqrt(sum(ICIdx & TMk == 1));
% MdDyn_IC = nanmedian(DTdyn(ICIdx & TMk == 1));
% SeDyn_IC = nanstd(DTdyn(ICIdx & TMk == 1))./sqrt(sum(ICIdx & TMk == 1));
% 
% MdExp_INC = nanmedian(DTexp(INCIdx & TMk == 1));
% SeExp_INC = nanstd(DTexp(INCIdx & TMk == 1))./sqrt(sum(INCIdx & TMk == 1));
% MdDyn_INC = nanmedian(DTdyn(INCIdx & TMk == 1));
% SeDyn_INC = nanstd(DTdyn(INCIdx & TMk == 1))./sqrt(sum(INCIdx & TMk == 1));
% 
% subplot(1, 2, 1);
% plot_barerrors([MdExp_IC MdDyn_IC], [SeExp_IC SeDyn_IC]);
% h(1) = gca;
% grid on;
% ylabel('Time [s]');
% set(gca, 'XTick', [1 2]);
% set(gca, 'XTickLabel', {'Exp', 'Dyn'});
% title('IC trials');
% 
% subplot(1, 2, 2);
% plot_barerrors([MdExp_INC MdDyn_INC], [SeExp_INC SeDyn_INC]);
% h(2) = gca;
% grid on;
% ylabel('Time [s]');
% set(gca, 'XTick', [1 2]);
% set(gca, 'XTickLabel', {'Exp', 'Dyn'});
% title('INC trials');
% 
% plot_set_limits(h, 'y', [0 5]);
% suptitle(['Subject ' subject]);


