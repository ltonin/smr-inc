clearvars; clc;

subject = 'c6';
alpha = 0.97;
datadir  = 'analysis/';
datapath = [datadir subject '.probabilities.raw.real.mat'];

%% Force function definition

sx = [0   0.1     0.2      0.4   0.5     0.6     0.8     0.9     1];
sy = [0  -0.03      0     0.01     0   -0.01       0    0.03     0];

degrees = 10;
scoeff  = polyfit(sx, sy, degrees);
F       = @(x, c) polyval(c, x);
U       = @(x, c) -cumsum(F(x, c));


%% Loading posterior probabilities
data = load(datapath);

pp = data.postprob;
nsamples = size(pp, 1);
Tk = data.Tk;
Ck = data.Ck;
Mk = data.Mk;

trialid = setdiff(unique(Tk), 0);
ntrials = max(unique(Tk));

Ydyn = 0.5*ones(nsamples, 1);
Yexp = 0.5*ones(nsamples, 1);

dt      = 0.0625;
for trId = 1:ntrials
    cindex    = Tk == trialid(trId);
    cnsamples = sum(cindex);

    cstart    = find(Tk == trialid(trId), 1, 'first');
    cstop     = find(Tk == trialid(trId), 1, 'last');
    
    
    for nId = cstart+1:cstop
        prevy_dyn = Ydyn(nId-1);                                               % Previous integrated values
        prevy_exp = Yexp(nId-1);
        
        currx_dyn = pp(nId, 1);                                                 % Current random input from distribution
        
        if(max(pp(nId, :)) < 0.55)
            currx_exp = Yexp(nId -1);
        else
            currx_exp  = pp(nId, 1);
        end
        
        curry_dyn = smrinc_dynamic_integrator(currx_dyn, prevy_dyn, scoeff, dt);    % Current integrated Y 
        curry_exp = alpha*prevy_exp + (1-alpha)*currx_exp;
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

% %% Apply integration
% 
% % Time initialization
% time    = [0 6];
% dt      = 0.0625;
% t       = time(1):dt:time(2);
% npoints = length(t);
% it      = 1:npoints;
% 
% % Simulated input distribution
% ldistr = {'Normal', 'BetaHigh', 'BetaLow', 'BetaHighLow'};
% ndistr = length(ldistr);
% input  = zeros(npoints, ndistr);
% 
% % Initial states initialization
% Y0      = 0:0.01:1;
% nstarts = length(Y0);
% Y = zeros(npoints, nstarts, ndistr);
% 
% for dId = 1:ndistr
%     
%     % Getting the current input distribution
%     cdistr = ldistr{dId};
%     input(:, dId) = smrinc_get_distribution(cdistr, npoints);
%     
%     for ss = 1:nstarts
%         Y(1, ss, dId) = Y0(ss);
% 
%         for n = 2:npoints
%             
%             prevy = Y(n-1, ss, dId);                                        % Previous integrated values
%             currx = input(n, dId);                                          % Current random input from distribution
%             curry = smrinc_dynamic_integrator(currx, prevy, scoeff, dt);    % Current integrated Y 
%             Y(n, ss, dId) = curry;
%         end
%     end
%     
% end
% 
% %% Plotting
% 
% % State space, force profile, potential profile
% state     = 0:0.01:1;
% force     = F(state, scoeff);
% potential = U(state, scoeff);
% 
% % Plots
% fig1 = figure;
% fig_set_position(fig1, 'All'); 
% 
% NumRows = 4;
% NumCols = 5;
% 
% % Force plot
% subplot(NumRows, NumCols, [1 NumCols+1]);
% hold on;
% plot(state, force);
% dforce = diff(force);
% plot(sx(sy==0), sy(sy ==0), 'or');
% plot(sx(sy~=0), sy(sy ~=0), 'og');
% xlim([-0.05 1.05])
% hold off;
% plot_hline(0, 'k-');
% plot_vline(0.5, 'k-');
% grid on;
% ylabel('F(x)')
% title('Force');
% 
% % Potential plot
% subplot(NumRows, NumCols, 1 + [2*NumCols 3*NumCols]);
% plot(state, potential);
% xlim([-0.05 1.05])
% plot_vline(0.5, 'k-');
% grid on;
% xlabel('x');
% ylabel('U(x) = - \int_{0}^{1} F(x) dx');
% title('Potential');
% 
% % Dynamic plots
% 
% for dId = 1:ndistr
%     subplot(NumRows, NumCols, [2 7 12] + dId-1);
%     plot(t, Y(:, :, dId));
%     xlabel('time [s]');
%     ylabel('x');
%     xlim([t(1) t(end)]);
%     ylim([-0.1 1.1]);
%     grid on;
%     plot_hline(0.5, 'k-');
%     plot_hline(0, 'k-');
%     plot_hline(1, 'k-');
%     title(ldistr{dId});
%     
%     subplot(NumRows, NumCols, 17 + dId -1);
%     edges = 0:0.05:1;
%     cnts  = histc(input(:, dId), 0:0.05:1);
%     bar(edges, 100*cnts./sum(cnts), 'histc');
%     xlim([0 1]);
%     ylim([0 30]);
%     grid on;
%     xlabel('x');
%     ylabel('%');
% end
% 
% suptitle('Dynamic integration study - random input');


