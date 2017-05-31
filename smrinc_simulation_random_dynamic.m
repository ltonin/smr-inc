clearvars; clc;


%% Force function definition
% sx = [0   0.1     0.2      0.3   0.5     0.7     0.8     0.9     1];
% sy = [0  -0.03      0     0.01     0   -0.01       0    0.03     0];
% degrees = 8;
% 
% th = 0.55;
% x1 = (1-th):0.01:th;
% y1 = -0.5*(x1-0.5);
% 
% sx = [0  (1-th)/2  x1 th+(1-th)/2  1];
% sy = [0 -0.05 y1 0.05 0];
% degrees = 10;
% 
% 
% % sx = [0   0.1     0.15     0.2  0.35  0.425  0.5   0.575    0.65   0.8     0.85     0.9     1];
% % sy = [0  -0.03      0      0.02  0.01   0.005  0   -0.005   -0.01  -0.02        0    0.03     0];
% % 
% % degrees = 10;
% scoeff  = polyfit(sx, sy, degrees);
inclim = 0.65;
nrpt   = 0.7;
[scoeff, support] = smrinc_get_forceprofile(inclim, nrpt);
F       = @(x, c) polyval(c, x);
U       = @(x, c) -cumsum(F(x, c));

%% Dynamic integrator parameters
phi  = 0.8; 


%% Apply integration

% Time initialization
time    = [0 6];
dt      = 0.0625;
t       = time(1):dt:time(2);
npoints = length(t);
it      = 1:npoints;

% Simulated input distribution
ldistr = {'Normal', 'BetaHigh', 'BetaLow', 'BetaHighLow'};
ndistr = length(ldistr);
input  = zeros(npoints, ndistr);

% Initial states initialization
Y0      = 0:0.01:1;
nstarts = length(Y0);
Y = zeros(npoints, nstarts, ndistr);
Yf = zeros(npoints, nstarts, ndistr);
Yb = zeros(npoints, nstarts, ndistr);

for dId = 1:ndistr
    
    % Getting the current input distribution
    cdistr = ldistr{dId};
    input(:, dId) = smrinc_get_distribution(cdistr, npoints);
    
    for ss = 1:nstarts
%         Yf(1, ss, dId) = Y0(ss);
%         Yb(1, ss, dId) = Y0(ss);
        Y(1, ss, dId) = Y0(ss);
        for n = 2:npoints
            
            prevy = Y(n-1, ss, dId);                                                % Previous integrated values
            currx = input(n, dId);                                                  % Current random input from distribution
            curry = smrinc_integrator_dynamic(currx, prevy, scoeff, phi, dt); % Current integrated Y 
            Y(n, ss, dId) = curry;
            
%             prevyf = Yf(n-1, ss, dId);  
%             prevyb = Yb(n-1, ss, dId);% Previous integrated values
%             currx = input(n, dId);                                                  % Current random input from distribution
%             [curryf, curryb] = smrinc_integrator_dynamic2(currx, prevyb, prevyf, scoeff, 1, 1, dt); % Current integrated Y 
%             Y(n, ss, dId) = (curryf+curryb)./2;
%             Yf(n, ss, dId) = curryf;
%             Yb(n, ss, dId) = curryb;
        end
    end
    
end

%% Plotting

% State space, force profile, potential profile
state     = 0:0.01:1;
force     = F(state, scoeff);
potential = U(state, scoeff);

% Plots
fig1 = figure;
fig_set_position(fig1, 'All'); 

NumRows = 4;
NumCols = 5;

% Force plot
subplot(NumRows, NumCols, [1 NumCols+1]);
hold on;
plot(state, force);
dforce = diff(force);
plot(support.anchor.x, support.anchor.y, 'or');
%plot(sx(sy~=0), sy(sy ~=0), 'og');
xlim([-0.05 1.05])
hold off;
plot_hline(0, 'k-');
plot_vline(0.5, 'k-');
grid on;
ylabel('F(x)')
title('Force');

% Potential plot
subplot(NumRows, NumCols, 1 + [2*NumCols 3*NumCols]);
plot(state, potential);
xlim([-0.05 1.05])
plot_vline(0.5, 'k-');
grid on;
xlabel('x');
ylabel('U(x) = - \int_{0}^{1} F(x) dx');
title('Potential');

% Dynamic plots

for dId = 1:ndistr
    subplot(NumRows, NumCols, [2 7 12] + dId-1);
    plot(t, Y(:, :, dId));
    xlabel('time [s]');
    ylabel('x');
    xlim([t(1) t(end)]);
    ylim([-0.1 1.1]);
    grid on;
    plot_hline(0.5, 'k-');
    plot_hline(0, 'k-');
    plot_hline(1, 'k-');
    title(ldistr{dId});
    
    subplot(NumRows, NumCols, 17 + dId -1);
    edges = 0:0.05:1;
    cnts  = histc(input(:, dId), 0:0.05:1);
    bar(edges, 100*cnts./sum(cnts), 'histc');
    xlim([0 1]);
    ylim([0 30]);
    grid on;
    xlabel('x');
    ylabel('%');
end

suptitle('Dynamic integration study - random input');


