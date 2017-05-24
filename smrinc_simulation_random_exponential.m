clearvars; clc;

%% Exponential defintion
alpha = 0.97;
rejection = 0.55;

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

for dId = 1:ndistr
    
    % Getting the current input distribution
    cdistr = ldistr{dId};
    input(:, dId) = smrinc_get_distribution(cdistr, npoints);
    
    for ss = 1:nstarts
        Y(1, ss, dId) = Y0(ss);

        for n = 2:npoints
            
            prevy = Y(n-1, ss, dId);                                                % Previous integrated values
            currx = input(n, dId);                                                  % Current random input from distribution
            curry = smrinc_integrator_exponential(currx, prevy, alpha, rejection);  % Current integrated Y 
            Y(n, ss, dId) = curry;
        end
    end
    
end

%% Plotting

% Plots
fig1 = figure;
fig_set_position(fig1, 'All'); 

NumRows = 4;
NumCols = 4;

% Exponential plots

for dId = 1:ndistr
    subplot(NumRows, NumCols, [1 5 9] + dId-1);
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
    
    subplot(NumRows, NumCols, 13 + dId -1);
    edges = 0:0.05:1;
    cnts  = histc(input(:, dId), 0:0.05:1);
    bar(edges, 100*cnts./sum(cnts), 'histc');
    xlim([0 1]);
    ylim([0 30]);
    grid on;
    xlabel('x');
    ylabel('%');
end

suptitle('Exponential integration study - random input');


