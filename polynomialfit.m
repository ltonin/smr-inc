clearvars; clc; close all;

x = [0   0.1     0.3      0.4   0.5     0.6     0.7     0.9     1];
y = [0  -0.03      0     0.01     0   -0.01       0    0.03     0];

% x = [0.5  0.6  0.7  0.9 1];
% y       = [0  -0.01    0 0.03 0];

p = polyfit(x, y, 10);

x1 = 0:0.01:1;
y1 = polyval(p, x1);

fig1 = figure;
fig_set_position(fig1, 'Top'); 
NumRows = 4;
NumCols = 5;

subplot(NumRows, NumCols, [1 NumCols+1]);
hold on;
plot(x1, y1);
dy = diff(y);
plot(x(y==0), y(y ==0), 'o');
xlim([-0.05 1.05])
hold off;
plot_hline(0, 'k-');
plot_vline(0.5, 'k-');
% plot_vline(x(y==0), 'k--');
grid on;


subplot(NumRows, NumCols, 1 + [2*NumCols 3*NumCols]);
plot(x1, -cumsum(y1));
xlim([-0.05 1.05])
plot_vline(0.5, 'k-');
% plot_vline(x(y==0), 'k--');
grid on;


itypes = {'Normal', 'BetaHigh', 'BetaLow', 'BetaHighLow'};
dt = 0.0625;
it = 1:1:100;
t = 0:dt:dt*(length(it)-1);
V0 = 0.0:0.01:1;
V  = zeros(length(it), length(V0), length(itypes));

c_value = 0.7;


for itype = 1:length(itypes)
    switch(itypes{itype})
        case 'Normal'
            distr = makedist('Beta', 'a', 1, 'b', 1);
            values = random(distr, length(it), 1);
        case 'BetaHigh'
            distr = makedist('Beta', 'a', 5, 'b', 1);
            values = random(distr, length(it), 1);
        case 'BetaLow'
            distr = makedist('Beta', 'a', 1, 'b', 5);
            values = random(distr, length(it), 1);
        case 'BetaHighLow'
            distr = makedist('Beta', 'a', 0.25, 'b', 0.25);
            values = random(distr, length(it), 1);
    end
% values = ones(length(it), 1);
% values(2:2:end) = 0;
% values = rand(length(it));
% beta = makedist('Beta', 'a', 5, 'b', 1);
% values = random(beta, length(it), 1);
    for v = 1:length(V0)
        V(1, v, itype) = V0(v);

        for i = 2:length(it)

            %V(i, v, itype) = V(i-1, v, itype) + dt*5*polyval(p, V(i-1, v, itype)) + dt*0.3*((values(i)-0.5)*exp((values(i)-0.5).^2));
            V(i, v, itype) = smrinc_dynamic_integrator(values(i), V(i-1, v, itype), p, dt);
        end
    end
    
    subplot(NumRows, NumCols, [2 7 12] + itype-1);
    plot(t, V(:, :, itype));
    xlabel('time [s]');
    ylabel('x (mV)');
    xlim([t(1) t(end)]);
    ylim([-0.1 1.1]);
    grid on;
    title(itypes{itype});
    
    subplot(NumRows, NumCols, 17 + itype -1);
    hist(values);
    xlim([0 1]);
end

