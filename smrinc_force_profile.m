clearvars; clc; 

% inc_l = 0.6;
% nr   = 0.8;
% 
% eq_x  = [0 1-nr 0.5 nr 1];
% eq_y  = zeros(1, length(eq_x));
% 
% inc_x = (1-inc_l):0.001:inc_l;
% inc_y = (0.5 - inc_x);
% 
% ic_x = [(1-nr)/2  (nr + 1)/2];
% ic_y = [1 1].*[-1 1];
% 
% xovl = [eq_x inc_x ic_x];
% yovl = [eq_y inc_y ic_y];
% [xforce, xId] = sort(xovl);
% yforce = yovl(xId);
% 
% degree = 8;
% coeff = polyfit(xforce, yforce, degree);

inclim = 0.6;
nrpt   = 0.8;
bias   = 0.5;
[coeff, support] = smrinc_get_forceprofile(inclim, nrpt, bias);%, degree);
s = 0:0.001:1;

subplot(2, 1, 1);
hold on;
plot(s, polyval(coeff, s)./max(abs(polyval(coeff, s))));
plot(support.anchor.x, support.anchor.y, 'og');
ylim([-1 1]);
plot_hline(0, 'k');
plot_vline(0.5, 'k');
plot_vline([1-nrpt nrpt], 'k--');
plot_vline([1-inclim inclim], 'k--');
grid on;
hold off;

subplot(2, 1, 2);
hold on;
plot(s, -cumsum(polyval(coeff, s)./max(abs(polyval(coeff, s)))));
plot(support.anchor.x, support.anchor.y, 'og')
plot_hline(0, 'k');
plot_vline(0.5, 'k');
plot_vline([1-nrpt nrpt], 'k--');
plot_vline([1-inclim inclim], 'k--');
grid on;
hold off;

% th = 0.55;
% x1 = (1-th):0.01:th;
% y1 = -0.5*(x1-0.5);
% 
% sx = [0  (1-th)/2 x1 th+(1-th)/2  1];
% sy = [0 -0.05  y1  0.05 0];
% degrees = 8;
% 
% scoeff  = polyfit(sx, sy, degrees);
% F       = @(x, c) polyval(c, x);
% U       = @(x, c) -cumsum(F(x, c));
% Ff      = @(x) polyval(scoeff, x);
% 
% state     = 0:0.01:1;
% force     = F(state, scoeff);
% potential = U(state, scoeff);
% 
% % Force plot
% subplot(2, 1, 1);
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
% subplot(2, 1, 2);
% plot(state, potential);
% xlim([-0.05 1.05])
% plot_vline(0.5, 'k-');
% grid on;
% xlabel('x');
% ylabel('U(x) = - \int_{0}^{1} F(x) dx');
% title('Potential');