clearvars; clc; close all;


C  = 10;    % [uF]

Gl = 19;    % [mS]
El = -67;   % [mV]

Gn = 74;    % [mS]
En = 60;    % [mV]
V12 = 1.5;  % [mV]
k   = 16;   % [mV]

leak = @(x, Gl, El) Gl*(x - El)./C;
m    = @(x, V12, k) 1./(1+exp((V12-x)./k));
fast = @(x, Gn, En, V12, k) Gn*(1./(1+exp((V12-x)./k))).*(x-En)./C; 
F    = @(leak, fast) -leak - fast;

dt = 0.0625;
it = 1:1:100;
t = 0:dt:dt*(length(it)-1);
V0 = -60:5:40;
V = zeros(length(it), length(V0));


for v = 1:length(V0)
    V(1, v) = V0(v);

    for i = 2:length(it)

        V(i, v) = V(i-1, v) + dt*F(leak(V(i-1, v), Gl, El), fast(V(i-1, v), Gn, En, V12, k));

    end
end

figure;
plot(t, V);
xlabel('time [s]');
ylabel('x (mV)');
xlim([t(1) t(end)]);


figure;
subplot(2, 2, 1);
xl = -100:0.1:0;
plot(xl, -leak(xl, Gl, El));
xlabel('x (mV)');
ylabel('F(x) = dx/dt');
title('Leak');
xlim([xl(1) xl(end)]);
plot_hline(0, 'k-');

subplot(2, 2, 2);
xf = -60:0.1:60;
plot(xf, F(leak(xf, Gl, El), fast(xf, Gn, En, V12, k)));
xlabel('x (mV)');
ylabel('F(x) = dx/dt');
title('Leak+Fast');
plot_hline(0, 'k-');
xlim([xf(1) xf(end)]);

subplot(2, 2, 3);
xl = -100:0.1:0;
plot(xl, -cumsum(-leak(xl, Gl, El)));
xlabel('x (mV)');
ylabel('U(x) = -|F(x)');
plot_hline(0, 'k-');
xlim([xl(1) xl(end)]);
title('Potential Leak');

subplot(2, 2, 4);
xf = -60:0.1:60;
plot(xf, -cumsum(F(leak(xf, Gl, El), fast(xf, Gn, En, V12, k))));
xlabel('x (mV)');
ylabel('U(x) = -|F(x)');
plot_hline(0, 'k-');
title('Potential Leak+Fast');
xlim([xf(1) xf(end)]);


%plot(xf, -cumsum(-leak(xf, Gl, El)-fast(xf, Gn, En,V12, k)));

