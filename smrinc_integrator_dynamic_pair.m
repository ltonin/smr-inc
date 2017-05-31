function y = smrinc_integrator_dynamic_pair(x, yp, p, phi, beta, dt)
    
    cx = 0:0.01:1;
    mf = max(abs(polyval(p, cx)));

%     plot(cx, polyval(p, cx)./mf);
%     keyboard


    Ff = phi*polyval(p, yp)./mf;
%     Fb = beta*(x-0.5)*exp(-(x-0.5).^2);
%     Fb = beta*(0.5-x)*exp(-(0.5-x).^2);
%      Fb = beta*(0.5-x);
    Fb = beta*polyval(p, x)./mf;

    y = yp + dt*(Ff + Fb);

end