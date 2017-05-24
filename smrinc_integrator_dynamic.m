function y = smrinc_integrator_dynamic(x, yp, p, dt)
    
    lambdaf = 5;
    lambdab = 0.3;
    
    y = yp + dt*lambdaf*polyval(p, yp) + dt*lambdab*(x-0.5)*exp((x-0.5).^2);


end