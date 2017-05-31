function y = smrinc_integrator_exponential(x, yp, alpha, rejection)

    if(x > 1 - rejection && x < rejection)
        x = yp;
    end
    
    y = alpha*yp + (1-alpha)*x;
    
end