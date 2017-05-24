function y = smrinc_integrator_exponential(x, yp, alpha, rejection)

    if(max(x) >= 1 - rejection && max(x) <= rejection)
        x = yp;
    end
    
    y = alpha*yp + (1-alpha)*x;
    
end