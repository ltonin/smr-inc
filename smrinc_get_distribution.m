function [y, d] = smrinc_get_distribution(type, npoints)
    
    switch(type)
        case 'Normal'
            d = makedist('Beta', 'a', 1, 'b', 1);
            y = random(d, npoints, 1);
        case 'BetaHigh'
            d = makedist('Beta', 'a', 5, 'b', 1);
            y = random(d, npoints, 1);
        case 'BetaLow'
            d = makedist('Beta', 'a', 1, 'b', 5);
            y = random(d, npoints, 1);
        case 'BetaHighLow'
            d = makedist('Beta', 'a', 0.25, 'b', 0.25);
            y = random(d, npoints, 1);
    end

end