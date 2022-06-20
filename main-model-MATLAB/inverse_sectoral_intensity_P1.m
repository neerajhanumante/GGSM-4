%  %  %  %  %  %  %  %  Functions - sectoral intensity  %  %  %  %  %  %  %
function [out] = inverse_sectoral_intensity_P1(x)
% this function returns change in state variable corresponding to the water
% demand change
    
    % x may be computed to be negative, hence, the preprocessing
    x = abs(x);

    if x <= 5.03381960534042
        out = 0.119293111*x + 0;
    end
    
    if x > 5.03381960534042
        if x <= 7.8087194753
            out = 0.0214723423*x + 0.4924121022;
        end
    end

    
    if x > 7.8087194753
        if x <= 15.6886359549
            out = 0.0025298856*x + 0.6403448331;
        end
    end
    
    if x > 15.6886359549
        if x <= 44.7839689743
            out = 0.0044683339*x + 0.6098979365;
        end
    end

    if x > 44.7839689743
        out = 0.0222968792*x -0.1885427457;
    end

end