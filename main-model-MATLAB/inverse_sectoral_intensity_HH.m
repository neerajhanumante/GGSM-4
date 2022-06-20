%  %  %  %  %  %  %  %  Functions - sectoral intensity  %  %  %  %  %  %  %
function [out] = inverse_sectoral_intensity_HH(x)
% this function returns change in state variable corresponding to the water
% demand change
    
    % x may be computed to be negative, hence, the preprocessing
    x = abs(x);

    if x <= 0.44944328155009916
        out = 77380.9584150069*x + 1273.37282012312;
    end
    
    if x > 0.44944328155009916
        if x <= 1.4533993599999862
            out = 8125.89575359953*x + 32399.5954466238;
        end
    end

    
    if x > 1.4533993599999862
        if x <= 8.055252854372617
            out = 2552.82178433151*x + 40499.4975867905;
        end
    end
    if x > 8.055252854372617
        out = 23694.7721372485*x - 129804.258340548;
    end

end
