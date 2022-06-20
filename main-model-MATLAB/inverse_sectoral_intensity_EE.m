%  %  %  %  %  %  %  %  Functions - sectoral intensity  %  %  %  %  %  %  %
function [out] = inverse_sectoral_intensity_EE(x)
% this function returns change in state variable corresponding to the water
% demand change
    
    % x may be computed to be negative, hence, the preprocessing
    x = abs(x);

    upper_limit = 0.2549338088;
    m = 0.0002632448;
    c = 0;
    if x <=  upper_limit
        out = m*x + c;
    end
    
    lower_limit = 0.2549338088;
    upper_limit = 0.4009877824;
    m = 0.0001362444;
    c = 0.0000323767;

    if x > lower_limit
        if x <= upper_limit
            out = m*x + c;
        end
    end

    lower_limit = 0.4009877824;
    upper_limit = 1.9835006780;
    m = 0.0000107384;
    c = 0.0000827040;

    if x > lower_limit
        if x <= upper_limit
            out = m*x + c;
        end
    end
    
    lower_limit = 1.9835006780;
    upper_limit = 8.8217314200;
    m = 0.0000033644;
    c = 0.0000973266;

    if x > lower_limit
        if x <= upper_limit
            out = m*x + c;
        end
    end

    lower_limit = 8.8217314200;
    m = 0.0001246449;
    c = -0.0009725840;

    if x > lower_limit
            out = m*x + c;
    end

end
