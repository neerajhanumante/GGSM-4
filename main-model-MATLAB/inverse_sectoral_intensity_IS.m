%  %  %  %  %  %  %  %  Functions - sectoral intensity  %  %  %  %  %  %  %
function [out] = inverse_sectoral_intensity_IS(x)
% this function returns change in state variable corresponding to the water
% demand change
    
    % x may be computed to be negative, hence, the preprocessing
    x = abs(x);

    upper_limit = 0.1529613116;
    m = 0.0000566222;
    c = 0;
    if x <=  upper_limit
        out = m*x + c;
    end
    
    lower_limit = 0.1529613116;
    upper_limit = 0.3035001313;
    m = 0.0000154701;
    c = 0.0000062947;

    if x > lower_limit
        if x <= upper_limit
            out = m*x + c;
        end
    end

    lower_limit = 0.3035001313;
    upper_limit = 1.2463762082;
    m = 0.0000015899;
    c = 0.0000105075;

    if x > lower_limit
        if x <= upper_limit
            out = m*x + c;
        end
    end
    
    lower_limit = 1.2463762082;
    upper_limit = 5.3236853827;
    m = 0.0000003186;
    c = 0.0000120929;

    if x > lower_limit
        if x <= upper_limit
            out = m*x + c;
        end
    end

    lower_limit = 5.3236853827;
    m = 0.0000089482;
    c = -0.0000338475;

    if x > lower_limit
            out = m*x + c;
    end

end
