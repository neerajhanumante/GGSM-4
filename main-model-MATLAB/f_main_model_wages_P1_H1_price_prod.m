
%  %  %  %  %  %  %  %  Functions - per capita demand, main model %  %  %  %  %  %  %
function [out] = f_main_model_wages_P1_H1_price_prod(i, aw, cw, ISbar,...
    ISmassdeficit, ISmass, theta, lambda, dw, numHH, MinConstant,...
    aP1, bP1, cP1, P1massdeficit, P1, P1bar,...
    aP1p, bP1p, cP1p,...
    aH1, bH1, cH1, H1massdeficit, H1, H1bar,...
    aH1p, bH1p, cH1p)

% this function returns the per capita demand for different sectors
% price setting mechanism is used for economic modelling


    W=max(aw+cw*(ISbar-(ISmassdeficit(i)+ISmass(i)))/(theta(i)+lambda)-dw*numHH(i),MinConstant);
    
    %
    %	Based on the Wage, industries set prices and their production (how
    %	much they would like to produce to maximize their profits based on
    %	their assumption as to what the demand for their products will be).
%	Here, a linear functional form is assumed for the supply.
    %
    if (P1(i)==0)
        pP1=0;
        P1production=0;
    else
        pP1=max(aP1+bP1*W-cP1*((P1massdeficit+P1(i))-P1bar),MinConstant);
        P1production=max(aP1p-bP1p*W-cP1p*((P1massdeficit+P1(i))-P1bar),MinConstant);
    end



    if (H1(i)==0)
        pH1=0;
        H1production=0;
    else
        pH1=max(aH1+bH1*W-cH1*((H1massdeficit(i)+H1(i))-H1bar),MinConstant);
        H1production=max(aH1p-bH1p*W-cH1p*((H1massdeficit(i)+H1(i))-H1bar),MinConstant);
    end
    out = [W; pP1; P1production; pH1; H1production];

end
