function [out] = f_energy_for_water(a)
% this function computes the energy requirement for water processing for
% each sector and stage. It returns total energy requirement for the water
% processing in kWh

% 50th percentile
% kWh/m3                P1+H1    IS+EE           HH
% Sourcing/conveyance 0.0967    0.0967         0.0967
% Treatment              -      0.1658          0.225
% Distribution           -       0              0.247
% Total               0.0967    0.2625          0.5687

% 75th percentile
% kWh/m3                P1+H1    IS+EE           HH
% Sourcing/conveyance 0.1333    0.1333         0.1333
% Treatment              -      0.3067          0.4
% Distribution           -       -              0.339
% Total               0.1333    0.44          0.8723

% Compute the total global water demands for each sector
d_W_P1 = a(1);
d_W_H1 = a(2);
d_W_IS = a(3);
d_W_EE = a(4);
d_W_HH = a(5);
d_W_ISEE = d_W_EE + d_W_IS;
d_W_P1H1 = d_W_P1 + d_W_H1;
% Compute the total Energy requirements for each sector
% Converting the demand from billion cu m to cu m


% 50th percentile
d_E_W_P1H1 = 0.0967 * d_W_P1H1 * 1e9;
d_E_W_ISEE = 0.2625 * d_W_ISEE * 1e9;
d_E_W_HH = 0.5687 * d_W_HH * 1e9;

% 75th percentile
d_E_W_P1H1 = 0.1333 * d_W_P1H1 * 1e9;
d_E_W_ISEE = 0.44 * d_W_ISEE * 1e9;
d_E_W_HH = 0.8723 * d_W_HH * 1e9;
% disp('75 percentile used.')
out = [d_E_W_P1H1; 
    d_E_W_ISEE;
    d_E_W_HH
    ];
end
