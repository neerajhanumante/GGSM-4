%  %  %  %  %  %  %  %  Function - agri inelasticity update %  %  %  %  %  %  %
function [array_agri_price_change, timestep_agri_price_limit_breach] = f_correction_agri_water_price_change_inelasticity( ...
    i, ...
    mat_percent_p_water_change, ...
    timestep_agri_price_limit_breach, ...
    inelastic_price_change_limit, ...
    mat_price_water)

% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
%   timestep_agri_price_limit_breach is a vector of the timestamps at
%   which regional water price change breaches the inelastic price change
%   limit. Updating elements of this vector need two conditions to be
%   satisfied: First, the regional water price change greater than the
%   inelastic price change limit, and second, presently the element value
%   is zero.
% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %

% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
% determine the index of the timestep_agri_price_limit_breach vector
% to be updated.
% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
    
    array_agri_price_change = mat_percent_p_water_change(:, 1, i);
    condition_limit = array_agri_price_change > inelastic_price_change_limit;
    condition_zero = timestep_agri_price_limit_breach == 0;
    index_to_be_updated = find(condition_limit & condition_zero);

% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
% updating timestep_agri_price_limit_breach vector.
% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
    timestep_agri_price_limit_breach(index_to_be_updated) = i;
    
% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
% Update the water price change matrix
% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %

% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
% Inelastic behaviour, set all values within the threshold to be zero. 
% This has to be carried out for all the ZERO index values.
% Necessary checks are already done while identifying the index in the 
% vector to be updated.
% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
    condition_limit = array_agri_price_change < inelastic_price_change_limit;
    array_agri_price_change(condition_limit) = 0;
% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
% Using the timesteps in the timestep_agri_price_limit_breach vector
% as base values compute the price change values for the agricultural
% sector. 
% This has to be carried out for all the NON-ZERO index values.
% Necessary checks are already done while identifying the index in the 
% vector to be updated.
% - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  - -  %
    condition_limit = array_agri_price_change > inelastic_price_change_limit;
    index_to_be_updated = find(condition_limit);
    if ~isempty(index_to_be_updated)    
        for counter_index = 1:length(index_to_be_updated)
            local_index = index_to_be_updated(counter_index);
            skip_initial_timesteps_agri = timestep_agri_price_limit_breach(local_index);
            numerator = mat_price_water(local_index, 1, i) - mat_price_water(local_index, 1, skip_initial_timesteps_agri);
            denomenator = mat_price_water(local_index, 1, skip_initial_timesteps_agri);
            array_agri_price_change(local_index) = numerator./denomenator * 100;
        end
    end
   
    
end                                       