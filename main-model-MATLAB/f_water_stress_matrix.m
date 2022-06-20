%  %  %  %  %  %  %  %  Functions - sectoral intensity  %  %  %  %  %  %  %
function [local_mat_water_stress, local_file_mat_water_stress] = f_water_stress_matrix(water_stress_Africa_total, ...
    water_stress_Asia_total, ...
    water_stress_Europe_total, ...
    water_stress_North_America_total, ...
    water_stress_Oceania_total, ...
    water_stress_South_America_total)

    % water stress experienced by all sectors in a regions is same
    delta_P1_Africa = water_stress_Africa_total;
    delta_P1_Asia = water_stress_Asia_total;
    delta_P1_Europe = water_stress_Europe_total;
    delta_P1_North_America = water_stress_North_America_total;
    delta_P1_Oceania = water_stress_Oceania_total;
    delta_P1_South_America = water_stress_South_America_total;
    
    delta_H1_Africa = water_stress_Africa_total;
    delta_H1_Asia = water_stress_Asia_total;
    delta_H1_Europe = water_stress_Europe_total;
    delta_H1_North_America = water_stress_North_America_total;
    delta_H1_Oceania = water_stress_Oceania_total;
    delta_H1_South_America = water_stress_South_America_total;
    
    delta_HH_Africa = water_stress_Africa_total;
    delta_HH_Asia = water_stress_Asia_total;
    delta_HH_Europe = water_stress_Europe_total;
    delta_HH_North_America = water_stress_North_America_total;
    delta_HH_Oceania = water_stress_Oceania_total;
    delta_HH_South_America = water_stress_South_America_total;
    
    delta_IS_Africa = water_stress_Africa_total;
    delta_IS_Asia = water_stress_Asia_total;
    delta_IS_Europe = water_stress_Europe_total;
    delta_IS_North_America = water_stress_North_America_total;
    delta_IS_Oceania = water_stress_Oceania_total;
    delta_IS_South_America = water_stress_South_America_total;
    
    delta_EE_Africa = water_stress_Africa_total;
    delta_EE_Asia = water_stress_Asia_total;
    delta_EE_Europe = water_stress_Europe_total;
    delta_EE_North_America = water_stress_North_America_total;
    delta_EE_Oceania = water_stress_Oceania_total;
    delta_EE_South_America = water_stress_South_America_total;
    
    local_mat_water_stress = [
    % Column headers - P1, H1, HH, IS, EE
                       delta_P1_Africa, delta_H1_Africa, delta_HH_Africa, delta_IS_Africa, delta_EE_Africa;      % Africa
                       delta_P1_Asia, delta_H1_Asia, delta_HH_Asia, delta_IS_Asia, delta_EE_Asia;      % Asia
                       delta_P1_Europe, delta_H1_Europe, delta_HH_Europe, delta_IS_Europe, delta_EE_Europe;      % Europe
                       delta_P1_North_America, delta_H1_North_America, delta_HH_North_America, delta_IS_North_America, delta_EE_North_America;      % North America
                       delta_P1_Oceania, delta_H1_Oceania, delta_HH_Oceania, delta_IS_Oceania, delta_EE_Oceania;      % Oceania
                       delta_P1_South_America, delta_H1_South_America, delta_HH_South_America, delta_IS_South_America, delta_EE_South_America;      % South America
                     ];
    
    local_file_mat_water_stress = [water_stress_Africa_total, water_stress_Asia_total, water_stress_Europe_total, water_stress_North_America_total, water_stress_Oceania_total, water_stress_South_America_total];

end