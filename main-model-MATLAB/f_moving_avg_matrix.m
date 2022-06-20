%  %  %  %  %  %  %  %  Functions %  %  %  %  %  %  %
function [out] = f_moving_avg_matrix(a, kb, timestep)
% this function returns the moving average of 3D matrix
% third dimension is the time dimension
% output gives moving average of each element of matrix a for the backward
% time kb
% local_matrix is a subset of a with last kb timesteps
local_matrix = a(:, :, timestep-kb:timestep);

for i = 1:size(local_matrix, 1)
    for j = 1:size(local_matrix, 2)
        a(i, j, timestep) = mean(local_matrix(i, j, :));
    end
end
out = a;
end