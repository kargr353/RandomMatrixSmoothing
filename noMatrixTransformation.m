function [M,M1,M2,I_x] = noMatrixTransformation(m)

% Assumes four states: 2D position, 2D velocity
% Assumes 2D extent

% States that have any effect on the function M
I_x = [];

% No matrix transformation
M = eye(2);

% First derivative
M1 = zeros(2);

% Second derivative
M2 = zeros(2);