function [Z] = generateMeasurements(Nt,p_D,xtrue,Xtrue,Nz,R,d)

% Allocate memory
Z = cell(1,Nt);
% Iteratve over time
for t = 1:Nt
    % Sample detection process
    if rand<p_D
        % Random measurements
        Z{t} = mvnrnd(xtrue(1:2,t)',Xtrue(:,:,t),Nz)' + mvnrnd([0 0],R,Nz)';
    else
        % No detection
        Z{t} = zeros(d,0);
    end
end