function [xtrue,Xtrue] = ...
    generateTrueTrack(Nt,x0,P0,Xbase,Q,T,maxTR)

% Allocate memory
xtrue = zeros(5,Nt);
Xtrue = zeros(2,2,Nt);
% Initial state
xtrue(:,1) = [mvnrnd([0 0],P0(1:2,1:2)).'; x0];
Xtrue(:,:,1) = Xbase;
for t=1:Nt-1
    % Process noise sample
    w = mvnrnd([0 0 0 0 0],Q)';
    % Simulate motion model, add noise
    xtrue(:,t+1) = LinearVelocityConstantTurnrate(xtrue(:,t),T) + w;
    % Heading
    heading = atan2(xtrue(4,t+1),xtrue(3,t+1));
    % Rotation matrix
    rm = [cos(heading) -sin(heading); sin(heading) cos(heading)];
    % Rotate extent
    rxr = rm*Xbase*rm';
    % Save
    Xtrue(:,:,t+1) = 0.5*(rxr+rxr.');
    
    % The following code can be used to make sure turn-rate is not 
    % unreasonably large
    if xtrue(end,t+1)>maxTR
        xtrue(end,t+1) = 2*maxTR-xtrue(end,t+1);
    elseif xtrue(end,t+1)<-maxTR
        xtrue(end,t+1) = -2*maxTR-xtrue(end,t+1);
    end
end