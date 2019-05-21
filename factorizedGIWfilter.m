function [mpred,Ppred,vpred,Vpred,mup,Pup,vup,Vup] = ...
    factorizedGIWfilter(m0,P0,v0,V0,Nt,Z,models)

%%%%%%%%%%%%%%%%%%%%%%
% Initialize filter
%%%%%%%%%%%%%%%%%%%%%%
m = m0;
P = P0;
v = v0;
V = V0;

% Allocate memory
nx = length(m);
d = size(V0,1);
mpred = zeros(nx,Nt);
Ppred = zeros(nx,nx,Nt);
vpred = zeros(1,Nt);
Vpred = zeros(d,d,Nt);
mup = zeros(nx,Nt);
Pup = zeros(nx,nx,Nt);
vup = zeros(1,Nt);
Vup = zeros(d,d,Nt);

for t = 1:Nt
    
    % Predict
    [m,P,v,V] = factorizedGIWprediction(m,P,v,V,models);
    % Store values
    mpred(:,t) = m;
    Ppred(:,:,t) = P;
    vpred(t) = v;
    Vpred(:,:,t) = V;

    % Update
    [m,P,v,V] = factorizedGIWupdate(Z{t},m,P,v,V,models);
    % Store values
    mup(:,t) = m;
    Pup(:,:,t) = P;
    vup(t) = v;
    Vup(:,:,t) = V;
    
end