function [GWD, KIN, EXT] = computeGWDmetric(xtrue,Xtrue,mpred,vpred,Vpred,mup,vup,Vup,msm,vsm,Vsm)

% Number of time steps
Nt = size(mpred,2);
% Dimension of extent
d = size(Vpred,1);
% Allocate memory
GWD = zeros(Nt,3);
KIN = zeros(Nt,3);
EXT = zeros(Nt,3);

for t = 1:Nt

    % Position
    xpred = mpred(:,t);
    xup = mup(:,t);
    xsm = msm(:,t);

    % Extent
    Xpred = Vpred(:,:,t)/(vpred(t)-2*d-2);
    Xup = Vup(:,:,t)/(vup(t)-2*d-2);
    Xsm = Vsm(:,:,t)/(vsm(t)-2*d-2);

    % Gaussian Wasserstein Distance
    [GWD(t,1), KIN(t,1), EXT(t,1)] = GaussianWassersteinDistanceMetric( xpred(1:2), Xpred, xtrue(1:2,t), Xtrue(:,:,t) );
    [GWD(t,2), KIN(t,2), EXT(t,2)] = GaussianWassersteinDistanceMetric( xup(1:2),   Xup,   xtrue(1:2,t), Xtrue(:,:,t) );
    [GWD(t,3), KIN(t,3), EXT(t,3)] = GaussianWassersteinDistanceMetric( xsm(1:2),   Xsm,   xtrue(1:2,t), Xtrue(:,:,t) );
    
end


function [gwd,kin_vec_norm,extent_norm] = GaussianWassersteinDistanceMetric(m1,P1,m2,P2)

% Computes the Wasserstein Distance between to Gaussian distributions with
% mean m1 and m2 and covariance P1 and P2.

% Kinematic vector norm
kin_vec_norm = norm(m1-m2,2);
% Extent norm
sqrtP1 = sqrtm(P1);
extent_norm = sqrt(trace(P1) + trace(P2) - 2*trace(sqrtm(sqrtP1*P2*sqrtP1)));
% Gaussian Wasserstein distance
gwd = sqrt(kin_vec_norm^2 + extent_norm^2);