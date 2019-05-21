function [mkk,Pkk,vkk,Vkk] = factorizedGIWupdate(W,mkk1,Pkk1,vkk1,Vkk1,models)

if size(W,2)>0

    % Models
    H = models.H;
    d = models.d;
    R = models.R;
    scaleFactor = models.scaleFactor;
    
    % Number of detections
    n = size(W,2);
    % Centroid detection
    zc = sum(W,2)/n;
    % Spread of detections
    Z = cov(W.')*(n-1);
    % Innovation
    innov = zc - H*mkk1;
    % Predicted expected extent matrix
    Xhat = Vkk1/(vkk1-2*d-2);
    sqrtXhat = sqrtm(Xhat);
    % Predicted measurement noise covariance
    Y = scaleFactor*Xhat + R;
    sqrtY = sqrtm(Y);
    % Innovation matrix
    S = H*Pkk1*H.' + (Y/n);
    sqrtS = sqrtm(S);
    % Gain
    K = (Pkk1*H.')/S;
    % Innovation spread
    Nhat = (sqrtXhat/sqrtS)*(innov*innov.')*(sqrtS\sqrtXhat);
    % Modified detection spread
    Zhat = (sqrtXhat/sqrtY)*Z*(sqrtY\sqrtXhat);
    
    % Update mean
    mkk = mkk1+K*innov;
    % Update covariance
    Pkk = Pkk1 - K*S*K.';
    % Update degrees of freedom
    vkk = vkk1+n;
    % Update shape matrix
    Vkk = Vkk1 + Nhat + Zhat;
else
    % No detections, posterior parameters are equal to predicted parameters
    mkk = mkk1;
    Pkk = Pkk1;
    vkk = vkk1;
    Vkk = Vkk1;
end

% Make sure P and V are symmetric
Pkk = 0.5*(Pkk + Pkk.');
Vkk = 0.5*(Vkk + Vkk.');