function [mkk,Pkk,vkk,Vkk] = conditionalGIWupdate(W,mkk1,Pkk1,vkk1,Vkk1,models)

if size(W,2)>0

    % Models
    H = models.H;
    Id = models.Id;
    
    % Number of detections
    n = size(W,2);
    % Centroid detection
    zc = sum(W,2)/n;
    % Spread of detections
    Z = cov(W.')*(n-1);
    % Innovation
    innov = zc - kron(H,Id)*mkk1;
    % Innovation matrix
    S = H*Pkk1*H.' + (1/n);
    % Gain
    K = (Pkk1*H.')/S;
    % Innovation spread
    N = S\(innov*innov.');
    
    % Update mean
    mkk = mkk1+kron(K,Id)*innov;
    % Update covariance
    Pkk = Pkk1 - K*S*K.';
    % Update degrees of freedom
    vkk = vkk1+n;
    % Update shape matrix
    Vkk = Vkk1 + N + Z;
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