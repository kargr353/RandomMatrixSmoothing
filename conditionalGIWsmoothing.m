function [mkT,PkT,vkT,VkT] = conditionalGIWsmoothing(mkk,Pkk,vkk,Vkk,mk1T,Pk1T,vk1T,Vk1T,vk1k,Vk1k,models)

% Models
F = models.F;
D = models.D;
Id = models.Id;
d = models.d;
n = models.n;

%% Gaussian

% Predicted mean
mk1k = kron(F,Id)*mkk;
% Predicted covariance
Pk1k = F*Pkk*F.' + D;
% Smoothing gain
G = Pkk*F.'/Pk1k;

% Smoothed mean
mkT = mkk + kron(G,Id)*(mk1T - mk1k);
% Smoothed covariance
PkT = Pkk - G*(Pk1k-Pk1T)*G.';

%% Inverse Wishart

% Division of IW k+1|T and IW k+1|k is porportional to IW(e,E)
e = vk1T-vk1k;
E = Vk1T-Vk1k;

% Make sure that n is large enough for the given e
if n<max([2*(d+1)^2/e+0.1 , 3*(d+1)-e+0.1])
    n = max([2*(d+1)^2/e+0.1 , 3*(d+1)-e+0.1]);
end

if n*e > (2*(d+1)^2) && n+e>3*(d+1)
    % Density approximations achieved by matching E[X] and E[X^-1] instead
    % of minimising the KL-div
    g = (e- (2*(d+1)^2)/n )/(1 +e/n - 3*(d+1)/n);
    G = E/(1+(e-3*(d+1))/n);
    
else
    % To avoid numerical issues
    g = 0;
    G = 0;   
end

% IW(Xk ; vkk, Vkk) IW(Xk ; g, G) is proportional to IW(Xk, vkk+g , Vkk+G)
vkT = vkk + g;
VkT = Vkk + G;

% Make sure P and V are symmetric
PkT = 0.5*(PkT + PkT.');
VkT = 0.5*(VkT + VkT.');

end



