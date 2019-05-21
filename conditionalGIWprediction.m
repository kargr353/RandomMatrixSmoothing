function [mk1k,Pk1k,vk1k,Vk1k] = conditionalGIWprediction(mkk,Pkk,vkk,Vkk,models)

% Models
F = models.F;
D = models.D;
Id = models.Id;
d = models.d;
n = models.n;

% Predicted mean
mk1k = kron(F,Id)*mkk;
% Predicted covariance
Pk1k = F*Pkk*F.'+D;

% Predicted degrees of freedom
vk1k = (vkk*(1+(d+1)/n)-(2*(d+1)^2)/n)/(1+vkk/n-2*(d+1)/n);
% Predicted shape matrix
Vk1k = Vkk*(1-(d+1)/n)/(1-(d+1)/n+(vkk-(d+1))/n);

% Make sure P and V are symmetric
Pk1k = 0.5*(Pk1k + Pk1k.');
Vk1k = 0.5*(Vk1k + Vk1k.');
