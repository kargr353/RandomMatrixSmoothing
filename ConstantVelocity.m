function [mp,Fm] = ConstantVelocity(m,Ts)

% Constant velocity model
Fm = [
    1 0 Ts 0;
    0 1 0  Ts;
    0 0 1  0;
    0 0 0  1
    ];
% Prediction
mp = Fm*m;