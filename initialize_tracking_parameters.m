function [models,models_conditional,Xbase] = initialize_tracking_parameters(T,siga,sigw,target_width,target_length)

% Process noise
G = [0.5*T^2*eye(2) zeros(2,1); T*eye(2) zeros(2,1); zeros(1,2) T];
Qa = blkdiag(eye(2)*0.5*siga^2,sigw^2);
Q = G*Qa*G';
% Process noise parameters for extent
tau = 10;
n_IW = inf;

siga_2 = siga;
tau_2 = 0;
n_IW_2 = 100;

% Measurement noise
sig_r = 0.0001;
R = diag([sig_r sig_r].^2);

% Target extent, aligned to x-axis. Extent is assumed to be aligned with
% the heading of the simulated target. The extent corresponds to Ns
% standard deviations
Ns = 2;
Xbase = diag((1*[target_length/2/Ns target_width/2/Ns]).^2);

% Factorized models
% Process model
models.predictionType = '';
models.motionModel = @LinearVelocityConstantTurnrate;
models.matrixTransformationFunction = @RotationMatrixLVCT;
models.inverseMatrixTransformationFunction = @inverseRotationMatrixLVCT;
models.Q = Q;
models.Ts = T;
models.tau = tau;
models.d = 2;
models.n = n_IW;
models.s_last = 100;
models.n_x = 5;
% Measurement model
models.updateType = '';
models.H = [eye(models.d) zeros(models.d,models.n_x-models.d)];
models.Id = eye(models.d);
models.R = R;
models.scaleFactor = 1;

% Conditional models
% Process model
models_conditional.predictionType = '';
models_conditional.d = 2;
models_conditional.F = [1 T; 0 1];
models_conditional.D = siga_2^2*(1-exp(-2*T/tau_2))*[0 0; 0 1];
models_conditional.Ts = T;
models_conditional.tau = tau;
models_conditional.n = n_IW_2;
models_conditional.A = eye(models_conditional.d);
% Measurement model
models_conditional.updateType = '';
models_conditional.H = [1 0];
models_conditional.Id = eye(models_conditional.d);
models_conditional.B = [];
models_conditional.scaleFactor = 1;
models_conditional.R = R;