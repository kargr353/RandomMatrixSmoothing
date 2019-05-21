function plotCovariance(m,P,col,nsig)

sqrtP = nsig*sqrtm(P);
phi = 0:pi/20:2*pi;
xy = sqrtP*[cos(phi); sin(phi)];
xy = xy + repmat(m,[1 length(phi)]);
plot(xy(1,:),xy(2,:),'-','color',col,'linewidth',3)