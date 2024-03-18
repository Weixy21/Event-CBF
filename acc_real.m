function dx = acc_real(t,x)
global u v0 sigma1 sigma2 sigma3;
m = 1650;
f0 = 0.1; f1 = 5; f2 = 0.25;
dx = zeros(3,1);
dx(1) = x(2);
dx(2) = (1/m)*(sigma3*u(1) - f0 - f1*x(2) - f2*x(2).^2) + sigma1;
dx(3) = v0 - x(2) + sigma2;