function dx = acc_nomi(t,x)
global u v0 h1 h2;
m = 1650;
g0 = 0.3; g1 = 10; g2 = 0.5;
%g0 = 5; g1 = 200; g2 = 10;
dx = zeros(3,1);
dx(1) = x(2);
dx(2) = (1/m)*(u(1) - g0 - g1*x(2) - g2*x(2).^2) + h1;
dx(3) = v0 - x(2) + h2;