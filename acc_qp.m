clc;
x = [0 20 100]; y = x; te = 600;
g = 9.81; m = 1650;
g0 = 0.3; g1 = 10; g2 = 0.5; %nominal dynamics
%g0 = 5; g1 = 200; g2 = 10;
f0 = 0.1; f1 = 5; f2 = 0.25; % real dynamics
vd = 24; eps = 10; psc = 1;
ca = 2; cd = 2;
result1 = zeros(te,5); time = [0,0];
c = 10;
set(0,'DefaultTextInterpreter','latex')
warning('off');
global u v0;
global sigma1 sigma2 sigma3 h1 h2
v0 = 13.89;


t_next = 0; dt = 0.05; h1 = 0; h2 = 0; u = [0;0];
s1 = 0.4; s2 = 0.5; w1 = 0.2; w2 = 1; dw21 = 0.5; dw22 = 0.2; %0.5, 0.2
%% main loop, event-driven
for i = 0:te
    i
% generate noise
sigma1 = 0.04*(0.5 - rand()); %0.4
sigma2 = 0.4*(0.5 - rand()); %4
sigma3 = 1 + 0.02*(0.5 - rand()); %0.2


if(dt*i == t_next)
    %% update nominal dynamics
    y = x;
    Fn = g0 + g1*y(2) + g2*y(2)^2;
    Fr = f0 + f1*x(2) + f2*x(2)^2;
    e2_dot = sigma2 + v0 - x(2) - (h2 + v0 - y(2));
    h2 = h2 + e2_dot;
    e2_ddot = -sigma1 + Fr/m - sigma3*u(1)/m - (-h1 + Fn/m - u(1)/m);
    h1 = h1 - e2_ddot;
    time = [time; [t_next, t_next - time(end,1)]];
    %% find the minimum terms in robust CBFs
%     b_dot = (h2 + v0 - v + e2_dot) + z + e2 - c > 0
%     z + e2 - c > 0
%      v - s1 < v < v + s1
%      z - s2 < z < z + s2
%      -w2 < e2 < w2
%      -dw21 < e2_dot < dw21
%      -dw22 < e2_ddot < dw22
     
     % v, z, e2, e2_dot, e2_ddot
     v_tk = y(2);
     z_tk = y(3);
     C2_a = [1, -1, -1, -1, 0];
     C2_b = h2 + v0 - c;
     C1_a = [0, -1, -1, 0, 0];
     C1_b = -c;
     v_a = [1, 0, 0, 0, 0; -1, 0, 0, 0, 0];
     v_b = [v_tk + s1; s1 - v_tk];
     z_a = [0, 1, 0, 0, 0; 0, -1, 0, 0, 0];
     z_b = [z_tk + s2; s2 - z_tk];
     e2_a = [0, 0, 1, 0, 0; 0, 0, -1, 0, 0];
     e2_b = [w2; w2];
     e2_dot_a = [0, 0, 0, 1, 0; 0, 0, 0, -1, 0];
     e2_dot_b = [dw21; dw21];
     e2_ddot_a = [0, 0, 0, 0, 1; 0, 0, 0, 0, -1];
     e2_ddot_b = [dw22; dw22];
     A = [C2_a; C1_a; v_a; z_a; e2_a; e2_dot_a; e2_ddot_a];
     b = [C2_b; C1_b; v_b; z_b; e2_b; e2_dot_b; e2_ddot_b];
     
     H = [g2/2, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0];
     F = [g1; 0; 0; 0; 0];
     
     options = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','off');
    [~,Fn_min,~,~,~] = ...
       quadprog(H,F,A,b,[],[],[],[],[],options); 
    Fn_min = g0 + Fn_min;
    
    f = [0, 0, 0, 0, 1];
    options = optimoptions('linprog','Algorithm','active-set','Display','off');
    [~,e2_ddot_min,~,~,~] = linprog(f,A,b,[],[],[],[],[],options);
    f = [-1, 0, 0, 0, 0];
    [~,v_min,~,~,~] = linprog(f,A,b,[],[],[],[],[],options);
    f = [0, 0, 0, 1, 0];
    [~,e2_dot_min,~,~,~] = linprog(f,A,b,[],[],[],[],[],options);
    f = [0, 1, 0, 0, 0];
    [~,z_min,~,~,~] = linprog(f,A,b,[],[],[],[],[],options);
    f = [0, 0, 1, 0, 0];
    [~,e2_min,~,~,~] = linprog(f,A,b,[],[],[],[],[],options);
    
    Lf_terms = -h1 + Fn_min/m + e2_ddot_min + 2*(h2 + v0 + v_min + e2_dot_min) + z_min + e2_min - c;
    Lg_term = 1/m;
    %% CLF for desired speed
    phi0 = -2*(x(2) - vd)*Fn/m + eps*(x(2) - vd)^2;
    phi1 = 2*(x(2) - vd)/m;
    %% QP for OC
    A = [phi1 -1;Lg_term 0;1 0;-1 0];
    b = [-phi0;Lf_terms;ca*m*g;cd*m*g];
    H = [2/(m^2) 0; 0 2*psc];
    F = [-2*Fr/(m^2); 0];
    options = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','off');
    [u,fval,exitflag,output,lambda] = ...
       quadprog(H,F,A,b,[],[],[],[],[],options);
end

%% update dynamics
C1 = x(3) - c;
C2 = C1 + sigma2 + v0 - x(2);

t = [0 dt];
[~,xx]=ode45('acc_real',t,x);
x = xx(end, 1:3);
result1(i+1,:) = [dt*i xx(end,2) u(1)/m C1 C2];

[~,xx]=ode45('acc_nomi',t,y);
y = xx(end, 1:3);

Fn = g0 + g1*y(2) + g2*y(2)^2;
Fr = f0 + f1*x(2) + f2*x(2)^2;

%% evaluate the next event time
%      v - s1 < v < v + s1
%      z - s2 < z < z + s2
%      -w2 < e2 < w2
%      -dw21 < e2_dot < dw21
%      -dw22 < e2_ddot < dw22
e2 = x(3) - y(3);
e2_dot = sigma2 + v0 - x(2) - (h2 + v0 - y(2));
e2_ddot = -sigma1 + Fr/m - sigma3*u(1)/m - (-h1 + Fn/m - u(1)/m);
if(y(2) > v_tk + s1 || y(2) < v_tk - s1 || y(3) > z_tk + s2 || y(3) < z_tk - s2 || abs(e2) > w2 || abs(e2_dot) > dw21 || abs(e2_ddot) > dw22)
    t_next = dt*i + dt;
end

end
% 
figure(1)
subplot(4, 1, 1)
plot(result1(:, 1),result1(:, 2),'b-',[0,30],[24,24], 'k--',[0,30],[13.89,13.89], 'k--')
grid on
subplot(4, 1, 2)
plot(result1(:, 1),result1(:, 3)/g, 'b-',[0,30],[ca,ca], 'k--',[0,30],[-cd,-cd], 'k--')
axis([0 30 -0.75 0.75]); 
grid on
subplot(4, 1, 3)
plot(result1(:, 1),result1(:, 4), 'b-',[0,30],[0,0], 'k--')
axis([0 30 -5 75]); 
grid on
subplot(4, 1, 4)
plot(result1(:, 1),result1(:, 5), 'b-')
grid on

figure(2)
plot(time(:, 1),time(:, 2), 'r-')
grid on
%% time-driven CBFs
x = [0 20 100]; y = x; dt = 0.1; te = 300; h1 = 0; h2 = 0;
result2 = zeros(te,5);
for i = 0:te
    i
sigma1 = 0.4*(0.5 - rand());
sigma2 = 4*(0.5 - rand());
sigma3 = 1 + 0.2*(0.5 - rand());
   y = x;
    
    Fn = g0 + g1*y(2) + g2*y(2)^2;
    Fr = f0 + f1*x(2) + f2*x(2)^2; 

    b = y(3) - c;
    C1 = x(3) - c;
    b_dot = v0 - y(2);
    LfB_e = 2*b_dot + b;
    C2 = C1 + sigma2 + v0 - x(2);
    LfB = LfB_e + Fn/m; % 2nd degree CBF
    LgB = 1/m;
    %% CLF
    phi0 = -2*(x(2) - vd)*Fn/m + eps*(x(2) - vd)^2;
    phi1 = 2*(x(2) - vd)/m;
    %% QP for OC
    A = [phi1 -1;LgB 0;1 0;-1 0];
    b = [-phi0;LfB;ca*m*g;cd*m*g];
    H = [2/(m^2) 0; 0 2*psc];
    F = [-2*Fr/(m^2); 0];
    options = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','off');
    [u,fval,exitflag,output,lambda] = ...
       quadprog(H,F,A,b,[],[],[],[],[],options);

%% update dynamics

t = [0 dt];
[~,xx]=ode45('acc_real',t,x);
x = xx(end, 1:3);
result2(i+1,:) = [dt*i xx(end,2) u(1)/m C1 C2];

[~,xx]=ode45('acc_nomi',t,y);
y = xx(end, 1:3);


end
% 
figure(3)
subplot(4, 1, 1)
plot(result2(:, 1),result2(:, 2),'b-',[0,30],[24,24], 'k--',[0,30],[13.89,13.89], 'k--')
grid on
subplot(4, 1, 2)
plot(result2(:, 1),result2(:, 3)/g, 'b-',[0,30],[ca,ca], 'k--',[0,30],[-cd,-cd], 'k--')
axis([0 30 -0.75 0.75]); 
grid on
subplot(4, 1, 3)
plot(result2(:, 1),result2(:, 4), 'b-',[0,30],[0,0], 'k--')
axis([0 30 -5 75]); 
grid on
subplot(4, 1, 4)
plot(result2(:, 1),result2(:, 5), 'b-')
grid on

