% An optimal control problem (OCP),
% solved with direct multiple-shooting.

%% ---- preliminaries ------
clearvars
close all
clc
%% ---- problem ------
disp('problem')

N=491;
% Size of the map
ct=sqrt(3)/2;
d = 5;
Hsize=9*ct*d;
Lsize=9*d;
sigma = 10; % Rayon moyen des obstacles
nsom = 12; % Nombre maximum de sommets
npol = 3; % Number of polygons

[P,D,L]=gen_polygons(Lsize,Hsize,sigma,npol,nsom,N);

N = 100/4; % number of control intervals, critical try with 5o & compare
T = 100;  % final time, 1000
speed = 0.4;
rmax = 1;
rmin = -1;
pos0 = [0;20];
yaw0 = 0;



import casadi.*


rad_g = d/2;
xrange = 0:2:Lsize;
yrange = 0:2:Hsize;

apf_fix = 1;
rise_fac = 100;
gamma_size = 20;
gamma_distsq = @(x) gamma_size*(pi/2-atan(rise_fac*((x/rad_g^2)-1)))/pi;

%% ---- Initialization ------
disp('Initialization')

% ---- Optimization initialization ------

opti = casadi.Opti(); 

% ---- decision variables ---------

Cellpos1 = repmat(xrange,length(yrange),1);
Cellpos2 = repmat((yrange)',1,length(xrange));
X = opti.variable(3,N+1); % state trajectory
pos  = X([1 2],:);
pos1 = X(1,:);
pos2 = X(2,:);
yaw = X(3,:);
U = opti.variable(N);   % control trajectory (throttle)

% ---- a-priori distribution ---------
apf_disc = zeros(length(yrange),length(xrange));
for i=1:length(yrange)
    for j=1:length(xrange)
        apf_disc(i,j) = apf_fix;
    end
end
for i=1:npol
    Pi1 = P(L==i,:);
    Pi = [Pi1;Pi1(1,:)];
    for j=1:length(xrange)
        for k=1:length(yrange)
            if inpolygon(xrange(j),yrange(k),Pi(:,1),Pi(:,2))
                apf_disc(k,j) = 0;
            end
        end
    end
end

% ---- objective          ---------
dt = T/N; % length of a control interval
cost = 0;
for i=1:length(yrange)
    for j=1:length(xrange)
        indic = 0;
        for k=1:(N+1)
            distsq=(pos1(k)-Cellpos1(i,j))^2+(pos2(k)-Cellpos2(i,j))^2;
            indic = indic+dt*gamma_distsq(distsq);
        end
        cost = cost + apf_disc(i,j)*exp(-indic);
    end
end

opti.minimize(cost);

% ---- dynamic constraints --------
f = @(x,u) [speed*cos(x(3));speed*sin(x(3));u]; % dx/dt = f(x,u)

for k=1:N % loop over control intervals
   % Runge-Kutta 4 integration
   k1 = f(X(:,k),         U(k));
   k2 = f(X(:,k)+dt/2*k1, U(k));
   k3 = f(X(:,k)+dt/2*k2, U(k));
   k4 = f(X(:,k)+dt*k3,   U(k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4); 
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

% ---- control constraints --------
opti.subject_to(rmin<=U<=rmax);% control is limited

% ---- boundary conditions --------
opti.subject_to(pos(:,1)==pos0);   % start at position 0 ...
opti.subject_to(yaw(1,1)==yaw0);   % start at yow 0 ...

%% ---- solve NLP              ------
disp('solve')
opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solve

% post-processing
t = linspace(0,sol.value(T),N+1);

pos1sol=sol.value(pos1);
pos2sol=sol.value(pos2);



%% ---- post-processing        ------
disp('post-processing')






Cellpos3 = Cellpos2;
for i=1:length(yrange)
    for j=1:length(xrange)
        indic = 0;
        for k=1:(N+1)
            distsq=(pos1sol(k)-Cellpos1(i,j))^2+(pos2sol(k)-Cellpos2(i,j))^2;
            indic = indic+gamma_distsq(distsq);
        end
        Cellpos3(i,j) = apf_disc(i,j)*exp(-indic);
    end
end

%% ---- plotting        ------
disp('plotting')

figure
plot(dist_sp,g_sp);
ylabel('gamma');
xlabel('dist (m)');
%print(['Figures' filesep 'detect_func'],'-dpng')

figure
hold on
plot(t,sol.value(pos1));
plot(t,sol.value(pos2));
xlabel('Time [s]');
legend('pos1','pos2','Location','northwest')
%print(['Figures' filesep 'OCP_sol'],'-dpng')

figure
hold on
contourf(Cellpos1,Cellpos2,Cellpos3)
plot(sol.value(pos1),sol.value(pos2),'r-');
xlabel('X');
ylabel('Y');
axis equal
grid on
legend('prob. of not detecting','traj','Location','northwest')
%print(['Figures' filesep 'OCP_traj1'],'-dpng')

figure
hold on
contourf(Cellpos1,Cellpos2,apf_disc)
plot(sol.value(pos1),sol.value(pos2),'r-');
xlabel('X');
ylabel('Y');
axis equal
grid on
legend('a-priori function','traj','Location','northwest')
%print(['Figures' filesep 'OCP_traj2'],'-dpng')

control=sol.value(U);
figure
hold on
plot(t(1:(end-1)),control,'r-')
legend('control','Location','northwest')
%print(['Figures' filesep 'OCP_speed'],'-dpng')

figure
spy(sol.value(jacobian(opti.g,opti.x)))
xlabel('decision variables')
ylabel('constraints')
%print(['Figures' filesep 'jac_sp'],'-dpng')