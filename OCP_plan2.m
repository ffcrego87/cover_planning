function [pos1sol,pos2sol,t] = OCP_plan2(Lsize,Hsize,d,P,~,L,speed,rmax,pos0,yaw0,N)
% An optimal control problem (OCP),
% solved with direct multiple-shooting.
% N - number of control intervals
% T - final time

import casadi.*

rmin = -rmax;

rad_g = d/2;
xrange = 0:(d/2):Lsize;
yrange = 0:(d/2):Hsize;

apf_fix = 1;
rise_fac = 10000;
gamma_size = 100;
gamma_distsq = @(x) gamma_size*(pi/2-atan(rise_fac*((x/rad_g^2)-1)))/pi;


xrange_penalty = 0:(d/16):Lsize;
yrange_penalty = 0:(d/16):Hsize;

% rad_penalty = d/4;
% rise_fac_penalty = 10000;
% gamma_size_penalty = 10;
% gamma_distsq_penalty = @(x) gamma_size_penalty*(pi/2-atan(rise_fac_penalty*((x/rad_penalty^2)-1)))/pi;

%% ---- Initialization ------
disp('Initialization')

% ---- Optimization initialization ------

opti = casadi.Opti();

% ---- decision variables ---------

Cellpos1 = repmat(xrange,length(yrange),1);
Cellpos2 = repmat((yrange)',1,length(xrange));
T = opti.variable(1,1);
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
for i=1:max(L)
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

% ---- penalty distribution ---------
penalty_disc = zeros(length(yrange_penalty),length(xrange_penalty));
for i=1:max(L)
    Pi1 = P(L==i,:);
    Pi = [Pi1;Pi1(1,:)];
    for j=1:length(xrange_penalty)
        for k=1:length(yrange_penalty)
            if inpolygon(xrange_penalty(j),yrange_penalty(k),Pi(:,1),Pi(:,2))
                penalty_disc(k,j) = apf_fix;
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

% ---- penalty          ---------
penalty = 0;
for i=1:length(yrange_penalty)
    for j=1:length(xrange_penalty)
        for k=1:(N+1)
            if penalty_disc(i,j)
                distsq=(pos1(k)-xrange_penalty(j))^2+(pos2(k)-yrange_penalty(i))^2;
                penalty = penalty + exp(-(distsq-(d/8)^2));
            end
        end
    end
end

opti.minimize(cost+penalty+10^-4*T);

% ---- dynamic constraints --------
f = @(x,u) [speed*cos(x(3));speed*sin(x(3));u]; % dx/dt = f(x,u)

% for k=1:N % loop over control intervals
%     % Runge-Kutta 4 integration
%     k1 = f(X(:,k),         U(k));
%     k2 = f(X(:,k)+dt/2*k1, U(k));
%     k3 = f(X(:,k)+dt/2*k2, U(k));
%     k4 = f(X(:,k)+dt*k3,   U(k));
%     x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
%     opti.subject_to(X(:,k+1)==x_next); % close the gaps
% end

for k=1:N % loop over control intervals
    % Euler integration
    x_next = X(:,k) + dt*f(X(:,k),U(k));
    opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

% ---- control constraints --------
opti.subject_to(rmin<=U<=rmax);% control is limited

% ---- boundary conditions --------
opti.subject_to(pos(:,1)==pos0);   % start at position 0 ...
opti.subject_to(yaw(1,1)==yaw0);   % start at yow 0 ...

% ---- state constraints --------
opti.subject_to(0<=pos(1,:)<=max(xrange))
opti.subject_to(0<=pos(2,:)<=max(yrange))

% ---- Initial guess

%% ---- solve NLP              ------
disp('solve')
opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solve

%% post-processing
t = linspace(0,sol.value(T),N+1);

pos1sol=sol.value(pos1);
pos2sol=sol.value(pos2);
end