function [pos1sol,pos2sol,t] = OCP_plan4(Lsize,Hsize,d,P,~,L,speed,pos1init,pos2init,Tinit)
% An optimal control problem (OCP),
% solved with direct multiple-shooting.
% N - number of control intervals
% T - final time

import casadi.*

N=length(pos1init)-1;

Xinit = [pos1init';pos2init'];

dtinit = Tinit/N;
Uinit = (Xinit(:,2:end)-Xinit(:,1:(end-1)))./dtinit;

rad_g = d/2;
xrange = 0:(d/5):Lsize;
yrange = 0:(d/5):Hsize;

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
X = opti.variable(2,N+1); % state trajectory
pos  = X;
pos1 = X(1,:);
pos2 = X(2,:);
U = opti.variable(2,N);   % control trajectory (throttle)

% ---- a-priori distribution ---------
apf_disc = zeros(length(yrange),length(xrange));
for i=1:length(yrange)
    for j=1:length(xrange)
        fprintf('apf_disk i:%d j:%d\n',i,j)
        apf_disc(i,j) = apf_fix;
    end
end
for i=1:max(L)
    Pi1 = P(L==i,:);
    Pi = [Pi1;Pi1(1,:)];
    for j=1:length(xrange)
        for k=1:length(yrange)
            fprintf('objective i:%d j:%d k:%d\n',i,j,k)
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
            fprintf('penalty i:%d j:%d k:%d\n',i,j,k)
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
            fprintf('objective i:%d j:%d k:%d\n',i,j,k)
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
            fprintf('penalty i:%d j:%d k:%d\n',i,j,k)
            if penalty_disc(i,j)
                distsq=(pos1(k)-xrange_penalty(j))^2+(pos2(k)-yrange_penalty(i))^2;
                penalty = penalty + exp(-(distsq-(d/8)^2));
            end
        end
    end
end

opti.minimize(10^-4*sum(sqrt(U(1,:).^2+U(2,:).^2))+cost+penalty);

% ---- dynamic constraints --------
f = @(x,u) u; % dx/dt = f(x,u)

for k=1:N % loop over control intervals
    % Euler integration
    x_next = X(:,k) + dt*f(X(:,k),U(:,k));
    opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

% ---- control constraints --------
opti.subject_to((U(1,:).^2+U(2,:).^2)<=speed^2);% control is limited

% ---- state constraints --------
opti.subject_to(pos(:,0)==[pos1init(1);pos2init(1)])
opti.subject_to(0<=pos(1,:)<=max(xrange))
opti.subject_to(0<=pos(2,:)<=max(yrange))

% ---- Initial guess ------
opti.set_initial(X,Xinit);
opti.set_initial(T,Tinit);
opti.set_initial(U,Uinit);

%% ---- solve NLP              ------
disp('solve')
p_opts = struct('expand',true);
s_opts = struct('max_iter',10^4,'hessian_approximation','limited-memory');
opti.solver('ipopt',p_opts,s_opts); % set numerical backend

try
    sol = opti.solve();   % actual solve
    %% post-processing
    pos1sol_prov=sol.value(pos1);
    pos2sol_prov=sol.value(pos2);
catch
    warning('Problem using casadi.');
    %% post-processing
    pos1sol_prov=opti.debug.value(pos1);
    pos2sol_prov=opti.debug.value(pos2);
end

waypoints = [];
for i=1:(length(pos1sol_prov)-1)
    path = pathfinder_poly([pos1sol_prov(i) pos2sol_prov(i)], [pos1sol_prov(i+1) pos2sol_prov(i+1)], P, L);
    waypoints = [waypoints;path(1:(end-1),:)];
end
waypoints = [waypoints;[pos1sol_prov(end) pos2sol_prov(end)]];

pos1sol = waypoints(:,1)';
pos2sol = waypoints(:,2)';

t = zeros(size(pos1sol));
for i=2:length(t)
    t(i)=t(i-1)+norm([pos1sol(i)-pos1sol(i-1);pos2sol(i)-pos2sol(i-1)])/speed;
end

end