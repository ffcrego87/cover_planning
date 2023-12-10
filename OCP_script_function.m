% An optimal control problem (OCP),
% solved with direct multiple-shooting.

%% ---- preliminaries ------
clearvars
close all
clc
%% ---- problem ------
disp('problem')


N = 493;
% Size of the map
Hsize=90;
Lsize=120;
d = 25;
sigma = 20;% Rayon moyen des obstacles
nsom = 12;% Nombre maximum de sommets
nk = 0;
npol = 2;

[P,D,L]=gen_polygons(Lsize,Hsize,sigma,npol,nsom,N);

%%

% Vehicle settings
speed = 0.4;
pos0 = [0;20];
yaw0 = 0;
T = 1500;  % final time

% OCP settings
N = 80; % number of control intervals

% Problem settings OCP
rmax = 5;

[pos1sol,pos2sol,t] = OCP_plan(Lsize,Hsize,d,P,D,L,speed,T,rmax,pos0,yaw0,N);

%% ---- post-processing        ------
disp('post-processing')

rad_g = d/2;
xrange = 0:(d/2):Lsize;
yrange = 0:(d/2):Hsize;

apf_fix = 1;
rise_fac = 100;
gamma_size = 20;
gamma_distsq = @(x) gamma_size*(pi/2-atan(rise_fac*((x/rad_g^2)-1)))/pi;

Cellpos1 = repmat(xrange,length(yrange),1);
Cellpos2 = repmat((yrange)',1,length(xrange));

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
hold on
plot(t,pos1sol);
plot(t,pos2sol);
xlabel('Time [s]');
legend('pos1','pos2','Location','northwest')
%print(['Figures' filesep 'OCP_sol'],'-dpng')

figure
hold on
contourf(Cellpos1,Cellpos2,Cellpos3)
plot(pos1sol,pos2sol,'r-');
xlabel('X');
ylabel('Y');
axis equal
grid on
legend('prob. of not detecting','traj','Location','northwest')
%print(['Figures' filesep 'OCP_traj1'],'-dpng')

figure
hold on
contourf(Cellpos1,Cellpos2,apf_disc)
plot(pos1sol,pos2sol,'r-');
xlabel('X');
ylabel('Y');
axis equal
grid on
legend('a-priori function','traj','Location','northwest')
%print(['Figures' filesep 'OCP_traj2'],'-dpng')