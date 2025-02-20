function [P,D,L]=gen_polygons(Lsize,Hsize,sigma,npol,nsom,N)


%npol Number of polygons
%sigma Rayon moyen des obstacles
%nsom Nombre maximum de sommets

rng(N)
P = zeros(1000,2);% Nodes of polygons
D = zeros(1000,2);% Edges
L = zeros(1000,1);% Liste des labels des polygones

nP = 0;% Nombre de points dans P
nk = 0;
for j=1:npol
    q = 0;
    nQ = 1;
    while q<=nQ
        while q == 0
            nk = nk+1;
            x=Lsize/4+(Lsize/2)*rand();
            y=Hsize/4+(Hsize/2)*rand();
            rho = sigma*rand(nsom,1);
            theta = 2*pi*rand(nsom,1);
            Q = [x+rho.*cos(theta),y+rho.*sin(theta)];
            Q0 = min(Q);
            Q1 = max(Q);
            if min(Q0)>0 && Q1(1)<Lsize && Q1(2)<Hsize
                [k,~] = convhull(Q);
                Q = Q(k,:);
                nQ = length(Q)-1;
                q = 1;
                p = 1;
                % Test if two polygons are included
                i=1;
                included = false;
                while i<j && not(included)
                    Pi1 = P(L==i,:);
                    Pi = [Pi1;Pi1(1,:)];
                    included = inpolygon(Q(1,1),Q(1,2),Pi(:,1),Pi(:,2));
                    included = included || inpolygon(Pi(1,1),Pi(1,2),Q(:,1),Q(:,2));
                    i = i+1;
                end
                if included
                    q = 0;
                end
            else
                q=0;
            end
        end
        if j>1
            q1 = mod(q,nQ)+1;
            R = GetIntersection(P,D(p,1),D(p,2),Q,q,q1);
            if isempty(R)
                p = p+1;
                if p>nP
                    p = 1;
                    q = q+1;
                end
            else
                q = 0;
            end
        else
            p = 1;
            q = q+1;
        end
    end
    P(nP+(1:nQ),:) = Q(1:nQ,:);% Only nodes
    D(nP+(1:nQ),:) = [nP+(1:nQ) ; nP+(2:nQ), nP+1 ]';
    L(nP+(1:nQ)) = j;
    nP = nP+nQ;
end
P=P(1:nP,:);
D=D(1:nP,:);
L=L(1:nP);

end

%%
function R = GetIntersection(P,i1,i2,Q,j1,j2)
b = (P(i1,1)-P(i2,1))*(Q(j2,2)-Q(j1,2))-(P(i1,2)-P(i2,2))*(Q(j2,1)-Q(j1,1));
a = (Q(j1,1)-P(i2,1))*(Q(j2,2)-Q(j1,2))-(Q(j1,2)-P(i2,2))*(Q(j2,1)-Q(j1,1));
t = a / b;
if t>0 && t<1
    R = [t*P(i1,1)+(1-t)*P(i2,1) , t*P(i1,2)+(1-t)*P(i2,2)];
    s = (R(1)-Q(j1,1))*(Q(j2,1)-Q(j1,1))+(R(2)-Q(j1,2))*(Q(j2,2)-Q(j1,2));
    s = s/((Q(j2,1)-Q(j1,1))^2+(Q(j2,2)-Q(j1,2))^2);
    if s>0 && s<1
        return
    end
end
R = [];

end
