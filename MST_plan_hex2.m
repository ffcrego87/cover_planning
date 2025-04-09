function [pos1sol,pos2sol,t] = MST_plan_hex2(MagFactor,Lsize,Hsize,d,nP,npol,P,D,L,pgons,speed)

ani=1;% 0 = isotrope = trajectoire erratiques => directions isotrope, en partie aléatoire
% 1 = anisotrope = trajetoire horizontale prvilégiée ("polarisation") => pénalisation des directions non privilégiées
[MSTMap] = generateMSTMapHex(P, L, D, ani, d, Lsize, Hsize, nP, npol, MagFactor);
path = extractSkeletonPath(MSTMap, pgons, d, MagFactor, Lsize, Hsize, npol);
pos1sol=path(:,1);
pos2sol=path(:,2);
t = zeros(size(pos1sol));
for i=2:length(t)
    t(i)=t(i-1)+norm([pos1sol(i)-pos1sol(i-1);pos2sol(i)-pos2sol(i-1)])/speed;
end

end