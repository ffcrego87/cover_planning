function uncovered_area = uncovered_area_comp(pos1sol,pos2sol,d,Lsize,Hsize,P,D,L)

disc = d/10;

xrange = 0:disc:Lsize;
yrange = 0:disc:Hsize;

Cellpos1 = kron(xrange,ones(size(yrange)));
Cellpos2 = kron(ones(size(xrange)),yrange);

idxs_to_remove = [];


figure
hold on

length(Cellpos1)


plot(pos1sol,pos2sol,'k-')

for i=1:max(L)
    Pi1 = P(L==i,:);
    Pi = [Pi1;Pi1(1,:)];
    plot(polyshape(P(L==i,1),P(L==i,2)));
    for k=1:length(Cellpos1)
        if inpolygon(Cellpos1(k),Cellpos2(k),Pi(:,1),Pi(:,2))
            idxs_to_remove = [idxs_to_remove;k];
        end
    end
end

Cellpos1(idxs_to_remove) = [];
Cellpos2(idxs_to_remove) = [];

idxs_to_remove = [];

for k=1:length(Cellpos1)
    for l = 1:length(pos1sol)
        if norm([Cellpos1(k);Cellpos2(k)]-[pos1sol(l);pos2sol(l)])<(d/2)
            idxs_to_remove = [idxs_to_remove;k];
            continue;
        end
    end
end

Cellpos1(idxs_to_remove) = [];
Cellpos2(idxs_to_remove) = [];

uncovered_area = length(Cellpos1)*disc^2;

plot(Cellpos1,Cellpos2,'k.')

end

