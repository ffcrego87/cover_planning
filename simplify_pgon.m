function Pout = simplify_pgon(Pin)

diff2 = Pin(2:end,:)-Pin(1:(end-1),:);
dists = sqrt(diff2(:,1).^2+diff2(:,2).^2);

perimeterorig = sum(dists);

perimetertarget = 0.95*perimeterorig;

number_of_points = size(Pin,1);

Pout = Pin;
perimeter = perimeterorig;
while(number_of_points>2 && perimeter>perimetertarget)
    disp(number_of_points)
    disp(perimeter)
    perimeter = 0;
    for i=2:(number_of_points-1)
        Paux = Pout;
        Paux(i,:) = [];
        diff2 = Paux(2:end,:)-Paux(1:(end-1),:);
        dists = sqrt(diff2(:,1).^2+diff2(:,2).^2);
        perimeter_aux = sum(dists);
        if(perimeter_aux>perimeter)
            perimeter = perimeter_aux;
            Pcand = Paux;
        end
    end
    Pout = Pcand;
    number_of_points = number_of_points-1;
end

end