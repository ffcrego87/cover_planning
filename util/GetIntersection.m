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