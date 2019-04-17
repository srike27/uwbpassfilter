function [a,b,c,d] = ytoabcd(y11,y12,y21,y22)
a = -y22./y21;
b = -1./y21;
c = (y12.*y21 - y11.*y22)./y21;
d = -y11./y21;
end