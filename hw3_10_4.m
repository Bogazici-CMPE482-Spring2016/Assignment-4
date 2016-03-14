%a------------

teta = pi/3;
F = [-cos(teta) sin(teta); sin(teta) cos(teta)];
J = [cos(teta) sin(teta); -sin(teta) cos(teta)];
v = [2;1];
fRes = F*v
jRes = J*v    

atan2(fRes(2),fRes(1))
atan2(v(2),v(1))
figure;
plotv(F*v,'-');
hold
plotv(v,'-');
plotv(J*v,'-');
 % J rotates plane by the angle teta in the clockwise direction
 % F takes the symmetric of v according to y axis, then rotates plane by teta
 % degree in the counter clockwise direction.
%-------------