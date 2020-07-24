% [1] S. Patala, C.A. Schuh, A continuous and one-to-one coloring scheme for misorientations, Acta Mater. 59 (2011) 554–562. doi:10.1016/j.actamat.2010.09.058.

function [h,s,v] = q2hsv(q)

%% convert to rodriguez parameters (Eq. 2.1)
x = q(:,2)./q(:,1);
y = q(:,3)./q(:,1);
z = q(:,4)./q(:,1);

%% Eq 2.2

x1 = x;
y1 = x.*(y+z)./(1-x);
z1 = x.*z.*(y+z)./(y.*(1-x));

doChange = x < 1/3 | mytan(z,y) < (1-2*x)./x;

x1(doChange) = x(doChange);
y1(doChange) = y(doChange);
z1(doChange) = z(doChange);

%% Eq. 2.3

g = q2gmat(rot2q(3*pi/8,pi/2,0));

x2y2z2 = g*[x1(:).'-tan(pi/8); y1(:).'; z1(:).'];

x2 = x2y2z2(1,:).';
y2 = x2y2z2(2,:).';
z2 = x2y2z2(3,:).';

%% Eq. 2.4

x3 = x2;
y3 = y2.*(1+((y2./z2)*tan(pi/8)));
z3 = z2+(y2*tan(pi/8));

%% Eq. 2.5

x4 = x3;
y4 = y3.*cos(pi/8)./tan(pi/8);
z4 = z3-(x3./cos(pi/8));

%% Eq. 2.6

% phi = atan(-x4./y4); % should this be atan2?
phi = atan2(-x4,y4); % should this be atan2?

x5 = x4.*(sin(phi)+abs(cos(phi)));
y5 = y4.*(sin(phi)+abs(cos(phi)));
z5 = z4;

%% Eq. 2.7

% phi1 = atan(-y5./x5);
phi1 = atan2(-y5,x5);

h = -sqrt((x5.^2)+(y5.^2)).*cos(2*phi1);
s = sqrt((x5.^2)+(y5.^2)).*sin(2*phi1);
v = z5;

end

% helper function that defines tan(0/0) = 0
function val = mytan(num,denom)

val = tan(num./denom);

isZeroNumerator = num == 0;
val(isZeroNumerator) = 0;

end