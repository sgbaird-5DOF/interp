% po2el_test
clear; close all
rng(10)
npts = 4;
ptstmp = rand(3,npts).'; %to preserve rng properties
ptstmp(:,3) = ptstmp(:,3)-0.5;
pts = normr(ptstmp);
t = n2c(pts);

figure
hold on
sphplot('axview',[1.350469711107535e+02,32.524570056836382])
quivplot()
plot3(t{:},'*')

[x,y,z] = t{:};
[az,el,r] = cart2sph(x,y,z);

text(x,y,z,n2c(1:npts))

disp(rad2deg(el))