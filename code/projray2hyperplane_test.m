%project ray to hyperplane test
clear; close all

d = 5; %dimensionality
nvec = rand(1,d);
nvec = nvec./norm(nvec);

pt = rand(1,d);

a = projray2hyperplane(nvec,pt);
disp(' ')
disp('intersection between ray and hyperplane: ')
disp(a)
