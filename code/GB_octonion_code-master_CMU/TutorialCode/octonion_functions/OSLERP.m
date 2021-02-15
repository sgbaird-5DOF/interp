function O2 = OSLERP(o1,o2,omega,nt)
arguments
   o1
   o2
   omega
   nt
%    extend(1,1) double = 0
end
%% INPUT DATA 
%
% o1, o2: symmetrized octonions: see example 3
% w: GBOM angle (see GBdist function, example 3)
% nt: is number of interpolated GB's to output (including end points)
%
%% OUTPUT
%
% O2: list of interpolated GB octonions

qA = o1(:,1:4); qB = o1(:,5:8);
% qC = o2(:,1:4); qD = o2(:,5:8);

%% OSLERP 

% nt = 10; %number of interpolated points
% start = 0-extend;
% finish = 1+extend;
% omega = (1+2*extend)*omega;
% t = linspace(start,finish,nt);
t = linspace(0,1,nt);

perms = [-qA qB; qA -qB; qA qB; -qA -qB];
nperms = length(perms(:,1));
tol = 1e-5;
check = false;
for k = 1:nperms
%     disp(k)
    o1i = perms(k,:);
    O2 = zeros(nt,8);
    
    for i = 1:nt
        ti = t(i);
        theta = omega/2;

        s = sin(theta);
        o = sin((1-ti)*theta)*(o1i)/s + sin(ti*theta)*(o2)/s;
        O2(i,:) = o;
    end
    
    testnorm = vecnorm(O2,2,2);
    checknorm = testnorm-sqrt(2);
    
    if abs(sum(checknorm)) < tol
        check = true;
        break
    end
   
end

if ~check
    disp('WARNING: normed OSLERP trajectory not found, are octonions symmetrized?')
end


end


