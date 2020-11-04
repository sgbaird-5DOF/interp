% CELL600  generates the vertices of a 600-cell (the 4D analog of the
% icosahedron), that is circumscribed within a unit hypersphere. This
% yields a regular grid on the unit hypersphere.
%-------------------------------------------------------------------------%
%Filename:  cell600.m
%Author:    Oliver Johnson
%Date:      8/8/2013
%
% Inputs:
%   s - A scalar indicating the scaling factor. Higher values of s
%       subdivide the grid and produce a denser set of points. The number 
%       of vertices for various values of s is provided below.
%
%                   s       # vertices      time (sec.)
%                  ---      ----------      -----------
%                   1          120             0.02
%                   2          2520            0.13
%                   3          34616           1.98
%
% Outputs:
%   w,t,p - Column vectors with #vertices elements, containing the
%           hyperspherical coordinates of the vertices.
%
% SEE ALSO: project4D.m, rot2q.m, q2rot.m
% http://en.wikipedia.org/wiki/600-cell
% http://en.wikipedia.org/wiki/Geodesic_grid
% http://en.wikipedia.org/wiki/Permutation_matrix
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/4558
%-------------------------------------------------------------------------%

function [w,t,p] = cell600(s)

if nargin < 1
    s = 1;
end
assert(isscalar(s) && s >= 1 && mod(s,1) == 0,'s must be a scalar integer greater than or equal to 1.')

%% Define standard 600-cell

phi = (1+sqrt(5))/2; %golden ratio
len = 1/phi; %edge length

%---first 16 points---%
a = zeros(16,4);
ind = 1;
for i = -1:2:1
    for j = -1:2:1
        for k = -1:2:1
            for l = -1:2:1
                a(ind,:) = [i j k l];
                ind = ind+1;
            end
        end
    end
end
a = a/2;

%---next 8 points---%
b = [1 0 0 0;...
    0 1 0 0;...
    0 0 1 0;...
    0 0 0 1];
b = [b;-b];

%---final 96 points---%
c0 = zeros(8,4);
ind = 1;
for i = [phi,-phi]/2
    for j = [1,-1]/2
        for k = [1/phi,-1/phi]/2
            c0(ind,:) = [i j k 0];
            ind = ind+1;
        end
    end
end

% get even permutations
p = perms(1:4);
I = eye(4);
iseven = false(size(p,1),1);
for i = 1:size(p,1)
    P = I(p(i,:),:); %permutation matrix of permutation p(i,:)
    if det(P) == 1
        iseven(i) = true;
    end
end
p = p(iseven,:);

% include permutations
c = zeros(0,4);
for i = 1:size(p,1)
    c = [c; c0(:,p(i,:))];
end

%---make unit quaternions---%
q = [a; b; c];
q = bsxfun(@rdivide,q,sqrt(sum(q.^2,2)));

%% Apply scale factor

if s > 1
    
    neighbs = nchoosek(1:4,3); %vertex neighbors in tetrahedra
    
    for i = 2:s
        
        %---get convex hull---%
        K = convhulln(q);
        
        for j = 1:size(K,1)
            %---subdivide each tetrahedron---%
            qtemp = q(K(j,:),:); %vertices of this cell
            q3 = cat(3,qtemp(neighbs(:,1),:),qtemp(neighbs(:,2),:),qtemp(neighbs(:,3),:));
            qnew = sum(q3,3)/3;
            
            %---project to hypersphere---%
            qnew = bsxfun(@rdivide,qnew,sqrt(sum(qnew.^2,2)));
            q = [q; qnew];
        end
    end
end

%% Convert to hyperspherical angles

[w,t,p] = q2rot(q);
