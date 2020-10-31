function s = vcell_solidangle(P, K, xyz)
% s = vcell_solidangle(P, K)
% s = vcell_solidangle(P, K, xyz)
%
% Compute the solid angles of All voronoi cell
%
% P is (3 x m) array of vertices, coordinates of the vertices of voronoi diagram
% K is (n x 1) cell, each K{j} contains the indices of the voronoi cell
% xyz is (3 x n) optional knot points to guide vcell_solidangle to compute
%   the solid angle of the "right" cell containing the node (and not the
%   complement cell)
%
% Restrictions:
% - P must be unit vectors
% - For each cell vertices must be counter-clockwise oriented when looking from outside
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 16/July/2014
%
% See also: voronoisphere

% Turn it to false to maximize speed
check_unit = true;

if check_unit
    u = sum(P.^2, 1);
    if any(abs(u-1) > 1e-6)
        error('vcell_solidangle:Pnotnormalized', ...
            'vcell_solidangle: P must be unit vectors');
    end
end

if nargin < 3
    s = cellfun(@(k) one_vcell_solidangle(P(:,k)), K);
else
    n = size(xyz,2);
    s = arrayfun(@(j) one_vcell_solidangle(P(:,K{j}),xyz(:,j)), (1:n).');
end


end % vcell_solidangle

%%
function omega = one_vcell_solidangle(v, center)
% omega = one_vcell_solidangle(v)
% Compute the solid angle of spherical polygonal defined by v
% v is (3 x n) matrix, each column is the coordinates of the vertex (unit vector, not checked)
% "correctly" oriented
%
% Ref: A. van Oosterom, J. Strackee: "A solid angle of a plane triangle."
%   IEEE Trans. Biomed. Eng. 30:2 (1983); 125–126. 
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 16/July/2014
%
% See also: vcell_solidangle

if nargin < 2
    n = size(v,2);
    s = zeros(1,n-2);
    for i = 2:n-1
        T = v(:,[1 i i+1]);
        num = det(T);
        denom = 1 + sum(sum(T .* T(:,[2 3 1]), 1), 2);
        s(i-1) = num / denom;
    end
else
    v(:,end+1) = v(:,1);
    n = size(v,2);
    s = zeros(1,n-1);
    for i = 1:n-1
        T = [center, v(:,[i i+1])];
        num = det(T);
        denom = 1 + sum(sum(T .* T(:,[2 3 1]), 1), 2);
        s(i) = num / denom;
    end
end

omega = atan(s);
omega = 2 * sum(omega);

end % one_vcell_solidangle
