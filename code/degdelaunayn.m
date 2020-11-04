function TRI = degdelaunayn(pts)
% DEGDELAUNAYN compute delaunay triangulation of d-1 hyperplane in d-dimensions
%  (i.e. degenerate convex hull) by projecting onto "thin"
%  dimension using singular value decomposition.
% 
% Input:
% 
% --- pts === n x d matrix of input points corresponding to d-1 hyperplane
% 
% Output:
% 
% --- K				=== indices of sub-facet vertices
% 
% --- vertices	=== n x d vertices
% 
% --- area			=== area of full facet (i.e. the degenerate convex hull)
%--------------------------------------------------------------------------
%}
%%
%project points to d-1
d = size(pts,2);
[U,S]=svd(bsxfun(@minus,pts,mean(pts)),0);
newpts = U*S(:,1:d-1);

%compute convex hull
TRI=delaunayn(newpts);

%convert back to d
%vertices=pts(K,:);

%trisurf(K,v(:,1),newpts(:,2),zeros(1,length(newpts)))
%axis equal


%https://www.mathworks.com/matlabcentral/answers/352830-convhull-convhulln-data-is-coplanar-data-is-degenerate-in-at-least-one-dimension-i-know-it-is-but

%{
%-------------------------------CODE GRAVEYARD-----------------------------
TESTING
%pts = -eye(d)+1;

% vec = {0:1};
% mat = repelem(vec,d);
% pts = allcomb(mat{:});
% deg = sum(pts,2);
% pts = pts(find(deg == d-1),:);
%}