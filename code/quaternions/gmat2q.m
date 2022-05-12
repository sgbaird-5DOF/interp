function q = gmat2q(g)
arguments
    g(3,3,:) double
end
% GMAT2Q converts rotation matrices to their equivalent quaternion
% representation.
%-------------------------------------------------------------------------%
%Filename:  gmat2q.m
%Author:    Oliver Johnson
%Date:      3/28/2020
%
% Inputs:
%   g - A 3-by-3-by-npts array of rotation matrices.
%
% Outputs:
%   q - An npts-by-4 array of equiavalent quaternions.
%
% CURRENT METHOD
% Shepperd's Method (as described in
% https://doi.org/10.1007/978-3-319-93188-3_5). Works correctly for both
% skew-symmetric and symmetric (e.g. rotations by pi) rotation matrices.
% 
% OLD METHOD
% https://doi.org/10.1007/978-3-319-93188-3_5 (Only works for
% skew-symmetric rotation matrices)
%
% OLD OLD METHOD (4/1/2013):
% Code adapted from http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
% This algorithm only works when the rotation matrix is stricly skew
% symmetric, in the case that it is symmetric the difference of the off
% diagonal terms is zero so the sign function returns zero. I need to get a
% better method to assign the signs consistently. Take a look at MTEX's
% mat2quat function.
%-------------------------------------------------------------------------%

%---pre-allocate---%
npts = size(g,3);
q = zeros(npts,4);

%---compute transformation---%
% compute the traces
t = g(1,1,:)+g(2,2,:)+g(3,3,:);

% figure out which term is largest (this will determine which case to use)
[~,id] = max([t,g(1,1,:),g(2,2,:),g(3,3,:)],[],2);

% case 1
id1 = id == 1;
val = sqrt(1+t(:,:,id1));
q(id1,1) = 0.5*squeeze(val);
q(id1,2) = 0.5*squeeze(( g(3,2,id1) - g(2,3,id1) ) ./ val);
q(id1,3) = 0.5*squeeze(( g(1,3,id1) - g(3,1,id1) ) ./ val);
q(id1,4) = 0.5*squeeze(( g(2,1,id1) - g(1,2,id1) ) ./ val);

% case 2
id2 = id == 2;
val = sqrt(1+g(1,1,id2)-g(2,2,id2)-g(3,3,id2));
q(id2,1) = 0.5*squeeze(( g(3,2,id2) - g(2,3,id2) ) ./ val);
q(id2,2) = 0.5*squeeze(val);
q(id2,3) = 0.5*squeeze(( g(1,2,id2) + g(2,1,id2) ) ./ val);
q(id2,4) = 0.5*squeeze(( g(3,1,id2) + g(1,3,id2) ) ./ val);

% case 3
id3 = id == 3;
val = sqrt(1-g(1,1,id3)+g(2,2,id3)-g(3,3,id3));
q(id3,1) = 0.5*squeeze(( g(1,3,id3) - g(3,1,id3) ) ./ val);
q(id3,2) = 0.5*squeeze(( g(1,2,id3) + g(2,1,id3) ) ./ val);
q(id3,3) = 0.5*squeeze(val);
q(id3,4) = 0.5*squeeze(( g(2,3,id3) + g(3,2,id3) ) ./ val);

% case 4
id4 = id == 4;
val = sqrt(1-g(1,1,id4)-g(2,2,id4)+g(3,3,id4));
q(id4,1) = 0.5*squeeze(( g(2,1,id4) - g(1,2,id4) ) ./ val);
q(id4,2) = 0.5*squeeze(( g(3,1,id4) + g(1,3,id4) ) ./ val);
q(id4,3) = 0.5*squeeze(( g(3,2,id4) + g(2,3,id4) ) ./ val);
q(id4,4) = 0.5*squeeze(val);

% choose overall sign that makes q0 positive
idChangeSign = q(:,1) < 0;
q(idChangeSign,:) = -q(idChangeSign,:);

% %% Old Method
% %---compute transformation---%
% 
% %q0
% t = g(1,1,:)+g(2,2,:)+g(3,3,:); %traces
% t = t(:);
% isPos = t > 0; 
% q(isPos,1) = 0.5*(sqrt(1+t(isPos)));
% q(~isPos,1) = 0.5*(sqrt(...
%     squeeze((g(3,2,~isPos)-g(2,3,~isPos)).^2+(g(1,3,~isPos)-g(3,1,~isPos)).^2+(g(2,1,~isPos)-g(1,2,~isPos)).^2)./...
%     squeeze(3-t(~isPos))));
% 
% %q1
% t = g(1,1,:)-g(2,2,:)-g(3,3,:);
% t = t(:);
% isPos = t > 0; 
% q(isPos,2) = 0.5*(sqrt(1+t(isPos)));
% q(~isPos,2) = 0.5*(sqrt(...
%     squeeze((g(3,2,~isPos)-g(2,3,~isPos)).^2+(g(1,2,~isPos)+g(2,1,~isPos)).^2+(g(3,1,~isPos)+g(1,3,~isPos)).^2)./...
%     squeeze(3-t(~isPos))));
% 
% %q2
% t = -g(1,1,:)+g(2,2,:)-g(3,3,:);
% t = t(:);
% isPos = t > 0; 
% q(isPos,3) = 0.5*(sqrt(1+t(isPos)));
% q(~isPos,3) = 0.5*(sqrt(...
%     squeeze((g(1,3,~isPos)-g(3,1,~isPos)).^2+(g(1,2,~isPos)+g(2,1,~isPos)).^2+(g(2,3,~isPos)+g(3,2,~isPos)).^2)./...
%     squeeze(3-t(~isPos))));
% 
% %q3
% t = -g(1,1,:)-g(2,2,:)+g(3,3,:);
% t = t(:);
% isPos = t > 0; 
% q(isPos,4) = 0.5*(sqrt(1+t(isPos)));
% q(~isPos,4) = 0.5*(sqrt(...
%     squeeze((g(2,1,~isPos)-g(1,2,~isPos)).^2+(g(3,1,~isPos)+g(1,3,~isPos)).^2+(g(3,2,~isPos)+g(2,3,~isPos)).^2)./...
%     squeeze(3-t(~isPos))));
% 
% %signs
% q(:,2) = copysign(q(:,2), squeeze(g(3,2,:)-g(2,3,:)));
% q(:,3) = copysign(q(:,3), squeeze(g(1,3,:)-g(3,1,:)));
% q(:,4) = copysign(q(:,4), squeeze(g(2,1,:)-g(1,2,:)));
% 
% % %% Old Old Method
% % %---compute transformation---%
% % t = g(1,1,:)+g(2,2,:)+g(3,3,:); %traces
% % t = t(:);
% % r = sqrt(1+t);
% % q(:,1) = 0.5*r;
% % q(:,2) = copysign(0.5*sqrt(1+g(1,1,:)-g(2,2,:)-g(3,3,:)), g(3,2,:)-g(2,3,:));
% % q(:,3) = copysign(0.5*sqrt(1-g(1,1,:)+g(2,2,:)-g(3,3,:)), g(1,3,:)-g(3,1,:));
% % q(:,4) = copysign(0.5*sqrt(1-g(1,1,:)-g(2,2,:)+g(3,3,:)), g(2,1,:)-g(1,2,:));
% 
% end
% 
% function x = copysign(x,y)
% 
% sgn = ones(size(y));
% sgn(y < 0) = -1;
% 
% x = sgn.*x; % note we don't need abs(x) because the equations above guarantee x >= 0 to begin with
% 
% end