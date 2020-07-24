function BPFZ_ids = inBPFZ(nAlist,vertices,R)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description:
% 
% Inputs:
%
% Outputs:
%
% Dependencies:
%
% Notes:
%
%--------------------------------------------------------------------------
%rotate points into lab frame
nArot = (R*nAlist.').'; %inverted rotation matrix
% nArot = nAlist;

%convert all query points to spherical
t = num2cell(nArot,1);
[az,el,~] = cart2sph(t{:});

az(az < 0) = az(az < 0) + pi;

%compute spherical triangle inequalities
% t = num2cell(vertices,1);
% [azvtx,elvtx,~] = cart2sph(t{:});
azvtx = vertices(:,1);
elvtx = vertices(:,2);

%get min & max spherical coordinates
azmin = min(azvtx);
elmin = min(elvtx);

azmax = max(azvtx);
elmax = max(elvtx);

tol = 1e-3;
r = @(a,b) a-b >= -tol;
%take only those in BP FZ
BPFZ_ids = r(az,azmin) & r(azmax,az) & r(el,elmin) & r(el,elmax); %alternative: inpolygon()

%---------------------------CODE GRAVEYARD---------------------------------
%{
BPFZ_ids = (abs(azmin-az) <= tol) & (abs(az-azmax) <= tol) & ...
	(abs(elmin-el) <= tol) & (abs(el-elmax) <= tol); %alternative: inpolygon()

% check1 = az >= azmin - tol;
% check2 = az <= azmax + tol;
% check3 = el >= elmin - tol;
% check4 = el <= elmax + tol;

%}
