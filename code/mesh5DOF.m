function [five,sept,o] = mesh5DOF(featureType,ctrcuspQ,varargin)
% MESH5DOF  generate "five" and "o" for a five degree-of-freedom fundamental zone (deprecated)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-02
%
% Description:
%
% Inputs:
%
% Outputs:
%
% Dependencies:
%		misFZfeatures.mat
%
%		meshFZ.m
%			-tetgen
%			neq.m
%
%		meshBP.m
%			findgeometry.m
%
%			symaxis.m
%
%			sphtri_subdiv.m
%				-facet_subdiv.m
%					--simplex_subdiv.m
%
%		neq.m
%--------------------------------------------------------------------------
%% setup
usual = 2; %usual # of variables

if nargin > usual
	% initialize new variables, add as input to var_names
	%----------------------
	res = [];
	nint = [];
	S = var_names(res,nint); %package into struct
	%----------------------
	vars = fields(S); %get fields (i.e. strings of variable names)
	
	% load varargin values into struct
	irange = 1:nargin-usual;
	for i = irange
		var = vars{i};
		S.(var) = varargin{i};
	end
	
	% unpack variables (for brevity)
	unpack_method = 'eval';
	switch unpack_method
		case 'eval'
			for i = 1:length(vars)
				var = vars{i};
				temp = S.(var); %#ok<NASGU> %temporary value of vName
				evalc([var '= temp']); %assign temp value to a short name
			end
			
		case 'manual'
			res = S.(vars{1});
			nint = S.(vars{2});
	end
end

%% mesh 5DOF

%mesh FZ if applicable
switch featureType
	case {'resolution','interior','exterior'}
		%mesh misorientation FZ using a constant angular resolution
		[FZmesh,meshTriangles,meshTetrahedra,idAboveAC,idBelowAC,idUnique] = meshFZ(res);
		rodpts = FZmesh.Points; %rodriguez vectors
end

%partition FZ if applicable
switch featureType
	case {'vertices','vtx_deleteO','misFZfeatures','vtx_deleteOz'}
		load('misFZfeatures.mat','dlist','qnames')
		
	case 'interior' %not sure if this is functioning
		extIDs = FZexteriorIDs(rodpts);
		rodpts(extIDs,:) = [];
		
	case {'exterior','ext_deleteO'}
		intIDs = ~FZexteriorIDs(rodpts);
		rodpts(intIDs,:) = [];
		disp(['FZ mesh points = ' int2str(size(rodpts,1))])
end

switch featureType
	case {'vtx_deleteO'}
		qnames = setdiff(qnames,'O');
	case 'ext_deleteO'
		O = [0,0,0];
		[row,~] = find(ismembertol(rodpts,O,1e-6,'ByRows',true));
		rodpts(row,:) = [];
	case 'misFZfeatures'
		qnames = fields(dlist);
end

switch featureType
	case {'vertices','vtx_deleteO','misFZfeatures','vtx_deleteOz'}
		npts = length(qnames);
		rodpts = zeros(npts,3); %initialize
		for i = 1:npts
			qname = qnames{i};
			rodpts(i,:) = dlist.(qname);
		end
end

%find and add folders to path
fnlist = {'GBfive2oct.m','rod2q.m','ax2qu.m'};
filepathgen = fullfile('**',fnlist);
for i = 1:length(filepathgen)
	octfile = dir(filepathgen{i});
	octfolder = octfile(1).folder;
	addpath(octfolder);
end

if isempty(rodpts)
	error('mesh is empty. Consider decreasing resDegrees')
end

%initialize
npts = size(rodpts,1);
qlist = zeros(npts,4); %quaternions
nApts = cell(npts,1); % BP normal vectors
A = nApts;
R = nApts;
TRI = nApts;
fivelist = nApts;
olist = nApts;

%construct octonions
for i = 1:npts
	d = rodpts(i,:);
	qlist(i,:) = rod2q(d);
	q = qlist(i,:);
	geometry = findgeometry(q);
	
	%get boundary plane mesh
	[nApts{i},A{i},R{i},TRI{i}] = meshBP(q,nint,ctrcuspQ);
	
	if strcmp(featureType,'vtx_deleteOz')
		qname = qnames{i};
		if strcmp(qname,'O')
			%assumes that nA = [0 0 1] is the first point
			nApts{i}.main(1,:) = [];
			nApts{i}.sub(1,:) = [];
		end
	end
	
	nAptsTemp = vertcat(nApts{i}.sub);
	npts2 = size(nAptsTemp,1);
	
	%convert from 5DOF to octonions
	fivelist{i}(npts2) = struct();
	fivelist{i}(1).q = [];
	fivelist{i}(1).nA = [];
	fivelist{i}(1).d = [];
	fivelist{i}(1).geometry = '';
	
	olist{i} = zeros(npts2,8);
	
	for j = 1:npts2
		nA = nAptsTemp(j,:);
		fivelist{i}(j).q = q;
		fivelist{i}(j).nA = nA;
		fivelist{i}(j).d = d;
		fivelist{i}(j).geometry = geometry;
		olist{i}(j,:) = GBfive2oct(q,nA);
	end
end

%catenate 5DOF data and octonions
five = horzcat(fivelist{:});
o = vertcat(olist{:});

if all(o(:,8) == 0)
	sept = o(:,1:7); %septonion?
else
	warning('last column of octonions ~= 0. Setting sept == o')
	sept = o;
end
pts = sept;

disp(['total of ' int2str(size(pts,1)) ' points in 5DOF mesh'])

end %mesh5DOF

%---------------------------HELPER FUNCTIONS-------------------------------
function extIDs = FZexteriorIDs(rodpts)
O = [0, 0, 0];                                 % point O
A = [sqrt(2)-1, 0, 0];                         % point A
B = [sqrt(2)-1, sqrt(2)-1, 0];                 % point B
C = [sqrt(2)-1, sqrt(2)-1, 3-2*sqrt(2)];       % point C
E = [1/3, 1/3, 1/3];                           % point E
D = [sqrt(2)-1, 1-(1/sqrt(2)), 1-(1/sqrt(2))]; % point D

%Plane CED
ED = D - E;ED = ED/norm(ED);
EC = C - E;EC = EC/norm(EC);
nCED = cross(ED,EC);
planeCED = @(x,y)(-(nCED(1)*(x - E(1)) + nCED(2)*(y - E(2)))/nCED(3) + E(3));
isOnCED = neq(rodpts(:,3),planeCED(rodpts(:,1),rodpts(:,2)),6);

%Plane OBA
isOnOBA = neq(rodpts(:,3),0,6);

%Plane ABCD
isOnABCD = neq(rodpts(:,1),sqrt(2)-1,6);

%Plane OECB (vertical, so use an equation for x)
OB = B - O;OB = OB/norm(OB);
OE = E - O;OE = OE/norm(OE);
nOECB = cross(OB,OE);
planeOECB = @(y,z)(-(nOECB(3)*(z - O(3)) + nOECB(2)*(y - O(2)))/nOECB(1) + O(1));
isOnOECB = neq(rodpts(:,1),planeOECB(rodpts(:,2),rodpts(:,3)),6);

%Plane OADE
OA = A - O;OA = OA/norm(OA);
OE = E - O;OE = OE/norm(OE);
nOADE = cross(OA,OE);
planeOADE = @(x,y)(-(nOADE(1)*(x - O(1)) + nOADE(2)*(y - O(2)))/nOADE(3) + O(3));
isOnOADE = neq(rodpts(:,3),planeOADE(rodpts(:,1),rodpts(:,2)),6);

extIDs = isOnCED | isOnOBA | isOnABCD | isOnOECB | isOnOADE;
end

%--------------------------------CODE GRAVEYARD----------------------------
%{
		load('misFZfeatures.mat','dlist','qnames')
		npts = length(qnames);
		rodptsDel = zeros(npts,3); %initialize, rodpts to delete
		for i = 1:npts
			qname = qnames{i};
			rodptsDel(i,:) = dlist.(qname);
		end
		keepIDs = find(~ismembertol(rodpts,rodptsDel,'ByRows',true));
		rodpts = rodpts(keepIDs,:);

%normalize octonions
o = normr(o);
%}
