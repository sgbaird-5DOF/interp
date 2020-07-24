function [FZmesh,meshTriangles,meshTetrahedra,idAboveAC,idBelowAC,idUnique] = meshFZ(resDegrees)
%--------------------------------------------------------------------------
% Author: Oliver Johnson
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
%		tetgen
%		
%--------------------------------------------------------------------------
% define resolution in rodriguez space
res = tan(deg2rad(resDegrees)/2);

% get current directory
my_dir = pwd;

% change to tetgen directory
file = dir(fullfile('**','tetgen.exe'));
tetgen_dir = fullfile(file(1).folder);
%tetgen_dir = 'tetgen1.5.1';
cd(tetgen_dir);

% % change to tetgen directory
% tetfile = dir('**/tetgen.cxx');
% tetfolder = tetfile.folder;
% 
% if isempty(tetfile)
% 	disp('need to download & install tetgen or change tetgen location or searchpath')
% else
% 	tetfile2 = dir('**/tetgen.exe');
% 	tetfolder = tetfile.folder;
% 	if isempty(tetfile2)
% 		disp('need to install tetgen - see instructions online')
% 	else
% 		addpath(tetfolder);
% 		cd(tetfolder);
% 	end
% end

%% Write FZ.poly

% define FZ borders
FZpoly.Vertices = [0, 0, 0;...                  % point O (1)
    sqrt(2)-1, 0, 0;...                         % point A (2)
    sqrt(2)-1, sqrt(2)-1, 0;...                 % point B (3)
    sqrt(2)-1, sqrt(2)-1, 3-2*sqrt(2);...       % point C (4)
    1/3, 1/3, 1/3;...                           % point E (5)
    sqrt(2)-1, 1-(1/sqrt(2)), 1-(1/sqrt(2))];   % point D (6)
FZpoly.Faces = [1 3 2 nan;... % OBA  (1)
    2 3 4 nan;...             % ABC  (2) 
    2 4 6 nan;...             % ACD  (3)
    1 5 4 3;...               % OECB (4)
    1 2 6 5;...               % OADE (5)
    4 5 6 nan];               % CED  (6)
FZpoly.Edges = [1 2;... % OA (1)
                1 3;... % OB (2)
                1 5;... % OE (3)
                2 3;... % AB (4)
                2 4;... % AC (5)
                2 6;... % AD (6)
                3 4;... % BC (7)
                4 5;... % CE (8)
                4 6;... % CD (9)
                5 6];   % ED (10)

% calculate length of each line
for i = 1:size(FZpoly.Edges,1)
    FZpoly.EdgeLength(i) = norm(FZpoly.Vertices(FZpoly.Edges(i,2),:)-FZpoly.Vertices(FZpoly.Edges(i,1),:));
end

% assign vertices
FZmesh.Points = FZpoly.Vertices;

% assign points along edges
for i = 1:size(FZpoly.Edges,1)
    t = linspace(0,1,ceil(FZpoly.EdgeLength(i)/res)+1).';
    t = t(2:(end-1));
    FZmesh.Edges{i} = [FZpoly.Edges(i,1), size(FZmesh.Points,1)+(1:numel(t)), FZpoly.Edges(i,2)];
    FZmesh.Points = [FZmesh.Points;...
        (1-t).*FZpoly.Vertices(FZpoly.Edges(i,1),:)+t*FZpoly.Vertices(FZpoly.Edges(i,2),:)];
end
FZmesh.Edges = FZmesh.Edges.';

% define faces
FZmesh.Faces = {[FZmesh.Edges{2},fliplr(FZmesh.Edges{4}(1:end-1)),fliplr(FZmesh.Edges{1}(2:end-1))];...
    [FZmesh.Edges{4},FZmesh.Edges{7}(2:end-1),fliplr(FZmesh.Edges{5}(2:end))];...
    [FZmesh.Edges{5},FZmesh.Edges{9}(2:end-1),fliplr(FZmesh.Edges{6}(2:end))];...
    [FZmesh.Edges{3},fliplr(FZmesh.Edges{8}(2:end-1)),fliplr(FZmesh.Edges{7}(2:end)),fliplr(FZmesh.Edges{2}(2:end))];...
    [FZmesh.Edges{1},FZmesh.Edges{6}(2:end-1),fliplr(FZmesh.Edges{10}(2:end)),fliplr(FZmesh.Edges{3}(2:end))];...
    [FZmesh.Edges{8},FZmesh.Edges{10}(2:end-1),fliplr(FZmesh.Edges{9}(2:end))]};

% write FZ.poly file
fid = fopen('FZ.poly','w');
fprintf(fid,'# vertices\r\n');
fprintf(fid,'# Part 1 - node list\r\n');
fprintf(fid,'%u  %u  0  0\r\n',size(FZpoly.Vertices));
fprintf(fid,'# Node index, node coordinates\r\n');
fprintf(fid,'%u  %.15f %.15f %.15f\r\n',[1:size(FZpoly.Vertices,1); FZpoly.Vertices.']);
fprintf(fid,'\r\n');
fprintf(fid,'# Part 2 - facet list\r\n');
fprintf(fid,'# facet count, no boundary marker\r\n');
fprintf(fid,'%u 0\r\n',size(FZpoly.Faces,1));
for i = 1:size(FZpoly.Faces,1)
    fprintf(fid,'1\r\n');
    str = ['%u  ',repmat('%u ',1,numel(FZpoly.Faces(i,~isnan(FZpoly.Faces(i,:)))))];
    str = [str(1:end-1),'\r\n'];
    num = [numel(FZpoly.Faces(i,~isnan(FZpoly.Faces(i,:)))), FZpoly.Faces(i,~isnan(FZpoly.Faces(i,:)))];
    fprintf(fid,str,num);
end
fprintf(fid,'\r\n');
fprintf(fid,'# Part 3 - hole list\r\n');
fprintf(fid,'0\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'# Part 4 - region list\r\n');
fprintf(fid,'0\r\n');
fclose(fid);

%% Write mesh sizing function

fid = fopen('FZ.mtr','w');
fprintf(fid,'%u 1\r\n',size(FZpoly.Vertices,1));
fprintf(fid,'%.15f\r\n',res*ones(1,size(FZpoly.Vertices,1)));
fclose(fid);

%% Write input nodes file

% write FZ.a.node file
fid = fopen('FZ.a.node','w');
fprintf(fid,'# vertices\r\n');
fprintf(fid,'# Part 1 - node list\r\n');
fprintf(fid,'%u  %u  0  0\r\n',size(FZmesh.Points(7:end,:)));
fprintf(fid,'# Node index, node coordinates\r\n');
fprintf(fid,'%u  %.15f %.15f %.15f\r\n',[1:size(FZmesh.Points(7:end,:),1); FZmesh.Points(7:end,:).']);
fprintf(fid,'\r\n');
fclose(fid);

%% Generate initial mesh

% options: 
% V (verbose)
% p (mesh the polyhedron in FZ.poly)  
% i (use the additional points specified in the FZ.a.node file)
% m (use the mesh sizing function in FZ.mtr)
% q1.1 (quality mesh with radius-edge ratio 1.1)

[status,result] = system('./tetgen -Vpimq1.1 FZ.poly');

%% Read in the surface mesh nodes

FZmesh.Points = readNODE('FZ.1.node');

%% Mirror Points on ACD to ACB

% identify points on ABCD
isOnABCD = neq(FZmesh.Points(:,1),sqrt(2)-1,6);

% identify points on ACD
isAboveAC = isOnABCD & ngt(FZmesh.Points(:,3)./FZmesh.Points(:,2),FZpoly.Vertices(4,3)./FZpoly.Vertices(4,2),6);
idAboveAC = find(isAboveAC);

% identify points on ACB
isBelowAC = isOnABCD & nlt(FZmesh.Points(:,3)./FZmesh.Points(:,2),FZpoly.Vertices(4,3)./FZpoly.Vertices(4,2),6);

% find barycentric coordinates of points on ACD
A = [FZpoly.Vertices(FZpoly.Faces(3,1:3),2).';...
    FZpoly.Vertices(FZpoly.Faces(3,1:3),3).';...
    ones(1,3)];
b = [FZmesh.Points(isAboveAC,2).'; FZmesh.Points(isAboveAC,3).'; ones(1,sum(isAboveAC))];
lambda = (A\b).';

% calculate new coordinates of points on ACB
newPoints = lambda*FZpoly.Vertices(FZpoly.Faces(2,[1 3 2]),:);
idBelowAC = size(FZmesh.Points,1)+(1:size(newPoints,1));

% remove old points from ACB
FZmesh.Points(isBelowAC,:) = [];

% add new points to ACB
FZmesh.Points = [FZmesh.Points; newPoints];

% update indices
idOld = (1:max(idBelowAC)).';
idOld([isBelowAC; false(size(newPoints,1),1)]) = [];
idNew = nan(max(idBelowAC),1);
idNew(idOld) = (1:size(FZmesh.Points,1)).';

idAboveAC = idNew(idAboveAC);
idBelowAC = idNew(idBelowAC);
idUnique = setdiff((1:size(FZmesh.Points,1)).',idBelowAC);

%% Recalculate delaunay tetrahedralization (removing any flat tetrahedra)

% change back to original directory
cd(my_dir)

% calculate delaunay tetrahedralization and remove flat tetrahedra
FZmesh = removeFZFaceTetra(FZmesh);

% extract tetrahedra and triangles
meshTetrahedra = FZmesh.ConnectivityList;
meshTriangles = zeros(0,3);
id = nchoosek(1:4,3);
for i = 1:size(id,1)
    meshTriangles = [meshTriangles; meshTetrahedra(:,id(i,:))];
end
meshTriangles = unique(sort(meshTriangles,2),'rows');

%% Make sure surface points are on the surface

warnID = 'MATLAB:structOnObject';
warning('off',warnID)
FZmesh = struct(FZmesh);
warning('on',warnID)
FZmesh.Points = moveIntoFZ(FZmesh.Points,0);
FZmesh = triangulation(FZmesh.ConnectivityList,FZmesh.Points);

end

%------------------------HELPER FUNCTIONS---------------------------------%

function DT = removeFZFaceTetra(FZmesh)
%--------------------------------------------------------------------------
%Author: Brandon Snow
%function DT = removeFZFaceTetra(FZmesh)
%finds points in the FZ mesh that are outside of the FZ and moves them
%back. Also eliminates any tetrahedra that are on the faces of the FZ

% % define FZ borders
% FZvertices = [0, 0, 0;...                       % point O
%     sqrt(2)-1, 0, 0;...                         % point A
%     sqrt(2)-1, sqrt(2)-1, 0;...                 % point B
%     sqrt(2)-1, sqrt(2)-1, 3-2*sqrt(2);...       % point C
%     1/3, 1/3, 1/3;...                           % point E
%     sqrt(2)-1, 1-(1/sqrt(2)), 1-(1/sqrt(2))];   % point D
% FZfaces = [1 3 2 nan;...    % OBA
%     2 3 4 6;...             % ABCD
%     1 5 4 3;...             % OECB
%     1 2 6 5;...             % OADE
%     4 5 6 nan];             % CED

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
isOnCED = neq(FZmesh.Points(:,3),planeCED(FZmesh.Points(:,1),FZmesh.Points(:,2)),6);

%Plane OBA
isOnOBA = neq(FZmesh.Points(:,3),0,6);

%Plane ABCD
isOnABCD = neq(FZmesh.Points(:,1),sqrt(2)-1,6);

%Plane OECB (vertical, so use an equation for x)
OB = B - O;OB = OB/norm(OB);
OE = E - O;OE = OE/norm(OE);
nOECB = cross(OB,OE);
planeOECB = @(y,z)(-(nOECB(3)*(z - O(3)) + nOECB(2)*(y - O(2)))/nOECB(1) + O(1));
isOnOECB = neq(FZmesh.Points(:,1),planeOECB(FZmesh.Points(:,2),FZmesh.Points(:,3)),6);

%Plane OADE
OA = A - O;OA = OA/norm(OA);
OE = E - O;OE = OE/norm(OE);
nOADE = cross(OA,OE);
planeOADE = @(x,y)(-(nOADE(1)*(x - O(1)) + nOADE(2)*(y - O(2)))/nOADE(3) + O(3));
isOnOADE = neq(FZmesh.Points(:,3),planeOADE(FZmesh.Points(:,1),FZmesh.Points(:,2)),6);

DT = delaunayTriangulation(FZmesh.Points);
DT = struct('Points',DT.Points,'ConnectivityList',DT.ConnectivityList);

%Find Tetrahedra on Planes
CEDtetra = sum(ismember(DT.ConnectivityList,find(isOnCED)),2)==4;
OBAtetra = sum(ismember(DT.ConnectivityList,find(isOnOBA)),2)==4;
ABCDtetra = sum(ismember(DT.ConnectivityList,find(isOnABCD)),2)==4;
OECBtetra = sum(ismember(DT.ConnectivityList,find(isOnOECB)),2)==4;
OADEtetra = sum(ismember(DT.ConnectivityList,find(isOnOADE)),2)==4;
planeTetra = find(CEDtetra+OBAtetra+ABCDtetra+OECBtetra+OADEtetra);

%Remove Plane Tetra
DT.ConnectivityList(planeTetra,:) = [];
DT = triangulation(DT.ConnectivityList,DT.Points);

%Report Removal
reportremovalQ = false;
if reportremovalQ
	a = 'Removed ';b = ' Tetrahedra from plane ';
	disp([a num2str(sum(CEDtetra)) b 'CED'])
	disp([a num2str(sum(OBAtetra)) b 'OBA'])
	disp([a num2str(sum(ABCDtetra)) b 'ABCD'])
	disp([a num2str(sum(OECBtetra)) b 'OECB'])
	disp([a num2str(sum(OADEtetra)) b 'OADE'])
end
end
%-----------------------END removeFZFaceTetra------------------------------

function pts = moveIntoFZ(pts,N,neqVar)
%--------------------------------------------------------------------------
%Move points within 10^-neqVar of the FZ faces to N*eps inside the FZ
O = [0, 0, 0];                                 % point O
A = [sqrt(2)-1, 0, 0];                         % point A
B = [sqrt(2)-1, sqrt(2)-1, 0];                 % point B
C = [sqrt(2)-1, sqrt(2)-1, 3-2*sqrt(2)];       % point C
E = [1/3, 1/3, 1/3];                           % point E
D = [sqrt(2)-1, 1-(1/sqrt(2)), 1-(1/sqrt(2))]; % point D

if nargin < 3
    neqVar = 5;
end

%Plane CED
ED = D - E;
ED = ED/norm(ED);
EC = C - E;
EC = EC/norm(EC);
nCED = cross(ED,EC);
planeCED = @(x,y)(-(nCED(1)*(x - E(1)) + nCED(2)*(y - E(2)))/nCED(3) + E(3));

%Plane OBA
planeOBA = 0;

%Plane ABCD
planeABCD = sqrt(2)-1;

%Plane OECB (vertical, so use an equation for y)
planeOECB = @(x,z)(x);

%Plane OADE
OA = A - O;
OA = OA/norm(OA);
OE = E - O;
OE = OE/norm(OE);
nOADE = cross(OA,OE);
planeOADE = @(x,y)(-(nOADE(1)*(x - O(1)) + nOADE(2)*(y - O(2)))/nOADE(3) + O(3));

%find points on the planes
isOnCED = neq(pts(:,3),planeCED(pts(:,1),pts(:,2)),neqVar);
isOnOBA = neq(pts(:,3),planeOBA,neqVar);
isOnABCD = neq(pts(:,1),sqrt(2)-1,neqVar);
isOnOECB = neq(pts(:,2),planeOECB(pts(:,1),pts(:,3)),neqVar);
isOnOADE = neq(pts(:,3),planeOADE(pts(:,1),pts(:,2)),neqVar);

%Check Edges
isOnOA = and(isOnOADE,isOnOBA);
isOnOB = and(isOnOECB,isOnOBA);
isOnAB = and(isOnABCD,isOnOBA);
isOnOE = and(isOnOADE,isOnOECB);
isOnAD = and(isOnOADE,isOnABCD);
isOnDE = and(isOnCED,isOnOADE);
isOnCD = and(isOnABCD,isOnCED);
isOnCE = and(isOnOECB,isOnCED);
isOnBC = and(isOnABCD,isOnOECB);

%Define lines to return the other two dimensions
lineOA = @(x)([zeros(size(x)),zeros(size(x))]);%y,z
lineOB = @(x)([x,zeros(size(x))]);
lineAB = @(y)([(sqrt(2)-1)*ones(size(y)),zeros(size(y))]);%x,z
lineOE = @(x)([(E(2)/E(1))*x,(E(3)/E(1))*x]);%y,z
lineAD = @(y)([(sqrt(2)-1)*ones(size(y)),A(3) + y*((D(3)-A(3))/(D(2)-A(2)))]);
lineDE = @(x)([D(2) + (x-D(1))*((E(2)-D(2))/(E(1)-D(1))),D(3) + (x-D(1))*((E(3)-D(3))/(E(1)-D(1)))]);
lineCD = @(y)([(sqrt(2)-1)*ones(size(y)),C(3) + (y-C(2))*((D(3)-C(3))/(D(2)-C(2)))]);
lineCE = @(x)([C(2) + (x-C(1))*((E(2)-C(2))/(E(1)-C(1))),C(3) + (x-C(1))*((E(3)-C(3))/(E(1)-C(1)))]);
lineBC = @(z)([(sqrt(2)-1)*ones(size(z)),(sqrt(2)-1)*ones(size(z))]);

%Check Vertices
isO = all(neq(pts,O,neqVar),2);
isA = all(neq(pts,A,neqVar),2);
isB = all(neq(pts,B,neqVar),2);
isC = all(neq(pts,C,neqVar),2);
isD = all(neq(pts,D,neqVar),2);
isE = all(neq(pts,E,neqVar),2);

%Remove Vertices from Edges
isOnOA = isOnOA - isO - isA > 0;
isOnOB = isOnOB - isO - isB > 0;
isOnAB = isOnAB - isA - isB > 0;
isOnOE = isOnOE - isO - isE > 0;
isOnAD = isOnAD - isA - isD > 0;
isOnDE = isOnDE - isD - isE > 0;
isOnCD = isOnCD - isC - isD > 0;
isOnCE = isOnCE - isC - isE > 0;
isOnBC = isOnBC - isB - isC > 0;

%Remove Edges from Planes
isOnCED = isOnCED - isOnCE - isOnDE - isOnCD > 0;
isOnOBA = isOnOBA - isOnOB - isOnOA - isOnAB > 0;
isOnABCD = isOnABCD - isOnAB - isOnBC - isOnCD - isOnAD > 0;
isOnOECB = isOnOECB - isOnOE - isOnCE - isOnBC - isOnOB > 0;
isOnOADE = isOnOADE - isOnOA - isOnAD - isOnDE - isOnOE > 0;

%Move Points
if nargin < 2
    N = 0;
end
%Move Plane Points
pts(isOnOBA,3) = planeOBA + N*eps;
pts(isOnCED,3) = planeCED(pts(isOnCED,1),pts(isOnCED,2)) - N*eps;
pts(isOnABCD,1) = planeABCD - N*eps;
pts(isOnOECB,2) = planeOECB(pts(isOnOECB,1),pts(isOnOECB,3)) - N*eps;
pts(isOnOADE,3) = planeOADE(pts(isOnOADE,1),pts(isOnOADE,2)) - N*eps;

%Edges
%Move Close points onto the Edges
pts(isOnOA,2:3) = lineOA(pts(isOnOA,1));
pts(isOnOB,2:3) = lineOB(pts(isOnOB,1));
pts(isOnAB,[1,3]) = lineAB(pts(isOnAB,2));
pts(isOnOE,2:3) = lineOE(pts(isOnOE,1));
pts(isOnAD,[1,3]) = lineAD(pts(isOnAD,2));
pts(isOnDE,2:3) = lineDE(pts(isOnDE,1));
pts(isOnCD,[1,3]) = lineCD(pts(isOnCD,2));
pts(isOnCE,2:3) = lineCE(pts(isOnCE,1));
pts(isOnBC,1:2) = lineBC(pts(isOnBC,3));

%Vertices
if any(isO)
    pts(isO,:) = O;
end
if any(isA)
    pts(isA,:) = A;
end
if any(isB)
    pts(isB,:) = B;
end
if any(isC)
    pts(isC,:) = C;
end
if any(isD)
    pts(isD,:) = D;
end
if any(isE)
    pts(isE,:) = E;
end

end

%----------------------END HELPER FUNCTIONS-------------------------------%

%% OLD STUFF
% %%
% 
% % define FZ borders
% FZpoly.Vertices = [0, 0, 0;...                  % point O
%     sqrt(2)-1, 0, 0;...                         % point A
%     sqrt(2)-1, sqrt(2)-1, 0;...                 % point B
%     sqrt(2)-1, sqrt(2)-1, 3-2*sqrt(2);...       % point C
%     1/3, 1/3, 1/3;...                           % point E
%     sqrt(2)-1, 1-(1/sqrt(2)), 1-(1/sqrt(2))];   % point D
% FZpoly.Faces = [1 3 2;... % OBA (1)
%     2 3 4;...             % ABC (2)
%     2 4 6;...             % ACD (3)
%     1 5 4;...             % OEC (4)
%     1 4 3;...             % OCB (5)
%     1 2 6;...             % OAD (6)
%     1 6 5;...             % ODE (7)
%     4 5 6];               % CED (8)
% FZpoly.Edges = unique(sort([FZpoly.Faces(:,1:2);...
%     FZpoly.Faces(:,2:3);...
%     FZpoly.Faces(:,[3,1])],2),'rows'); % AC -> FZpoly.Edges(7,:)
% 
% % mesh points
% FZmesh.Points = FZmesh.Vertices;
% 
% % define incidence matrices
% FZmesh.PointFace = false(size(FZmesh.Points,1),size(FZmesh.Faces,1));
% FZmesh.PointFace(sub2ind(size(FZmesh.PointFace),FZmesh.Faces(:),reshape(repmat((1:size(FZmesh.Faces,1)).',1,3),[],1))) = true;
% 
% FZmesh.PointEdge = false(size(FZmesh.Points,1),size(FZmesh.Edges,1));
% FZmesh.PointEdge(sub2ind(size(FZmesh.PointEdge),FZmesh.Edges(:),reshape(repmat((1:size(FZmesh.Edges,1)).',1,2),[],1))) = true;
% 
% FZmesh.EdgeFace = false(size(FZmesh.Edges,1),size(FZmesh.Faces,1));
% for i = 1:size(FZmesh.Edges,1)
%     FZmesh.EdgeFace(i,:) = (sum(ismember(FZmesh.Faces,FZmesh.Edges(i,:)),2) == 2).';
% end
% 
% % calculate length of each line
% for i = 1:size(FZmesh.Edges,1)
%     FZmesh.EdgeLength(i) = norm(FZmesh.Vertices(FZmesh.Edges(i,2),:)-FZmesh.Vertices(FZmesh.Edges(i,1),:));
% end
% 
% % assign points along edges
% for i = 1:size(FZmesh.Edges,1)
%     t = linspace(0,1,ceil(FZmesh.EdgeLength(i)/res)+1).';
%     t = t(2:(end-1));
%     FZmesh.Points = [FZmesh.Points;...
%         t.*FZmesh.Vertices(FZmesh.Edges(i,1),:)+(1-t)*FZmesh.Vertices(FZmesh.Edges(i,2),:)];
%     FZmesh.PointFace = [FZmesh.PointFace;...
%         repmat(all(FZmesh.PointFace(FZmesh.Edges(i,:),:),1),numel(t),1)];
%     FZmesh.PointEdge = [FZmesh.PointEdge;...
%         repmat(all(FZmesh.PointEdge(FZmesh.Edges(i,:),:),1),numel(t),1)];
% end
% 
% % calculate area of each triangular face
% for i = 1:size(FZmesh.Faces,1)
%     a = FZmesh.Points(FZmesh.Faces(i,1),:).';
%     b = FZmesh.Points(FZmesh.Faces(i,2),:).';
%     c = FZmesh.Points(FZmesh.Faces(i,3),:).';
%     FZmesh.FaceArea(i) = (1/2)*norm(cross(b-a,c-a));
% end
% 
% % assign points on faces
% targetTriangleArea = (sqrt(3)/4)*(res^2);
% for i = 1:size(FZmesh.Faces,1)
%     if i == 3 % we are doing ACD, so reuse same barycentric coordinates from ABC, but permuted
%         pts = simplexCoordinates*FZmesh.Points(FZmesh.Faces(i,[1 3 2]),:);
%     else
%         nSamples = max(1,ceil(((FZmesh.FaceArea(i)/targetTriangleArea)-sum(FZmesh.PointFace(:,i))+2)/2)); % add at least 1 point
%         simplexCoordinates = simplexsamples(nSamples,3);
%         pts = simplexCoordinates*FZmesh.Points(FZmesh.Faces(i,:),:);
%     end
%     FZmesh.Points = [FZmesh.Points; pts];
%     FZmesh.PointFace = [FZmesh.PointFace;...
%         repmat(all(FZmesh.PointFace(FZmesh.Faces(i,:),:),1),nSamples,1)];
%     FZmesh.PointEdge = [FZmesh.PointEdge;... 
%         false(nSamples,size(FZmesh.PointEdge,2))]; % they can't be on an edge since they are on face interiors
% end
% 
% % equilibrate points on faces
% nIter = 100;
% for i = 1:size(FZmesh.Faces,1)
%     disp(['Meshing Face: ',num2str(i)])
%     % get 2D coordinates (rotate if necessary)
%     switch i
%         case 1 % OBA
%             x = FZmesh.Points(FZmesh.PointFace(:,i),1);
%             y = FZmesh.Points(FZmesh.PointFace(:,i),2);
%         case 2 % ABC
%             x = FZmesh.Points(FZmesh.PointFace(:,i),2);
%             y = FZmesh.Points(FZmesh.PointFace(:,i),3);
%         case 3 % ACD
%             x = FZmesh.Points(FZmesh.PointFace(:,i),2);
%             y = FZmesh.Points(FZmesh.PointFace(:,i),3);
%         case 4 % OEC
%             v = (q2gmat(rot2q(deg2rad(45),0,0))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the y-z plane for 2D delaunay triangulation
%             x = v(:,2);
%             y = v(:,3);
%         case 5 % OCB
%             v = (q2gmat(rot2q(deg2rad(45),0,0))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the y-z plane for 2D delaunay triangulation
%             x = v(:,2);
%             y = v(:,3);
%         case 6 % OAD
%             v = (q2gmat(rot2q(deg2rad(-45),deg2rad(90),0))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the x-y plane for 2D delaunay triangulation
%             x = v(:,1);
%             y = v(:,2);
%         case 7 % ODE
%             v = (q2gmat(rot2q(deg2rad(-45),deg2rad(90),0))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the x-y plane for 2D delaunay triangulation
%             x = v(:,1);
%             y = v(:,2);
%         case 8 % CED
%             v = (q2gmat(rot2q(-acos(1/sqrt(3)),deg2rad(90),deg2rad(135)))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the x-y plane for 2D delaunay triangulation
%             x = v(:,1);
%             y = v(:,2);
%     end
%     
%     % get indices
%     id = find(FZmesh.PointFace(:,i)); % indices of these points relative to the full set
%     idMoveLocal = find(~any(FZmesh.PointEdge(FZmesh.PointFace(:,i),:),2)); % indices to move relative to points on this face
%     idMove = id(idMoveLocal); % indices to move relative to the full set
%     vertexIDLocal = find(sum(FZmesh.PointEdge(FZmesh.PointFace(:,i),:),2) > 1); % indices of vertices of this face relative to points on this face
%     vertexID = id(vertexIDLocal); % indices of vertices of this face relative to full set
%     
%     % identify boundary points of this face
%     isbnd = any(FZmesh.PointEdge(FZmesh.PointFace(:,i),:),2);
%     
%     for iter = 1:nIter % lloyds algorithm
%         
%         % get voronoi cells
%         [~,vorvx] = polybnd_voronoi([x,y],[x(isbnd),y(isbnd)]);
%             
%         % move points to centroid of neighbors (if they are interior to face)
%         for j = 1:numel(idMove)
%             
%             % compute centroid
%             xi = vorvx{idMoveLocal(j)}(:,1);
%             xi1 = circshift(vorvx{idMoveLocal(j)}(:,1),-1);
%             yi = vorvx{idMoveLocal(j)}(:,2);
%             yi1 = circshift(vorvx{idMoveLocal(j)}(:,2),-1);
%             A = 0.5*sum((xi1+xi).*(yi1-yi));
%             cx = (1/(6*A))*sum((xi+xi1).*(xi.*yi1-xi1.*yi));
%             cy = (1/(6*A))*sum((yi+yi1).*(xi.*yi1-xi1.*yi));
% 
%             x(idMoveLocal(j)) = cx;
%             y(idMoveLocal(j)) = cy;
% 
%         end
%     end
%     
%     % convert plane coordinates to simplex weights
%     A = [x(vertexIDLocal).'; y(vertexIDLocal).'; ones(1,3)];
%     b = [x(idMoveLocal).'; y(idMoveLocal).'; ones(1,numel(idMoveLocal))];
%     lambda = (A\b).';
%     
%     % convert simplex weights to 3D coordinates
%     FZmesh.Points(idMove,:) = lambda*FZmesh.Points(vertexID,:);
% end
% 
% figure; 
% patch('vertices',FZmesh.Points(1:6,:),'faces',FZmesh.Faces,'facecolor','none','edgecolor','k','marker','.'); 
% view(3); 
% axis equal tight vis3d; 
% hold on;
% plot3(FZmesh.Points(:,1),FZmesh.Points(:,2),FZmesh.Points(:,3),'r.');
% 
% % assign points in interior
% aFZ = sum(FZmesh.FaceArea);
% [~,vFZ] = convhull(FZmesh.Vertices);
% nSamples = ceil((vFZ/aFZ)*(1/4)-size(FZmesh.Points,1)); % heuristic
% nSamples = 20:10:1000;
% for i = 1:numel(nSamples)
%     workbar(i/numel(nSamples));
% [A,b] = vert2lcon(FZmesh.Vertices); %get linear constraints representation of polygon
% temp_samples = cprnd(nSamples(i),A,b);
% pts = [FZmesh.Points;temp_samples];
% Edges = edges(delaunayTriangulation(pts));
% isInteriorEdge = ~all(ismember(Edges,1:size(FZmesh.Points,1)),2);
% EdgeLength = sqrt(sum((pts(Edges(:,2),:)-pts(Edges(:,1),:)).^2,2));
% meanEdge(i) = mean(EdgeLength(isInteriorEdge));
% medianEdge(i) = median(EdgeLength(isInteriorEdge));
% minEdge(i) = min(EdgeLength(isInteriorEdge));
% maxEdge(i) = max(EdgeLength(isInteriorEdge));
% end
% 
% % equilibrate points in interior
% 
% 
% end
% 
% % 
% % Edges = zeros(0,2);
% % for i = 1:size(FZmesh.Faces,1)
% %     % get 2D coordinates (rotate if necessary)
% %     switch i
% %         case 1 % OBA
% %             x = FZmesh.Points(FZmesh.PointFace(:,i),1);
% %             y = FZmesh.Points(FZmesh.PointFace(:,i),2);
% %         case 2 % ABC
% %             x = FZmesh.Points(FZmesh.PointFace(:,i),2);
% %             y = FZmesh.Points(FZmesh.PointFace(:,i),3);
% %         case 3 % ACD
% %             x = FZmesh.Points(FZmesh.PointFace(:,i),2);
% %             y = FZmesh.Points(FZmesh.PointFace(:,i),3);
% %         case 4 % OEC
% %             v = (q2gmat(rot2q(deg2rad(45),0,0))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the y-z plane for 2D delaunay triangulation
% %             x = v(:,2);
% %             y = v(:,3);
% %         case 5 % OCB
% %             v = (q2gmat(rot2q(deg2rad(45),0,0))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the y-z plane for 2D delaunay triangulation
% %             x = v(:,2);
% %             y = v(:,3);
% %         case 6 % OAD
% %             v = (q2gmat(rot2q(deg2rad(-45),deg2rad(90),0))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the x-y plane for 2D delaunay triangulation
% %             x = v(:,1);
% %             y = v(:,2);
% %         case 7 % ODE
% %             v = (q2gmat(rot2q(deg2rad(-45),deg2rad(90),0))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the x-y plane for 2D delaunay triangulation
% %             x = v(:,1);
% %             y = v(:,2);
% %         case 8 % CED
% %             v = (q2gmat(rot2q(-acos(1/sqrt(3)),deg2rad(90),deg2rad(135)))*(FZmesh.Points(FZmesh.PointFace(:,i),:)).').'; % rotate points to the x-y plane for 2D delaunay triangulation
% %             x = v(:,1);
% %             y = v(:,2);
% %     end
% %     
% %     % get delaunay
% %     DT = delaunayTriangulation([x,y]);
% %     
% %     % get edges
% %     e = edges(DT);
% %     
% %     % map back to original indices
% %     id = find(FZmesh.PointFace(:,i)); % indices of these points relative to the full set
% %     Edges = [Edges; id(e)]; % indices to move relative to the full set
% % 
% % end
% % Edges = unique(sort(Edges,2),'rows');
% % EdgeLength = sqrt(sum((FZmesh.Points(Edges(:,2),:)-FZmesh.Points(Edges(:,1),:)).^2,2));
