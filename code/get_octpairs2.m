function [octvtx,usv,five,omega23_pair,omega23_GBdist] = get_octpairs2(pts,five,savename,o2addQ,varargin)
%  GET_OCTPAIRS2  deprecated version of get_octpairs()
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description:
%
% Inputs:
%
%		o2addQ ===	logical, whether to add o2 (in this case, point 'O' with
%						nA = [0 0 1]) to "five", and if not, then get rid of
%						pts(1,:) which corresponds to the same.
%
% Outputs:
%
% Dependencies:
%		sphconvhulln.m
%
%		5DOF_vtx_deleteOz_octsubdiv1.mat
%
%		misFZfeatures.mat
%
%		GBdist.m
%
%		GBdist2.m
%--------------------------------------------------------------------------
usual = 4;
if nargin - usual == 1
	plotQ = varargin{1};
else
	plotQ = false;
end

fnames = {'PGnames.mat','olist.mat','misFZfeatures.mat'};
for i = 1:length(fnames)
	fname = fnames{i};
	fpath = fullfile('**',fname);
	file = dir(fpath);
	if ~isempty(file)
		addpath(file.folder);
	end
end

load('misFZfeatures.mat','qlist','dlist')

%% correct octonions if necessary
if size(pts,2) == 7
	pts = [pts zeros(size(pts,1),1)]; % add column of zeros
end

% if norm(o) == 1 within tolerance, multiply by sqrt(2)
if abs(norm(pts(1,:)) - 1) < 1e-6
	pts = pts*sqrt(2);
elseif abs(norm(pts(1,:)) - sqrt(2)) > 1e-6
	error('norm of octonions ~= 1 || sqrt(2)')
end

%% get two reference octonions (o1 and o2)

%unpack interior point
name1 = 'interior';
qA = qlist.(name1);
% qA = normr(qA+0.1*rand(1,4));

%unpack no boundary point
name2 = 'O';
qB = normr(qlist.(name2));
% qB = normr(qB+0.1*rand(1,4));

%load normals (both should just be [0 0 1])
[~,RA] = symaxis(qA,name1);
nA = (RA*[0 0 1].').';
% nA = normr(nA+0.1*rand(1,3));

[~,RB] = symaxis(qB,name2);
nB = normr((RB*[0 0 1].').');
% nB = normr(nB+0.1*rand(1,3));

%package some "five" output for saving
fiveref1.q =qA;
fiveref1.nA = nA;
fiveref1.d = q2rod(qA);
fiveref1.geometry = name1;

%convert to octonions
o1 = GBfive2oct(qA,nA); %input
o2 = GBfive2oct(qB,nB);
% o2 = [-1 0 0 0 1 0 0 0]; %input

[omega0,oct_sym0,zeta0] = GBdist2([o1 o2],32,false);

%take the symmetrized versions for comparison
% o1 = oct_sym0(1,1:8);
o2 = oct_sym0(1,9:16);

octvtx(1,:) = o2;

%% get minimized distance octonions relative to oct pairs
% o3 = pts(1,:);
% [octvtx(1,:),~,omega23_pair(1),omega23_GBdist(1)] = GBpair(o1,o2,o3);

%loop through pairs relative to interior point. Each pair contains (+z) origin point
npts = size(pts,1);
for i = 1:npts
	%unpack other octonion in pair
	%(o2 and o3 form a pair, each is compared to o1)
	o3 = pts(i,:); %input
	[~,octvtx(i+1,:),omega23_pair(i+1),omega23_GBdist(i+1)] = GBpair(o1,o2,o3);
end

five = GBoct2five(octvtx);

if ~o2addQ
	five(1) = [];
	octvtx(1,:) = [];
end

if plotQ
	figure
	plotFZrodriguez_vtx();
	hold on
	disQ = true;
	if disQ
		t = num2cell(q2rod(disorientation(vertcat(five.q),'cubic')),1);
	else
		t=num2cell(vertcat(five.d),1);
	end
	plot3(t{:},'*')
	title(['disQ == ' int2str(disQ)])
end

%compute spherical convex hull
tol = 1e-6;
[octvtx2,usv] = proj_down(octvtx,tol);
if ~isempty(octvtx2)
	disp('null dimension was found')
	sphK = sphconvhulln(octvtx2);
else
	maxnormQ = true;
	sphK = sphconvhulln(octvtx,maxnormQ);
end

%package output
oref1 = o1;
oref2 = o2;

pts = octvtx2;
disp(savename)
save(savename,'pts','usv','oref1','oref2','five','fiveref1','sphK','octvtx')

end

%--------------------------HELPER FUNCTIONS--------------------------------
function [o12_out,o13_out,mymin,omega23] = GBpair(o1,o2,o3)
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
%--------------------------------------------------------------------------
%calculate distances
[omega12,oct_sym12,zeta12,wveclist12,octonion_pair_sym_list12] = GBdist2([o1 o2],32,false);
[omega13,oct_sym13,zeta13,wveclist13,octonion_pair_sym_list13] = GBdist2([o1 o3],32,false);

%take away imaginary part if any
omega12 = abs(omega12);
wveclist12 = abs(wveclist12);
omega13 = abs(omega13);
wveclist13 = abs(wveclist13);

%unpack symmetrized octonions
o12_sym1 = oct_sym12(1:8);
o12_sym2 = oct_sym12(9:16);

o13_sym1 = oct_sym13(1:8);
o13_sym2 = oct_sym13(9:16);

%find all symmetrized octonions with same omega
minIDs12 = find(ismembertol(wveclist12,omega12,1e-6,'DataScale',1));
octpairsymlist12 = octonion_pair_sym_list12(minIDs12,:);

minIDs13 = find(ismembertol(wveclist13,omega13,1e-6,'DataScale',1));
octpairsymlist13 = octonion_pair_sym_list13(minIDs13,:);

%remove duplicate rows (low tol OK b.c. matching against 16 numbers)
[~,minIA] = uniquetol(round(octpairsymlist12,12),1e-3,'ByRows',true,'DataScale',1);
[~,minIA2] = uniquetol(round(octpairsymlist13,12),1e-3,'ByRows',true,'DataScale',1);

min12 = octpairsymlist12(minIA,:);
min13 = octpairsymlist13(minIA2,:);

%get omega values of combinations of second octonions from sets
%initialize
wveclist23 = zeros(1,size(min12,1)*size(min13,1));
k = 0;
for i = 1:size(min12,1)
	for j = 1:size(min13,1)
		k = k+1;
		wveclist23(k) = get_omega(min12(i,9:16),min13(j,9:16));
		min12list(k,:) = min12(i,9:16);
		min13list(k,:) = min13(j,9:16);
	end
end

%get minimum omega value within precision
prec = 6; % number of decimals
tol = 1e-6;
[mymin,~] = min(round(wveclist23,prec));

%get corresponding octonions
myminIDs = find(abs(round(wveclist23 - mymin,prec)) < tol);
o12 = min12list(myminIDs,:); %output
o13 = min13list(myminIDs,:); %output

o12 = uniquetol(round(o12,prec),tol,'ByRows',true,'DataScale',1);
o13 = uniquetol(round(o13,prec),tol,'ByRows',true,'DataScale',1);

%arbitrarily take first octonion
if size(o12,1) > 2
	disp('')
end

% qA_0 > qB_0 convention added based on discussion with Toby Francis
if (o12(1,1) > o12(1,5)) || (size(o12,1) == 1)
	o12_out = o12(1,:);
else
	o12_out = o12(2,:);
end

if (o13(1,1) > o13(1,5)) || (size(o13,1) == 1)
	o13_out = o13(1,:);
else
	o13_out = o13(2,:);
end

%calculate distance again using GBdist (for comparison)
[omega23,oct_sym23,zeta23] = GBdist([o12_out o13_out],32,false);

%%display results
mat = [omega12;omega13;mymin;omega23];
T = array2table(mat,'VariableNames',{'values'},...
	'RowName',{'Omega12','Omega13','Omega23_pair','Omega23_GBdist'});
% 	disp([name1 '-->' name2 ', ' name1 '-->' name3])
disp(T);


end
%---------------------------END GBpair()-----------------------------------


%--------------------------CODE GRAVEYARD----------------------------------
%{
fivetemp.q = qB;
fivetemp.nA = nB;
fivetemp.d = q2rod(qB);
fivetemp.geometry = name2;

if plotQ
	figure
	
end

if o2addQ
	five = [fivetemp five];
else
	octvtx(1,:) = [];
end

	if plotQ
		plotFZrodriguez_vtx();
		hold on
		%plot rodrigues point
		fivetmp = GBoct2five(octvtx(i+1,:));
		disp(fivetmp)
		
% 		t=num2cell(vertcat(fivetmp.d),1);
		t = num2cell(q2rod(disorientation(vertcat(fivetmp.q),'cubic')),1);
		
		plot3(t{:},'k*')
		hold off
	end


%compute spherical convex hull
tol = 1e-6;
% [octvtx2,usv] = proj_down(octvtx,tol);
% if ~isempty(octvtx2)
% 	sphK = sphconvhulln(octvtx2);
% else
maxnormQ = true;
sphK = sphconvhulln(octvtx,maxnormQ);
% end

%calculate distances
[omega12,oct_sym12,zeta12,wveclist12,octonion_pair_sym_list12] = GBdist2([o1 o2],32,false);
[omega13,oct_sym13,zeta13,wveclist13,octonion_pair_sym_list13] = GBdist2([o1 o3],32,false);
%}
