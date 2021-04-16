function [octvtx,oref,fiveref,ids] = get_octpairs(pts,epsijk,nv)
arguments
	pts(:,8) double {mustBeSqrt2Norm}
    epsijk(1,1) double = 1
	nv.o2addQ(1,1) logical = false
	nv.pgnum(1,1) double = 32
	nv.wtol double = []
    nv.fiveref = []
    nv.oref(1,8) double = get_ocubo(1,'random',[],10)
    nv.dispQ = []
    nv.nNN(1,1) double = 1 %number of NNs
    nv.IncludeTies(1,1) {mustBeLogical} = true
end
% GET_OCTPAIRS  Get a set of octonions that are symmetrized with respect to a fixed reference GB (default rng seed == 10)
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Inputs:
%  pts - rows of octonions
%
%  o2addQ -	logical, whether to add o2 (in this case, point 'O' with
%		nA = [0 0 1]) to "five", and if not, then get rid
%       of pts(1,:) which corresponds to the same.
%
% Outputs:
%  octvtx - rows of octonions that form a mesh
%
%  usv - struct to use with proj_down.m and proj_up.m
%
%  five - struct containing misorientation quaternions (q), BP normals (nA,
%  grain A reference frame), rodrigues vectors (d), and misorientation
%  fundamental zone feature type (geometry)

% Dependencies:
%  misFZfeatures.mat
%  GBdist4.m
%  mustBeSqrt2Norm.m (argument validation function)
%--------------------------------------------------------------------------
dispQ = nv.dispQ;
nNN = nv.nNN;
IncludeTies = nv.IncludeTies;
if isempty(dispQ)
    if size(pts,1) <= 1000
        dispQ = false;
    else
        dispQ = true;
    end
end

fnames = {'PGnames.mat'};
addpathdir(fnames)

%% Unpack 5DOF reference (empty is OK)
fiveref = nv.fiveref;

%% get reference octonion
if isempty(fiveref)
    oref = nv.oref;
else
    oref = five2oct(fiveref,epsijk);
end

% if isempty(savename)
%     savename = 'temp.mat';
% end

%% get minimized distance octonions relative to oct pairs
if dispQ
    disp('get_octpairs ')
end
npts = size(pts,1);
% orefrep = repmat(oref,npts,1);
[dmin,octvtx] = GBdist4(oref,pts,nv.pgnum,'norm',nv.wtol,dispQ,epsijk,'IncludeTies',IncludeTies,'nNN',nNN);

len = cellfun(@(x) size(x,1),octvtx);
ids = arrayfun(@(x) repelem(x,len(x)),1:npts,'UniformOutput',false);
ids = [ids{:}];
% ids = cumsum(cellfun(@(x) size(x,1),octvtx));
% %check if multiple octonions found (rare, otherwise might indicate an error)
% idstmp = cellfun(@(oct) size(oct,1),octvtx) > 1;
% nids = sum(idstmp);
% if nids > 0
%     disp(['nids: ' int2str(nids)])
%     %display the id since it's a rare occurrence
%     disp(find(idstmp))
%     %replace octonions with first octonion
%     ids = find(idstmp);
%     for i = 1:length(ids)
%         id = ids(i);
%         octvtx{id} = octvtx{id}(1,:);
%     end
% end
%catenate
octvtx = vertcat(octvtx{:});

if nv.o2addQ
    %add reference octonion
	octvtx = [oref; octvtx];
end

% %save data
% if ~isempty(savename)
%     if exist('./data','dir') == 7
%         savepath = fullfile('data',savename);
%     else
%         savepath = savename;
%     end
%     disp(savepath)
%     save(savepath,'pts','oref','octvtx')
% end
end

%--------------------------HELPER FUNCTIONS--------------------------------
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



%unpack symmetrized octonions
o12_sym1 = oct_sym12(1:8);
o12_sym2 = oct_sym12(9:16);

o13_sym1 = oct_sym13(1:8);
o13_sym2 = oct_sym13(9:16);





prec = 6;
tol = 1e-6;

method = 2;
switch method
	case 1
		%both with respect to o3
		
		%calculate distances
		[omega12,oct_sym12,zeta12,wveclist12,octonion_pair_sym_list12] = GBdist2([o1 o3],32,false);
		[omega13,oct_sym13,zeta13,wveclist13,octonion_pair_sym_list13] = GBdist2([o2 o3],32,false);
		
		wveclist3 = zeros(size(wveclist12,1),size(wveclist13,1));
		for i = 1:size(wveclist12)
			for j = 1:size(wveclist13)
				wveclist3(i,j) = wveclist12(i)+wveclist13(2);
			end
		end
		
		[omega3,minID3] = min(wveclist3);
		%find all symmetrized octonions with same omega
		minIDs12 = find(ismembertol(wveclist3,omega3,1e-6,'DataScale',1));
		
	case 2
		%both with respect to o1
		
		%calculate distances
		[omega12,oct_sym12,zeta12,wveclist12,octonion_pair_sym_list12] = GBdist2([o1 o2],32,false);
		[omega13,oct_sym13,zeta13,wveclist13,octonion_pair_sym_list13] = GBdist2([o1 o3],32,false);
end

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
		k = k +1;
		wveclist23(k) = get_omega(min12(i,9:16),min13(j,9:16));
		min12list(k,:) = min12(i,9:16);
		min13list(k,:) = min13(j,9:16);
	end
end

%get minimum omega value within precision
[mymin,~] = min(round(wveclist23,prec));

%get corresponding octonions
myminIDs = find(abs(round(wveclist23 - mymin,prec)) < tol);
o12 = min12list(myminIDs,:); %output
o13 = min13list(myminIDs,:); %output

o12 = uniquetol(round(o12,prec),tol,'ByRows',true);
o13 = uniquetol(round(o13,prec),tol,'ByRows',true);

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
% mat = [omega12;omega13;mymin;omega23];
% T = array2table(mat,'VariableNames',{'values'},...
% 	'RowName',{'Omega12','Omega13','Omega23_pair','Omega23_GBdist'});
% % 	disp([name1 '-->' name2 ', ' name1 '-->' name3])
% disp(T);

% if method == 1
	wveclist3 = zeros(1,size(min1,1),size(min2,1));
% else
% 	wveclist3 = zeros(1,size(min1,1)*size(min2,1));
% 	min1list
% 	k = 0;
% end

% 				k = k +1;
% 				wveclist3(k) = get_omega(min1(i,9:16),min2(j,9:16));
% 				min1list(k,:) = min1(i,9:16);
% 				min2list(k,:) = min2(j,9:16);

o12 = min1list(row,:); %output
o13 = min2list(col,:); %output

% 		wveclist3 = zeros(size(wveclist1,1),size(wveclist2,1));
% 		for i = 1:size(wveclist1)
% 			for j = 1:size(wveclist2)
% 				wveclist3(i,j) = wveclist1(i)+wveclist2(2);
% 			end
% 		end
%
% 		[omega3,minID3] = min(wveclist3);
% 		%find all symmetrized octonions with same omega
% 		minIDs1 = find(ismembertol(wveclist3,omega3,1e-6,'DataScale',1));

		% myminIDs = find(abs(round(wveclist3 - mymin,prec)) < tol);

% function [o2_out,o3_out,omega3,omega3_GBdist] = GBpair(o1,o2,o3)


for i = 1:npts
	%unpack other octonion in pair
	%(o2 and o3 form a pair, each is compared to o1)
	o3 = pts(i,:); %input
	[~,octvtx(i+1,:),omega23_pair(i+1),omega23_GBdist(i+1)] = GBpair(o1,o2,o3);
end


% if ~isempty(octvtx2)
% 	disp('null dimension was found')
% 	sphK = sphconvhulln(octvtx2);
% else
% 	maxnormQ = true;
% 	sphK = sphconvhulln(octvtx,maxnormQ);
% end


% if size(pts,2) == 7
% 	pts = [pts zeros(size(pts,1),1)]; % add column of zeros
% end

%% correct octonions if necessary

% if norm(o) == 1 within tolerance, multiply by sqrt(2)
if abs(norm(pts(1,:)) - 1) < 1e-6
	pts = pts*sqrt(2);
elseif abs(norm(pts(1,:)) - sqrt(2)) > 1e-6
	error('norm of octonions ~= 1 || sqrt(2)')
end
pts = sqrt2norm(pts);


load_type = 'evalc'; %'evalc', 'manual'
switch load_type
	case 'evalc'
		vars = fields(opts);
		for i = 1:length(vars)
			var = vars{i};
			temp = opts.(var); %#ok<NASGU> %temporary value of vName
			evalc([var '= temp']); %assign temp value to the field name
		end
	case 'manual'
		o2addQ = opts.o2addQ;
		plotQ = opts.plotQ;
		method = opts.method;
end

% default_o2addQ = true;
% defaultplotQ = false;
% defaultmethod = 2;
%
% P = inputParser;
% addRequired(P,'pts',@isnumeric);
% addRequired(P,'five',@isstruct);
% addRequired(P,'savename',@ischar);
% addParameter(P,'plotQ',defaultplotQ,@islogical);
% addParameter(P,'method',defaultmethod,@isscalar);
% addParameter(P,'o2addQ',default_o2addQ,@islogical);
% parse(P,pts,five,savename,varargin{:});
%
% plotQ = P.Results.plotQ;
% method = P.Results.method;
% o2addQ = P.Results.o2addQ;



% [~,oct_sym0] = GBdist4(o1,o2,32,'norm',nv.wtol);

%unpack no boundary point
% name2 = 'O';
% disp(['name2 = ' name2])
% qB = normr(qlist.(name2));
% qB = normr(qB+0.05*rand(1,4));

% [~,RB] = symaxis(qB,name2);
% nB = normr((RB*[0 0 1].').');
% nB = normr(nB+0.05*rand(1,3));
% o2 = GBfive2oct(qB,nB);
% o2 = [-1 0 0 0 1 0 0 0]; %input
% o2 = sqrt(2)*[1 0 0 0 0 0 0 0]; %certainly seems to speed things up

% [omega0,oct_sym0,zeta0] = GBdist2([o1 o2],32,false);
% [omega0,oct_sym0,zeta0] = GBdist([o1 o2],32,false);

%take the symmetrized versions for comparison
% o2 = oct_sym0{1}(1,:);

% o1 = oct_sym0(1:8);
% o2 = oct_sym0(9:16);


% octvtx{1} = oct_sym0(9:16);

%textwaitbar setup
% D = parallel.pool.DataQueue;
% afterEach(D, @nUpdateProgress);
% N=npts;
% p=1;
% reverseStr = '';
% nreps2 = floor(N/20);
% nreps = nreps2;

% 	function nUpdateProgress(~)
% 		percentDone = 100*p/N;
% 		msg = sprintf('%3.0f', percentDone); %Don't forget this semicolon
% 		fprintf([reverseStr, msg]);
% 		reverseStr = repmat(sprintf('\b'), 1, length(msg));
% 		p = p + nreps;
% 	end


% parfor i = 1:npts %parfor compatible
% 	%text waitbar
% 	if mod(i,nreps2) == 0
% 		send(D,i);
% 	end
% 	
% 	%unpack other octonion in pair
% 	o3 = pts(i,:); %input
% 	%symmetrized pairs
% 	[octvtx{i+1},omega3(i+1)] = GBpair(o1,o2,o3,nv.pgnum,nv.method,nv.wtol);
% end

% o3 = pts(1,:);
% [octvtx(1,:),~,omega23_pair(1),omega23_GBdist(1)] = GBpair(o1,o2,o3);

%loop through pairs relative to interior point. Each pair contains (+z) origin point


% npts = size(pts,1);

% octvtx = cell(1,npts);
% octvtx{1} = o1;
% t = num2cell(o2,2);
% [octvtx{2:end}] = t{:};
% 
% if ~nv.o2addQ
% 	octvtx{1} = [];
% end
% octvtx = vertcat(octvtx{:});


	nv.method char {mustBeMember(nv.method,{'standard','pairwise'})} = 'pairwise'


name1 = 'random';
switch name1
    case 'interior'
        qA = qlist.(name1);
        
        %load normals (both are arbitrary set to [0 0 1])
        [~,RA] = symaxis(qA,name1);
        nA = (RA*[0 0 1].').';
        
        %package some "five" output for saving
        fiveref1.q =qA;
        fiveref1.nA = nA;
        fiveref1.d = q2rod(qA);
        fiveref1.geometry = name1;
        
        %convert to octonions
        oref = GBfive2oct(qA,nA);
        
    case 'random'
        % o1 = get_ocubo(1,'random',[],10);
        oref = get_ocubo;
end

% tol = 1e-3;
% [octvtx2,usv] = proj_down(octvtx,tol,'zeroQ',true);


% pts = octvtx2;



% load('misFZfeatures.mat','qlist')

fnames = {'PGnames.mat','olist.mat','misFZfeatures.mat'};



	nv.plotQ(1,1) logical = false

if nv.plotQ
    % compute 5DOF representation
    five = GBoct2five(octvtx,true);
    figure
	plotFZrodriguez_vtx();
	hold on
    t = num2cell(q2rod(disorientation(vertcat(five.q),'cubic')),1);
	plot3(t{:},'*')
	title(['disQ == ' int2str(disQ)])
end

% fnames = {'PGnames.mat','olist.mat'};
% addpathdir(fnames)

%}
