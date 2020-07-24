function [Ktr,K_out,newpts] = hypersphere_subdiv(pts,K,nint,varargin)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-03
%
% Description: Subdivide a "spherical" convex hull and collapse the
%		triangulation.
%
% Inputs:
%
% Outputs:
%
% Dependencies:
%		facet_subdiv.m
%
%		tricollapse.m

%		normr.m
%--------------------------------------------------------------------------

if nargin == 4
	tricollapseQ = varargin{1};
else
	tricollapseQ = true;
end

%% compute top-level convex hull
if isempty(K) % && tricollapseQ
	maxnormQ = true;
	K = sphconvhulln(pts,maxnormQ);
	% 	K = convhulln(pts);
% elseif isempty(K)
% 	K = sphconvhulln(pts);
% 	if isempty(K)
% 		K = convhulln(pts);
% 		disp('computing regular convex hull instead.')
% 	end
end

%% Ktree top-level
npts = length(pts);
nfacets = size(K,1);
Ktr.main.K = K;
Ktr.main.pts = pts;
Ktr.sub = cell(1,nfacets);

KtrtempPts = cell(1,nfacets);
KtrtempK = KtrtempPts;

%% Ktree level 2
%initialize counter variables
% nptsrm = 0;
tic

% sphbaryQ = false;

%textwaitbar setup
waitbarQ = true;
if waitbarQ
	%comment these lines & the function "nUpdateProgress", and set waitbarQ
	%== false if you want to add variables during debugging
	
	D = parallel.pool.DataQueue;
	afterEach(D, @nUpdateProgress);
	N=nfacets;
	p=1;
	reverseStr = '';
	if N > 100
		nreps = floor(N/100);
	else
		nreps = 1;
	end
	nreps2 = nreps;
end

K2temp = num2cell(K,2);
ptstemp = cellfun(@(K2) pts(K2,:),K2temp,'UniformOutput',false); 

parfor i = 1:nfacets %parfor compatible, uncomment "send" lines if using parfor (enables text waitbar)
	%extract facet IDs and pts
% 	K2 = K(i,:);
% 	mpts2 = pts(K2,:);
	
	mpts2 = ptstemp{i};
	
	% 	if sphbaryQ
	% 		mpts2 = projfacet2hyperplane(mean(mpts2),mpts2);
	% 	end
	
	%subdivide facet
	if nint > 1
		if i == 1
			delaunayQ = true;
			[mpts2a,K2a] = facet_subdiv(mpts2,nint,delaunayQ);
		else
			delaunayQ = false;
			mpts2a = facet_subdiv(mpts2,nint,delaunayQ); %assumes K2a can apply to next set of pts, probably only valid if data lies on hyperhemisphere
			K2a = [];
		end
	else
		mpts2a = mpts2;
		K2a = 1:size(mpts2a,2);
	end
	
	if ~isempty(mpts2a)
		if sum(myismember(mpts2,mpts2a)) ~= size(mpts2,1)
			1+1;
		end
	end
	
% 	%renormalize to unit hypersphere
%  	mpts2a = normr(mpts2a);
	
	%add subdivision to K-tree
	KtrtempK{i} = K2a;
	KtrtempPts{i} = mpts2a;

	if waitbarQ
		if mod(i,nreps2) == 0
			send(D,i);
		end
	end
end

	function nUpdateProgress(~)
		percentDone = 100*p/N;
		msg = sprintf('hypersphere_subdiv percent done: %3.1f ', percentDone); %Don't forget this semicolon
		fprintf([reverseStr, msg]);
		reverseStr = repmat(sprintf('\b'), 1, length(msg));
		p = p + nreps;
	end

disp(' ')
toc
disp(' ')

K2a = KtrtempK{1};

%update counter for index
nptstemp = size(KtrtempPts{1},1);
for i = 1:nfacets
	%unpackage pts
	mpts2a = KtrtempPts{i};
	
	%update K2a triangulation values (i.e. add constant) so that IDs of new
	%vertices are unique
	if i ~= 1
		K2a = K2a+nptstemp;
	end
	
	%package back into K tree (Ktr)
	Ktr.sub{i}.main.K = K2a;
	Ktr.sub{i}.main.pts = mpts2a;
end

%catenate pts and hull
lvltwo = vertcat(Ktr.sub{:});
lvltwo = vertcat(lvltwo.main);
lvltwoK = vertcat(lvltwo.K);
lvltwoPts = vertcat(lvltwo.pts);

%collapse to single convex hull
if tricollapseQ
 	disp('tricollapse')
	[K_out, newpts] = tricollapse(lvltwoK,lvltwoPts);
else
	disp('uniquetol')
	[~,ia] = uniquetol(round(lvltwoPts,12),'ByRows',true); %careful to not output the rounded values directly
	newpts = lvltwoPts(ia,:);
	K_out = [];
end

assert(sum(ismembertol(pts,newpts,'ByRows',true)) == npts,...
	'points lost during subdiv. check facet_subdiv input degeneracy')

end %hypersphere_subdiv

%----------------------------HELPER FUNCTIONS------------------------------
% 	function o = octnorm(o)
%
%
% 	end


%---------------------------------NOTES------------------------------------
%{
Doing two projections might be faster than collapsing the convex hull. Even
so, I would still want the # of unique vertices from unique(), and I would
need to find connecting facets from a much larger list of points. I'd
probably need irepsets anyway, but maybe only two calls to ismember() per
point.

I can use a cell scheme K = cell(1,2), where K{1} is the main
facet (array) with two arrays? - i.e. level 1, and K{2} is contains
information about the second level. K{2}{1}
%}


%-------------------------------CODE GRAVEYARD-----------------------------
%{

%loop through facets of main hull
for i = 1:nint
	for j = 1:length(K{i})
		K{i,j}
	end
end
	
facet_subdiv(meshpts,nint)




Kcell = num2cell(K,2); %cell array containing rows from K

% rmIDs = find(ismembertol(pts,rmpts,1e-6,'ByRows',true)); % find IDs of remove points

nptsrm = 0;



% 	nfacetstot = nfacetstot + nfacets2;
% 	nfacets2 = length(K2a);


% 	nptstot = nptstot + nptsrm;

nptstot = 0;


% this ended up messing with finding intersections, actually shifted points outside of hull I think
	%renormalize points (individual quaternions)
% 	if size(mpts2a,2) == 7
% 		mpts2a(:,1:4) = normr(mpts2a(:,1:4));
% 		mpts2a(:,5:7) = normr(mpts2a(:,5:7));
% 		mpts2a = 1/sqrt(2)*mpts2a;
% 	end


reverseStr = '';

	if mod(i,floor(nfacets/100)) == 0
		percentDone = 100*i/nfacets;
		msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
		fprintf([reverseStr, msg]);
		reverseStr = repmat(sprintf('\b'), 1, length(msg));
	end


	%update K2a triangulation values (i.e. add constant) so that IDs of new
	%vertices are unique
% 	K2a = K2a + nptsrm;
	
	%update counter for index
% 	nptsrm = size(mpts2a,1);


		KtrtempK{i} = K2a;
	% 	Ktr.sub{i}.main.K = K2a;
	% 	Ktr.sub{i}.main.pts = mpts2a;

KtrtempPts = KtrtempK;


msg = ['looping through ' int2str(nfacets) ' facets..'];


%% compute top-level convex hull
if isempty(K) && tricollapseQ
	maxnormQ = true;
	K = sphconvhulln(pts,maxnormQ);
	% 	K = convhulln(pts);
elseif isempty(K)
	K = sphconvhulln(pts);
% 	if isempty(K)
% 		K = convhulln(pts);
% 		disp('computing regular convex hull instead.')
% 	end
end

%}