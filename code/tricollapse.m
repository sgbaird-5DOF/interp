function [uniqueK, uniquePts] = tricollapse(K,pts)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-03
%
% Description: collapse triangulation of points (i.e. a triangulation
% involving repeat points).
%
% Inputs:
%		K				=== convex hull triangulation of pts
%
%		pts			=== rows of points (allowed to have repeat points)
%
% Outputs:
%		uniqueK		=== collapsed triangulation
%
%		uniquePts	=== rows of unique pts
%
% Dependencies:
%		get_repsets.m
%
%--------------------------------------------------------------------------

%initalize
nptstot = length(pts);

%get sets of degenerate & non-degenerate values
irepsets = get_repsets(pts);

nsets = length(irepsets);

%take first point ID from each repset
disp('take first point ID of each repset')
icuList = zeros(1,nsets);
parfor i = 1:nsets
	icuList(i) = irepsets{i}(1);
end

%icuList = cellfun(@(irep) irep(1),irepsets); %unique ic values


%make copy of K
uniqueK = K;

%% find fixQs
tic
fixQ = cell(1,nsets);
for_type = 'parfor';
switch for_type
	case 'for'
		for i = 1:nsets
			fixQ{i} = find(ismember(K,irepsets{i}));
			K(fixQ{i}) = 0; %might make sorting faster for next repetition (but doesn't work with parallelization)
		end
	case 'parfor'
		%textwaitbar setup
		D = parallel.pool.DataQueue;
		afterEach(D, @nUpdateProgress);
		N=nsets;
		p=1;
		reverseStr = '';
		if nsets > 100
			nreps = floor(nsets/100);
		    nreps2 = floor(nsets/100);
		else
			nreps = 1;
         nreps2 = 1;
		end
		
		disp('loop through sets')
		parfor i = 1:nsets
			fixQ{i} = find(ismember(K,irepsets{i}));
			
			if mod(i,nreps2) == 0
				send(D,i);
			end
		end
		
		disp(' ')
end

%% reformat uniqueK
for i = 1:nsets
	uniqueK(fixQ{i}) = i+nptstot; % replace degenerate locations with a unique ID
end

uniqueK = uniqueK - nptstot; % shift IDs to go from 1:nsets (in same order as icuList)

toc
disp(' ')

%% reformat pts
uniquePts = pts(icuList,:);

	function nUpdateProgress(~)
		percentDone = 100*p/N;
		msg = sprintf('Tricollapse: %3.1f ', percentDone); %Don't forget this semicolon
		fprintf([reverseStr, msg]);
		reverseStr = repmat(sprintf('\b'), 1, length(msg));
		p = p + nreps;
	end

end %tricollapse


%--------------------------------CODE GRAVEYARD----------------------------
%{

for i = 1:nsets
	uniqueK(fixQ{i}) = i-nptstot; %give it a non-positive index to keep IDs unique during recursion
end

uniqueK = uniqueK + nptstot; %shift IDs to go from 1:nsets


for i = 1:nsets
	irepset = irepsets{i}; %unpack
	uID = irepsets(1); % "unique ID", take first ID from a set of degenerate IDs
	uniqueK(fixQ{i}) = icu+nptstot; % replace degenerate locations with unique ID
end

for i = 1:nsets
	icu = icuList(i);
	uniqueK(fixQ{i}) = icu-nptstot; % replace degenerate locations with unique ID
end	


% 		msg2 = ['looping through ' int2str(nsets) ' ireps'];
% 		f = waitbar(0,msg2);
		for i = 1:nsets
			fixQ{i} = find(ismember(K,irepsets{i}));
			K(fixQ{i}) = 0; %might make sorting faster for next repetition (but doesn't work with parallelization)
% 			if mod(i,100) == 0
% 				waitbar(i/nsets,f);
% 			end
		end
% 		close(f)
%}
