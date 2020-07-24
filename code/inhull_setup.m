function oint = inhull_setup(datapts,usv,xyz,tess,five,savename)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-15
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
ndatapts = size(datapts,1);

%symmetrically equivalent octonions for each datapoint
% pgnum = 32;
% olist = osymsets(datapts,pgnum,usv);
% olist2 = vertcat(olist{:});


o2addQ = false;
[olist,usv2] = get_octpairs(datapts,five,savename,o2addQ);
olist2 = proj_up(olist,usv2);

% if size(datapts,2) == 7
% 	datapts = [datapts zeros(ndatapts,1)];
% end
% load('octvtx_pairmin.mat','oref1')
% for i = 1:ndatapts
% 	[~, ~, ~,wveclist{i},octonion_pair_sym_list{i}] = GBdist2([oref1 datapts(i,:)],32,false);
% end
% olist2 = vertcat(octonion_pair_sym_list{:});
% olist2 = olist2(:,9:16);

if size(xyz,2) == 7
	xyzup = proj_up(xyz,usv);
else
	xyzup = xyz;
end

d = size(xyzup,2);
pts2 = [zeros(1,d); xyzup];
f = waitbar(0);
k = 0;
tol = 1e-3;
for i = 1:size(olist2,1)
% for i = 1:10000
	pts2(1,:) = olist2(i,:);
	newpts = proj_down(pts2,tol,usv);
	
	if ~isempty(newpts)
		k = k+1;
		newpts2(k,:) = newpts(1,:);
	end
	
	if mod(i,100) == 0
		waitbar(i/size(olist2,1),f)
	end
end
close(f)
xyzdown = proj_down(xyz,tol,usv);
tic
intIDs = intersect_facet(xyzdown,tess,newpts2,1e-1);
toc
intIDs2 = vertcat(intIDs{:});

disp('')
% 
% o2addQ = true;
% olist3 = get_octpairs(pts,five,savename,o2addQ);
% 
% intersect_facet
% sphconvhulln
% 
% [orthoPts,orthoK] = orthoplex(8);
% intfacetIDs = intersect_facet(orthoPts,orthoK,olist{1});
% intfacetIDs2 = intersect_facet(orthoPts,orthoK,xyz);
% intfacetIDs3 = intersect_facet(xyz,tess,olist{1});
% intfacetIDs3 = intersect_facet(normr(xyz),convhulln(normr(xyz)),olist{1}(:,1:7));
% 
% intfacetIDs3 = intersect_facet(normr(xyz),tess,normr(olist3));
% 
% xyz = proj_up(xyz,usv.V,usv.avg);
% [olistproj,usv] = proj_down(olist{1});

% in = inhull(olist2,[ones(1,8); xyz],[],1e-12);
% % in = inhull(olist2,xyz,tess,1e-12);
% inIDs = find(in);

%correlate back with olist cell array
%setup
nocts = cellfun(@(omat) size(omat,1),olist); % # of octonions in each cell
octbins = cumsum(horzcat(nocts{:}))+1; % +1 to change from 0 to 1 indexing
Y = discretize(inIDs,[1 octbins]);
%loop through cells
npts = size(datapts,1);
ointlist = cell(1,npts);
oint = zeros(npts,8);
for i = 1:npts
	ointlist{i} = olistcat(Y == i); %intersecting octonions
	
	if size(ointlist{i},1) > 1
		disp(['octonions found for datapoint ' int2str(i) ':' int2str(size(ointlist{i},1))])
	end
	oint(i,:) = ointlist{i}(1,:); %arbitrarily take first one
end

end


%-------------------------------CODE GRAVEYARD-----------------------------
%{

%check if in convex hull


% o2addQ = false;
% olist = get_octpairs(pts,five,savename,o2addQ);
%}