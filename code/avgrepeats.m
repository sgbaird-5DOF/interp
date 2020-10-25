function [meshList,propList] = avgrepeats(meshList,propList)
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-07-27
%
% Description: %average property values for degenerate mesh points and then
% remove all but one mesh point from the degenerate set.
% 
% Inputs:
%  meshList - rows of mesh points (e.g. coordinates)
%  propList - rows of property values (can have more than one, but untested)
%
% Outputs:
%  meshList - rows of unique mesh points propList - rows of unique property
%  propList - property values of degenerate mesh points were averaged
%
% Usage:
%  [meshList,propList] = avgrepeats(meshList,propList);
%
% Dependencies:
%  get_repsets.m
%
% Notes:
%  *
%--------------------------------------------------------------------------
nmeshpts = size(meshList,1);
nprops = size(propList,1);
errmsg = ['nmeshpts: ' int2str(nmeshpts) ', nprops: ' int2str(nprops) ', should be same'];
assert(nmeshpts == nprops,errmsg)

repsets = get_repsets(meshList,2);
nrepsets = length(repsets);
rmIDs = cell(1,nrepsets);
for i = 1:nrepsets
	%unpack
	repset = repsets{i};
	
	%octonion to keep
	keepID = repset(1);
	
	%degenerate octonions to remove (after loop)
	rmIDs{i} = repset(2:end);
	
	%average properties
	propList(keepID,:) = mean(propList(repset,:),1);
end

%remove repeat octonions
rmIDlist = horzcat(rmIDs{:});
meshList(rmIDlist,:) = [];
propList(rmIDlist,:) = [];

end
