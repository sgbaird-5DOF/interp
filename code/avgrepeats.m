function [meshList,propList,rmIDs,keepIDs,rmIDcell] = avgrepeats(meshList,propList,avgtype)
arguments
   meshList
   propList
   avgtype char {mustBeMember(avgtype,{'mean','min'})} = 'mean'
end
% AVGREPEATS  average values for duplicate points, remove all but one from degenerate set
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-07-27
% 
% Inputs:
%  meshList - rows of mesh points (e.g. coordinates)
%  propList - rows of property values (can have more than one, but untested)
%
% Outputs:
%  meshList - rows of unique mesh points propList - rows of unique property
%  propList - property values of degenerate mesh points were averaged
%  rmIDlist - list of IDs that got removed from meshList and propList
%
% Usage:
%  [meshList,propList,rmIDlist] = avgrepeats(meshList,propList);
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

if (nprops == 1) && (size(propList,2) ~= 1)
    error('inputs should be row vectors')
end

switch avgtype
    case 'mean'
        avgfn = @(x) mean(x,1);
    case 'min'
        avgfn = @(x) min(x,[],1);
end

repsets = get_repsets(meshList,2);
keepIDs = cellfun(@(repset) repset(1),repsets);
nrepsets = length(repsets);
rmIDcell = cell(1,nrepsets);
for i = 1:nrepsets
	%unpack
	repset = repsets{i};
	
	%degenerate octonions to remove (after loop)
	rmIDcell{i} = repset(2:end);
    
    %octonion to keep
    keepID = repset(1);
    
	%average properties
	propList(keepID,:) = avgfn(propList(repset,:));
end

%remove repeat octonions
rmIDs = horzcat(rmIDcell{:});
meshList(rmIDs,:) = [];
propList(rmIDs,:) = [];

end
