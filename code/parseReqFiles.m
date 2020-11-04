function [fpathshort,nameExt,plist] = parseReqFiles(fname,expression,NV)
arguments
    fname = 'interp5DOF.m'
    expression = 'interp*'
    NV.sep char {mustBeMember(NV.sep,{'/','\'})} = filesep()
    NV.removefname(1,1) logical = true
end
% PARSEREQFILES  parse required files & products using fileseparator sep
% (default is native separator of computer running MATLAB) for files that
% are accessible via MATLAB path.
[flist,plist] = matlab.codetools.requiredFilesAndProducts(fname);
nfiles = length(flist);

%initialize outputs
[fpathshort,nameExt] = deal(cell(nfiles,1));

%initialize intermediates
[fpath,names,ext] = deal(cell(nfiles,1));

%get names
for i = 1:nfiles
    f = flist{i};
    [fpath{i},names{i},ext{i}] = fileparts(f);
    fpath{i} = strrep(fpath{i},'/',NV.sep);
    fpath{i} = strrep(fpath{i},'\',NV.sep);
    [~,endidx] = regexp(fpath{i},expression);
    fpathshort{i} = [fpath{i}(endidx+2:end) NV.sep];
    nameExt{i} = [names{i},ext{i}];
end

if NV.removefname
    %remove top-level file
    topIDtmp = strfind(nameExt,fname);
    topID = cellfun(@(x)~isempty(x),topIDtmp);
    fpathshort(topID) = [];
    nameExt(topID) = [];
end
end