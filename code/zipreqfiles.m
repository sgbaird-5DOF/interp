function [fpathshort,nameExt,plist,fnames] = zipreqfiles(fname,expression,zipfname,NV)
arguments
    fname = 'interp5DOF.m'
    expression = 'interp*'
    zipfname = 'temp.zip'
    NV.sep char {mustBeMember(NV.sep,{'/','\'})} = filesep()
    NV.removefname(1,1) logical = true
end
% ZIPREQFILES  zip required files of a function using parseReqFiles() and zip()
NVpairs = namedargs2cell(NV);
[fpathshort,nameExt,plist]=parseReqFiles(fname,expression,NVpairs{:});
fnames = cellfun(@(fpath,name) [fpath,name],fpathshort,nameExt,'UniformOutput',false);
zip(zipfname,fnames)
end