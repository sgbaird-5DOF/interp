%parse required files & products
function [fpathshort,nameExt,plist] = parseReqFiles(fname)
[flist,plist] = matlab.codetools.requiredFilesAndProducts(fname);
nfiles = length(flist);
%initialization
init1 = cell(nfiles,1);
fpath = init1;
names = init1;
ext = init1;
fpathshort = init1;
nameExt = init1;
%get names
for i = 1:nfiles
    f = flist{i}; %file
    [fpath{i},names{i},ext{i}] = fileparts(f);
    [~,endidx] = regexp(fpath{1},'interp-5DOF*');
    fpathshort{i} = [fpath{i}(endidx+2:end) filesep()];
    nameExt{i} = [names{i},ext{i}];
end
end
