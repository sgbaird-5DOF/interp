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
    startidx = regexp(fpath{1},'octonion-inference*');
    fpathshort{i} = fpath{i}(startidx:end);
    nameExt{i} = [names{i},ext{i}];
end
end