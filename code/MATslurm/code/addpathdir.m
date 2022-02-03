function addpathdir(fnames)
%ADDPATHDIR add folders of filenames using addpath() and dir() 
% search through subfolders and add first folder found for each filename
% (if any) to path.
%   Since it only adds the first folder found, then you avoid a disaster that would
%   otherwise come about by e.g. navigating to 'C:\' or '~/' and running
%   addpathdir('*'). However this does limit the functionality somewhat if you
%   wanted to add all folders with files matching the form myname*.m. Instead,
%   this will only add the first file found.

%check if cell input or single name
if iscell(fnames)
    n = length(fnames);
else
    n = 1;
end

for i = 1:n
    %unpack name
    if iscell(fnames)
        fname = fnames{i};
    elseif i == 1
        fname = fnames;
    end
    fpath = fullfile('**',fname);
    %find file(s)
    file = dir(fpath);
    if ~isempty(file)
        %add first folder to path
        addpath(file(1).folder);
    else
        disp(['DNE: ' fname '. Is file name/extension correct? Is file contained in subfolder of working directory?'])
    end
end
end