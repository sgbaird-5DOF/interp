function addpathdir(fnames)

if iscell(fnames)
	n = length(fnames);
else
	n = 1;
end

for i = 1:n
	if iscell(fnames)
		fname = fnames{i};
	elseif i == 1
		fname = fnames;
	end
	fpath = fullfile('**',fname);
	file = dir(fpath);
	if ~isempty(file)
		addpath(file.folder);
	end
end
end