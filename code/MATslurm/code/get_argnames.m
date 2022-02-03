function argnames = get_argnames(fn)
% GET_ARGNAMES get anonymous function argument names
fntxt = func2str(fn);
startID = 3;
endID = strfind(fntxt,')')-1;
argnames = strsplit(fntxt(startID:endID),',');
end