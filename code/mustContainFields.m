function mustContainFields(S,checknames)
% MUSTCONTAINFIELDS  check fieldnames(S) and make sure every checkname exists
	varnames = fieldnames(S);
	errmsg = ['input needs minimum fields: ' strjoin(checknames) ...
		' but contains fields: ' strjoin(varnames)];
	assert(all(ismember(checknames,varnames)),errmsg)
end