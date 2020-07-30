function mustContainFields(S,checknames)
	varnames = fieldnames(S);
	errmsg = ['input needs minimum fields: ' strjoin(checknames) ...
		' but contains fields: ' strjoin(varnames)];
	assert(all(ismember(checknames,varnames)),errmsg)
end