function tbl = tblfilt(tbl,pars)
% TBLFILT  filter a table based on a struct of parameters (deprecated, see built-in findgroups() and splitapply())
names = fieldnames(pars);
checkQ = false;
for i = 1:length(names)
    name = names{i};
    if ismember(name,tbl.Properties.VariableNames)
        ids = any(tbl.(name) == pars.(name),2);
        tbl = tbl(ids,:);
        if sum(ids) == 0 && ~checkQ
            warning([name ' and possibly others resulted in no matches'])
            checkQ = true;
        end
    end
end