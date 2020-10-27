function tbl = parfilter(tbl,pars)
% PARFILTER  filter a table based on a struct of parameters
names = fieldnames(pars);
for i = 1:length(names)
    name = names{i};
    if ismember(name,tbl.Properties.VariableNames)
        ids = any(tbl.(name) == pars.(name),2);
        tbl = tbl(ids,:);
    end
end