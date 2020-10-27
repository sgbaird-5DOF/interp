%parfilter_test

A = categorical({'a','b','c'});
B = 1:3;
pars = struct('A',A(1),'B',B(1));
tbl = struct2table(struct('A',A.','B',B.'));
parfilter(tbl,pars)