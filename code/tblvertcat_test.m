%tblvertcat_test
clear; close all

nrows = 3;
%define table columns
doubles = rand(nrows,1);
chars = repelem('a',nrows,1);
cells = repelem({rand(10)},nrows,1);

%make two tables
tbl1 = table(doubles,chars);
tbl2 = table(chars,cells);

%catenate the tables
tblout = tblvertcat(tbl1,tbl2)

%simple outerjoin version
tbl1.ids = (1:3).';
tbl2.ids = (4:6).';
tbloutjoin = outerjoin(tbl1,tbl2,'Key',{'ids','chars'},'MergeKeys',true);
removevars(tbloutjoin,'ids')