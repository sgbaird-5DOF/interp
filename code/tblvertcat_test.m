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