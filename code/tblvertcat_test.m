%tblvertcat_test
clear; close all

%define table columns
doubles = rand(10,1);
chars = repelem('a',10,1);
cells = repelem(cell(1,1),10,1);
structs = repelem(struct('var',1),10,1);

%make two tables
tbl1 = table(doubles,chars,cells);
tbl2 = table(chars,cells,structs);

%catenate the tables
tblout = tblvertcat(tbl1,tbl2)