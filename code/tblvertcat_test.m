%tblvertcat_test
clear; close all

nrows = 3;
%define table columns
doubles = rand(nrows,1);
chars = repelem('a',nrows,1);
cells = repelem({rand(10)},nrows,1);
structs = repelem({struct()},nrows,1);

%make two tables
tbl1 = table((1:3).',doubles,chars,structs);
tbl2 = table((4:6).',chars,cells);
cells = repelem({rand(11)},nrows,1);
tbl3 = table((7:9).',cells);

%catenate the tables
%chars will get converted to strings
tblout = tblvertcat(tbl1,tbl2,tbl3)

%chars will stay chars
tblout2 = outerjoin(tbl1,tbl2,'Key',{'Var1','chars'},'MergeKeys',true)
%error due to cells with non-char vector values
tblout3 = outerjoin(tbl2,tbl3,'Key',{'Var1','cells'},'MergeKeys',true)

%with outerjoin, I can't merge structs with fields that have different
%dimensions it would seem (something about how sort isn't defined for
%structs). Outerjoin also can't deal with cells that have non-character
%vector values.