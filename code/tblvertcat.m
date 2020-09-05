function tblout = tblvertcat(tbl)
arguments (Repeating)
    tbl table
end
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-09-04
%
% Description: vertically catenate tables with different variables, filling
% in data where necessary.
%
% Inputs:
%  tbl - table, where each table can have a different number of rows and
%  same and/or different variables*
%
% Outputs:
%  tblout - vertically catenated table
%
% Usage:
%  tblout = tblvertcat(tbl1,tbl2);
%
% Notes:
%  See https://www.mathworks.com/matlabcentral/answers/179290-merge-tables-with-different-dimensions
%
%  types 'cell' and 'struct' are not supported by missing, and so are
%  implemented manually.
%
%  *variables of the same name must also be of the same datatype
%--------------------------------------------------------------------------
%number of tables
ntbls = length(tbl);
for n = 1:ntbls
    for p = n:ntbls
        if p ~= n
            %% unpack table pair
            t1 = tbl{n};
            t2 = tbl{p};
            
            %% find missing variables
            %get variable names from t2 that are not in t1
            [missingtmp1,ia1] = setdiff(t2.Properties.VariableNames,t1.Properties.VariableNames);
            %get variable names from t1 that are not in t2
            [missingtmp2,ia2] = setdiff(t1.Properties.VariableNames,t2.Properties.VariableNames);
            
            %% make struct and cell tables
            % struct (technically a table of structs with a 'dummy' field and NaN Value)
            [structtbl1,sreplaceNames1] = replacevartbl(t2,ia1,struct());
            [structtbl2,sreplaceNames2] = replacevartbl(t1,ia2,struct());
            
            % cell (cell with 0x0 double inside)
            [celltbl1,creplaceNames1] = replacevartbl(t2,ia1,cell(1));
            [celltbl2,creplaceNames2] = replacevartbl(t1,ia2,cell(1));
            
            % remove values that are represented in cell and struct tables
            missing1 = setdiff(missingtmp1,[sreplaceNames1 creplaceNames1],'stable');
            missing2 = setdiff(missingtmp2,[sreplaceNames2 creplaceNames2],'stable');
            
            %% splice the missing variable tables into original tables
            % matrices of missing elements to splice into original
            missingmat1 = repelem(missing,height(t1),numel(missing1));
            missingmat2 = repelem(missing,height(t2),numel(missing2));
            
            %tables to splice into original tables
            missingtbl1 = array2table(missingmat1,'VariableNames',missing1);
            missingtbl2 = array2table(missingmat2,'VariableNames',missing2);
            
            %perform the splice
            tbl{n} = [t1, missingtbl1, structtbl1, celltbl1];
            tbl{p} = [t2 missingtbl2, structtbl2, celltbl2];
        end
    end
end
%catenate all tables
tblout = vertcat(tbl{:});
end

%% Replace Variable Table
function [replacetbl,replaceNames] = replacevartbl(t,ia,replaceval)
%replace type
replacetype = class(replaceval);

%% missing variable IDs and names
%variable names
varnames = t.Properties.VariableNames;

%variable types
vartypes=varfun(@class,t,'OutputFormat','cell');

%variable IDs of some type
IDtmp = find(strcmp(replacetype,vartypes));

%missing variable IDs of different types
ID = intersect(ia,IDtmp);

%missing variable names of different types
replaceNames = varnames(ID);

%% construct table with replacement values 
%table dimensions
nrows = height(t);
nvars = length(ID);

if isstruct(replaceval) && isempty(replaceval)
    error('if type struct, cannot be empty. Instead supply struct with no fields via struct()')
end

replacemat = repelem(replaceval,nrows,nvars);
replacetbl = array2table(replacemat);
replacetbl.Properties.VariableNames = replaceNames;

end


%% CODE GRAVEYARD
%{
            %get variable types for each
            vartypes1=varfun(@class,t1,'OutputFormat','cell');
            vartypes2=varfun(@class,t2,'OutputFormat','cell');
            
            %% find variable IDs of different types
            %struct
            structIDtmp1 = find(strcmp('struct',vartypes1));
            structIDtmp2 = find(strcmp('struct',vartypes2));
            
            %cell
            cellIDtmp1 = find(strcmp('cell',vartypes1));
            cellIDtmp2 = find(strcmp('cell',vartypes2));
            
            %% find missing variable IDs of different types
            %struct
            structID1 = union(ia1,structIDtmp1);
            structID2 = union(ia2,structIDtmp2);
            
            %cell
            cellID1 = union(ia1,cellIDtmp1);
            cellID2 = union(ia2,cellIDtmp2);


% for i = 1:n
%     strcmp('struct',varfun(@class,tbl{i},'OutputFormat','cell')
% end

%             ia1([structID1,cellID1]) = [];
%             ia2([structID2,cellID2]) = [];

if isstruct(replaceval)
    for i = ID
        varname = varnames{ID};
        sfields = fields(t(varname));
        for j = 1:length(sfields)
            sfield = sfields{j};
            sfieldtype = 
            replaceval.(sfield) = 
%}