%%
clear; close all
% User input character vector:
str = 'stores{2}.name';
% Create a 2xN cell array {index type ; field/index} from that char vector:
tkn = regexp(['.',str],'(\W)(\w+)','tokens'); % assumes the first part is always a field.
tkn = vertcat(tkn{:}).';
tkn(1,:) = strrep(strrep(tkn(1,:),'{','{}'),'(','()');
% Convert linear indices to numeric:
vec = str2double(tkn(2,:));
idx = ~isnan(vec);
tkn(2,idx) = num2cell(num2cell(vec(idx)));
% Create indexing structure from the cell array:
sbs = substruct(tkn{:});

val = subsref(config,sbs)

%%
K = cell(1,2); K{2} = cell(1,2);
S.type = '{}';
S.subs = {2};
val = subsref(K,S)

%%
a2 = {{'asdf',{'asdf','asdf',{'asdf'},{'asdf'},{'asdf',{'asdf'}},'asdf'},{'asdf',{'asdf'},{{'asdf'}}}},'asdf'}

while any(cellfun(@iscell,a2))
    a2 = [a2{cellfun(@iscell,a2)} a2(~cellfun(@iscell,a2))];
end


%%
%{
How to apply a tree-based indexing scheme for a triangulation?

%}