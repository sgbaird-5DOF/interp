function lbls = lblcat(numlist,mkr)
arguments
    numlist = [1,10,100]
    mkr = '-'
end
% LBLCAT  catenate an array of numbers as a string with marker separations (mkr)
lbls = cellfun(@char,num2cell(string(numlist)),'UniformOutput',false);
lbls = strcat(lbls,mkr);
lbls{end}(end) = []; %remove last catenation marker
lbls = [lbls{:}];

end