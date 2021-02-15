function captionlist = get_captionlist(descriptions)
arguments
    descriptions = {'one','two','three'} %[1,2,3]
end

ndescript = length(descriptions);
if (~iscell(descriptions) && (ndescript ~= 1)) || (~iscell(descriptions) && isnumeric(descriptions(1)))
    descriptions = num2cell(descriptions);
elseif ~iscell(descriptions) && ndescript == 1
    warning('descriptions should be a cell array of chars or strings, or a numeric array')
end

if isnumeric(descriptions{1})
    for i = 1:ndescript
        descriptions{i} = char(string(descriptions{i}));
    end
end
charlbl = get_charlbl();

switch ndescript
    case 1
        captionlist = [charlbl{1} ' ' descriptions{1}];
    case 2
        captionlist = [charlbl{1} ' ' descriptions{1} ', and ' charlbl{2} ' ' descriptions{2}];
    otherwise
        if ndescript >= 3
            captionlist = '';
            for i = 1:ndescript-1
                captionlist = [captionlist, charlbl{i} ' ' descriptions{i} ', '];
            end
        captionlist = [captionlist, ' and ' charlbl{ndescript} ' ' descriptions{ndescript}];
        else
            warning('description should be a non-negative integer')
        end
end
end