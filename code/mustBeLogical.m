function mustBeLogical(value)
% MUSTBELOGICAL  must be a logical value, error otherwise
%NOTE: do not use "arguments" syntax here since this is a validation fn
%---------------------HELPER VALIDATION FUNCTION---------------------------
mustBeA(value,'logical')
end