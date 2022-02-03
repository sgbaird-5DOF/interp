function Sout = structvertcat(S)
arguments (Repeating)
    S struct
end
% STRUCTVERTCAT  "vertically" concatenate structs with different variables, filling in dummy values as needed
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-09-05
%
% Inputs:
%  S - struct, where each struct corresponds to a set of e.g.
%  conditions/parameters/results for a unique experiment, and each struct
%  can have different numbers of rows and same and/or different variables*
%
% Outputs:
%  Sout - "vertically" concatenated structure
%
% Usage:
%  Sout = tblvertcat(S1,S2);
%  Sout = tblvertcat(S1,S2,S3);
%
% Dependencies:
%  tblvertcat.m (FEX)
%
% Notes:
%--------------------------------------------------------------------------

% convert structures to tables
T = cellfun(@(S) struct2table(S,'AsArray',true),S,'UniformOutput',false);
% concatenate tables
Tcat = tblvertcat(T{:});
% convert back to structure
Sout = table2struct(Tcat);
end
