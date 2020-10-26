function Sout = structhorzcat(S)
arguments (Repeating)
    S(1,1) struct
end
% STRUCTHORZCAT  "horizontally" concatenate any number of structures with
% different variables (no overlapping variables, scalar struct).
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
%  Sout - "horizontally" concatenated structure
%
% Usage:
%  Sout = tblvertcat(S1,S2);
%  Sout = tblvertcat(S1,S2,S3);
%
% Dependencies:
%
% Notes:
%--------------------------------------------------------------------------

% convert structures to tables
T = cellfun(@(S) struct2table(S,'AsArray',true),S,'UniformOutput',false);
% concatenate tables
Tcat = horzcat(T{:});
% convert back to structure
Sout = table2struct(Tcat);
end
