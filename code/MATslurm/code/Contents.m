% CODE
%
% Files
%   addpathdir         - add folders of filenames using addpath() and dir() 
%   exec_argfn         - execute function_handle by automatically retrieving argument values from a structure
%   exec_combs         - execute parameter combinations (task-level on SLURM)
%   get_argnames       - get anonymous function argument names
%   get_argnames_test  - get_argnames_test
%   get_args           - get argument names and values for a function from a structure
%   get_argvals        - get argument values of a function by parsing values from a struct
%   get_cmd            - get slurm command to submit an sbatch script
%   get_cmd_test       - get_cmd test
%   get_gitcommit      - get git commit version (or return empty '' if error)
%   get_uuid           - get unique ID (8 characters, mixture of numbers and letters) via java.util.UUID.randomUUID
%   get_walltimes      - get walltimes based on a function of parameters and
%   sortbins           - sort parameter combinations into bins based on their walltimes
%   sortbins_test      - sortbins_test
%   structhorzcat      - "horizontally" concatenate structures with different variables
%   structhorzcat_test - test function for struct horizontal concatenation
%   structvertcat      - "vertically" concatenate structs with different variables, filling in dummy values as needed
%   structvertcat_test - structvertcat_test
%   submit_sbatch      - submit multiple sbatch jobs that depend on a parameter file via writeparfile.m
%   var_names          - take variables and output a combined struct with each of the variable names
%   writeparfile       - generate a parameters file for an sbatch submission
%   writeparfile_test  - writeparfile_test





