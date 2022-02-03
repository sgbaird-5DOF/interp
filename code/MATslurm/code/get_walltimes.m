function walltimes = get_walltimes(parcombs,walltimefn)
% GET_WALLTIMES  get walltimes based on a function of parameters and
% parameter combinations

%% walltime setup
%number of combinations
ncombs = size(parcombs,1);

%extract arguments for walltime function
argnames = get_argnames(walltimefn); %argument names
nargs = length(argnames); %number of arguments

%% walltimes
if ~isempty(argnames{1})
    walltimes = zeros(ncombs,1);
    argvals = cell(ncombs,nargs);
    
    %loop through walltime variables
    for i = 1:nargs
        %unpack argument name
        argname = argnames{i};
        %concatenate and convert to cell
        if ischar(parcombs(1).(argname))
            argvals(:,i) = {parcombs.(argname)};
        else
            argvals(:,i) = num2cell(vertcat(parcombs.(argname)));
        end
    end
    
    %loop through combinations
    for i = 1:ncombs
        %plug into walltime function
        walltimes(i) = walltimefn(argvals{i,:});
    end
else
    %assign a constant value
    walltime = walltimefn();
    walltimes = repelem(walltime,ncombs,1);
end

%add small random value
walltimes = walltimes + 0.01*rand(ncombs,1); %helps with binning in histcounts

end %get_walltimes

%% CODE GRAVEYARD
%{

        %         argval = cellfun(@(arg) arg(i,:),argvals,'UniformOutput',false);

%     argIDs = cellfun(@(argname) find(strcmp(parnames,argname)),argnames); %argument positions
%     for i = 1:length(argnames)
%         argname = argnames{i};
%         argvals{i} = [parcombs.(argname)];
%     end
%     argvals = parcombs(argIDs); %argument values
%}