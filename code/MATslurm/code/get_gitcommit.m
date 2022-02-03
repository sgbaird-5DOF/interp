function [gitcommit, comment, warnedQ] = get_gitcommit()
% GET_GITCOMMIT  get git commit version (or return empty '' if error)
[status,cmdout] = system('git rev-parse HEAD');
if status == 0
    gitcommit = cmdout(1:7);
    comment = cmdout(9:end);
else
    gitcommit = '';
    comment = '';
end

[status,cmdout] = system('git status --untracked-files=no --porcelain');
if status == 0
    if ~isempty(cmdout)
        warning('Working directory not clean (i.e. uncommitted/unpushed) files exist. Use !git commit -am "<message>", then !git push')
        warnedQ = true;
    else
        warnedQ = false;
    end
end