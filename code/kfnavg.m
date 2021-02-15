function cov = kfnavg(kfns,X1,X2)
arguments
    kfns function_handle
    X1(:,:) double
    X2(:,:) double
end
% KFNAVG  compute the average covariance matrix based on several covariance functions
% example: kfn = @(X1,X2) kfnavg(kfns,X1,X2);

n = length(kfns);

cov = [];
for i = 1:n
    cov = cov+kfns{i}(X1,X2);
end
cov = cov/n;

end

%% CODE GRAVEYARD
%{
% kfns = cellfun(@(fn) fn(X1,X2),kfns,'UniformOutput',false);
%
% switch n
%     case 2
%         kfn = @(X1,X2) (kfns{1}(X1,X2)+kfns{2}(X1,X2))/2;
%     case 3
%         kfn = @(X1,X2) (kfns{1}(X1,X2)+kfns{2}(X1,X2)+kfns{3}(X1,X2))/3;
%     case 4
%         kfn = @(X1,X2) (kfns{1}(X1,X2)+kfns{2}(X1,X2)+kfns{3}(X1,X2)+kfns{4}(X1,X2))/4;
%     case 5
%         kfn = @(X1,X2) (kfns{1}(X1,X2)+kfns{2}(X1,X2)+kfns{3}(X1,X2)+kfns{4}(X1,X2)+...
%             kfns{5}(X1,X2))/5;
%     case 6
%         kfn = @(X1,X2) (kfns{1}(X1,X2)+kfns{2}(X1,X2)+kfns{3}(X1,X2)+kfns{4}(X1,X2)+...
%             kfns{5}(X1,X2)+kfns{6}(X1,X2))/6;
%     case 7
%         kfn = @(X1,X2) (kfns{1}(X1,X2)+kfns{2}(X1,X2)+kfns{3}(X1,X2)+kfns{4}(X1,X2)+...
%             kfns{5}(X1,X2)+kfns{6}(X1,X2))/6;
%
%     otherwise
%         error('kfnavg was manually programmed for 2 to 10 kfns')
% end
%
% kfn = @(X1,X2) mean(cellfun(@(fn) fn(X1,X2),kfns),3);
%}