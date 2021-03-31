function pd_fix = get_pd_fix(o)
arguments
    o(:,8) double
end
% get_pd_fix get a pairwise distance matrix with the first octonion always fixed

%% setup
n = size(o,1);
% pd_fix = zeros(n,n);

C = nchoosek(1:n,2);

id1 = C(:,1);
id2 = C(:,2);

o1 = o(id1,:);
o2 = o(id2,:);

%% distance matrix
%calculate
pd_fix = GBdist4(o1,o2,32,'omega',1e-12,true);
pd_fix = squareform(pd_fix);

end

%% CODE GRAVEYARD
%{
% for i = 1:n-1
%     o1 = repmat(o(i,:),n-i,1);
%     o2 = o(i+1:n,:);
%     pd_fix(i,i+1:n) = GBdist4(o1,o2,32,'omega');
%     %text waitbar
%     if mod(i,nreps2) == 0
%         if waitbarQ
%             send(D,i);
%         end
%     end
% end


% waitbarQ = true;
% 
% %% textwaitbar setup
% D = parallel.pool.DataQueue;
% afterEach(D, @nUpdateProgress);
% nsets = n;
% ninterval = 100;
% N=nsets;
% p=1;
% reverseStr = '';
% if nsets > ninterval
% 	nreps2 = floor(nsets/ninterval);
% 	nreps = nreps2;
% else
% 	nreps2 = 1;
% 	nreps = nreps2;
% end
% 
% function nUpdateProgress(~)
% 	percentDone = 100*p/N;
% 	msg = sprintf('%3.0f', percentDone); %Don't forget this semicolon
% 	fprintf([reverseStr, msg]);
% 	reverseStr = repmat(sprintf('\b'), 1, length(msg));
% 	p = p + nreps;
% end

[row,col] = find(triu(ones(n))-diag(ones(1,n)));

%}