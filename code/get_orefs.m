function orefs = get_orefs(K)
arguments
    K(1,1) double = 10
end
% GET_OREFS  get K reference octonions
orefs = zeros(K,8);
orefs(1,:) = get_ocubo(1,'random',[],10);
orefs(2:K,:) = get_ocubo(K-1);

end

% npts = 10000;
% % o = get_ocubo(npts);
% % o = get_octpairs(o);
% S = load('oct50000.mat','pts');
% pts = S.pts;
% o = pts(1:npts,:);
% pd = squareform(pdist(o));
% ncheck = 10;
% checkvecs = randi(npts,ncheck,K-1);
% rmse = zeros(ncheck,1);
% 
% ind = cell(1,ncheck);
% tic
% parfor i = 1:ncheck
%     checkvec = checkvecs(i,:);
%     checks = nchoosek(checkvec,2);
%     t = n2c(checks);
%     ind{i} = sub2ind(size(pd),t{:});
%     
%     %     pdtmp = pd(ind{i});
%     o = pts(checks(:,1),:);
%     o2 = pts(checks(:,2),:);
%     pdavg = mean(pdtmp)*ones(size(pdtmp));
%     %error metrics
%     e = pdavg-pdtmp; %error
%     ae = abs(e); %absolute error
%     mae = mean(ae,'all'); %mean absolute error
%     se = e.^2; %square error
%     mse = mean(se,'all'); %mean square error
%     rmse(i) = sqrt(mse); %root mean square error
% end
% toc
% [minrmse,id] = min(rmse);
% [orefID1,orefID2] = ind2sub(size(pd),ind{id});
% orefIDs = unique([orefID1,orefID2]);
% orefs = pts(orefIDs,:);
% orefs = [get_ocubo(1,'random',[],10);orefs]; %prepend the original oref

% looks like I'll have to ensemble the ensemble to get a good set of reference octonions...

%that are spaced far from each other