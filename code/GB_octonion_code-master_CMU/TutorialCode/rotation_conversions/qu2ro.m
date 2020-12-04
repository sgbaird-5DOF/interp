% from quaternions to Rodrigues vector

function q = qu2ro(qq,epsijk)
arguments
    qq(:,4) double
    epsijk(1,1) double = 1
end

% global epsijk

thr = 1e-10;

q(:,1:3) = qq(:,2:4);
q(:,4) = 0.0;

ia = qq(:,1) < thr;
q(ia,4) = Inf;

iatemp = find(ia);
ib = setdiff(1:size(qq,1),iatemp);

if ~isempty(ib)
    
    s = sqrt(sum(q(ib,1:3).*q(ib,1:3),2));
    
    ic = find(s < thr);
    if ~isempty(ic)
        q(ic,:) = repmat([0.0, 0.0, epsijk, 0.0],length(ic),1);
    end
    
    id = find(s >= thr);
    if ~isempty(id)
        t = tan(acos(qq(id,1)));
        q(id,:) = [q(id,1)./s(id), q(id,2)./s(id), q(id,3)./s(id), t];
    end
    
    q(abs(q(ib,1)) < thr,:) = 0;
    
%     q(:,4) = [];
    
else
%     q(:,4) = [];
end

%------------------------------CODE GRAVEYARD------------------------------
%{
% if (qq(1)<thr)
%   q(4) = Inf;
%   return
% end

% % set values very close to 0 as 0
% if (abs(q(1))-0)<thr
%     q(1)=0;
% elseif (abs(q(2))-0)<thr
%     q(2)=0;
% elseif (abs(q(3))-0)<thr
%     q(3)=0;
% elseif (abs(q(4))-0)<thr
%     q(4)=0;
% end

    for i = 1:4
        ie = find(abs(q(ib,1)-0) < thr);
        if ~isempty(ie)
            q(ie,:) = 0;
        end
    end

%}