function q = qu2ax(qq)
arguments
    qq(:,4) double
end
% QU2AX Convert from quaternions to axis-angle pair
%--------------------------------------------------------------------------
% Date: 2020-08-21
% 
% Inputs:
%  qq - rows of quaternions
%
% Outputs:
%  q - rows of axis-angle pairs (first 3 columns, axis, last column angle)
%
% Usage:
%  q = qu2ax(qq);
%
% Notes:
%  Vectorized by SGB (2020-08-21)
%--------------------------------------------------------------------------
thr = 1e-10;
omega = 2*acos(qq(:,1));

ids = omega < thr;
if any(ids)
    nids = sum(ids);
    q(ids,:) = repmat([0, 0, 1, 0],nids,1);
end

ids2 = ~ids & (qq(:,1) < thr);
if any(ids2)
    nids2 = sum(ids2);
    q(ids2,:) = [qq(ids2,2), qq(ids2,3), qq(ids2,4), repelem(pi,nids2,1)];
end

ids3 = ~ids & (qq(:,1) >= thr);
if any(ids3)
    s = (qq(ids3,1)./abs(qq(ids3,1)))./rssq(qq(ids3,2:4),2);
    q(ids3,:) = [qq(ids3,2:4).*[s s s], omega(ids3,:)];
end

% set values very close to 0 as 0
q(abs(q) < thr) = 0; %was missing abs() 2020-09-01 SGB

end



%---------------------------CODE GRAVEYARD---------------------------------
%{
% s =  (qq(~ids2,1)/abs(qq(~ids2,1)))/sqrt(qq(~ids2,2).^2+qq(~ids2,3).^2+qq(~ids2,4).^2);

%------Original CMU group function------
function q = qu2ax(qq)

thr = 1e-8;
omega = 2.0 * acos(qq(1));

if ((omega-0.0)<thr)
  q = [0.0, 0.0, 1.0, 0.0];
else
    
  if ((qq(1)-0.0)<thr)
      q = [qq(2), qq(3), qq(4), pi];
  else
    s =  (qq(1)/abs(qq(1)))/sqrt(qq(2)^2+qq(3)^2+qq(4)^2);
    q = [ qq(2)*s, qq(3)*s, qq(4)*s, omega];
  end
 
end

% set values very close to 0 as 0
if (abs(q(1))-0)<thr
    q(1)=0;
elseif (abs(q(2))-0)<thr
    q(2)=0;
elseif (abs(q(3))-0)<thr
    q(3)=0;
end

    s = (qq(~ids2,1)./abs(qq(~ids2,1)))./sqrt(qq(~ids2,2).^2+qq(~ids2,3).^2+qq(~ids2,4).^2);
%}