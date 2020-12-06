function eu = qu2eu(q,epsijk)
arguments
    q(:,4) double
    epsijk(1,1) double = 1
end
% QU2EU  from quaternions to Euler angles, vectorized by SGB 2020-12-05

q0 = q(:,1);
q1 = q(:,2);
q2 = q(:,3);
q3 = q(:,4);

q03 = q0.^2+q3.^2;
q12 = q1.^2+q2.^2;
chi = sqrt(q03.*q12);

%% case 1 (if)
ids = (chi==0) & (q12==0);

npts = nnz(ids);
Zero = zeros(npts,1);
eu(ids,:) = [atan2(-2*epsijk*q0(ids,:).*q3(ids,:),q0(ids,:).^2-q3(ids,:).^2),Zero,Zero];

%% case 2 (elseif)
ids2 = ~ids & (chi==0) & (q03==0);

npts2 = nnz(ids2);
Pi = repelem(pi,npts2,1);
Zero = zeros(npts2,1);
eu(ids2,:) = [atan2(2*q1(ids2,:).*q2(ids2,:),q1(ids2,:).^2-q2(ids2,:).^2),Pi,Zero];

%% case 3 (else)
ids3 = ~ids & ~ids2;

npts3 = nnz(ids3); %#ok<NASGU>

eu(ids3,:) = [atan2((q1(ids3,:).*q3(ids3,:)-epsijk*q0(ids3,:).*q2(ids3,:)) ...
    ./chi(ids3,:),(-epsijk*q0(ids3,:).*q1(ids3,:)-q2(ids3,:).*q3(ids3,:))./chi(ids3,:)), ...
    ...
    atan2(2*chi(ids3,:),q03(ids3,:)-q12(ids3,:)),...
    ...
    atan2((epsijk*q0(ids3,:).*q2(ids3,:)+q1(ids3,:).*q3(ids3,:)) ...
    ./chi(ids3,:),(-epsijk*q0(ids3,:).*q1(ids3,:)+q2(ids3,:).*q3(ids3,:))./chi(ids3,:))];

% reduce Euler angles to definition ranges (and positive values only)

npts4 = numel(eu(eu < 0));
Pi = repelem(pi,npts4,1);
eu(eu < 0) = mod(eu(eu < 0)+100*Pi,2*Pi);

end

%% CODE GRAVEYARD
%{

% if (eu(1)<0.0)
%     eu(1) = mod(eu(1)+100.0*pi,2*pi);
% end
% if (eu(2)<0.0)
%     eu(2) = mod(eu(2)+100.0*pi,pi);
% end
% if (eu(3)<0.0)
%     eu(3) = mod(eu(3)+100.0*pi,2*pi);
% end


Original CMU Function:

q0 = qq(1);
q1 = qq(2);
q2 = qq(3);
q3 = qq(4);

q03 = q0^2+q3^2;
q12 = q1^2+q2^2;
chi = sqrt(q03*q12);

if chi==0 && q12==0
    q = [atan2(-2*epsijk*q0*q3,q0^2-q3^2),0,0];
   
elseif chi==0 && q03==0
    q = [atan2(2*q1*q2,q1^2-q2^2),pi,0];
    
else
    q = [atan2((q1*q3-epsijk*q0*q2)/chi,(-epsijk*q0*q1-q2*q3)/chi), ...
         atan2(2*chi,q03-q12),...
         atan2((epsijk*q0*q2+q1*q3)/chi,(-epsijk*q0*q1+q2*q3)/chi)];
end

% reduce Euler angles to definition ranges (and positive values only)
if (q(1)<0.0)
    q(1) = mod(q(1)+100.0*pi,2*pi);
end
if (q(2)<0.0)
    q(2) = mod(q(2)+100.0*pi,pi);
end
if (q(3)<0.0)
    q(3) = mod(q(3)+100.0*pi,2*pi);
end


%commented code from original CMU function
% if (chi==0.0)
%   if (q12==0.0)
%    if (epsijk==1)
%     Phi = 0.0;
%     phi2 = 0.0;                 % arbitrarily due to degeneracy
%     phi1 = atan2(-2.0*q0*q3,q2^2-q3^2);
%    else
%     Phi = 0.0;
%     phi2 = 0.0;                 % arbitrarily due to degeneracy
%     phi1 = atan2( 2.0*q0*q3,q2^2-q3^2);
%    end
%   else
%    if (epsijk==1)
%     Phi = pi;
%     phi2 = 0.0;                 % arbitrarily due to degeneracy
%     phi1 = atan2(2.0*q1*q2,q1^2-q2^2);
%    else
%     Phi = pi;
%     phi2 = 0.0;                 % arbitrarily due to degeneracy
%     phi1 = atan2(2.0*q1*q2,q1^2-q2^2);
%    end
%   end
% else            % this is not a special degenerate case
%   if (epsijk==1)
%     Phi = atan2( 2.0*chi, q03-q12 );
%     chi = 1/chi;
%     phi1 = atan2( (-q0*q2+q1*q3)*chi, (-q0*q1-q2*q3)*chi );
%     phi2 = atan2( (q0*q2+q1*q3)*chi, (-q0*q1+q2*q3)*chi );
%   else
%     Phi = atan2( 2.0*chi, q03-q12 );
%     chi = 1/chi;
%     phi1 = atan2( (q0*q2+q1*q3)*chi, (q0*q1-q2*q3)*chi );
%     phi2 = atan2( (-q0*q2+q1*q3)*chi, (q0*q1+q2*q3)*chi );
%   end
% end
%
% q = [ phi1, Phi, phi2 ];
%}