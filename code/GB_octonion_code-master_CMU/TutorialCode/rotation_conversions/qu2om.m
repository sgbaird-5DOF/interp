function om = qu2om(qq,epsijk,method)
arguments
    qq(:,4) double
    epsijk(1,1) double = 1
    method = 'new'
end
% QU2OM  from quaternions to rotation matrix

npts = size(qq,1);
% global epsijk
switch method
    case 'new'
        qbar = qq(:,1).*qq(:,1)-(qq(:,2).*qq(:,2)+qq(:,3).*qq(:,3)+qq(:,4).*qq(:,4));
        
        om = zeros(3,3,npts);
        
        fn = @(x) reshape(x,1,1,npts);
        om(1,1,:) = fn(qbar + 2.0*qq(:,2).*qq(:,2));
        om(1,2,:) = fn(2.0*(qq(:,2).*qq(:,3)-qq(:,1).*qq(:,4)));
        om(1,3,:) = fn(2.0*(qq(:,2).*qq(:,4)+qq(:,1).*qq(:,3)));
        om(2,1,:) = fn(2.0*(qq(:,3).*qq(:,2)+qq(:,1).*qq(:,4)));
        om(2,2,:) = fn(qbar + 2.0*qq(:,3).*qq(:,3));
        om(2,3,:) = fn(2.0*(qq(:,3).*qq(:,4)-qq(:,1).*qq(:,2)));
        om(3,1,:) = fn(2.0*(qq(:,4).*qq(:,2)-qq(:,1).*qq(:,3)));
        om(3,2,:) = fn(2.0*(qq(:,4).*qq(:,3)+qq(:,1).*qq(:,2)));
        om(3,3,:) = fn(qbar + 2.0*qq(:,4).*qq(:,4));
        
        thr = eps;
        om(abs(om) < thr) = 0;
        
        if epsijk ~= 1
            om = permute(om,[2 1 3]);
        end
        
    case 'original'
        qbar = qq(1)*qq(1)-(qq(2)*qq(2)+qq(3)*qq(3)+qq(4)*qq(4));
            
        om(1,1) = qbar + 2.0*qq(2)*qq(2); %om1
        om(2,2) = qbar + 2.0*qq(3)*qq(3); %om5
        om(3,3) = qbar + 2.0*qq(4)*qq(4); %om9
        
        om(1,2) = 2.0*(qq(2)*qq(3)-qq(1)*qq(4)); %om2
        om(2,3) = 2.0*(qq(3)*qq(4)-qq(1)*qq(2)); %om6
        om(3,1) = 2.0*(qq(4)*qq(2)-qq(1)*qq(3)); %om7
        om(2,1) = 2.0*(qq(3)*qq(2)+qq(1)*qq(4)); %om4
        om(3,2) = 2.0*(qq(4)*qq(3)+qq(1)*qq(2)); %om8
        om(1,3) = 2.0*(qq(2)*qq(4)+qq(1)*qq(3)); %om3
        
        if (epsijk ~= 1)
            om = transpose(om);
        end
        
        thr = 1e-8;
        for i=1:3
            for j=1:3
                if (abs(om(i,j))< thr)
                    om(i,j) = 0.0;
                end
            end
        end
end


%% CODE GRAVEYARD
%{
%
% qbar = qq(1)*qq(1)-(qq(2)*qq(2)+qq(3)*qq(3)+qq(4)*qq(4));
%
%
% q(1,1) = qbar + 2.0*qq(2)*qq(2);
% q(2,2) = qbar + 2.0*qq(3)*qq(3);
% q(3,3) = qbar + 2.0*qq(4)*qq(4);
%
% q(1,2) = 2.0*(qq(2)*qq(3)-qq(1)*qq(4));
% q(2,3) = 2.0*(qq(3)*qq(4)-qq(1)*qq(2));
% q(3,1) = 2.0*(qq(4)*qq(2)-qq(1)*qq(3));
% q(2,1) = 2.0*(qq(3)*qq(2)+qq(1)*qq(4));
% q(3,2) = 2.0*(qq(4)*qq(3)+qq(1)*qq(2));
% q(1,3) = 2.0*(qq(2)*qq(4)+qq(1)*qq(3));
%
% if (epsijk ~= 1)
%     q = transpose(q);
% end
%
% thr = 1e-8;
% for i=1:3
%   for j=1:3
%     if (abs(q(i,j))< thr)
%         q(i,j) = 0.0;
%     end
%   end
% end



% a cell-based version
        qbar = qq(:,1).*qq(:,1)-(qq(:,2).*qq(:,2)+qq(:,3).*qq(:,3)+qq(:,4).*qq(:,4));
        
        om1 = qbar + 2.0*qq(:,2).*qq(:,2);
        om2 = 2.0*(qq(:,2).*qq(:,3)-qq(:,1).*qq(:,4));
        om3 = 2.0*(qq(:,2).*qq(:,4)+qq(:,1).*qq(:,3));
        om4 = 2.0*(qq(:,3).*qq(:,2)+qq(:,1).*qq(:,4));
        om5 = qbar + 2.0*qq(:,3).*qq(:,3);
        om6 = 2.0*(qq(:,3).*qq(:,4)-qq(:,1).*qq(:,2));
        om7 = 2.0*(qq(:,4).*qq(:,2)-qq(:,1).*qq(:,3));
        om8 = 2.0*(qq(:,4).*qq(:,3)+qq(:,1).*qq(:,2));
        om9 = qbar + 2.0*qq(:,4).*qq(:,4);
        
        thr = eps;
        
        om1(abs(om1) < thr) = 0;
        om2(abs(om2) < thr) = 0;
        om3(abs(om3) < thr) = 0;
        om4(abs(om4) < thr) = 0;
        om5(abs(om5) < thr) = 0;
        om6(abs(om6) < thr) = 0;
        om7(abs(om7) < thr) = 0;
        om8(abs(om8) < thr) = 0;
        om9(abs(om9) < thr) = 0;
        
        npts = size(qq,1);
        om = cell(npts,1);
        for i = 1:npts
            om{i} = [...
                om1(i),om2(i),om3(i);
                om4(i),om5(i),om6(i);
                om7(i),om8(i),om9(i)]; %I'm sure there's a better way than a for loop
            if epsijk ~= 1
                om{i} = transpose(om{i});
            end
        end

%}