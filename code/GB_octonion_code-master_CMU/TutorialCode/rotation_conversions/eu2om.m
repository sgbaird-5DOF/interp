function om = eu2om(eu)
arguments
    eu(:,3)
end
% EU2OM  from Euler angles to rotation matrix
thr = 1e-10;
npts = size(eu,1);

c1 = cos(eu(:,1));
c3 = cos(eu(:,3));
c2  = cos(eu(:,2));

s1 = sin(eu(:,1));
s3 = sin(eu(:,3));
s2  = sin(eu(:,2));

om1 = c1.*c3-s1.*c2.*s3;
om2 = s1.*c3+c1.*c2.*s3;
om3 = s2.*s3;
om4 = -c1.*s3-s1.*c2.*c3;
om5 = -s1.*s3+c1.*c2.*c3;
om6 = s2.*c3;
om7 = s1.*s2;
om8 = -c1.*s2;
om9 = c2;

om1(abs(om1) < thr) = 0;
om2(abs(om2) < thr) = 0;
om3(abs(om3) < thr) = 0;
om4(abs(om4) < thr) = 0;
om5(abs(om5) < thr) = 0;
om6(abs(om6) < thr) = 0;
om7(abs(om7) < thr) = 0;
om8(abs(om8) < thr) = 0;
om9(abs(om9) < thr) = 0;

om = cell(npts,1);
for i = 1:npts
    om{i} = [...
        om1(i),om2(i),om3(i);
        om4(i),om5(i),om6(i);
        om7(i),om8(i),om9(i)]; %I'm sure there's a better way than a for loop
end

end

%% CODE GRAVEYARD
%{
%-------- original ---------
% from Euler angles to rotation matrix

function q = eu2om(eu)

thr = 1e-10;

c1 = cos(eu(1));
c3 = cos(eu(3));
c2  = cos(eu(2));

s1 = sin(eu(1));
s3 = sin(eu(3));
s2  = sin(eu(2));

q = [ c1*c3-s1*c2*s3,  s1*c3+c1*c2*s3, s2*s3; ...
     -c1*s3-s1*c2*c3, -s1*s3+c1*c2*c3, s2*c3; ...
           s1*s2    ,      -c1*s2    ,  c2   ];

for i=1:3
  for j=1:3
    if (abs(q(i,j))< thr)
        q(i,j) = 0.0;
    end
  end
end
%------------


om = [...
    om1,om2,om3;
    om4,om5,om6;
    om7,om8,om9];
%}