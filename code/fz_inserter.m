function [aa3, disorientation, OM_sst, r_sst, th_sst, sizes] = fz_inserter(super1,super2,input_str)

%Compute misorientations in the fundamental zone for two lists of
%orientations
%INPUT: two arbitrary list of orientations expressed as
% super1 and super 2 :: N x 9 arrays corresponding to lists of orientations. 
% Each 1x9 row of array is flattened orientation matrix: see example below 
%EX: orientation matrix A = [1 1 1; 1 -1 0; 1 1 -2] should be expressed by normalizing rows of
% matrix (unnormalized rows are directions [h k l] in miller indices)
% An = normr(A)
% Rotation matrix should then be reshaped using reshape(An,[1 9])
% valid input: [1/sqrt(3) 1/sqrt(3) 1/sqrt(3) 1/sqrt(2) -1/sqrt(2) 0 1/sqrt(6) 1/sqrt(6)
% -2 /sqrt(6)]
% 
%input_str:: optional argument. If input_str = 'on', plot will be generated
%of disorientation angle distribution and SST projection of axes

%OUTPUT: Rotational D.O.F of misorientations represented in fundamental zone: 
%disorientation angle corresponds to angle minimizing rotation accounting for cubic symmetry operators
%rotation axis placed in sst (choice of sst is arbitrary)

% aa3:: rotations in fundamental zone as axis angle pair. 
% Fifth index is a size used to represent frequency of the rotation, but
% can be ignored. 
% disorientation:: list of disorientation angles (can be used to generate a
% mackenzie distribution, e.g.)

% omtest :: N x 9 list of orientation matrices in flattened format
% described above
% rtest :: polar coordinate of axes for stereographic projection
% thtest :: polar coordinate of axes for stereographic projection 
% sizes :: relative sizes for plotting purposes (ignore this unless you
% need to represent different frequencies)


[omlist1, sz1] = omcluster(super1); 
[omlist2, sz2] = omcluster(super2);

if nargin == 2
    input_str = 'off';
end 
%% Define Symmetries & Variables 
% qu = super_cu{120,7}(:,1:4);
% input_str = 'on';

% ct = repmat([1/sqrt(3) 1/sqrt(3) 1/sqrt(3) 60],[4 1]);
% ct_om = vrrotvec2mat([1/sqrt(3) 1/sqrt(3) 1/sqrt(3) pi/3]);
% ct_omlist = reshape(ct_om,[1 9]);
% om_ex1 = vrrotvec2mat([1/sqrt(2) 1/sqrt(2) 0 pi/2]);
% om_ex2 = ct_om*om_ex1;
% % omlist1 = reshape(om_ex1,[1 9]);
% % omlist2 = reshape(om_ex2,[1 9]);
% 
% omlist1 = super_cu{185,5}(:,5:13);
% omlist2 = super_nb{185,5}(:,5:13);
% qu = ax2qulist(fzqu3(om2qulist2(om)));
% qu2 = ax2qulist(fzqu3(om2qulist2(om2)));
% om = reshape(super_cu{2,1}(:,5:13),[3 3])';
% om2 = reshape(super_nb{2,1}(:,5:13),[3 3])';
% qu = om2qulist2(super_cu{164,5}(:,5:13));
% qu2 = om2qulist2(super_nb{164,5}(:,5:13));


% qu = ax2qulist(fzqu3(om2qulist2(super_cu{185,5}(:,5:13))));
% qu2 = ax2qulist(fzqu3(om2qulist2(super_cu{185,5}(:,5:13))));

n = length(omlist1(:,1));
n2 = length(omlist2(:,1));
disorientation = zeros(1,n*n2);
axes_sst = zeros(n*n2,3);
sizes = zeros(1,n*n2);

%cubic symmetry operators as quaternions.
%this representation differs from commonly used representation,
%in that the angular component comes first, follow by the axes components. 
q_cubic_sym_flip = ... 
[1 0 0 0; ...
 sqrt(2)/2 0 0 sqrt(2)/2; ... 
 sqrt(2)/2 sqrt(2)/2 0 0; ... 
 0 0 0 1; ... 
 0.5 0.5 -0.5 0.5; ... 
 sqrt(2)/2 0 0 -sqrt(2)/2; ... 
 0 0 -sqrt(2)/2 sqrt(2)/2; ... 
 0.5 0.5 0.5 -0.5; ... 
 0 1 0 0;
 0 0 sqrt(2)/2 sqrt(2)/2; ... 
 0 sqrt(2)/2 0 sqrt(2)/2; ... 
 0.5 0.5 -0.5 -0.5; ... 
 0 -sqrt(2)/2 0 sqrt(2)/2; ... 
 sqrt(2)/2 -sqrt(2)/2 0 0; ...
 0 0 1 0; ... 
 0.5 -0.5 -0.5 -0.5; ... 
 0 -sqrt(2)/2 sqrt(2)/2 0; ... 
 0.5 -0.5 0.5 0.5; ... 
 sqrt(2)/2 0 -sqrt(2)/2 0; ... 
 0.5 -0.5 0.5 -0.5; ... 
 sqrt(2)/2 0 sqrt(2)/2 0; ... 
 0.5 -0.5 -0.5 0.5; ... 
 0 sqrt(2)/2 sqrt(2)/2 0; ... 
 0.5 0.5 0.5 0.5];

%putting in format [qx qy qz q0]
q_cubic_sym = zeros(size(q_cubic_sym_flip));
q_cubic_sym(:,1:3) = q_cubic_sym_flip(:,2:4);
q_cubic_sym(:,4) = q_cubic_sym_flip(:,1);

%convert nx4 quaternion list to nx9 om list
%note that the function qu2omlist is based on Degraef's set of conversion
%functions
om_cubic_sym_list = qu2omlist(q_cubic_sym); %9x24 matrix
om_cubic_sym = om_reconstruct(om_cubic_sym_list); %cell with 24 3x3 matrices

%% Placing in FZ
k = 0;
for j = 1:n
    k = k+1;
    om_1 = omlist1(j,:);
    O1 = reshape(om_1,[3 3]);
    sz_1 = sz1(j);
    for m = 1:n2
        k = k+1;
        om_2 = omlist2(m,:);
        O2 = reshape(om_2,[3 3]);
        sz_2 = sz2(m);
        % Define misorientation g as rotation required to bring ori 1 into ori 2: O2 = misorientation*O1

        mis_init = O2*O1';
        mis_aa_init = vrrotmat2vec(mis_init);
%         mis_init = qmult((qu_1),qinv((qu_2)));
        curr_min = 2*pi;
        ax_choose = (mis_aa_init(1:3));
        for i = 1:24
            O1_sym_var = om_cubic_sym{1,i}*O1;
%             O1_sym_variants{1,i} = O1_sym_var;

            %calculate misorientations
            g_sym_var_om = O2*O1_sym_var';
            g_sym_var_aa = vrrotmat2vec(g_sym_var_om);
            
            g_4 = abs(g_sym_var_aa(4)); %we seek to maximize fourth component of quaternion

            if g_4 < curr_min
                curr_min = g_4;
                ax_choose = (g_sym_var_aa(1:3));
            end  
        end
     %       g_sym_variants_qu(i,:) = qg_sym_var;
    

    %postom variants control for conversion math (just using orientation
    %matrices)
%     postom_variants_qulist = om2qulist(g_sym_variants_omlist);
%     postom_variants_alist = qu2alist(postom_variants_qulist);

    disorientation_angle = rad2deg((curr_min));
    sortedax = sort(abs(ax_choose));
    curr_sz = sz_1+sz_2;
% 
%     [disorientation_angle, dis_index] = min(rad2deg(2*acos((abs(g_sym_variants_qu(:,4))))));%min(postom_variants_alist(:,4));
    disorientation(k) = disorientation_angle;
    sizes(k) = curr_sz;
%     
    %angle placed into the SST convenient for my plotting routine
%     sortedax = sort(abs(g_sym_variants_qu(dis_index,1:3)));
    axes_ordered = sortedax(:,[2 1 3]);
    axes_sst(k,:) = axes_ordered./sqrt(sum(abs(axes_ordered).^2,2));
    [th_dis,r_dis] = stereo(axes_sst);
    th_sst = th_dis;
    r_sst = r_dis;
    end
    
end


d2 = disorientation(disorientation ~= 0);
a2 = axes_sst(disorientation ~= 0,:);
s2 = sizes(disorientation ~= 0);

%ignore above section for now!
aa3 = zeros(length(d2),5);
aa3(:,5) = s2;
aa3(:,4) = d2;
aa3(:,1:3) = a2;
aa3 = abs(sortrows(-aa3,5));

sizes_sum = sum(aa3(:,5));
aa3(:,5) = aa3(:,5)/10000;%/sizes_sum;

aa4 = zeros(size(aa3(:,1:4)));
aa4(:,1:3) = aa3(:,1:3);
aa4(:,4) = deg2rad(aa3(:,4));
OM_sst = zeros(length(aa3(:,4)),10);
for i = 1:length(aa4(:,4))
    OM_sst(i,1:9) = reshape(vrrotvec2mat(aa4(i,:)),[1 9]);
    OM_sst(i,10) = aa3(i,5);
end


%% Plotting

if ~strcmp(input_str,'off')
%     figure
%     histogram(disorientation(disorientation > 0.5),100);
%     title(['Disorientation Statistics for Randomly Sampled Cubic Orientations, sample size = ',...
%         num2str(n*n2)])
%     xlabel('Disorientation Angle (°)')
%     ylabel('Count')

    figure
    pax = polaraxes;
    fz_edge_p = fzedge_find_m(linspace(0,pi/4));
    polarplot(pax,fz_edge_p(:,1),fz_edge_p(:,2),'k')
    hold on;
    p = polarscatter(th_sst(disorientation > 0.5),r_sst(disorientation > 0.5),sizes(disorientation > 0.5)/100, ...
disorientation(disorientation > 0.5),'filled');
    hold on;
    p_ks = polarscatter(0.7854, 0.1276, 1000, 42.848,'p','filled');
    hold off
    pax.ThetaLim = [0 45];
    pax.RLim = [0 fzedge_find(pi/4)];
%     title(['Axis Statistics for Randomly Sampled Cubic Orientations, sample size = ',...
%         num2str(n*n2)])
    title('OR statistics')
    alpha(p,0.5);
    alpha(p_ks,0.3);
    c = colorbar;
    c.Label.String = 'Disorientation Angle (°)';
    caxis([0 62])

end


end

function g = axis_insert(g_un,s)
%Insert axis into FZ

if nargin == 1
    s = 'off';
end 

g_n = abs(g_un)./sqrt(sum(abs(g_un.^2),2)); %normalize rows, take absolute value
g1 = zeros(length(g_n(:,1)),3);

for i = 1:length(g_n(:,1))
    g1(i,:) = (sort((g_n(i,:))));
end

g = g1(:,[2 1 3]); %this is the permutation that works to put axes in FZ!

if strcmp(s,'on')
    [th_g,r_g] = stereo(g);
    figure
    pax = polaraxes;
    fz_edge_p = fzedge_find_m(linspace(0,pi/4));
    polarplot(pax,fz_edge_p(:,1),fz_edge_p(:,2),'k')
    hold on;
    polarscatter(th_g,r_g,'filled')
    hold off
    pax.ThetaLim = [0 45];
    pax.RLim = [0 fzedge_find(pi/4)];
end


end
% take similar elements of an omlist and only considers them once. 
function [omout, szout, th_out, r_out] = omcluster(supertest, input_str)
% supertest = super_cu{185,6};
% qutest = supertest(:,1:4);

if nargin == 1
    input_str = 'off';
end
if length(supertest(1,:)) == 14
    omtest = supertest(:,5:13); 
    sztest = supertest(:,14);
elseif length(supertest(1,:)) == 10
    omtest = supertest(:,1:9);
    sztest = supertest(:,10);
else 
    omtest = supertest;
    sztest = 50000;
end

n = length(omtest(:,1));
omnew = zeros(size(omtest));

%sort each basis vector 
for i = 1:n
    omline = abs(omtest(i,:));
    newomline = [sort(omline(1:3)) sort(omline(4:6)) sort(omline(7:9))];
    omnew(i,:) = newomline;
end 

%for each orientation line, compute difference between other sorted orientation lines
%choose some cutoff below which we call an orientation "similar enough" to
%consider the same. Output unique OM's only to new list.
omold = omnew;
szold = sztest;
indold = 1:length(sztest);

omuniq = zeros(size(omnew));
szuniq = zeros(size(sztest));
induniq = zeros(size(sztest));

k = 0;
while ~isempty(szold)
    
    k = k+1;
    % take peak with most atoms
    [~,ind_max] = max(szold);
    
    %find similar peaks
    cutoff = 0.15;
    omdiff = sum(abs(omold - omold(ind_max,:)),2);
    ind_to_combine = find(omdiff < cutoff);
    ind_to_keep = setdiff(1:length(szold),ind_to_combine);
    %adding sizes of peaks that will be collapsed
    sz_sum = sum(szold(ind_to_combine));

    %output unique orientation and size sum, and index in original model
    omuniq(k,:) = omold(ind_max,:);
    szuniq(k,:) = sz_sum;
    induniq(k,:) = indold(ind_max);
    

    %delete current equivalent orientations from structure
    omold(ind_to_combine,:) = [];
    szold(ind_to_combine) = [];
    indold(ind_to_combine) = [];
end 

induniq(szuniq == 0,:) = [];
omuniq(szuniq == 0,:) = [];
szuniq(szuniq == 0) = [];

omuniq(:,end+1) = szuniq;
omuniq(:,end+1) = induniq;

omuniq = abs(sortrows(-omuniq,10));

% omout = omuniq(:,1:9);
szout = omuniq(:,10);
indout = omuniq(:,11);
omout = omtest(indout,:);

[th_out, r_out] = triple_stereo(omout);

if ~strcmp(input_str,'off')
ratio_sz = szout/sum(szout);
    for i = 1:length(omout(:,1))
        disp(['Peak ',num2str(i),' : ', num2str(round(ratio_sz(i)*100,2)), '% of atoms in layer'])
        disp(reshape(omout(i,:),[3 3])')
    end
end

end

function omlist = qu2omlist(qq)
%quaternion matrix to orientation matrix cell



qbar = qq(:,4).*qq(:,4)-(qq(:,1).*qq(:,1)+qq(:,2).*qq(:,2)+qq(:,3).*qq(:,3));

q = cell(3);

q{1,1} = qbar + 2.0*qq(:,1).*qq(:,1);
q{2,2} = qbar + 2.0*qq(:,2).*qq(:,2);
q{3,3} = qbar + 2.0*qq(:,3).*qq(:,3);

q{1,2} = 2.0*(qq(:,1).*qq(:,2)-qq(:,4).*qq(:,3));
q{2,3} = 2.0*(qq(:,2).*qq(:,3)-qq(:,4).*qq(:,1));
q{3,1} = 2.0*(qq(:,3).*qq(:,1)-qq(:,4).*qq(:,2));
q{2,1} = 2.0*(qq(:,2).*qq(:,1)+qq(:,4).*qq(:,3));
q{3,2} = 2.0*(qq(:,3).*qq(:,2)+qq(:,4).*qq(:,1));
q{1,3} = 2.0*(qq(:,1).*qq(:,3)+qq(:,4).*qq(:,2));

omlist = [q{1,:} q{2,:} q{3,:}];

thr = 1e-7;
for i = 1:length(omlist(:,1))
    omi = omlist(i,:);
    omi(abs(omi) < thr) = 0.0;
    omlist(i,:) = omi;
end 


end

function omcell = om_reconstruct(omlist)
%takes n x 9 list, reconstructs rows as 9 linear indices of matrix 

dim = size(omlist);
if dim(1) == 1 && dim(2) == 9
    omcell = zeros(3);
    omcell(1,:) = omlist(1:3);
    omcell(2,:) = omlist(4:6);
    omcell(3,:) = omlist(7:9);
else
    n = length(omlist(:,1));
    omcell = cell(1,n);

    for i = 1:n
        omtemp = zeros(3);
        omtemp(1,:) = omlist(i,1:3);
        omtemp(2,:) = omlist(i,4:6);
        omtemp(3,:) = omlist(i,7:9);
        omcell{1,i} = omtemp;
    end 

end

end 

function [th,r] = triple_stereo(test_omp1)

if length(test_omp1(1,:)) >= 9
    o_x1 = test_omp1(:,1:3);
    o_y1 = test_omp1(:,4:6);
    o_z1 = test_omp1(:,7:9);

    th = zeros(size(o_x1));
    r = zeros(size(o_x1));

    [thx,rx] = stereo(axis_insert(o_x1));
    [thy,ry] = stereo(axis_insert(o_y1));
    [thz,rz] = stereo(axis_insert(o_z1));

    th(:,1) = thx;
    th(:,2) = thy;
    th(:,3) = thz;

    r(:,1) = rx;
    r(:,2) = ry;
    r(:,3) = rz;
    
elseif length(test_omp1(1,:)) < 9 %quaternion
    o_x1 = test_omp1(:,1:3);
    [th,r] = stereo(axis_insert(o_x1));

end

end

function [Theta,R] = stereo(aa)
%Stereographic projection of rotation axes
% n x 3 matrix, columns are unit vectors [ax ay az]
n = length(aa(:,1));
% r = zeros(1,n);
R = zeros(1,n);
Theta = zeros(1,n);
%Phi = zeros(1,n);

for i = 1:n
    a = aa(i,1:3);
    r = norm(a);
    if a(1) < 0
        Theta(i) = atan(a(2)/a(1)) + pi;
    else 
        Theta(i) = atan(a(2)/a(1));
    end 
    phi = acos(a(3)/r);
    R(i) = sin(pi-phi)/(1-cos(pi-phi));
end 

R(isnan(R)) = 0;
Theta(isnan(Theta)) = 0;
% 
% figure
% pax = polaraxes;
% polarscatter(Theta,R)
% title('stereographic projection of rotation axes')
% pax.ThetaLim = [0 90];
% pax.RLim = [0 2];
% 
% figure
% polarhistogram(Theta)
% 
% figure
% histogram(R)
% figure
% histogram2(Theta,R)
% xlabel('Theta')
% ylabel('R')

end

function fz_edge = fzedge_find_m(th_list)
%Finds edge of standard triangle for specified quadruplet of r / theta
%pairs. 
%

fz_edge = zeros(length(th_list),2);
fz_shape = stereomat([0 0 1; 1 1 1; 1 0 1; -1 0 -1]);

r_new = (fz_shape(4,2) + fz_shape(3,2))/2.;
r_old = fz_shape(2,2);

theta1_span = fz_shape(2,1);
theta2_span = pi-acos((r_old^2 - r_new^2 - (r_new-fz_shape(3,2))^2)/(2*r_new*(r_new-fz_shape(3,2)))); %law of cosines, draw picture
%value above is insensitive to large angles 

th_ratio = theta1_span/theta2_span; 
% now we express the new curve in terms of polar coordinates of the first
fz_edge_pre = th_list/th_ratio;
fz_edge_x = r_new*cos(fz_edge_pre)-(r_new-fz_shape(3,2));%(r_new-r_old*cos(theta1_span));
fz_edge_y = r_new*sin(fz_edge_pre);
fz_edge_r = sqrt(fz_edge_x.^2 + fz_edge_y.^2);
fz_edge_th = atan(fz_edge_y ./ fz_edge_x);

fz_edge(:,1) = fz_edge_th;
fz_edge(:,2) = fz_edge_r;

end

function thr = stereomat(aa)
%Stereographic projection of rotation axes, theta column and r column
% n x 3 matrix, columns are unit vectors [ax ay az]
% assumes only positive axes values
n = length(aa(:,1));
thr = zeros(n,2);

ax = aa(:,1);
ay = aa(:,2);
az = aa(:,3);
r = sqrt(ax.^2 + ay.^2 + az.^2);
thr(:,1) = atan(ay./ax); %theta
phi = acos(az./r);
thr(:,2) = sin(pi*ones(length(phi),1)-phi)./(ones(length(phi),1)-cos(pi*ones(length(phi),1)-phi)); %R

% set [100] projection to origin
thr(isnan(thr)) = 0;

% 
% Theta = thr(:,1);
% R = thr(:,2);
% 
% figure
% pax = polaraxes;
% polarscatter(Theta,R)
% title('stereographic projection of rotation axes')
% pax.ThetaLim = [0 90];
% pax.RLim = [0 2];
% 
% figure
% polarhistogram(Theta)
% 
% figure
% histogram(R)
% figure
% histogram2(Theta,R)
% xlabel('Theta')
% ylabel('R')

end
function fz_edge_r = fzedge_find(th_list)
%Finds edge of standard triangle for specified quadruplet of r / theta
%pairs. 
%


fz_shape = stereomat([0 0 1; 1 1 1; 1 0 1; -1 0 -1]);

r_new = (fz_shape(4,2) + fz_shape(3,2))/2.;
r_old = fz_shape(2,2);

theta1_span = fz_shape(2,1);
theta2_span = pi-acos((r_old^2 - r_new^2 - (r_new-fz_shape(3,2))^2)/(2*r_new*(r_new-fz_shape(3,2)))); %law of cosines, draw picture
%value above is insensitive to large angles 

th_ratio = theta1_span/theta2_span; 
% now we express the new curve in terms of polar coordinates of the first
fz_edge_pre = th_list/th_ratio;
fz_edge_x = r_new*cos(fz_edge_pre)-(r_new-fz_shape(3,2));%(r_new-r_old*cos(theta1_span));
fz_edge_y = r_new*sin(fz_edge_pre);
fz_edge_r = sqrt(fz_edge_x.^2 + fz_edge_y.^2);
%fz_edge_th = atan(fz_edge_y ./ fz_edge_x);


end
