
function [omega_out] = GBpd(olm_pairs,olm_oct,oi,pgnum,printbool,fname)
%This function adds to a pairwise distance matrix that has already been
%created, or creates one from scratch

% INPUT: 

% 1/2. olm_pairs, olm_oct
% Existing pairwise distance matrix along with list of octonions (unsymmetrized or symmetrized is fine)

% 3. oi
% List of octonions to add to pairwise distance matrix

% 4. pgnum 
% point group number 

% 5. printbool, fname
% true if you want to output pairwise distance matrix to a text file named
% fname (include .txt in name)

% IMPORTANT LIMITING CASE: construct pd matrix from scratch. 
% olm_pairs = [], olm_oct = [], oi = relevant octonion list

% olm_pairs = []; olm_oct = []; 
% printbool = true; fname = 'little_test.txt';
% pgnum = 30; 
% 
% test = importdata('../Data/olm_octonion_list.txt',' ',1); %list of GB octonions with number of octonions in file at top
% data0 = test.data;
% oi = data0(1:3,:);

% OUTPUT 

% 1. omega_out - pairwise distance matrix 



% % Import total octonion list 
% olm_oct_path = 'oct_plus_42_339_385.txt';
% olm_oct0 = importdata(olm_oct_path);
% olm_oct = olm_oct0.data;
% 
% % Import pairwise distance matrix
% olm_pairwise_path = 'additive_42_339_385.txt'; %CHANGE
% olm_pairs0 = importdata(olm_pairwise_path);
% olm_pairs = olm_pairs0.data;
% 
% % Import octonion list to add 
% interp_path = 'interp_1_3_cubic.txt';
% oi0 = importdata(interp_path);
% oi1 = oi0.data; %first and last octonion are already in dataset above? If so, use next line.
% oi = oi1(2:end-1,:); %interpolated list of octonions to add to analysis


%% Concatenate old + new list of octonions
nold = length(olm_oct);
nnew = length(oi(:,1));
n = nold + nnew; 

%octonions go into a "data" matrix
data = zeros(n,8);
data(1:nold,:) = olm_oct;
data(nold+1:n,:) = oi;

%instantiate data structures to track pairwise distances
oct_new = zeros(n,n,16); %symmetrized octonion container
omega_new = zeros(n,n); %pairwise distance matrix
zeta_new = zeros(n,n); %U(1) angle matrix

%fill in pairwise distances which already exist
omega_new(1:nold,1:nold) = olm_pairs;

%% load crystal symmetry

symnames = load('PGnames.mat');
symops = load('PGsymops.mat');

pgname = symnames.PG_names{pgnum};
disp('loading point group:')
disp(pgname)

qpt = symops.Q{pgnum}; %choose point group symmetry 
npt = length(qpt(:,1));

%% MAIN ROUTINE: GBOM distance calculator 


k = 1;

for n1 = 1:n
    for n2 = (nold+1):n
        if n2 > n1
%         for pair_id = pair_list
%         %     disp(k)
%             
        %   
        
%             disp('entry:')
%             disp([n1 n2])
            if mod(k,100)==0
                disp(['k ',num2str(k)])
            end
%             GBO_super = data(pair_id,:);

            o1 = data(n1,1:8); %octonion 1, unnormalized
            o2 = data(n2,1:8); %octonion 2, unnormalized

            qA0 = o1(1:4); qB0 = o1(5:8); 
            qC0 = o2(1:4); qD0 = o2(5:8); 
            
            qA = qA0./norm(qA0); qB = qB0./norm(qB0); 
            qC = qC0./norm(qC0); qD = qD0./norm(qD0); 
            
            GBOM_curr = 2000; %current GBOM angle. We want to lower this value via crystal symmetry!
            omega_keep = GBOM_curr;

            min_rep_count = 0; %keep track of angles as they are minimized.
            diff = 1000; %large value
            min_rep_oct = [];
            min_rep_GBOM = [];
            min_rep_zeta = [];

            for i = 1:npt
                for j = 1:npt 
                    for m = 1:npt 
                        for l = 1:npt
                            %define symmetry operators 
                            Si = qpt(i,:);
                            Sj = qpt(j,:);
                            Sm = qpt(m,:);
                            Sl = qpt(l,:);

                            qSA = qmult(Si,qA);
                            qSB = qmult(Sj,qB);
                            qSC = qmult(Sm,qC);
                            qSD = qmult(Sl,qD);

                            %now we implement U(1) symmetry 

                            %1. (A B C'(zeta) D'(zeta))
                            zm1 = zeta_min(qSA,qSB,qSC,qSD);
                            qzm1 = [cos(zm1/2) 0 0 sin(zm1/2)];
                            qCz1 = qmult(qSC,qzm1);
                            qDz1 = qmult(qSD,qzm1);
                            w1 = 2*acos(abs(sum(qSA.*qCz1)-sum(qSB.*qDz1))/2);
                            w5 = 2*acos(abs(sum(-qSA.*qCz1)-sum(qSB.*qDz1))/2);

                            %2. (B A C'(sigma) D'(sigma))

                            sm1 = zeta_min(qSB,qSA,qSC,qSD);
                            qsm1 = [cos(sm1/2) 0 0 sin(sm1/2)];
                            qCs1 = qmult(qSC,qsm1);
                            qDs1 = qmult(qSD,qsm1);
                            w2 = 2*acos(abs(sum(qSB.*qCs1)-sum(qSA.*qDs1))/2);
                            w6 = 2*acos(abs(sum(-qSB.*qCs1)-sum(qSA.*qDs1))/2);

                            %3. (A -B C'(zeta') D'(zeta'))

                            zm2 = zeta_min(qSA,-qSB,qSC,qSD);
                            qzm2 = [cos(zm2/2) 0 0 sin(zm2/2)];
                            qCz2 = qmult(qSC,qzm2);
                            qDz2 = qmult(qSD,qzm2);
                            w3 = 2*acos(abs(sum(qSA.*qCz2)-sum(-qSB.*qDz2))/2);
                            w7 = 2*acos(abs(sum(-qSA.*qCz2)-sum(-qSB.*qDz2))/2);


                            %4. (B -A C'(sigma') D'(sigma'))

                            sm2 = zeta_min(qSB,-qSA,qSC,qSD);
                            qsm2 = [cos(sm2/2) 0 0 sin(sm2/2)];
                            qCs2 = qmult(qSC,qsm2);
                            qDs2 = qmult(qSD,qsm2);
                            w4 = 2*acos(abs(sum(qSB.*qCs2)-sum(-qSA.*qDs2))/2);
                            w8 = 2*acos(abs(sum(-qSB.*qCs2)-sum(-qSA.*qDs2))/2);
                            %store Omega values 
                            wvec = [w1 w2 w3 w4 w5 w6 w7 w8];

        %                     rad2deg(min(wvec))

        %                     omega_test1 = GBOM(-qSa,qSb,qSc,qSd);
        %                     omega_test2 = GBOM(qSa,-qSb,qSc,qSd);
        %                     omega_test3 = GBOM(qSa,qSb,qSc,qSd);
        %                     omega_test4 = GBOM(qSb,qSa,qSc,qSd);

                            octonion_pair_sym = [qSA qSB qCz1 qDz1;
                                  qSB qSA qCs1 qDs1;
                                  qSA -qSB qCz2 qDz2;
                                  qSB -qSA qCs2 qDs2;
                                  -qSA qSB qCz1 qDz1;
                                  -qSB qSA qCs1 qDs1;
                                  -qSA -qSB qCz2 qDz2;
                                  -qSB -qSA qCs2 qDs2]; %symmetrized octonion pairs;

        %                     
                            [omega_test,iwmin] = min(wvec);
                            zeta_sym = [zm1 sm1 zm2 sm2 zm1 sm1 zm2 sm2];  

                            if (omega_test) <= omega_keep+1e-5 

                                omega_keep = omega_test;
                                oct_keep = octonion_pair_sym(iwmin,:);
                                zeta_keep = zeta_sym(iwmin);
                                min_rep_count = min_rep_count+1; 
                                min_rep_GBOM(min_rep_count) = omega_keep;
        
                                if abs(omega_keep - 0) < 1e-5
                                    disp(['warning, epsilon close found for n1,n2:',num2str(n1),',',num2str(n2)])
                                end

                                if min_rep_count > 10
                                    diff = abs(omega_keep - min_rep_GBOM(min_rep_count-1));
                                    if diff < 1e-5 
        %                                 disp('break')
                                        break
                                    end
                                end


                            end

                            if diff < 1e-5
        %                         disp('break')
                                break
                            end

                        end

                        if diff < 1e-5
        %                         disp('break')
                                break
                        end
                    end

                    if diff < 1e-5
        %                         disp('break')
                                break
                    end
                end
            end
        min_rep_filter = min_rep_GBOM(min_rep_GBOM > 1e-5);
        oct_new(n1,n2,:) = oct_keep; %keep octonion
        omega_new(n1,n2) = min(min_rep_filter); %omega_keep; %keep 
        zeta_new(n1,n2) = zeta_keep;

    %     repoct_cell{k} = min_rep_oct;
    %     repomega_cell{k} = min_rep_GBOM;
    %     repzeta_cell{k} = min_rep_zeta;

        k = k+1;
        end
    end
end

%% Saving Pairwise distance matrix
% have to deal with matrix in four parts
% omega_out = omega_new;

omega_out = omega_new;
omega2 = omega_new(1:nold,nold+1:n);
omega4 = omega_new(nold+1:n,nold+1:n);

omega_out(nold+1:n,nold+1:n) = omega4 + transpose(omega4); %add last square
omega_out(nold+1:n,1:nold) = transpose(omega2);


if printbool

    fmt_cell = cell(1,length(omega_out));
    for i = 1:length(omega_out)
        fmt_cell{i} = '%6.8f ';
    end

    fID = fopen(fname,'w');
    % fprintf(fID,'oct \n');
    fprintf(fID,['#pairwise octonion distance matrix',' \n']);
    fprintf(fID,[cell2mat(fmt_cell),' \n'], omega_out');
end

end

%% ZETA MIN 

function zm = zeta_min(qA,qB,qC,qD)
%%% zeta is twist angle of U(1) symmetry (6 --> 5 DOF)
%%% GBOM angle can be analytically minimized w.r.t. zeta (EQN 25)

% [cA,sA,aA,~] = q2ax(qA);
% [cB,sB,aB,~] = q2ax(qB);
% [cC,sC,aC,~] = q2ax(qC);
% [cD,sD,aD,~] = q2ax(qD);

%quaternion dot products = cos(omega/2), omega is misorientation angle

qdot_AC = sum(qA.*qC); % dot(qA,qC);%qdot(cA,cC,sA,sC,aA,aC); 
qdot_BD = sum(qB.*qD); %dot(qB,qD);%qdot(cB,cD,sB,sD,aB,aD); 

mu_num1 = qA(4)*qC(1)-qC(4)*qA(1)+qB(4)*qD(1)-qD(4)*qB(1);
crossAC = crossp(qA(2:4),qC(2:4));
crossBD = crossp(qB(2:4),qD(2:4));

mu_arg = (mu_num1 + crossAC(3) + crossBD(3))/(qdot_AC+qdot_BD);
mu = 2*atan(mu_arg);

if mu >= 0
    zm = mu;
else
    zm = mu + 2*pi;
end

end

function Omega = GBOM(qA,qB,qC,qD)
%%%octonion dot product = cos(Omega/2), Omega is GBOM angle
%%%be careful about normalization here: see eqns (21,22)

[cA,sA,aA,~] = q2ax(qA);
[cB,sB,aB,~] = q2ax(qB);
[cC,sC,aC,~] = q2ax(qC);
[cD,sD,aD,~] = q2ax(qD);

%quaternion dot products = cos(omega/2), omega is misorientation angle

qdot_AC = qdot(cA,cC,sA,sC,aA,aC); 
qdot_BD = qdot(cB,cD,sB,sD,aB,aD); 

Omega = 2*acos((qdot_AC-qdot_BD)/2); %normalized octonions account for factor of 1/2 in arg

end


function out = qdot(cA,cB,sA,sB,aA,aB)
%%% dot product between two quaternions, expressed in terms of c,s,a
%%% out = cos(omega/2), where omega is misorientation angle 
out = cA*cB + sA*sB*sum(aA.*aB); %dot(aA,aB);

end

function out = qmult(pp,qq)
%%% multiply two quaternions p*q
p = pp(2:4); q = qq(2:4);

qr = pp(1)*qq(1)-sum(p.*q); %dot(p,q);
qi = pp(1)*q + qq(1)*p + crossp(p,q);

out = [qr qi];
end


% function out = odot(q1,q2)
% end

function [c,s,a,theta] = q2ax(qq)
%%% input: qq, quaternion with real component first
%%% output: c = cos(theta/2), s = sin(theta/2), axis a, angle theta
%%% qq = [c,s*a], where theta is angle and a is axis in axis angle pair

thr = 1e-8;
theta = 2.0 * acos(qq(1));

if ((qq(1)-0.0)<thr)
    c = 0; s = 1; 
    theta = pi;
    a = [qq(2) qq(3) qq(4)];
else 
    c = qq(1);
    ss =  sign(qq(1))/sqrt(qq(2)^2+qq(3)^2+qq(4)^2);
    a = [qq(2)*ss, qq(3)*ss, qq(4)*ss];
    s = sin(theta/2);
end

end

function c = crossp(a,b)
%%% faster implementation of matlab cross product 

c = [a(2).*b(3)-a(3).*b(2) a(3).*b(1)-a(1).*b(3) a(1).*b(2)-a(2).*b(1)];
end