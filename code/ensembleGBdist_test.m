% generate some random GBs
nGB = 10;
K=10;
coords = get_ocubo(nGB); % norm sqrt(2)
coords = get_octpairs(coords); %symmetrize (using default reference GBO)
coords = normr(coords); % unit norm

% K=50 ensemble using 'omega' distance
pdO = ensembleGBdist(sqrt(2)*coords,[],K,'omega');

% K=50 ensemble using 'euclidean' distance
pdE = ensembleGBdist(sqrt(2)*coords,[],K,'euclidean');

% unsymmetrized scaled euclidean distance
distanceScalingFactor = 2;
pdInputData = distanceScalingFactor*pdist(coords).'; % sub-diagonal lower-triangular values

% true symmetrized toby distance
k = reshape((1:(nGB^2)).',nGB,nGB);
[i,j] = ind2sub(size(k),k(tril(true(size(k)),-1)));
pdTrue = GBdist4(sqrt(2)*coords(i,:),sqrt(2)*coords(j,:),32,'omega'); % should be norm sqrt(2)

% Plot
figure
plot(rad2deg(pdTrue),rad2deg(pdO),'k.'); 
axis equal tight; 
xlabel('GBdist4(''omega'') [$d_{\Omega}^{\circ}$]','interpreter','latex'); 
ylabel('ensembleGBdist($K=50$,''omega'') [$d_{\Omega}^{\circ}$]','interpreter','latex');

figure
plot(pdTrue,pdE,'k.'); 
axis equal tight; 
xlabel('GBdist4(''omega'') [$d_{\Omega}^{\circ}$]','interpreter','latex'); 
ylabel('ensembleGBdist($K=50$,''euclidean'') [$d_{\Omega}^{\circ}$]','interpreter','latex');

figure
plot(rad2deg(pdTrue),rad2deg(pdInputData),'k.'); 
axis equal tight; 
xlabel('GBdist4(''omega'') [$d_{\Omega}^{\circ}$]','interpreter','latex'); 
ylabel('2*pdist [$d_{\Omega}^{\circ}$]','interpreter','latex');