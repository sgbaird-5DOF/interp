function meshpts = meshgen(meshMethod)
% MESHGEN  generate octonions randomly or from literature (deprecated)
switch meshMethod
	case 'Rohrer2009'
		step = 10; %degrees
		phi1_list = 0:step:90;
		cap_phi_list = cos(phi1_list);
		phi2_list = phi1_list;
		theta_list = cap_phi_list;
		Phi_list = phi1_list;
		
		meshpts = allcomb(phi1_list,cap_phi_list,phi2_list,theta_list,Phi_list);
		
	case 'Kim2011'
		meshTable = readtable('Kim2011_FeGBEnergy.txt','HeaderLines',9,'ReadVariableNames',true);
		varNames = meshTable.Properties.VariableNames; %Euler misorientation angles (deg), polar & azimuth inclinations (deg), GB Energy (mJ/m^2)
		meshptsTemp = table2array(meshTable);
		meshpts = meshptsTemp(:,1:5);
		
	case 'hypersphereGrid'
		%% quaternion setup
		
		n = 3; %sidelength of the n-dimensional grid, odd is preferred
		d = 7; %dimension of hypersphere
		subdivType = 'hypersphere'; %'orthant', 'hypersphere'
		
		meshpts = hypersphereSetup(n,d,subdivType);
		
	case '5DOF'
		featureType = 'resolution';
		
	case '5DOF_vtx'
		featureType = 'vertices';
		
	case 'Olmsted2004'
		%add folders to path
		fnlist = {'olm_octonion_list.txt','olm_properties.txt'};
		filepathgen = fullfile('**',fnlist);
		for i = 1:length(filepathgen)
			file = dir(filepathgen{i});
			name = file(1).name;
			folder = file(1).folder;
			filepath{i} = fullfile(folder,name);
		end
		
		meshpts = readmatrix(filepath{1},'NumHeaderLines',1,'Delimiter',' ');
		meshpts = meshpts(:,1:8); %octonion representation, gets rid of a random NaN column..
		meshprops = readmatrix(filepath{2},'NumHeaderLines',1,'Delimiter',' '); %properties
		
		meshvals = meshprops(:,1);
		
end

%5DOF cases
switch meshMethod
	case {'5DOF','5DOF_vtx'}
		resdegrees = 10; %does nothing for featureType == '5DOF_vtx'
		nint = 1;
		[five,sept,o] = mesh5DOF(resdegrees,nint,featureType);
		meshpts = sept;
end

end
