% DATAGEN_TEST  test generation of octonion data
clear; close all force

%'random', 'Kim2011', 'Rohrer2009', '5DOF', '5DOF_interior',
%'5DOF_exterior', '5DOF_vtx', '5DOF_vtx_deleteO'

% sampleMethodList = {'5DOF_vtx_deleteO','5DOF_vtx','5DOF_exterior',...
% 	'5DOF_misFZfeatures','5DOF_interior','5DOF_vtx_deleteOz','oct_vtx'};

% sampleMethodList = {'5DOF_vtx_deleteOz','5DOF_vtx_deleteO','5DOF_vtx'};

% sampleMethodList = {'5DOF_exterior_pseudo'};

sampleMethodList = {'5DOF_interior'};

% octsubdivlist = 1:3;
octsubdivlist = 1:2;

res = 10;
nint = 2;

loadQ = true;

for i = 1:length(sampleMethodList)
	for octsubdiv = octsubdivlist
		sampleMethod = sampleMethodList{i};
		sampleType = 'data';
		
		savename = get_fname(sampleMethod,res,nint,octsubdiv);
		
		%generate data
		
		if any(cellfun(@(pattern) contains(sampleMethod,pattern),{'resolution','interior','exterior'}))
			
			loadname = get_fname(sampleMethod,res,nint,1);
			
			if (octsubdiv > 1) && (exist(loadname,'file') ~= 0) && loadQ
				%load spherical convex hull for octsubdiv == 1 if available
				load(loadname,'sphK')
				[pts,props,sphK,five,usv,Ktr] = ...
					datagen(sampleMethod,octsubdiv,sampleType,res,nint,sphK);
				
			else %generate data normally
				[pts,props,sphK,five,usv,Ktr] = ...
					datagen(sampleMethod,octsubdiv,sampleType,res,nint);
			end
		else
			[pts,props,sphK,five,usv,Ktr] = datagen(sampleMethod,octsubdiv,sampleType,res,nint);
		end
		disp(savename)
		save(savename,'pts','props','sphK','five','usv','Ktr')
		disp(' ')
	end
end


%----------------------------CODE GRAVEYARD----------------------------------
%{
	switch sampleMethod
			case {'5DOF_resolution','5DOF_interior','5DOF_exterior'}
				
				loadname = get_fname(sampleMethod,res,nint,1);
				
				if (octsubdiv > 1) && (exist(loadname,'file') ~= 0)
					%load spherical convex hull for octsubdiv == 1 if available
					load(loadname,'sphK')
					[pts,props,sphK,five,usv,Ktr] = ...
						datagen(sampleMethod,octsubdiv,sampleType,res,nint,sphK);
					
				else %generate data normally
					[pts,props,sphK,five,usv,Ktr] = ...
						datagen(sampleMethod,octsubdiv,sampleType,res,nint);
				end
			otherwise
				[pts,props,sphK,five,usv,Ktr] = datagen(sampleMethod,octsubdiv,sampleType,res,nint);
		end
		
%}