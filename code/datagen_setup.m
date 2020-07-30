function S = datagen_setup(sampleMethod,opts,loadQ)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description:
% 
% Inputs:
%
% Outputs:
%
% Dependencies:
%		get_fname.m
%
%		datagen.m
%
%		var_names.m
%--------------------------------------------------------------------------
% res = opts.res;
% nint = opts.nint;
% octsubdiv = opts.octsubdiv;
% ocuboOpts = opts.ocuboOpts;
fname = get_fname(sampleMethod,opts);

if loadQ && (exist(fname,'file') ~= 0)
	disp(fname)
	S = load(fname,'pts','props','sphK','fname','five','Ktr','opts','sampleMethod','usv');
	
	load_type = 'evalc'; %'evalc', 'manual'
	switch load_type
		case 'evalc'
			vars = fields(S);
			for i = 1:length(vars)
				var = vars{i};
				temp = S.(var); %#ok<NASGU> %temporary value of vName
				evalc([var '= temp']); %assign temp value to the field name
			end
		case 'manual'
			pts = S.pts;
			sphK = S.sphK;
			props = S.props;
			five = S.five;
			usv = S.usv;
			Ktr = S.Ktr;
	end
	if ~contains(sampleMethod,'pseudo') && any(cellfun(@isempty,{pts,sphK,props}))
		computeQ = true;
	else
		computeQ = false;
	end
else
	computeQ = true;
end

if computeQ
	%assumption is that meshprops is a 1D array 2020-07-09
	[pts,props,sphK,five,usv,Ktr] = datagen(sampleMethod,'data',opts); %property values for validating sph. bary interp.
	if contains(sampleMethod,'pseudo')
		Ktr = [];
	end
	fpath = fullfile('data',fname);
   save(fpath,'pts','props','sphK','five','usv','Ktr','opts','sampleMethod','fname')
	disp(fname)
end

d = size(pts,2);
nmeshpts = length(pts);
disp(['# vertices: ',int2str(nmeshpts)])

if isempty(usv)
	S = var_names(pts,props,sphK,fname,five,Ktr,opts,sampleMethod);
else
	S = var_names(pts,props,sphK,fname,five,Ktr,opts,sampleMethod,usv);
end

end
