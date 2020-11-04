function n2cplot(pts,varargin)
% N2CPLOT  convert rows of points to cells of points, convenient for plotting.
% e.g. plot(n2c([0 0; 1 1; 2 2),'*')

d = size(pts,2);
t = num2cell(pts,1);

switch nargin
	case 2
		if ischar(varargin{1})
			opts = varargin{1};
			plotfn = [];
		else
			plotfn = varargin{1};
			opts = [];
		end
		
		if isempty(plotfn)
			switch d
				case 2
					opts = '*';
					
				case 3
					plotfn = @plot3;
					opts = '*';
			end
		end
		
	case 3
		plotfn = varargin{1};
		
		if strcmp(func2str(plotfn),'trisurf')
			K = varargin{2};
		else
			opts = varargin{2};
		end
		
	otherwise
		switch d
			case 2
				plotfn = @plot;
				opts = '*';
				
			case 3
				plotfn = @plot3;
				opts = '*';
		end
end

plotfn(t{:},opts{:})