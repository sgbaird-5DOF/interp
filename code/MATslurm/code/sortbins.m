function [Ntrim, parcombsets, jobwalltimes] = ...
    sortbins(parstruct, walltimes, NV)
arguments
    parstruct struct
    walltimes double
    NV.maxwalltime(1,1) double = 59
    NV.maxtasks(1,1) double = 500
    NV.maxtaskstot(1,1) double = 5000
    NV.walltimecushion(1,1) double = 20
end
% SORTBINS  sort parameter combinations into bins based on their walltimes
% for SLURM sbatch commands
%--------------------------------------------------------------------------
% Inputs:
%  a - a
%
% Outputs:
%  b - b
%
% Usage:
%  a = b(a);
%
% Dependencies:
%  *
%
% Notes:
%  Consider adding functionality for first binning by different memory
%  types (e.g. 'memsortQ', and maybe just a for loop through the unique
%  memory values, catenating Ntrim, parcombsets, and jobwalltimes along the
%  way. Could also try implementing the for-loop at a higher-level, e.g. in
%  WRITEPARFILE
%
% see also WRITEPARFILE
%
% Author(s): Sterling Baird
%
% Date: 2020-09-07
%--------------------------------------------------------------------------
%% unpack
maxwalltime = NV.maxwalltime;
maxtasks = NV.maxtasks;
maxtaskstot = NV.maxtaskstot;
walltimecushion = NV.walltimecushion;

%% sort parstruct by walltimes
%find sorted indices
[sortedwalltimes,walltimeidx] = sortrows(walltimes);
%rearrange parameter structure
sortedparstruct = parstruct(walltimeidx);

%% 
maxwalltime = maxwalltime/1.1;
idlist = zeros(1,maxtaskstot+1);
%----sort into max_walltime blocks----
while length(idlist) > maxtaskstot
    i = 0;
    j = 0;
    idlist = [];
    timetot_list = [];
    maxwalltime = maxwalltime*1.1;
    disp(['max_walltime = ',num2str(maxwalltime), ' minutes']);
    while i < length(sortedwalltimes)
        j = j+1;
        timetot_temp = 0;
        idlist_temp = [];
        while timetot_temp < maxwalltime
            i = i+1;
            if i == length(sortedwalltimes)+1
                break
            end
            timetot_temp = timetot_temp+sortedwalltimes(i);
            idlist_temp = [idlist_temp i];
        end
        idlist{j} = idlist_temp; %#ok<*AGROW>
        timetot_list(j) = timetot_temp;
    end
end
[timetotlist_sorted, id_temp] = sort(timetot_list);
idlist_sorted = idlist(id_temp);
bin_limits = [timetotlist_sorted(1),timetotlist_sorted(end)];
%^^^^sort into max_walltime blocks^^^^

%-----------divide into jobs----------
%{
	note: I can either favor more jobs and more equal
	calculated walltimes within a job or fewer jobs
	and more tasks assigned to each job, but with more
	spread out walltimes. The first approach would involve
	just doing an N = histcounts() followed by one or two
	while loops with the max_tasks and max_taskstot criteria.
%}
if length(timetot_list) < 1001
    method = 'more_jobs'; %'more_jobs', 'fewer_jobs'
else
    method = 'fewer_jobs';
end
switch method
    case 'more_jobs'
        N = histcounts(timetotlist_sorted,'BinLimits',bin_limits);
        if length(N) == 1
            nbins = 2;
        else
            nbins = length(N);
        end
        while max(N) > maxtasks
            nbins = nbins*2;
            N = histcounts(timetotlist_sorted,nbins,'BinLimits',bin_limits);
            disp(['max(N) = ',int2str(max(N))]);
        end
        [N, edges] = histcounts(timetotlist_sorted,nbins,'BinLimits',bin_limits);
        Nidx = find(N~=0);
        Nidxlength = length(Nidx);
        
        Ntrim = N(Nidx);
        
        %----separate for all but last bin------
        for i = 1:Nidxlength-1
            edge1 = edges(Nidx(i));
            edge2 = edges(Nidx(i)+1);
            temp_idx = (edge1 <= timetotlist_sorted) & (timetotlist_sorted < edge2);
            idx = idlist(temp_idx);
            for j = 1:length(idx)
                parcombsets{i,j} = sortedparstruct(idx{j});
            end
            if ~isempty(parcombsets{i,j})
                jobwalltimes(i) = ceil(max(timetotlist_sorted(temp_idx))+walltimecushion);
            end
        end
        
        %-------separate for last bin-----------
        i = Nidxlength;
        edge1 = edges(Nidx(i));
        edge2 = edges(Nidx(i)+1);
        temp_idx = (edge1 <= timetotlist_sorted) & (timetotlist_sorted <= edge2);
        idx = idlist(temp_idx);
        for j = 1:length(idx)
            parcombsets{i,j} = sortedparstruct(idx{j});
        end
        if ~isempty(parcombsets{i,j})
            jobwalltimes(i) = ceil(max(timetotlist_sorted(temp_idx))+walltimecushion);
        end
        
    case 'fewer_jobs'
        l = length(idlist_sorted);
        remainder = rem(l,maxtasks);
        
        first_bins = reshape(idlist_sorted(1:l-remainder),maxtasks,[]);
        last_bin = {idlist_sorted(l+1-remainder:l)};
        
        nfirst_bins = size(first_bins,2);
        bins_temp = mat2cell(first_bins',ones(1,nfirst_bins),maxtasks);
        bins = vertcat(bins_temp,last_bin);
        
        Nidxlength = length(bins);
        Ntrim = cellfun(@(list) length(list),bins);
        jtot = 0;
        for i = 1:Nidxlength
            lbins_temp = length(bins{i});
            for j = 1:lbins_temp
                parcombsets{i,j} = sortedparstruct(bins{i}{j});
            end
            jtot = jtot+j;
            jobwalltimes(i) = ceil(timetotlist_sorted(jtot)+walltimecushion);
        end
        %unfinished
end
%^^^^^^^^^^^divide into jobs^^^^^^^^^^
end