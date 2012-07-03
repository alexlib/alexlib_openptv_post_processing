% function data = target_files_to_pos_lst(directory,filelist)
% to be added - READ_RT_IS_FILES(DIRECTORY,FILELIST)
% Example:

directory = '/Users/alex/Documents/PTV/ptv_alex/img_35_try';
newdirectory = [directory,'_new'];
if ~isdir(newdirectory), mkdir(newdirectory); end

verbose = 0;

for ncam = 1:4
    basename = ['cam',int2str(ncam)];
    filelist = dir(fullfile(directory,[basename,'*_targets']));
    
    % assuming filelist is the structure obtained by dir function:
    tmp = textread(fullfile(directory,filelist(1).name));
    numRows = tmp(1);
    
    data(1).xr = tmp(2:end,2);
    data(1).yr = tmp(2:end,3);
    
    
    extpos = find(filelist(1).name == '.',1,'last');
    extpos2 = find(filelist(1).name == '_',1,'last');
    extension = filelist(1).name(extpos+1:extpos2-1);
    data(1).t = str2double(extension);
    % end
    
    
    data(length(filelist)) = data(1);
    
    for k = 2:length(filelist)
        % for i = 1:40000
        tmp = textread(fullfile(directory,filelist(k).name));
        numRows = tmp(1);
        
        data(k).xr = tmp(2:end,2);
        data(k).yr = tmp(2:end,3);
        
        data(k).t = data(k-1).t + 1;
        clear tmp
        %
    end
    
    for i = 1:length(rtis)
        data(i).t = repmat(data(i).t,length(data(i).xr),1);
    end
    pos_lst = cat(2,cat(1,data.xr),cat(1,data.yr),cat(1,data.t));
    
    param.mem = 5;
    param.dim = 2;
    param.good = 5;
    param.quiet = 0;
    
    
    result = track(pos_lst, 12, param);
    nTraj = result(end);
    
    
    traj = repmat(struct('xf',[],'yf',[],'t',[],'trajid',[]),nTraj,1);
    
    minLength = param.good;
    
    i = 0;
    for k = 1:nTraj
        id = find(result(:,end) == k);
        i = i + 1;
        trajLen(i) = length(id);
        traj(i).xf = result(id,1);
        traj(i).yf = result(id,2);
        traj(i).t = result(id,3);
        traj(i).trajid = result(id,end);
    end
    if verbose
        hf = figure;
        plot(pos_lst(:,1),pos_lst(:,2),'c.');
        hold on
        plot_2d_trajectories(traj,param.good,hf);
    end
    
    % Smooth the trajectories
    newtraj = traj;
    for i = 1:length(traj)
        newtraj(i) = link_2d_trajectories_smoothn(traj(i),10);
    end
    if verbose
        plot_2d_trajectories(newtraj,param.good);
        hold on
        plot(pos_lst(:,1),pos_lst(:,2),'rx','MarkerSize',2);
        hold off
    end
    
    
    % Reconstruct the position list
    pos_lst_1 = cat(2,cat(1,newtraj.xf),cat(1,newtraj.yf),cat(1,newtraj.t));
    pos_lst_1 = sortrows(pos_lst_1,3);
    
    param.mem = 5;
    param.dim = 2;
    param.good = 15;
    param.quiet = 0;
    
    
    result = track(pos_lst_1, 12, param);
    nTraj = result(end);
    
    traj_1 = repmat(struct('xf',[],'yf',[],'t',[],'trajid',[]),nTraj,1);
    
    minLength = param.good;
    
    i = 0;
    for k = 1:nTraj
        id = find(result(:,end) == k);
        i = i + 1;
        trajLen(i) = length(id);
        traj_1(i).xf = result(id,1);
        traj_1(i).yf = result(id,2);
        traj_1(i).t = result(id,3);
        traj_1(i).trajid = result(id,end);
    end
    if verbose
        hf = figure;
        plot(pos_lst(:,1),pos_lst(:,2),'c.','MarkerSize',2);
        hold on
        plot_2d_trajectories(traj_1,param.good,hf);
        hold off
    end
    
    
    
    %     % write targets into new directory
    result = sortrows(result,3); % sort by the frame number
    frames = unique(result(:,3));
    
    for i = 1:length(frames)
        ind = (result(:,3) == frames(i));
        x = result(ind,1);
        y = result(ind,2);
        n   = ones(size(x));
        sumg = n;
        nx = n;
        ny = n;
        np = length(n);
        pnr = 0:np-1;
        tnr = -1*ones(np,1);
        filename = fullfile(newdirectory,sprintf('%s.%d_targets',basename,frames(i)));
        fid = fopen(filename,'w');
        fprintf(1,sprintf('Writing %s\n',filename))
        fprintf(fid,'%d\n',np);
        fprintf(fid,'%4d %9.4f %9.4f %5d %5d %5d %5d %5d\n',[pnr',x,y,n,nx,ny,sumg,tnr]');
        fclose(fid);
    end
    
end






