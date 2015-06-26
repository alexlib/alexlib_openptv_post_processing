function varargout = TrajPoint2HDF_v8(hfile,varargin)
% TRAJPOINT2HDF_v8 Reads trajPoint.* files into Matlab and saves them 
% into HDF format.
%           TRAJPOINT2HDF_v35(HFILE) will open GUI to read in specific
% set of trajPoint.* files and convert them into multiple file series
% the the base name HFILE. The first and the most important file will be
% called HFILE.HDF, and it will be also distrubuted in many files, titled
% HFILE.HDF.1, HDFILE.HDF.2 ... HFILE.HDF.NUM_OF_FRAMES
%
%         TRAJPOINT2HDF(HFILE,FILEPATH,FIRSTFILE,LASTFILE)
% Will read trajPoint.* files from trajPoint.FIRSTFILE to
% trajPoint.LASTFILE in the directory FILEPATH
%
% Example:
%
%           trajPoint2HDF_v35('traj.hdf','e:\ptv\mar2004\res\',10000,12000);
%
% See also: HELP READTRAJHDF HDF HDFTOOL
%
% 

% Modified at: 28-Nov-2003
% Combination of writeTrajHDF_v1, mat2sdsid, and TrajPoint2HDF_v0
% Each ASCII file trajPoint.* is opened, read in, directly written to HDF
% file and closed.
% Modified at: 03-Dec-2003
% on IHW_Liberzon
% - the structure is different now - due to the limit of SD bases in Matlab
% I will try to save it as a Vdata objects, each one will be of class
% 'trajectory', and its name will be the number of this trajectory,
% fields are 32 fields in traj{} structure, see below and 
% every column will be the variable column. It is almost the same structure
% as Beat's trajPoint.* files
%
% Modified at: 03-Dec-2003
% Version 3 works fine, but produces very large file
% this version would try to use externalfile option.
% name of the file is now an input
%
% Modified at: 04-Dec-2003 on ihw-liberzon
% Version 4 works fine, but some small question arised when I 
% checked the order of the written data. The user manual for HDF
% recommends to combine the fields into multi-component fields, like
% 'position' with order 3 (x,y,z) insteald of 3 fields 'x','y','z'
% In addition, the full interlace is recommended for the writing efficiency
% 
% Release version: trajPoint2HDF_v5.m
% 
% New release v.6
% this is due to the some trajPoint. file change, it seems to be more
% columns than before.
% in addition, directories with an enormous number of files takes for
% Windows too much to list it. The GUI approach is not that smart, 
% we use just command line input
%
% version 7
% - it is a limit of number of objects,
% therefore we need to have one object per file
% not one per trajectory.
% 
% version 8, March 25, 2004
% - one file instead of multi-files

error(nargchk(1,5,nargin))


% ---------------------------------------------------------------------------------------------------------        

wd = cd;
baseName = 'trajPoint.';
MINLENGTH = 0;

% quantities = varargin{1};       % 'vel', or 'der', or 'all'


traj = struct('x',[],'y',[],'z',[],...
    'u',[],'v',[],'w',[],...
    'ax',[],'ay',[],'az',[],...
    'w1',[], 'w2',[],'w3',[],...
    's11',[],'s12',[],'s13',[],...
    's22',[],'s23',[],'s33',[],...
    'ww1',[],'ww2',[],'ww3',[],...
    'wws',[],'sss',[],'R',[],'Q',[],...
    'diss',[],'div',[],'trajnum',0,'t',[],...
    'alx',[],'aly',[],'alz',[]);

fieldname_list = 'x,y,z,u,v,w,ax,ay,az,w1,w2,w3,s11,s12,s13,s22,s23,s33,ww1,ww2,ww3,wws,sss,R,Q,diss,div,trajnum,t,alx,aly,alz';

units = {'m','m','m',...                % x,y,z
        'm/s','m/s','m/s',...           % u,v,w
        'm/s^2','m/s^2','m/s^2',...     % ax,ay,az
        '1/s','1/s','1/s',...           % w1,w2,w3
        '1/s','1/s','1/s',...           % s11,...
        '1/s','1/s','1/s',...           % s22,...
        '1','1','1',...                 % wW1,...
        '1/s^3','1/s^3','1','1',...     % wws,sss,R,Q
        'm^2/s^3','1','1','1',...       % dissipation, divergence, traj, t, 
        'm/s^2','m/s^2','m/s^2'};         % alx,aly,alz
% vclass = {'position','position','position',...                % x,y,z
%         'velocity','velocity','velocity',...           % u,v,w
%         'acceleration','acceleration','acceleration',...     % ax,ay,az
%         'vorticity','vorticity','vorticity',...           % w1,w2,w3
%         'strain','strain','strain',...           % s11,...
%         'strain','strain','strain',...           % s22,...
%         'stretching','stretching','stretching',...                 % wW1,...
%         'wws','sss','R','Q',...     % wws,sss,R,Q
%         'dissipation','divergence','trajectory','time',...       % dissipation, divergence, t0, ii, 
%         'acceleration','acceleration','acceleration'};         % alx,aly,alz

dtype = repmat({'double'},32,1);
% dtype{28} = 'int16';
% dtype{29} = 'int16';



% ---------- HDF part --------------
try
    % Check if file exists:
    if hdf('H','ishdf', hfile)
        access = 'write';
        %     answer = input('Do you want to overwrite ... (y/n)  ','s');
        %     if strcmp(answer,'y')
        %       hdfml('closeall'), delete(hfile), access = 'create';
        %       feval('TrajPoint2HDF_v7',hfile,varargin{:});
        %     else
        error('Do not write to the same file again') % this option should be changed, if somebody wants to append 28.11.03
        % end
    else
        access = 'create';
    end
    
    % open new or existing hfile
    file_id = hdf('H','open', hfile, access, 0);
    if file_id == -1
        error('HDF hopen failed');
    end
    
    % initialize the V interface
    status = hdf('V','start',file_id);
    if status == -1
        error('HDF vstart failed');
    end
    
    
    % initialize the SD interface
    sd_id = hdf('SD','start', hfile, 'write');
    if sd_id == -1
        error('HDF sdstart failed');
    end
    
    % assign optional top level SD attributes
    astr.Date       =  date;
    astr.Author     = 'Alex Liberzon';
    astr.Comment    = sprintf('By trajPoint2HDF_v%d ', 8);
    astr.Directory        = cd;
    
    afields = fieldnames(astr);
    for j = 1 : length(afields)
        status = hdf('SD','setattr', sd_id, afields{j}, astr.(afields{j}));
        if status == -1
            error('HDF sdsetattr failed')
        end
    end
    
    % ---------------------------------------------------------------------------------------------------------        
    
    if nargin < 4
        
        [filename1, pathname] = uigetfile('trajPoint.*', 'Pick a FIRST trajectory file');
        
        if isequal(filename1,0) | isequal(pathname,0)
            cd(wd); error('Wrong selection')
        end
        
        cd(pathname)
        [filename2, pathname] = uigetfile('trajPoint.*', 'Pick a LAST trajectory file');
        if isequal(filename2,0) | isequal(pathname,0)
            cd(wd); error('Wrong selection')
        end
        
        % First and last index of the files
        firstFileIndx =  eval(filename1(findstr(filename1,'.')+1:end));
        lastFileIndx =  eval(filename2(findstr(filename2,'.')+1:end));
        
    else  % given HDF name, directory, first file, last file
        pathname = varargin{1};
        firstFileIndx = varargin{2};
        lastFileIndx =  varargin{3};
    end  
    
    fprintf(1,'Wait please ...  ')
    disp('');
    trajIndx = 1;
    %
    
    
    
    % Loop through all files
    for n = firstFileIndx:lastFileIndx
        
        runIndx = n - firstFileIndx + 1; % running counter
        
        fprintf(1,'\b\b\b\b\b');
        fprintf(1,'%5d',n); 
        
        % f = load(sprintf('%s%d',baseName,n),'ASCII');
        fid = fopen(sprintf('%s%s%s%d',pathname,filesep,baseName,n));
        if fid == -1, continue, end
        f = fscanf(fid,'%g',[35 inf]); % It has two rows now.
        % f = fscanf(fid,'%g',[32 inf]); % It has two rows now.
        if isempty(f), continue , end; % if the file is empty, skip to next file
        f = f(1:32,:)';
        fclose(fid);
        
        
        
        
        [startTraj, endTraj, lenTraj] = singleTraj(f(:,29)); % prepare indices of the single trajectories
        
        % pick only long trajectories, longer than 'len'
        longTraj = find(lenTraj >= MINLENGTH);
        if isempty(longTraj), fprintf(1,'!'), continue, end; % if there are no long trajectories, move to the next file
        
        
        numLongTraj = length(longTraj);     % number of long trajectories
        startTraj = startTraj(longTraj);    % starting  
        endTraj = endTraj(longTraj);        % ending points
        lenTraj = lenTraj(longTraj);        % lengthes
        %     maxLenTraj = max(lenTraj);          % maximum length
        
        validIndx = zeros(size(f,1),1);
        for i = 1:numLongTraj
            validIndx(startTraj(i):endTraj(i)) = 1;
            f(startTraj(i):endTraj(i),28) = trajIndx;
            trajIndx = trajIndx + 1;
        end
        validIndx = find(validIndx);
        
        
        
        
        % ------------- HDF part -------------------
        % create a new vs data 
        access = 'w';
        vdata_ref = -1;                                         % flag to create
        vdata_id = hdf('VS','attach', file_id, vdata_ref, access);
        %         vdata_id = hdf('V',s('attach', file_id, vdata_ref, access);
        if vdata_id == -1
            error('HDF vdata attach failed')
        end
        
        %  offset = 0;
        
        
        % assign the group the name "i"
        % sname = sprintf('%d', lastTrajIndx);
        sname = sprintf('%d', runIndx);
        status = hdf('VS','setname', vdata_id, sname);
        if status == -1
            error('HDF vsetname failed')
        end
        
        % assign the group the class "Trajectory"
        sclass = 'Trajectory';
        status = hdf('VS','setclass', vdata_id, sclass);
        if status == -1
            error('HDF vsetclass failed')
        end
        
        dfields = fieldnames(traj);
        
        for j = 1 : length(dfields)                 % 32 fields
            status = hdf('VS','fdefine',vdata_id, dfields{j},'double',1);
        end
        
        status = hdf('VS','setfields',vdata_id,fieldname_list);
        if status == -1
            error('HDF VS setfields failed')
        end
        
        % Interlace mode is set to full, meaning that data is written row
        % by row of the record, rather than by field (not all x, then all
        % y, but all x,y,z, ... for the record 1, then for 2, etc)
        status = hdf('VS','setinterlace',vdata_id,'full');
        if status == -1
            error('HDF VS setinterlace failed')
        end
        
        f(validIndx,1:3) = 1000 * f(validIndx,1:3);
        
        % write the buffered data into the first vdata with full interlace mode
        num_of_records = hdf('VS','write',vdata_id,{f(validIndx,1)',f(validIndx,2)', f(validIndx,3)', ...
                f(validIndx,4)',f(validIndx,5)', f(validIndx,6)', ...
                f(validIndx,7)',f(validIndx,8)', f(validIndx,9)', ...
                f(validIndx,10)',f(validIndx,11)', f(validIndx,12)', ...
                f(validIndx,13)',f(validIndx,14)', f(validIndx,15)', ...
                f(validIndx,16)',f(validIndx,17)', f(validIndx,18)', ...
                f(validIndx,19)',f(validIndx,20)', f(validIndx,21)', ...
                f(validIndx,22)',f(validIndx,23)', f(validIndx,24)', ...
                f(validIndx,25)',f(validIndx,26)', f(validIndx,27)', ...
                f(validIndx,28)',f(validIndx,29)'+runIndx, ...
                f(validIndx,30)', f(validIndx,31)',f(validIndx,32)'});
        
        % nbytes = hdf('VS','sizeof',vdata_id,fieldname_list)
        
        % offset = offset + nbytes*num_of_records % *128;
        
        %     status = hdf('VS','setattr',vdata_id,'vdata','MaxLength',max(lenTraj));
        %     if status == -1
        %         error('HDF vsetattr(1) failed')
        %     end
        %     status = hdf('VS','setattr',vdata_id,'vdata','MaxLength',max(lenTraj));
        %     if status == -1
        %         error('HDF vsetattr(1) failed')
        %     end
        %     status = hdfvs('setattr', vdata_id, 27 , 'lengths' , lenTraj );
        %     if status == -1
        %         error('HDF vsetattr(2) failed')
        %     end
        
        % assign the group the class "Trajectory"
        %     sclass = 'Trajectory';
        sclass = sprintf('%d',max(lenTraj));
        status = hdf('VS','setclass', vdata_id, sclass);
        if status == -1
            error('HDF vsetclass failed')
        end
        
        
        % %     status = hdf('VS','setexternalfile',vdata_id,sprintf('%s.%d',hfile,runIndx), 0);
        % %     if status == -1
        % %         error('HDF v external failed')
        % %     end
        
        
        % detach from the first vdata
        status = hdf('VS','detach',vdata_id);
        if status == -1, error('HDF vsdetach failed'), end
        %         lastTrajIndx = lastTrajIndx + 1; 
        %     end  % of all trajectories in the file
        
    end % FILES
    % ---------------------------------------------------------------------------------------------------------        
    
    
    % end vgroup interface access
    status = hdf('V','end',file_id);
    if status == -1
        error('HDF v end failed')
    end
    
    % end vgroup interface access
    status = hdf('SD','end',sd_id); 
    if status == -1
        error('HDF sd end failed')
    end
    
    % close the HDF file
    status = hdf('H','close',file_id);
    if status == -1
        error('HDF hclose failed')
    end
    % ---------------------------------------------------------------------------------------------------------
    
    
    cd(wd);
    
catch
    hdfml('closeall')
    cd(wd)
end

return
% ---------------------------------------------------------------------------------------------------------


function [startTraj, endTraj, lenTraj] = singleTraj(age)
% Split the file into single trajectories, according to the
% AGE vector, column 29 of the trajPoint.* files
startTraj = find(age == 0);      % start points of the trajectories
endTraj = startTraj(2:end)-1;    % end points of the trajectories
endTraj = [endTraj;length(age)];
lenTraj = endTraj - startTraj + 1;
%    numTraj = length(startTraj);
return



