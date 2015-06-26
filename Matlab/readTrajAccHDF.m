function [dstr, astr] = readTrajAccHDF(hfile,varargin)
%
% [TRAJECTORIES,ATTRIBUTES] = READTRAJACCHDF(HDF_FILE,FIELD,VALUE,FIELD,VALUE)
%
% ReadTrajAccHDF reads the particle and flow quantities along the trajectories 
% from the HDF files (new Acc.HDF files, HDF version release 4.2, new Vdata
% format, July 2004)
% 
% Inputs:
%           HDF_FILE       - File name, including (!) the extension .hdf
%           (i.e. 'tmp.hdf')
%
%           FIELDS      - Field names, predefined in lower-case, given as a
%           string, e.g., 'x','u','ww1', etc. If not even one variable will
%           be provided, the whole set of 32 variables will be read in.
%           Three special FIELDS are 'minlength' and 'maxlength' that limit
%           the length of the trajectories, and the 'frames' that limits
%           the number of the frames to read.
%
% struct('x',[],'y',[],'z',[],...
%     'u',[],'v',[],'w',[],...
%     'dudx',[],'dudy',[],'dudz',[],...
%     'dvdx',[],'dvdy',[],'dvdz',[],...
%     'dwdx',[],'dwdy',[],'dwdz',[],...
%     'dudt',[],'dvdt',[],'dwdt',[],...
%     'ax',[],'ay',[],'az',[],...
%     'daxdx',[],'daxdy',[],'daxdz',[],...
%     'daydx',[],'daydy',[],'daydz',[],...
%     'dazdx',[],'dazdy',[],'dazdz',[],...
%     'trajnum',0,'t',[]);
%
%           VALUE - Limiting values of the FIELDS
%
%
%
% Outputs:
%           TRAJECTORIES    - Matlab structure of arrays, like traj.x,
%           traj.y, etc. Each array of different size, equal to the length
%           of that trajectory, if the condition is length, then the
%           structure will include all the variables in the file
%           (x,y,z,u,v,w ...) but of the certain length
%           If the condition is 'variable', than all trajectories will be
%           loaded, but only with one field, equal to the condition value
%           (e.g. 'x');
%
%           ATTRIBUTES - Structure of attributes to show the original
%           directory of the run, date, author, etc. 
%
% Example: 
%           [traj,attr] = readTrajAccHDF('test.hdf','u',[],'dudx',[],'daxdx',[],'minlength',20,'frames',5);
%
% See also: HELP TRAJPOINT2HDF_V32 TRAJPOINT2HDF_V35 READTRAJHDF_V9
%



% Parse inputs

if nargin == 0 % temproral solution for the debugging process
%     hfile = 'polymer.hdf';
%     fields.y = [];
%     fields.x = [-Inf Inf];
%     fields.z = [-100 0];
%     fields.minlength = 10;
%     fields.maxlength = 100;
%     fields.frames = [1 6000];
error('No inputs, see help for usage example');
elseif mod(nargin,2)                     % even number of varargins
    for i = 1:2:nargin-1
% fields.(lower(varargin{i})) = varargin{i+1}; % field names,
        fields.(varargin{i}) = varargin{i+1}; % field names, 
    end
end

% Length restrictions
if isfield(fields,'minlength'),     
    minLength   = fields.minlength; 
    fields = rmfield(fields,'minlength');
else 
    minLength = -Inf;
end

if isfield(fields,'maxlength'),     
    maxLength   = fields.maxlength; 
    fields = rmfield(fields,'maxlength');
else 
    maxLength = Inf;
end
if isfield(fields,'frames'),     
    frames   = fields.frames; 
    fields = rmfield(fields,'frames');
else 
    frames = [1 Inf];
end

% Preparing the list of variables to read from the data
% Always read traj
fieldNames = fieldnames(fields);
% If it is empty, read all the fields
if isempty(fieldNames)
%     fieldNames  = {'x', 'y', 'z', ...       % 1 - 3 
%             'u', 'v', 'w', ...          % 4 - 6
%             'ax', 'ay', 'az', ...       % 7 - 9
%             'w1',  'w2', 'w3', ...      % 10 - 12
%             's11', 's12', 's13', ...    % 13 - 15
%             's22', 's23', 's33', ...    % 16 - 18
%             'ww1', 'ww2', 'ww3', ...    % 19 - 21
%             'wws', 'sss', 'R', 'Q', ... % 22 - 25
%             'diss', 'div', 'trajnum','t', ... % 26 - 29
%             'alx', 'aly', 'alz'};       % 30 - 32 <-- New columnt, 30.01.04
fieldNames = {'x','y','z',...
    'u','v','w',...
    'dudx','dudy','dudz',...
    'dvdx','dvdy','dvdz',...
    'dwdx','dwdy','dwdz',...
    'dudt','dvdt','dwdt',...
    'ax','ay','az',...
    'daxdx','daxdy','daxdz',...
    'daydx','daydy','daydz',...
    'dazdx','dazdy','dazdz',...
    'trajnum','t'};

    for i = 1:length(fieldNames)
%         fields.(lower(fieldNames{i})) = []; % field names, 
                fields.(fieldNames{i}) = []; % field names, 

    end
end
i =  find(strcmp(fieldNames,'trajnum'));
fieldNames{end+1} = 'trajnum';
if ~isempty(i) 
    fieldNames(i) = [];
end
% Reformat into long string, comma-separated list of variables
% for HDFVS
fields2read = fieldNames{1};
for i = 2:length(fieldNames)
    fields2read = strcat(fields2read,',',deblank(fieldNames{i}));
end   


% ------------------------------------------------------
% open the HDF file
file_id = hdfh('open', hfile, 'read', 0);
if file_id == -1
    error('HDF hopen failed');
end

% initialize the SD interface
sd_id = hdfsd('start', hfile, 'read');
if sd_id == -1
    error('HDF sdstart failed');
end

% general info about SD contents
[ndatasets,nglobal_attr,status] = hdfsd('fileinfo',sd_id);
if status == -1
    error('HDF sdfileinfo failed');
end


% read global attributes

% loop on global attributes
for j = 0 : nglobal_attr - 1
    [aname,data_type,attr_count,status] = hdfsd('attrinfo', sd_id, j);
    if status == -1
        error('HDF sdattrinfo failed');
    end
    
    [adata,status] = hdfsd('readattr', sd_id, j);
    if status == -1
        error('HDF sdreadattr failed');
    end
    
    % copy attribute name and value
    % astr.<aname> = adata;
    eval(sprintf('astr.%s = adata;', aname));
    
end



% initialize the V interface
status = hdfv('start',file_id);
if status == -1
    error('HDF vstart failed');
end

% % get number of top-level v data groups
[vrefs, count] = hdfvs('lone', file_id, 1);
[vrefs, count] = hdfvs('lone', file_id, count);


for i = 1:size(fieldNames,1)
    dstr.(fieldNames{i}) = repmat(0, 250, 1); % assume 500 values in each trajectory
end
dstr = repmat(dstr,1,min(count,diff(frames))*250); % assume 500 trajectories per frame

% loop on vdatas / structure records
k = 0;
disp(' ');

for i = max(1,frames(1)):min(count,frames(2)) % for all Vdatas or for all selected frames
    
    fprintf(2,'\b\b\b\b\b\b')
    fprintf(2,'%6d',i)
    
    
    vdata_id = hdfvs('attach', file_id, vrefs(i), 'r');
    if vdata_id == -1
        error('HDF vattach failed')
    end
    
    [vdata_class, status] = hdfvs('getclass', vdata_id);
    if status == -1
        error('HDF vgetname failed')
    end
    
    maxTrajLength = eval(vdata_class);
    
    % if there is a minumum length restriction
    if maxTrajLength < minLength, continue, end
    
    % Inquire the information about the vdata, now it is a number of points
    % of all trajectories in this file
    numrecords = hdfvs('elts',vdata_id);
    
    % Read only the selected variables, fields2read
    status = hdfvs('setfields',vdata_id,fields2read);
    if status == -1
        error('HDF setfields failed');
    end
    
    [data,status] = hdfvs('read',vdata_id,numrecords);
    if status == -1
        error('HDF Vsread failed');
    end
    
    tmp = zeros(numrecords,length(fieldNames));
    % Conditioning on all properties
    for j = 1:length(fieldNames)
        tmp(:,j) = double([data{j}])';
    end
    

    
    
    ind = [1:numrecords]';
    for j = 1:length(fieldNames)
%          if ~isempty(fields.(lower(fieldNames{j})))
       if ~isempty(fields.(fieldNames{j}))
            ind = intersect(ind, find(tmp(:,j) >= fields.(fieldNames{j})(1) & tmp(:,j) <= fields.(fieldNames{j})(2)));
            % new development of cutting of trajectories
            ind1 = diff([ind;Inf]);
            ind2 = [0; find(ind1 > 1)];
            ind3 = diff(ind2);
            [i,j] = max(ind3);
            ind = ind(ind2(j)+1:ind2(j+1));
        end
    end      
    if isempty(ind), continue, end
    tmp = tmp(ind,:);
    
    
        
    % Separate the data into trajectories, look for long enough
    % trajectories, if there is a minimum length restriction:
    %
    lenTraj     = diff( [0;find(diff( [tmp(:,end);tmp(end,end)+1] )) ] );                      % all lengthes
    startTraj   = [1;find(diff( [tmp(:,end)] ))+1];                                        % all starting points
    
    % Find all the trajectories with an appropriate length
    ind  = find( lenTraj >= minLength & lenTraj <= maxLength );

    for m = 1:length(ind)
        k = k + 1;
        dstr(k).trajLen  = lenTraj(ind(m));
        for j = 1:length(fieldNames)
            dstr(k).(fieldNames{j}) = tmp(startTraj(ind(m)):startTraj(ind(m))+lenTraj(ind(m))-1,j);            
        end
        
    end
    clear tmp data firstRecord lenRecord ind
    
    % end access to the current vgroup
    status = hdfvs('detach', vdata_id);
    if status == -1
        error('HDF vdetach failed')
    end
end % loop on vgroups


% end SD interface access
status = hdfsd('end',sd_id);
if status == -1
    error('HDF sdend failed')
end

% end vgroup interface access
status = hdfv('end',file_id);
if status == -1
    error('HDF vend failed')
end

% close the HDF file
status = hdfh('close',file_id);
if status == -1
    hdfml('closeall')
    warning('HDF hclose failed')
end

dstr = dstr(1:k);

return