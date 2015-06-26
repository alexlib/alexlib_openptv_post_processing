function [dstr, astr] = readStatisticsHDF_4Beat(hfile,varargin)
%
% [DATA] = READStatisticsHDF(HDF_FILE,FIELD,VALUE,FIELD,VALUE)
%
% ReadStatisticsHDF reads the statistics from Eulerean data, from the HDF files 
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
% dfields = {'lambda1','lambda2','lambda3','sss','wws','s2','u2','w2','cosss','cosww','coswl1','coswl2',...
%             'coswl3','reldiv','ii','cosulambda1m','cosulambda2m','cosulambda3m'};
%
%
%           VALUE - Limiting values of the FIELDS
%
%
%
% Outputs:
%
% Example: 
%           [data] = readStatisticsHDF('meanCorField_1_5.hdf','lambda1',[],'lambda2',[],'lambda3',[],'ii',[]);
%
% See also: HELP STATISTICSFROMGRIDDIV2HDF
%

% Created: 18-Aug-2004
% Author: Alex Liberzon 
% E-Mail : liberzon@ihw.baug.ethz.ch 
% Phone : +41 (0)1 633 3754 
% Copyright (c) 2004 IHW, ETH, Zuerich 
%
% Modified at: 18-Aug-2004
% $Revision: 1.0 $  $Date: 18-Aug-2004 09:45:11$ 



% Parse inputs

if nargin == 0 % temproral solution for the debugging process
    %     hfile = 'polymer.hdf';
    %     fields.y = [];
    %     fields.x = [-Inf Inf];
    %     fields.z = [-100 0];
    %     fields.minlength = 10;
    %     fields.maxlength = 100;
    %     fields.frames = [1 6000];
    error('No inputs, see >> help readStatisticsHDF');
elseif mod(nargin,2)                     % even number of varargins
    for i = 1:2:nargin-1
        % fields.(lower(varargin{i})) = varargin{i+1}; % field names,
        if isempty(varargin{i+1})
            eval(['fields.',varargin{i},' = [];']); % field names, 
        else
            eval(['fields.',varargin{i},' = ',varargin{i+1},';']); % field names, 
            %         fields.(varargin{i}) = varargin{i+1}; % field names, 
        end
    end
end

% Preparing the list of variables to read from the data
% Always read traj

%fieldNames = fieldnames(fields)
% If it is empty, read all the fields
% if isempty(fieldNames)

    fieldNames = {'lambda1','lambda2','lambda3',...
            'sss','wws','s2','u2','w2','cosss',...
            'cosww','coswl1','coswl2','coswl3',...
            'reldiv','ii','cosulambda1m',...
            'cosulambda2m','cosulambda3m'};
    for i = 1:length(fieldNames)
        %fields.(fieldNames{i}) = []; % field names, 
        eval(['fields.',fieldNames{i},' = [];']); % field names, 
        
    end
% end

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

% initialize the V interface
status = hdfv('start',file_id);
if status == -1
    error('HDF vstart failed');
end

% % get number of top-level v data groups
[vrefs, count] = hdfvs('lone', file_id, 1);
[vrefs, count] = hdfvs('lone', file_id, count);


for i = 1:size(fieldNames,1)
    % dstr.(fieldNames{i}) = repmat(0, 1000, 1); % assume 500 values in each FILE
    eval(['dstr.',fieldNames{i},' = repmat(0,1000,1);']); % field names, 
    
end
dstr = repmat(dstr,1,count*100); % 

% loop on vdatas / structure records
k = 0;
disp(' ');

for i = 1:count % for all Vdatas or for all selected frames
    
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
    
    % Inquire the information about the vdata, now it is a number of points
    % of all trajectories in this file
    numrecords = hdfvs('elts',vdata_id);
    
    if numrecords == 0, continue, end;
    
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
        eval(['tmpField = fields.',fieldNames{j},';']);
        
        if ~isempty(tmpField) % fields.(fieldNames{j}))
            %             ind = intersect(ind, find(tmp(:,j) >= fields.(fieldNames{j})(1) & tmp(:,j) <= fields.(fieldNames{j})(2)));
            
            ind = intersect(ind, find(tmp(:,j) >= tmpField(1) & tmp(:,j) <= tmpField(2)));
            %             % new development of cutting of trajectories
            %             ind1 = diff([ind;Inf]);
            %             ind2 = [0; find(ind1 > 1)];
            %             ind3 = diff(ind2);
            %             [i,j] = max(ind3);
            %             ind = ind(ind2(j)+1:ind2(j+1));
        end
    end      
    if isempty(ind), continue, end
    tmp = tmp(ind,:);
    
    
    
    %     for m = 1:length(ind)
    k = k + 1;
    for j = 1:length(fieldNames)
        eval(['dstr(k).',fieldNames{j},' = tmp(:,j);']);            
    end
    
    %     end
    clear tmp data firstRecord lenRecord ind
    
    % end access to the current vgroup
    status = hdfvs('detach', vdata_id);
    if status == -1
        error('HDF vdetach failed')
    end
end % loop on vgroups


% end vgroup interface access
status = hdfv('end',file_id);
if status == -1
    error('HDF vend failed')
end

% close the HDF file
status = hdfh('close',file_id);
if status == -1
    hdfml('closeall')
    %     warning('HDF hclose failed')
end

dstr = dstr(1:k);

return