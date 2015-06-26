function writeTrajHDF_v1
% First version of writing trajectories 
% from the Matlab structure into HDF file
% Based on the idea of example 3 in HDF manual

% Load the test file
load traj2cm 

% Define file name
hfile = 'traj2cm.hdf';

% Attributes of the file
file_attr.comment   = 'Experiment: 2cm from the rotating disk, 75 rpm';
file_attr.date      = '24.11.03';
file_attr.author    = 'Alex Liberzon';



% Very short way, use the h4sgwrite, as it is used in sdtest1.m
traj_h4sgwrite(hfile, traj2cm, file_attr);

return
% -----------------------------------------------------------------------


% -------------------------------------------------------------------------



DFTAG_NDG   = 720;
pxy         = rand(30,2);
plist       = reshape([1:90],30,3);
tmp         = 10*(1:30);

% check if hfile already exists
if hdfh('ishdf', hfile)
  access = 'write'
else
  access = 'create'
                            error('dont write twice to the same file')
end

% try 
    
% open new or existing hfile
file_id = hdfh('open', hfile, access, 0);
if file_id == -1
  error('HDF hopen failed');
else
  fprintf(2, 'opened %s mode %s\n', hfile, access);
end

% initialize the V interface
status = hdfv('start',file_id);
if status == -1
  error('HDF vstart failed');
end

% create the vgroup
  access = 'w';
  vgroup_ref = -1; % flag to cr