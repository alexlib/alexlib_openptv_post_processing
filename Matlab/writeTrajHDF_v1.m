function writeTrajHDF_v1(hfile, dstr, astr)
% WRITETRAJHDF_V1(HFILE, TRAJ, ATTR) Writes the matlab structure TRAJ 
% (structure of arrays) to the HDF file (ver. 4)
% In addition, the file attributes are set to ATTR cell array of strings
%
% Example: within the al_readTraj()
%
%       writeTrajHDF_v1('traj.hdf',traj,struct('comment','2 cm from disk, 75 rpm','date',date,'author','Alex'); 
%          
% 
% Based completely on H4SGWRITE, only the name of the 
% structure was modified from record # to number, and
% the class of structure is Trajectory.
%
% NAME
%
%   h4sgwrite -- write a matlab structure array as a vgroups of SD's
%
% SYNOPSIS
%
%   h4sgwrite(hfile, dstr, astr)
%
% INPUTS
%
%   hfile  - hdf file name
%   dstr   - 1-d array of structures containing data to write
%   astr   - structure of global attributes 
%
% OUTPUT
%
%   an HDF file with a vgroup of SDS's for each structure index
%
% BUGS
%
%   - still "in progress", passes some basic tests
%   - HDF seems to have a hard limit of 5000 SDS's total
%   - reopen and write to an existing big file is slow
%   - needs the capability to update single fields
%   - probably some performance tuning would help, maybe
%     calculate num_dds_block from size of inp