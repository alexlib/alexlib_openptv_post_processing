function [varargout] = anisotropymap(varargin)
%ANISOTROPYMAP *Insert a one line summary here*
%   [VARARGOUT] = ANISOTROPYMAP(VARARGIN)
%
%   inputs:
%      varargin  - *Insert description of input variable here*
%
%   outputs:
%      varargout  - *Insert description of output variable here*
%
%   example:
%
%   notes:
%
%   See also HELP 
%

% Created: 24-Jan-2005
% Author: Alex Liberzon 
% E-Mail : liberzon@ihw.baug.ethz.ch 
% Phone : +41 (0)1 633 3754 
% Copyright (c) 2005 IHW, ETH Zurich 
%
%
% $Revision: 1.0 $  $Date: 24-Jan-2005 22:52$ 

hf = figure;
hold on;

a = [1 0 0; 
    0 0 0; 
    0 0 0];
[II,III] = anisotropy(a);
plot(III,-II,'o')

a = [1 0 0; 
    0 1 0; 
    0 0 0];
[II,III] = anisotropy(a);
plot(III,-II,'o')

a = [1 0 0; 
    0 1 0; 
    0 0 1];
[II,III] = anisotropy(a);
plot(III,-II,'o')

a = [1 1 0; 
    1 0 0; 
    0 0 0];
[II,III] = anisotropy(a);
plot(III,-II,'o')

a = [1 1 0; 
    1 1 0; 
    0 0 0];
[II,III] = anisotropy(a);
plot(III,-II,'o')

a = [1 1 0; 
    1 1 0; 
    0 0 1];
[II,III] = anisotropy(a);
plot(III,-II,'o')



a = [1 1 1; 
    1 0 0; 
    1 0 0];
[II,III] = anisotropy(a);
plot(III,-II,'.')

a = [1 1 1; 
    1 1 0; 
    1 0 0];
[II,III] = anisotropy(a);
plot(III,-II,'.')

a = [1 1 1; 
    1 1 0; 
    1 0 1];
[II,III] = anisotropy(a);
plot(III,-II,'.')


a = [1 1 1; 
    1 0 1; 
    1 1 0];
[II,III] = anisotropy(a);
plot(III,-II,'.')

a = [1 1 1; 
    1 1 1; 
    1 1 0];
[II,III] = anisotropy(a);
plot(III,-II,'.')

a = [1 1 1; 
    1 1 1; 
    1 1 1];
[II,III] = anisotropy(a);
plot(III,-II,'.')


varargout{1} = hf;
hold off;