function [varargout] = anisotropy(X,varargin)
%ANISOTROPY Calculates deviator, or tensor of anisotropy of square matrix X
%of size 3 x 3, and then it's two modified invariants, II and III
% For more details, see Lumley 1973, "Back to isotropy in homogeneous
% turbulence".
%
%   [II,III,AX] = ANISOTROPY(X)
%
%   inputs:
%      X  - 3 x 3 matrix
%
%   outputs:
%      II,III  - Second and third 'modified' invariants of the deviator AX
%      AX       - Deviator of the square matrix is simply = X./trace(X) - 1/3*eye(3)
%
%   example:
%           [II,III] = anisotropy(B);
%           scatter(II,III);
%
%   notes:
%
%   See also HELP INVARIANTS, DEVIATOR
%

% Created: 15-Nov-2004
% Author: Alex Liberzon & Michele Guala
% E-Mail : liberzon@ihw.baug.ethz.ch
% Phone : +41 (0)1 633 3754
% Copyright (c) 2004 IHW, ETH Zurich
%
%
% $Revision: 1.0 $  $Date: 15-Nov-2004 22:11$

% $Revision: 3.0 $  $Date: 15-Nov-2004 22:11$
% All is replaced by calling functions:

%% Functional approach
[I,II,III] = invariants(deviator(X));
III = -III; % note that in Lumley and others, the third invariant is defined as 
% minus of the invariant of Perry and Chong. for some reason ???

%% Detailed approach 
% dev = X./sum(diag(X)) - 1/3*eye(3);
% 
% % dev = uvtensor./sum(diag(uvtensor)) - 1/3*eye(3);
% 
% a111=dev(1,1)*dev(1,1)*dev(1,1);
% a222=dev(2,2)*dev(2,2)*dev(2,2);
% a333=dev(3,3)*dev(3,3)*dev(3,3);
% a112=dev(1,1)*dev(1,2)*dev(1,2);
% a113=dev(1,1)*dev(1,3)*dev(1,3);
% a221=dev(2,2)*dev(1,2)*dev(1,2);
% a223=dev(2,2)*dev(2,3)*dev(2,3);
% a331=dev(3,3)*dev(1,3)*dev(1,3);
% a332=dev(3,3)*dev(2,3)*dev(2,3);
% a123=dev(1,2)*dev(2,3)*dev(1,3);
% 
% I = trace(dev);
% 
% if nargin > 2 & ~isempty(findstr(varargin{1},'not'))
%     % not-modified invariants
%     II = 0.5*(I^2 - trace(dev^2)); 
%     III = -det(dev);
% else % modified invariants
%     % ref: Turbulent stress invariant analysis: classification of existing
%     % terminology by Simonsen A.J. and Krogstad P.-A., 2004
%     
%     II = -0.5*(dev(1,1)^2+dev(2,2)^2+dev(3,3)^2 ...
%         +dev(1,2)^2+dev(2,3)^2+dev(1,3)^2 ...
%         +dev(2,1)^2+dev(3,2)^2+dev(3,1)^2); 
% 
%     III = a111+a222+a333+ ...
%         3*(a112+a113+a221+a223+a331+a332)+ ...
%         6*a123;
% end

%%
if nargout == 1
    varargout{1} = II;
elseif nargout == 2
    varargout{1} = II;
    varargout{2} = III;
elseif nargout == 3
    varargout{1} = I;
    varargout{2} = II;
    varargout{3} = III;
else
    varargout{1} = [II,III];
end