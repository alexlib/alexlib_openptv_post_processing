
function [xd,yd,p] =  create_fit_set_of_points (figureHandle,N,varargin)
% create_fit_set_of_points(figureHandle) solves the 
% fit problem for the figure which is a set of separate  points
% rather than a single vector. 
% Default fit is linear, optioinal N can be used:
% >> create_set_of_points(9,2);
% 
% Default marker is just a patch or a line, but if specific selection 
% is needed, provide an additional input:
% >> chreate_set_of_points(gcf,1,'o');
%
% Example:
% 
% hf = figure;
% hold on
% for i = 1:5
%   plot(rand(1),rand(1),'o');
% end
% hold off
% % try linear fit
% create_fit_set_of_points(hf,1,'o');
% 

if nargin < 2
    N = 1; % default is linear
end

if nargin > 2
    ch = findobj(figureHandle,'type','line','marker',varargin{1});
else
    ch = findobj(figureHandle,'type','line');
    if isempty(ch)
        ch = findobj(figureHandle,'type','patch');
    end
end
    
if ~isempty(ch)
    xd = zeros(length(ch),1);
    yd = zeros(length(ch),1);
    for i = 1:length(ch)
        xd(i) = get(ch(i),'XData');
        yd(i) = get(ch(i),'YData');
    end
end
[xd,i] = sort(xd);
yd = yd(i);

figure
 p = polyfit(xd,yd,N);
 plot(xd,yd,'o',xd,polyval(p,xd),'-');
