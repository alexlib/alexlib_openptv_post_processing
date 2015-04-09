
function [xd,yd,p] = create_fit_set_of_points(figureHandle,N)
% create_fit_set_of_points(figureHandle) solves the 
% fit problem for the figure which is a set of separate  points
% rather than a single vector. 
% Default fit is linear, optioinal N can be used:
% >> create_set_of_points(9,2);
%

if nargin < 2
    N = 1; % default is linear
end
ch = findobj(figureHandle,'type','patch');
if ~isempty(ch)
    xd = zeros(length(ch),1);
    yd = zeros(length(ch),1);
    for i = 1:length(ch)
        xd(i) = get(ch(i),'XData');
        yd(i) = get(ch(i),'YData');
    end
end
figure
 p = polyfit(xd,yd,N);
 plot(xd,yd,'o',xd,polyval(p,xd),'-');
