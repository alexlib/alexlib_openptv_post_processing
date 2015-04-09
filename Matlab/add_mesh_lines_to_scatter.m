openfig('C:\Users\hadar\Dropbox\resuspension\tau_thesis\matlab_thesis\scatter\scatter_u_v_normalized_110.fig')
legend off
ch = get(gca,'children');
xd = get(ch(1),'xdata'); yd = get(ch(1),'ydata');
for i = 2:length(ch)
    xd = cat(2,xd,get(ch(i),'xdata')); 
    yd = cat(2,yd,get(ch(i),'ydata')); 
end

hold on


edges = {linspace(min(xd),max(xd),15),linspace(min(yd),max(yd),15)};
%
hist3([xd(:),yd(:)],'edges',edges);
% [c,h] = contour(edges{1},edges{2},smoothn(cnt'),'linewidth',2);
% clabel(c,h),
% colorbar
