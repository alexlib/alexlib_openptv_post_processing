clear,clc

uiopen('Reynolds_stress_scaled.fig',1);
legend(gca,'hide')
lines = get(gca,'children');

N = 100;

yd1= get(lines(4),'ydata') + N;

for i = 1:length(lines)-1
    xd = get(lines(i),'xdata');
    yd = get(lines(i),'ydata');
    xdata = xd;
    ydata = (yd+N)./yd1;
    set(lines(i),'xdata',xdata);
    set(lines(i),'ydata',(ydata - 1)*N+1);
    set(get(gca,'ylabel'),'string','magnitude/magnitude(1.5Hz)') 
    set(get(gca,'xlabel'),'string','x/L')
end

legend(gca,'show')

