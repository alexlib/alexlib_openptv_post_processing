% 
[XD,xD,yD,zD,count] = deal(zeros(300,1));
for i = 1:length(newtraj)
        
    x2 = newtraj(i).xf(2:end) - newtraj(i).xf(1);
    
    y2 = newtraj(i).yf(2:end) - newtraj(i).yf(1);
    
    z2 = newtraj(i).zf(2:end) - newtraj(i).zf(1);
    
    X2 = sum([x2,y2,z2].^2,2);
    x2 = x2.^2; y2 = y2.^2;z2 = z2.^2;
    
    len = length(x2);
    
    XD(1:len) = XD(1:len) + X2;
    xD(1:len) = xD(1:len) + x2;
    yD(1:len) = yD(1:len) + y2;
    zD(1:len) = zD(1:len) + z2;
    count(1:len) = count(1:len) + ones(length(x2),1);
    
end



xD = xD./count;
yD = yD./count;
zD = zD./count;
XD = XD./count;
figure, plot(xD), hold on, plot(yD,'g--'); plot(zD,'m-.'); plot(XD,'r-o')


% velocity 

[XD,xD,yD,zD,count] = deal(zeros(300,1));
for i = 1:length(newtraj)
        
    x2 = newtraj(i).uf(2:end) - newtraj(i).uf(1);
    
    y2 = newtraj(i).vf(2:end) - newtraj(i).vf(1);
    
    z2 = newtraj(i).wf(2:end) - newtraj(i).wf(1);
    
    X2 = sum([x2,y2,z2].^2,2);
    x2 = x2.^2; y2 = y2.^2;z2 = z2.^2;
    
    len = length(x2);
    
    XD(1:len) = XD(1:len) + X2;
    xD(1:len) = xD(1:len) + x2;
    yD(1:len) = yD(1:len) + y2;
    zD(1:len) = zD(1:len) + z2;
    count(1:len) = count(1:len) + ones(length(x2),1);
    
end



xD = xD./count;
yD = yD./count;
zD = zD./count;
XD = XD./count;
figure, plot(xD), hold on, plot(yD,'g--'); plot(zD,'m-.'); plot(XD,'r-o')