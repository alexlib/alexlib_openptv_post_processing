%PDF calculations

%reading the data
tmp2= readXUAPFiles('res_scene32');
tmp1=readXUAPFiles('reut2');
x1=[];x2=[];
y1=[];y2=[];
z1=[];z2=[];
vf_vel1=[];vf_vel2=[];
% 
for i=2:length(tmp1)  %the first file is zero file so we start at the second one
    ind= (find(abs(tmp1(i).uf) > 10e-5 & abs(tmp1(i).vf) > 10e-5 & abs(tmp1(i).wf) > 10e-5));
    if length(ind) > 10
        x1=[x1,(tmp1(i).xf(ind))']; % joining the x arrays into one in order to craete pdf
        y1=[y1,(tmp1(i).yf(ind))'];
        z1=[z1,(tmp1(i).zf(ind))'];
     % in order to check settling velocity we also create y velocity vector
        vf_vel1=[vf_vel1,(tmp1(i).vf(ind))'];
    end
end
for i=2:length(tmp2)  %the first file is zero file so we start at the second one
    ind= (find(abs(tmp2(i).uf) > 10e-5 & abs(tmp2(i).vf) > 10e-5 & abs(tmp2(i).wf) > 10e-5));
    if length(ind) > 10
        x2=[x2,(tmp2(i).xf(ind))']; % joining the x arrays into one in order to craete pdf
        y2=[y2,(tmp2(i).yf(ind))'];
        z2=[z2,(tmp2(i).zf(ind))'];
     % in order to check settling velocity we also create y velocity vector
        vf_vel2=[vf_vel2,(tmp2(i).vf(ind))'];
    end
end
% histogram
hist(z1,100)
figure
hist(z2,100)



hist(x,100)
title('histogram  of x axis scattering of the particles')
xlabel('location [m] from up to downsteam wall')
ylabel('number of occurence')

a=find(y2>=-0.04);
y2=y2(a);
hist(y_new,100)
title('histogram  of y axis scattering of the particles')
xlabel('location [m] from floor to the driving band')
ylabel('number of occurence')

hist(z,100)
title('histogram of z axis scattering of the particles')
xlabel('location in [m] from distant to the near wall regard to cameras 1,2')
ylabel('number of occurence')

[f_x,x_x] = ksdensity(x); 
plot(x_x,f_x); 
title('Density estimation for location x') 
xlabel('location x')
ylabel(' pdf ')
grid

[f_y,x_y] = ksdensity(y_new); 
plot(x_y,f_y); 
title('Density estimation for location y') 
xlabel('location y')
ylabel(' pdf ')
grid

[f_z,x_z] = ksdensity(z); 
plot(x_z,f_z); 
title('Density estimation for location z') 
xlabel('location z')
ylabel(' pdf ')
grid

% plotting the w velocity pdf and sqewness
[f_vf,x_vf] = ksdensity(vf_vel);
plot(x_vf,f_vf)
title('histogram of  y velocity values')
xlabel('y velocity value')
ylabel('number of occurence')

mean(vf_vel)
std(vf_vel)
var(vf_vel)
skewness(vf_vel)
kurtosis(vf_vel)