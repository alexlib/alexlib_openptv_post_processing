% acceleration pdf %
% reading the data %
tmp = readXUAPFiles('reut2');
ax=[];
ay=[];
az=[];
for i=2:length(tmp)  %the first file is zero file so we start at the second one
    ind= (find(abs(tmp(i).uf) > 10e-5 & abs(tmp(i).vf) > 10e-5 & abs(tmp(i).wf) > 10e-5));
    if length(ind) > 10
        ax=[ax,(tmp(i).axf(ind))']; % joining the ax arrays into one
        ay=[ay,(tmp(i).ayf(ind))'];
        az=[az,(tmp(i).azf(ind))'];
    end
end

%normalization by a_rms
ax_rms=(mean(ax.^2))^0.5;
ax=ax./ax_rms;
ay=ay./(mean(ay.^2))^(0.5);
az=az./(mean(az.^2))^(0.5);

%plotting

[f_x,x_x] = ksdensity(ax); 
[f_y,x_y] = ksdensity(ay);
[f_z,x_z] = ksdensity(az); 
figure;
hold on
plot(x_x,f_x,'r'); 
plot(x_y,f_y,'g'); 
plot(x_z,f_z,'b');
hold off
title('Density estimation ') 
xlabel('acceleration values')
ylabel(' pdf ')
grid

figure;
nhist(ax)
hold on
nhist(ay)
nhist(az)
