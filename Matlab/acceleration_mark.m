% acceleration pdf %
% reading the data %
tmp = readXUAPFiles('I:\Experiments\PTV\test\res');
%%
ax=[];
ay=[];
az=[];
u=[];
v=[];
w=[];

% xmin=-0.05;
% xmax=0.05;
% ymin=-0.05;
% ymax=0.05;
% zmin=-0.05;
% zmax=0.05;

xmin=-0.04606;
xmax=0.02344;
ymin=-0.0327;
ymax=0.03803;
zmin=-0.05;
zmax=0.05;

for i=2:length(tmp)  %the first file is zero file so we start at the second one
    ind= (find(abs(tmp(i).uf) > 10e-5 & abs(tmp(i).vf) > 10e-5 & abs(tmp(i).wf) > 10e-5));
    if (size(ind,2)>0)
        for j=1:size(ind,1)
            if (tmp(i).xf(ind(j))>xmin) && (tmp(i).xf(ind(j))<xmax) && (tmp(i).yf(ind(j))>ymin) && (tmp(i).yf(ind(j))<ymax) && (tmp(i).zf(ind(j))>zmin) && (tmp(i).zf(ind(j))<zmax)
                ax=[ax,(tmp(i).axf(ind(j)))']; % joining the ax arrays into one
                ay=[ay,(tmp(i).ayf(ind(j)))'];
                az=[az,(tmp(i).azf(ind(j)))'];
                u=[u,(tmp(i).uf(ind(j)))'];
                v=[v,(tmp(i).vf(ind(j)))'];
                w=[w,(tmp(i).wf(ind(j)))'];
            end
        end
    end
end

%%

u_rms=u.^2;
u_mean=mean(u_rms);
u_rms=u_mean^0.5;
u_norm=u/u_rms;

v_rms=v.^2;
v_mean=mean(v_rms);
v_rms=v_mean^0.5;
v_norm=v/v_rms;

w_rms=w.^2;
w_mean=mean(w_rms);
w_rms=w_mean^0.5;
w_norm=w/w_rms;

ax_rms=ax.^2;
ax_mean=mean(ax_rms);
ax_rms=ax_mean^0.5;
ax_norm=ax/ax_rms;

ay_rms=ay.^2;
ay_mean=mean(ay_rms);
ay_rms=ay_mean^0.5;
ay_norm=ay/ay_rms;

az_rms=az.^2;
az_mean=mean(az_rms);
az_rms=az_mean^0.5;
az_norm=az/az_rms;

%%
%plotting

[f_x,x_x] = ksdensity(ax); 
[f_y,x_y] = ksdensity(ay);
[f_z,x_z] = ksdensity(az); 
figure;
hold on
plot(x_x,log10(f_x),'r'); 
plot(x_y,log10(f_y),'g'); 
plot(x_z,log10(f_z),'b');
hold off
title('Density estimation for location x') 
xlabel('acceleration values')
ylabel(' log(pdf) ')
grid

[r_x,x_x] = ksdensity(u); 
[r_y,x_y] = ksdensity(v);
[r_z,x_z] = ksdensity(w); 
figure;
hold on
plot(x_x,r_x,'r'); 
plot(x_y,r_y,'g'); 
plot(x_z,r_z,'b');
hold off
title('Density estimation for location x') 
xlabel('velocity values')
ylabel(' pdf ')
grid

%%
%plotting normalized graphs

[f_x,x_x] = ksdensity(ax_norm); 
[f_y,x_y] = ksdensity(ay_norm);
[f_z,x_z] = ksdensity(az_norm); 
figure;
hold on
plot(x_x,log10(f_x),'r'); 
plot(x_y,log10(f_y),'g'); 
plot(x_z,log10(f_z),'b');
hold off
title('Density estimation for location x') 
xlabel('acceleration values')
ylabel(' log(pdf) ')
grid

[r_x,x_x] = ksdensity(u_norm); 
[r_y,x_y] = ksdensity(v_norm);
[r_z,x_z] = ksdensity(w_norm); 
figure;
hold on
plot(x_x,r_x,'r'); 
plot(x_y,r_y,'g'); 
plot(x_z,r_z,'b');
hold off
title('Density estimation for location x') 
xlabel('velocity values')
ylabel(' pdf ')
grid

%%
% Plotting all particles
figure
hold on
for i=1:5000
    for j=1:length(tmp(i).xf)
        plot3(tmp(i).xf(j),tmp(i).yf(j),tmp(i).zf(j),'--o');
    end
end
hold off
grid

% Plotting all linked particles
figure
hold on
for i=1:10000
    for k=1:length(tmp(i).next)
        if tmp(i).next(k)>0
            plot(tmp(i).xf(k),tmp(i).yf(k),'--o');
        end
    end
end
hold off

% % Plotting all particles
% figure
% scatter(tmp.xf,tmp.yf);
%  
% % Plotting all linked particles
% figure
% for i=1:10000
%     for k=1:length(tmp(i).next)
%         if tmp(i).next(k)>0
%             scatter(tmp(i).xf(k),tmp(i).yf(k));
%         end
%     end
% end
