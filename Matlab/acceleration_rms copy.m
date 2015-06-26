%  LAgrangian time acceleration correlation 
data=readXUAPFiles('res_scene32'); % write here the directory name
data=building_trajectories(data);
trajectories=show_length_of_trajectories(data);

%creating arrays for the data needed to process
xf = cat(1,data.xf);
yf = cat(1,data.yf);
zf = cat(1,data.zf);
axf = cat(1,data.axf);
ayf = cat(1,data.ayf);
azf = cat(1,data.azf);

%choosing  the trajectories to correlate
trajID = cat(2,data.traj);

%choosing trajectory to work on
% example trajID2Corr=177723;
trajID2Corr=find(trajectories>=20);

ax=[];ay=[];az=[];
% Lagrangian  time correlation 
Traj_t=zeros(max(trajectories),8); % Creating maximum length array
% creating Traj array to process the data into more proper form
for l=1:length(trajID2Corr)
    traj = find(trajID == trajID2Corr(l));
    Traj=zeros(length(traj)-2,6);
    for i=1:(length(traj)-2) %inserting the data to the trajectorie array
        Traj(i,1)=xf(traj(i+1));Traj(i,2)=yf(traj(i+1));Traj(i,3)=zf(traj(i+1));
        Traj(i,4)=axf(traj(i+1));Traj(i,5)=ayf(traj(i+1));Traj(i,6)=azf(traj(i+1));
    end
    % processing  for time correlation   
    ax=[ax;Traj(:,4)];
    ay=[ay;Traj(:,5)];
    az=[az;Traj(:,6)];
end

ax_rms=(mean(ax.^2))^0.5;
ay_rms=(mean(ay.^2))^0.5;
az_rms=(mean(az.^2))^0.5;
ax=ax./ax_rms;
ay=ay./ay_rms;
az=az./az_rms;

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