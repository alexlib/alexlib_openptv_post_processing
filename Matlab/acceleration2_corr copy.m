%  LAgrangian time acceleration correlation 
data=readXUAPFiles('reut2'); % write here the directory name
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


% Lagrangian  time correlation 
Traj_at=zeros(max(trajectories),3); % Creating maximum length array
% creating Traj array to process the data into more proper form
for l=1:length(trajID2Corr)
    traj = find(trajID == trajID2Corr(l));
    Traj=zeros(length(traj)-2,6);
    for i=1:(length(traj)-2) %inserting the data to the trajectorie array
        Traj(i,1)=xf(traj(i+1));Traj(i,2)=yf(traj(i+1));Traj(i,3)=zf(traj(i+1));
        Traj(i,4)=axf(traj(i+1));Traj(i,5)=ayf(traj(i+1));Traj(i,6)=azf(traj(i+1));
    end
    % processing  for time correlation   
    for i=1:length(Traj)+1 % index of jump
        for j=1:(length(Traj)-i) % index that run through the Traj array
            r_t=[Traj(i+j-1,1)-Traj(j,1),Traj(i+j-1,2)-Traj(j,2),Traj(i+j-1,3)-Traj(j,3)]; % connecting vector
            r_t=r_t/norm(r_t);
            % going on acceleration in direction of connecting vector
            aaf=[Traj(j,4),Traj(j,5),Traj(j,6)]; % full aceleration of the particle at first point
            abf=[Traj(i+j-1,4),Traj(i+j-1,5),Traj(i+j-1,6)];% full aceleration of the particle at second point
            aa=aaf*r_t'; %projection on connection vector of the first point
            ab=abf*r_t';%projection on connection vector of the second point
            % inserting the data into the array
            Traj_at(i,1)=aaf*abf'+Traj_at(i,1); % counting array
            Traj_at(i,2)=aa*ab'+Traj_at(i,2); % counting array
            Traj_at(i,3)=1+Traj_at(i,3); % counting array
        end
    end
end


% normalization of the Traj_at vector
for i=1:length(Traj_at)
    Traj_at(i,1)=Traj_at(i,1)./Traj_at(i,3);
    Traj_at(i,2)=Traj_at(i,2)./Traj_at(i,3);
end

% Rafaf correlation
Rafaf=Traj_at(:,1)./Traj_at(1,1);
del_t=1/150; % here is the frame rate 
x=0:del_t:del_t*(length(Traj_at)-1);
plot(x,Rafaf,'MarkerSize',3)
xlabel('time [s]')
ylabel('Ra_fa_f correlation')

%Raxax correlation
Raxax=Traj_at(:,2)./Traj_at(3,2);
Raxax=Raxax(3:length(Raxax));
del_t=1/150; % here is the frame rate 
x=0:del_t:del_t*(length(Traj_at)-3);
plot(x,Raxax,'MarkerSize',3)
xlabel('time [s]')
ylabel('Ra_xa_x correlation')



tou_0_ax=del_t*10;
in_ax=Raxax(1:10);
tou_int_ax=trapz(0:del_t:(tou_0_ax-del_t),in_ax);

tou_0_af=del_t*10;
in_af=Rafaf(1:10);
tou_int_af=trapz(0:del_t:(tou_0_af-del_t),in_af);