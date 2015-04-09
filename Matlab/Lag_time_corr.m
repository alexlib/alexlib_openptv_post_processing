% a program to Alex of LAgrangian time correlation run all
data=readXUAPFiles('Directory_Name'); % write here the directory name
data=building_trajectories(data);
trajectories=show_length_of_trajectories(data);

%creating arrays for the data needed to process
xf = cat(1,data.xf);
yf = cat(1,data.yf);
zf = cat(1,data.zf);
uf = cat(1,data.uf);
vf = cat(1,data.vf);
wf = cat(1,data.wf);

%choosing  the trajectories to correlate
trajID = cat(2,data.traj);

%choosing trajectory to work on
% example trajID2Corr=177723;
trajID2Corr=find(trajectories>=10);


% Lagrangian  time correlation 
Traj_t=zeros(max(trajectories),8); % Creating maximum length array
% creating Traj array to process the data into more proper form
for l=1:length(trajID2Corr)
    traj = find(trajID == trajID2Corr(l));
    Traj=zeros(length(traj)-2,6);
    for i=1:(length(traj)-2) %inserting the data to the trajectorie array
        Traj(i,1)=xf(traj(i+1));Traj(i,2)=yf(traj(i+1));Traj(i,3)=zf(traj(i+1));
        Traj(i,4)=uf(traj(i+1));Traj(i,5)=vf(traj(i+1));Traj(i,6)=wf(traj(i+1));
    end
    % processing  for time correlation   
    for i=1:length(Traj) % index of jump
        for j=1:(length(Traj)-i) % index that run through the Traj array
            r_t=[Traj(i+j,1)-Traj(j,1),Traj(i+j,2)-Traj(j,2),Traj(i+j,3)-Traj(j,3)]; % connecting vector
            r_t=r_t/norm(r_t);
            % going on velocity in direction of connecting vector
            ua=[Traj(j,4),Traj(j,5),Traj(j,6)]*r_t'; %velocity normal to the conecting vector of the first point
            ub=[Traj(i+j,4),Traj(i+j,5),Traj(i+j,6)]*r_t';%velocity normal to the conecting vector of the second point
            % going on velocity in direction perpendicular to the connecting vector
            r_n=cross([Traj(j,4),Traj(j,5),Traj(j,6)],r_t);
            r_n=r_n/norm(r_n);
            va=[Traj(j,4),Traj(j,5),Traj(j,6)]*r_n';
            vb=[Traj(i+j,4),Traj(i+j,5),Traj(i+j,6)]*r_n';
            % inserting the data into the array
            Traj_t(i,1)=ua*ub+Traj_t(i,1);  %  ua*ub
            Traj_t(i,2)=ua*ua+Traj_t(i,2);  % ua*ua
            Traj_t(i,3)=ub*ub+Traj_t(i,3);  % ub*ub
            Traj_t(i,4)=va*vb+Traj_t(i,4);  % va*vb
            Traj_t(i,5)=va*va+Traj_t(i,5);  % va*va
            Traj_t(i,6)=vb*vb+Traj_t(i,6);  % vb*vb
            Traj_t(i,7)=ua*vb+Traj_t(i,7);  % ua*vb
            Traj_t(i,8)=1+Traj_t(i,8); % counting array
        end
    end
end


Traj_t_1=zeros(length(Traj_t)-3,7);
% normalization of the Traj_t vector
for i=1:length(Traj_t_1)
    Traj_t_1(i,1)=Traj_t(i,1)./Traj_t(i,8);
    Traj_t_1(i,2)=Traj_t(i,2)./Traj_t(i,8);
    Traj_t_1(i,3)=Traj_t(i,3)./Traj_t(i,8);
    Traj_t_1(i,4)=Traj_t(i,4)./Traj_t(i,8);
    Traj_t_1(i,5)=Traj_t(i,5)./Traj_t(i,8);
    Traj_t_1(i,6)=Traj_t(i,6)./Traj_t(i,8);
    Traj_t_1(i,7)=Traj_t(i,7)./Traj_t(i,8);
end

% Ruu correlation
Ruu=zeros(length(Traj_t_1),1);
for i=1:length(Traj_t_1)
    Ruu(i,1)=Traj_t_1(i,1)/ (Traj_t_1(i,2)* Traj_t_1(i,3))^0.5;
end
del_t=1/150; % here is the frame rate 
x=del_t:del_t:del_t*length(Traj_t_1);
plot(x,Ruu,'rs','MarkerSize',1)
title('Ruu versus time (full velocity of the particle) ')
xlabel('time ( delta_t=6.66*10^-^3) [s]')
ylabel('Ruu correlation')

%Rvv correlation
Rvv=zeros(length(Traj_t_1),1);
for i=1:length(Traj_t_1)
    Rvv(i,1)=Traj_t_1(i,4)/ (Traj_t_1(i,5)* Traj_t_1(i,6))^0.5;
end
del_t=1/150;
x=del_t:del_t:del_t*length(Traj_t_1);
plot(x,Rvv,'MarkerSize',1)
title('Rvv versus time (full velocity of the particle) ')
xlabel('time ( delta_t=6.66*10^-^3) [s]')
ylabel('Rvv correlation')