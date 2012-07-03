clear all; 
close all; 


%plot velocity relations 

 load('H:\PTV\exp_031110\trajectories\traj92.mat');
 load('H:\PTV\exp_031110\trajectories\traj93.mat');
 load('H:\PTV\exp_171010\trajectories\scene63_traj.mat');
 
 dt63=250;%fps
 dt93=125;
 dt92=100;
 
 
 ux93=(-1)*cat(1,traj93.uf)*dt93/10;% units conversion : mm/frames->cm/s
 uy93=(-1)*cat(1,traj93.vf)*dt93/10;

 
 ux92=(-1)*cat(1,traj92.uf)*dt92/10;
 uy92=(-1)*cat(1,traj92.vf)*dt92/10;
 
 ux63=(-1)*cat(1,traj.uf)*dt63/10;
 uy63=(-1)*cat(1,traj.vf)*dt63/10;
 
 figure;
 scatter(ux93/max(ux93),uy93/max(ux93),'b','.');hold on;                                             
 scatter(ux63/max(ux63),uy63/max(ux63),'r','.');
 
 figure;
 subplot(1,2,1);
 scatter(ux93/max(ux93),uy93/max(ux93),'b','.');
 xlabel('u_x/u_xmax');
 ylabel('u_y/u_xmax');
 title('f=0.25 Hz, q_av=1.54');
 axis([-0.6 0.6 -0.4 0.4]);
 
 subplot(1,2,2);
 scatter(ux92/max(ux92),uy92/max(ux92),'r','.');
 xlabel('u_x/u_xmax');
 ylabel('u_y/u_xmax');
 title('f=0.5 Hz, q_av=1.82');
  axis([-0.6 0.6 -0.4 0.4]);
 