clear all; 
close all; 


%plot velocity relations 

 load('H:\PTV\exp_031110\trajectories\traj92.mat');
 load('H:\PTV\exp_031110\trajectories\traj93.mat');
 load('H:\PTV\exp_171010\trajectories\scene63_traj.mat');
 
 dt63=250;%fps
 dt93=125;
 dt92=100;
 
 
 ax93=(-1)*cat(1,traj93.axf)*(dt93)^2/10;% units conversion : mm/frames->cm/s2
 ay93=(-1)*cat(1,traj93.ayf)*(dt93)^2/10;

 
 ax92=(-1)*cat(1,traj92.axf)*(dt92)^2/10;
 ay92=(-1)*cat(1,traj92.ayf)*(dt92)^2/10;
 
 ax63=(-1)*cat(1,traj.axf)*(dt63)^2/10;
 ay63=(-1)*cat(1,traj.ayf)*(dt63)^2/10;
 
 figure;
 scatter(ax93/max(ax93),ay93/max(ax93),'b','.');hold on;                                             
 scatter(ax63/max(ax63),ay63/max(ax63),'r','.');
 xlabel('a_x/a_xmax');
 ylabel('a_y/a_xmax');
 legend('f=0.25,q_av=4.31,no reverse flow','f=0.25,q_av=1.54, revrse flow');
axis([-0.6 0.6 -0.4 0.4]);
 
 figure;
 subplot(1,2,1);
 scatter(ax93/max(ax93),ay93/max(ax93),'b','.');
 xlabel('a_x/a_xmax');
 ylabel('a_y/a_xmax');
 title('f=0.25 Hz, q_av=1.54');
axis([-0.6 0.6 -0.4 0.4]);
 
 subplot(1,2,2);
 scatter(ax92/max(ax92),ay92/max(ax92),'r','.');
 xlabel('a_x/a_xmax');
 ylabel('a_y/a_xmax');
 title('f=0.5 Hz, q_av=1.82');
 axis([-0.6 0.6 -0.4 0.4]);
 
 
 
 