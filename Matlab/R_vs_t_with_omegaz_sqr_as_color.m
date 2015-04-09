clear all
%load('/home/mark/Mark_Ron_Alex_22_02_15_spherical_grid/PTV_data/MAT_files/traj_water_freq_0_430_frames_10000_11858_vorticity.mat')
load('/home/mark/Mark_Ron_Alex_22_02_15_spherical_grid/PTV_data/MAT_files/traj_5ppm_freq_0_430_frames_10000_12584_vorticity.mat')

%Agitation origin 
x0 = -39.13; %[mm]
y0 = 50;%[mm]
z0 = 0;%[mm]

Fps=100;

x=cat(2,Traj.xf);
y=cat(2,Traj.yf);
z=cat(2,Traj.zf);

omegax=cat(2,Traj.omegax);
omegay=cat(2,Traj.omegay);
omegaz=cat(2,Traj.omegaz);

omegaz_sqr=omegaz.^2;

t=cat(2,Traj.t)+10000-4;

[~,~,r]=cart2sph(x-x0,y-y0,z-z0);

Idx=find(omegaz_sqr>=0.01);

scatter(r(Idx),t(Idx),5,log(omegaz_sqr(Idx)))
