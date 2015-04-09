clear all
%load('/home/mark/Mark_Ron_Alex_22_02_15_spherical_grid/PTV_data/MAT_files/traj_water_freq_0_430_frames_10000_11858_vorticity.mat')
load('/home/mark/Mark_Ron_Alex_22_02_15_spherical_grid/PTV_data/MAT_files/traj_5ppm_freq_0_430_frames_10000_12584_vorticity.mat')
TrajLengthTreshold=20; %Minimum trajectory length to show
%Agitation origin 
x0 = -39.13; %[mm]
y0 = 50;%[mm]
z0 = 0;%[mm]
%Fps dt
Fps=100;

%Cat all 
x=cat(2,Traj.xf);
y=cat(2,Traj.yf);
z=cat(2,Traj.zf);

% u=cat(2,Traj.uf);
% v=cat(2,Traj.vf);
% w=cat(2,Traj.wf);
% 
% ax=cat(2,Traj.axf);
% ay=cat(2,Traj.ayf);
% az=cat(2,Traj.azf);

omegax=cat(2,Traj.omegax);
omegay=cat(2,Traj.omegay);
omegaz=cat(2,Traj.omegaz);

omegaz_sqr=omegaz.^2;

t=cat(2,Traj.t)+10000-4;
[~,~,r]=cart2sph(x-x0,y-y0,z-z0);
Idx=find(omegaz_sqr>=0.01);

%% 
close all
set(0,'defaulttextinterpreter','latex')
window=100;
z_depth=30;%[mm]
y_depth=-48;%[mm]
Theta_min=-35;%[degrees]
Theta_max=0;%[degrees] 

f5=figure(5);
set(f5,'Position',[0,0,1900,1000])

ProfileIdx=1;
R_Tnti=zeros(1,length(10000:window:11858));
clear R_for_profile Mean_Scalar
for FrameIdx=10000:window:11858
    
    clear fastScalarP slowScalarP
    
    frame1 = intersect(find(t <= FrameIdx+window),find(t >= FrameIdx));
    frame2 = intersect(find(z <= z_depth/2),find(z >= -z_depth/2));
    frame4 = intersect(frame1,frame2);      
  
    %Theta average profiles through the whole field
    R_Tolerance=4; %[mm]
    [~,~,R] = cart2sph(x-x0,y-y0,z-z0);
    R_PartitionNum=floor((max(R)-min(R))/R_Tolerance);
    R_QueryPoints=linspace(min(R),max(R),R_PartitionNum);
    FreeIdx=1;
    for RIdx=1:R_PartitionNum
        R_QueryIdxes1=find( R<( R_QueryPoints(RIdx)+(R_Tolerance/2) ) );
        R_QueryIdxes2=find( R>( R_QueryPoints(RIdx)-(R_Tolerance/2) ) );
        R_QueryIdxes=intersect(R_QueryIdxes1,R_QueryIdxes2);
        if ~isempty(R_QueryIdxes)
            
            Scalar_query_idxes=intersect(frame4,R_QueryIdxes); %Theta Phi and R are already from frame4
            %Scalar_per_radius=omegax(Scalar_query_idxes).^2+omegay(Scalar_query_idxes).^2+omegaz(Scalar_query_idxes).^2; %All scalar values on current radius
            Scalar_per_radius=omegaz(Scalar_query_idxes).^2; %All scalar values on current radius
            Mean_Scalar{ProfileIdx}(FreeIdx)=nanmean(Scalar_per_radius); %Average scalar value for this radius
          %  Mean_Scalar{ProfileIdx}=smooth(Mean_Scalar{ProfileIdx},3);
            R_for_profile{ProfileIdx}(FreeIdx)=R_QueryPoints(RIdx);
            clear Scalar_query_idxes
            FreeIdx=FreeIdx+1;
            
        end
        
        for InProfileIdx=1:length(Mean_Scalar{ProfileIdx})
            if (Mean_Scalar{ProfileIdx}(length(Mean_Scalar{ProfileIdx})-InProfileIdx+1)>=0.01) %Threshold 0.01 [1/sec^2]
                R_Tnti(ProfileIdx)=R_for_profile{ProfileIdx}(length(Mean_Scalar{ProfileIdx})-InProfileIdx+1);
                R_Tnti_time(ProfileIdx)=(FrameIdx-10000+(window/2))/100;
                break
            end
        end
        
        clear EnsZ
    end
    
    semilogy(R_for_profile{ProfileIdx},Mean_Scalar{ProfileIdx})
    hold on
    if exist('R_Tnti','var') && ~isempty(R_Tnti(ProfileIdx))
        y_line=linspace(0,1,10000);
        line(R_Tnti(ProfileIdx),y_line,'LineWidth',3);
    end
    xlim([0,150])
    xlabel('R[mm]','FontSize',20)
    ylabel('$\omega^2[1/sec^2]$','FontSize',20)
    title(['Frames=',num2str(FrameIdx),':',num2str(FrameIdx+window)]);
    drawnow
    clear Theta Phi R
   
    
    ProfileIdx=ProfileIdx+1;
end

figure(6)
%plot(R_Tnti_time,R_Tnti);
t2=(t-10000)/100;
scatter(r(Idx),t2(Idx),5,log(omegaz_sqr(Idx)))
hold on
    errorbar(smooth(R_Tnti),R_Tnti_time,repmat(0.5,1,length(R_Tnti)),'-o','LineWidth',3,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r','Color','r')
    herrorbar(smooth(R_Tnti),R_Tnti_time,repmat(2,1,length(R_Tnti)),'r')
xlabel('r[mm]','FontSize',20)
ylabel('t[sec]','FontSize',20)
title('Water','FontSize',20)
colorbar
hold off

