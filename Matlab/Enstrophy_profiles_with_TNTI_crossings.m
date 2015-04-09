clear all

% Solution='Water';
% load('D:/Mark_ptv/MAT_files/traj_water_freq_0_430_frames_10000_11858_vorticity.mat')

Solution='Polymer';
load('../../MAT_files/traj_5ppm_freq_0_430_frames_10000_12584_vorticity.mat')


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

trajid=cat(2,Traj.trajid);

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

t=cat(2,Traj.t)+10000-4; %In frames
t2=(t-10000)/100; %In seconds
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

% f1=figure(1);
% set(f1,'Position',[0,0,1900,1000])

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
    
%     semilogy(R_for_profile{ProfileIdx},Mean_Scalar{ProfileIdx})
%     hold on
%     if exist('R_Tnti','var') && ~isempty(R_Tnti(ProfileIdx))
%         y_line=linspace(0,1,10000);
%         line(R_Tnti(ProfileIdx),y_line,'LineWidth',3);
%     end
%     xlim([0,150])
%     xlabel('R[mm]','FontSize',20)
%     ylabel('$\omega^2[1/sec^2]$','FontSize',20)
%     title(['Frames=',num2str(FrameIdx),':',num2str(FrameIdx+window)]);
%     drawnow

    clear Theta Phi R
    ProfileIdx=ProfileIdx+1;
end

InterfaceTimeError=(window/100)/2; %plus-minus this value
InterfaceRadiusError=R_Tolerance/2; %plus-minus this value

Smooth_R_Tnti=smooth(R_Tnti);
% f2=figure(2);
% set(f2,'Position',[0,0,1900,1000])
% %plot(R_Tnti_time,R_Tnti);
% scatter(r(Idx),t2(Idx),5,log(omegaz_sqr(Idx)))
% hold on
% clear R_Tnti
% errorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceTimeError,1,length(Smooth_R_Tnti)),'-o','LineWidth',3,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r','Color','r')
% herrorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceRadiusError,1,length(Smooth_R_Tnti)),'r')
% xlabel('r[mm]','FontSize',20)
% ylabel('t[sec]','FontSize',20)
% title(Solution,'FontSize',20)
% colorbar
% hold off

%Find particles which are close to interface. And their omegaz_sqr>0.01
for R_Tnti_Idx=1:length(Smooth_R_Tnti)
    Cut1=find(r > (Smooth_R_Tnti(R_Tnti_Idx)-InterfaceRadiusError));
    Cut2=find(r < (Smooth_R_Tnti(R_Tnti_Idx)+InterfaceRadiusError));
    CloseToTNTI_r=intersect(Cut1,Cut2);
    
    Cut3=find(t2 > (R_Tnti_time(R_Tnti_Idx)-InterfaceTimeError));
    Cut4=find(t2 < (R_Tnti_time(R_Tnti_Idx)+InterfaceTimeError));
    CloseToTNTI_t2=intersect(Cut3,Cut4);
    
    CloseToTNTI{R_Tnti_Idx}=intersect(CloseToTNTI_r,CloseToTNTI_t2);
    CloseToTNTI_w_thresh{R_Tnti_Idx}=intersect(CloseToTNTI{R_Tnti_Idx},Idx);
end

% f3=figure(3);
% set(f3,'Position',[0,0,1900,1000])
% %Plot particles which are close to interface. And their omegaz_sqr>0.01
% for R_Tnti_Idx=1:length(Smooth_R_Tnti)
%     scatter(r(CloseToTNTI{R_Tnti_Idx}),t2(CloseToTNTI{R_Tnti_Idx}),10,log(omegaz_sqr(CloseToTNTI{R_Tnti_Idx})))
%     hold on
%     errorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceTimeError,1,length(Smooth_R_Tnti)),'-o','LineWidth',3,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r','Color','r')
%     herrorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceRadiusError,1,length(Smooth_R_Tnti)),'r')
% end
% xlabel('r[mm]','FontSize',20)
% ylabel('t[sec]','FontSize',20)
% title([Solution,', all trajectories that cross the TNTI.'],'FontSize',20)
% c=colorbar;
% xlabel(c,'$log((\omega_{z})^2)$','FontSize',20)
% hold off

%Find Trajectory IDs which are close to interface. And their omegaz_sqr>0.01
for R_Tnti_Idx=1:length(Smooth_R_Tnti)
    TrajID_close_to_TNTI{R_Tnti_Idx}=trajid(CloseToTNTI{R_Tnti_Idx});
end


% %Plot particles which are close to interface.
% f4=figure(4);
% set(f4,'Position',[0,0,1900,1000])
% InterfaceTimeError=(window/100)/2; %plus-minus this value
% InterfaceRadiusError=R_Tolerance/2; %plus-minus this value
% for R_Tnti_Idx=1:length(Smooth_R_Tnti)
%     for TrajIdx=1:length(TrajID_close_to_TNTI{R_Tnti_Idx});
%         Traj_r=sqrt((Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).xf-x0).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).yf-y0).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).zf-z0).^2);
%         Traj_t=(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).t-4)/100;
%         Traj_omegaz_sqr=(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).omegaz).^2;
%         %Kinetic energy per unit mass 3D
%         Traj_k=((Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).uf).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).vf).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).wf).^2)/2;
%         scatter(Traj_r,Traj_t,10,log(Traj_k))
%         hold on
%         clear Traj_t Traj_r Traj_omegaz_sqr
%     end
%     errorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceTimeError,1,length(Smooth_R_Tnti)),'-o','LineWidth',3,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r','Color','r')
%     herrorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceRadiusError,1,length(Smooth_R_Tnti)),'r')
%     
% end
% xlabel('r[mm]','FontSize',20)
% ylabel('t[sec]','FontSize',20)
% title([Solution,', all trajectories that cross the TNTI.'],'FontSize',20)
% c=colorbar;
% xlabel(c,'$log((\omega_{z})^2)$','FontSize',20)
% hold off
%saveas(f4,[Solution,'Whole_trajectories_all_color_log(KE)'],'fig');
%close f4

% %Plot only those that exit the patch
% f5=figure(5);
% set(f5,'Position',[0,0,1900,1000])
% for R_Tnti_Idx=1:length(Smooth_R_Tnti)
%     for TrajIdx=1:length(TrajID_close_to_TNTI{R_Tnti_Idx});
%         Traj_r=sqrt((Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).xf-x0).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).yf-y0).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).zf-z0).^2);
%         Traj_t=(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).t-4)/100;
%         if mean(diff(Traj_r)/diff(Traj_t))>0
%             
%             %Enstrophy in z direction
%             %Traj_omegaz_sqr=(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).omegaz).^2;
%             %scatter(Traj_r,Traj_t,5,Traj_omegaz_sqr)
%             
%             %Kinetic energy per unit mass 3D
%             Traj_k=((Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).uf).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).vf).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).wf).^2)/2;
%             scatter(Traj_r,Traj_t,10,log(Traj_k))
%             hold on
%             clear Traj_t Traj_r Traj_omegaz_sqr Traj_k
%         end
%     end
%     errorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceTimeError,1,length(Smooth_R_Tnti)),'-o','LineWidth',3,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r','Color','r')
%     herrorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceRadiusError,1,length(Smooth_R_Tnti)),'r')
% end
% xlabel('r[mm]','FontSize',20)
% ylabel('t[sec]','FontSize',20)
% title([Solution,', only those exiting the patch.'],'FontSize',20);
% c=colorbar;
% xlabel(c,'$log((u^{2}+v^{2}+w^{2})/2))$','FontSize',20)
% hold off
%saveas(f5,[Solution,'_whole_trajectories_only_exiting_color_log(KE)'],'fig');
%close f5


%Plot only those that enter the patch
f6=figure(6);
set(f6,'Position',[0,0,1900,1000])
for R_Tnti_Idx=1:length(Smooth_R_Tnti)
    for TrajIdx=1:length(TrajID_close_to_TNTI{R_Tnti_Idx});
        Traj_r=sqrt((Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).xf-x0).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).yf-y0).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).zf-z0).^2);
        Traj_t=(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).t-4)/100;
        if mean(diff(Traj_r)/diff(Traj_t))<0
            
            %Enstrophy in z direction
            %Traj_omegaz_sqr=(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).omegaz).^2;
            %scatter(Traj_r,Traj_t,5,Traj_omegaz_sqr)
            
            %Kinetic energy per unit mass 3D
            Traj_k=((Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).uf).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).vf).^2+(Traj(TrajID_close_to_TNTI{R_Tnti_Idx}(TrajIdx)).wf).^2)/2;
            scatter(Traj_r,Traj_t,10,log(Traj_k))
            hold on
            clear Traj_t Traj_r Traj_omegaz_sqr Traj_k
        end
    end
    errorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceTimeError,1,length(Smooth_R_Tnti)),'-o','LineWidth',3,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r','Color','r')
    herrorbar(Smooth_R_Tnti,R_Tnti_time,repmat(InterfaceRadiusError,1,length(Smooth_R_Tnti)),'r')
end
xlabel('r[mm]','FontSize',20)
ylabel('t[sec]','FontSize',20)
title([Solution,', only those entering the patch.'],'FontSize',20);
c=colorbar;
xlabel(c,'$log((u^{2}+v^{2}+w^{2})/2))$','FontSize',20)
hold off
% saveas(f6,[Solution,'_whole_trajectories_only_entering_color_log(KE)'],'fig');
% close f6

% clear all


