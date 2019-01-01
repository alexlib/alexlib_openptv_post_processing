function traj_bin = plot_traj(traj, traj_bin, ptv, fig_id, time)

% global traj;
% global traj_bin;
% global ptv;

figure(fig_id); hold on;
for i=1:length(traj)
    x(i)=ptv(traj(i,1),traj(i,2),23);
    y(i)=ptv(traj(i,1),traj(i,2),24);
    z(i)=ptv(traj(i,1),traj(i,2),25);
    
    totalpix(i)=ptv(traj(i,1),traj(i,2),6);
    xpx(i)=ptv(traj(i,1),traj(i,2),7);
    ypx(i)=ptv(traj(i,1),traj(i,2),8);
    sumgrv(i)=ptv(traj(i,1),traj(i,2),9);
    
    m_u(i)= ptv(traj(i,1),traj(i,2),10);
    m_v(i)= ptv(traj(i,1),traj(i,2),11);
    m_w(i)= ptv(traj(i,1),traj(i,2),12);
    m_vel(i)= ptv(traj(i,1),traj(i,2),13);
    ux(i)= ptv(traj(i,1),traj(i,2),14);
    uy(i)= ptv(traj(i,1),traj(i,2),15);
    uz(i)= ptv(traj(i,1),traj(i,2),16);
    vx(i)= ptv(traj(i,1),traj(i,2),17);
    vy(i)= ptv(traj(i,1),traj(i,2),18);
    vz(i)= ptv(traj(i,1),traj(i,2),19);
    wx(i)= ptv(traj(i,1),traj(i,2),20);
    wy(i)= ptv(traj(i,1),traj(i,2),21);
    wz(i)= ptv(traj(i,1),traj(i,2),22);
    strain(i)= ptv(traj(i,1),traj(i,2),26);
    lambda1(i)= ptv(traj(i,1),traj(i,2),27);
    lambda2(i)= ptv(traj(i,1),traj(i,2),28);
    lambda3(i)= ptv(traj(i,1),traj(i,2),29);
    grad_meas(i)= ptv(traj(i,1),traj(i,2),30);
    size_meas1(i)= ptv(traj(i,1),traj(i,2),31);
    size_meas2(i)= ptv(traj(i,1),traj(i,2),32);
    
    
    id(i)=i;
    frame_i(i)=traj(i,1);
    j(i)=traj(i,2);
end
ok=1;

if 1<2 %min(x)<0.01
    plot3(x,y,totalpix,'b');
    %scatter3(x,y,z,5,totalpix);
    traj_bin=[traj_bin; [frame_i' j' id' x' y' z' totalpix' m_u' m_v' m_w' ux' uy' uz' vx' vy' vz' wx' wy' wz' strain' lambda1' lambda2' lambda3' grad_meas',xpx',ypx',size_meas1',size_meas2',sumgrv']];%%%%%%%%%%%%%%%%%%%%%
end

% figure(fig_id+1); hold on;
% if length(traj)>3 & min(x)<0.01
%     scatter3(x,y,z,5,xpx);
% end
% figure(fig_id+2); hold on;
% if length(traj)>3 & min(x)<0.01
%     scatter3(x,y,z,5,ypx);
% end
% figure(fig_id+3); hold on;
% if length(traj)>3 & min(x)<0.01
%     scatter3(x,y,z,5,sumgrv);
% end