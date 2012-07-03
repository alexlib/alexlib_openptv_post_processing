function [output,newtraj] = link_trajectories_spline(traj,n)

Px = cat(1,traj.xf);
Py = cat(1,traj.yf);
Pz = cat(1,traj.zf);

figure
hold on
plot3(Px,Py,Pz,'ro','linewidth',1);
Tension=0.001;
output = [];
for k=1:length(Px)-3
    [MatOut3]=crdatnplusoneval([Px(k),Py(k),Pz(k)],[Px(k+1),Py(k+1),Pz(k+1)],[Px(k+2),Py(k+2),Pz(k+2)],[Px(k+3),Py(k+3),Pz(k+3)],Tension,n);
    output = cat(2,output,MatOut3);
end
t = cat(1,traj.t);
t = linspace(t(1),t(end),length(output));
% for i = 1:3
%     output(i,:) = smooth(t,output(i,:),0.1,'loess');
% end

plot3(output(1,:),output(2,:),output(3,:),'b','linewidth',2)

newtraj.xf = output(1,:);
newtraj.yf = output(2,:);
newtraj.zf = output(3,:);
newtraj.t = t;
newtraj.uf = gradient(newtraj.xf,newtraj.t/50);
newtraj.vf = gradient(newtraj.yf,newtraj.t/50);
newtraj.wf = gradient(newtraj.zf,newtraj.t/50);

newtraj.axf = gradient(newtraj.uf,newtraj.t/50);
newtraj.ayf = gradient(newtraj.vf,newtraj.t/50);
newtraj.azf = gradient(newtraj.wf,newtraj.t/50);






title('\bf 3D Cardinal Spline')
view(3);
box;
