%Conditional average

trajLen=traj_length(traj);
traj=traj(trajLen>5);

dt=250;%fps

ux=(-1)*cat(1,traj.uf)*dt/10;%units conversion mm/frames->cm/s
uy=(-1)*cat(1,traj.vf)*dt/10;
 


p=(-0.5:0.15:1)*max(ux);


clear ux_av;
clear uy_av;

j=1;
for i=1:length(p)

ind=find(ux<=p(i)+0.01 & ux>=p(i)-0.01 );
if ~ isempty(ind);
ux_av(j)=p(i)/max(ux);
uy_av(j)=mean(abs(uy(ind)))/max(ux);
j=j+1;
  
end
end

figure;
plot(ux_av,uy_av,'LineStyle','none','Marker','o','Color','g');
legend('scene 93');
xlabel('Ux/uxmax');
ylabel('Uy/Ux_max')
