%reading tracers data
% data1=readXUAPFiles('reut2');
% data1=building_trajectories_2(data1);
% newtraj_trac = xuap2traj(data1);
% save newtraj_trac

%reading particles data
% data2=readXUAPFiles('res_scene32');
% data2=building_trajectories_2(data2);
% newtraj_par = xuap2traj(data2);
% save newtraj_par


%load of data
% load newtraj_trac
% or load newtraj_par; newtraj_trac = newtraj_par;
trajlen = zeros(length(newtraj_trac),1);
for i = 1:length(trajlen)
    trajlen(i) = length(newtraj_trac(i).uf);
end;

longtraj = find(trajlen > 20);

%processing
taubins = 0:3:300;
duda = 0*taubins;
count = duda;
% cos_trac=[];cos_par=[];
for k = 1:length(longtraj); % temp solution ->length(newtraj_trac)
   traj = newtraj_trac(longtraj(k));
   for i = 1:length(traj.uf)-1
       for j = i:length(traj.uf)
           u = traj.uf;
           v = traj.vf;
           w = traj.wf;
           ax = traj.axf;
           ay = traj.ayf;
           az = traj.azf;
           cos_trac=[cos_trac,cosine([u(j) - u(i), v(j) - v(i), w(j) - w(i)],[ax(j) - ax(i), ay(j) - ay(i), az(j) - az(i)])];
%             tmp = (u(j)-u(i))*(ax(j)-ax(i))+...
%                 (v(j)-v(i))*(ay(j)-ay(i))+...
%                 (w(j)-w(i))*(az(j)-az(i));
%             tau = j - i;
%             ind = bindex(tau,taubins);
%             duda(ind) = duda(ind) + tmp;
%             count(ind) = count(ind) + 1;
       end
   end
end
figure, plot(taubins/150,duda./count);
xlabel('$\tau$, [sec]','Interpreter','Latex');
ylabel('$\langle \Delta u \cdot \Delta a \rangle (\tau) $','Interpreter','Latex')

% cosine([u(j) - u(i), v(j) - v(i), w(j) - w(i)],[ax(j) - ax(i), ay(j) - ay(i), az(j) - az(i)])
[trac_f]=ksdensity(cos_trac,[-1:0.001:1]);
[par_f]=ksdensity(cos_par,[-1:0.001:1]);
plot([-1:0.001:1],trac_f,'b');
figure;
plot([-1:0.001:1],par_f,'r');
