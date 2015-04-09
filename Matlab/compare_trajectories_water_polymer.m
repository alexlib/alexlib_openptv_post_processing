%%  plot trajectories that show strong omega
ids = zeros(length(traj),1);

for id1 = 1:length(traj)
    omega1 = traj(id1).omegax.^2 + traj(id1).omegay.^2 + traj(id1).omegaz.^2;
    if mean(omega1) > .25 % std(omega1) > .2
        ids(id1) = 1;
    end
end

ids = find(ids);
plot_long_trajectories(traj(ids),20);

%%  compare two trajectories, 
id1 = ids(742); % 11421; % turbulent
id2 = ids(772); %13508; % non-turbulent

omega1 = traj(id1).omegax.^2 + traj(id1).omegay.^2 + traj(id1).omegaz.^2;
omega2 = traj(id2).omegax.^2 + traj(id2).omegay.^2 + traj(id2).omegaz.^2;

figure, plot(omega1); hold on; plot(omega2,'--');


%%  plot trajectories that show strong Reynolds stress 
ids = zeros(length(traj),1);

for id1 = 1:length(traj)
    % omega1 = traj(id1).omegax.^2 + traj(id1).omegay.^2 + traj(id1).omegaz.^2;
    Rs = traj(id1).uf.*traj(id1).vf + traj(id1).uf.*traj(id1).wf + traj(id1).wf.*traj(id1).vf;
    % Rs = Rs./traj(id1).uf.^2;
    Rs = Rs./std(traj(id1).uf)./std(traj(id1).vf);
    if mean(Rs) > 1e-5 % std(omega1) > .2
        ids(id1) = 1;
    end
end

ids = find(ids);
plot_long_trajectories(traj(ids),20);

