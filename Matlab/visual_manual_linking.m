function tmptraj = manual_linking(newtraj,listoftraj)
if nargin < 2
    dt = 1/50;
end

tmptraj = newtraj(listoftraj);

    plot_trajectory_properties_vs_t(tmptraj,0);
    for i = 2:length(listoftraj)
        % newtraj(first) = combine_trajectories(newtraj(first),newtraj(listoftraj(i)));
        tmptraj(1) = link_trajectories_cubicspline(tmptraj([1,i]));
    end
    % newtraj(first) = link_trajectories_rbf(newtraj(first),smoothness,dt);
    
    tmptraj(2:end) = [];
    newtraj = ensure_order_traj(tmptraj);
    plot_trajectory_properties_vs_t(tmptraj(1),0)
end