function traj = graphical_cleaning_traj(traj,mode)
% GRAPHICAL_CLEANING_TRAJ allows for the graphical input of a polygon.
% the points found in the polygon will be removed from the trajectory
% database. Works in 2D only, provide 'xy' or 'xz' or 'yz' input. Default
% is 'xy';
% 'xyz' is added using Delete pop-menu (right-button of your mouse)

if nargin < 2
    help graphical_cleaning_traj
end

save tmp traj

switch mode
    case{'xy','yx'}
        x = cat(1,traj.xf);
        y = cat(1,traj.yf);
    case{'xz','zx'}
        x = cat(1,traj.xf);
        y = cat(1,traj.zf);
    case{'yz','zy'}
        x = cat(1,traj.yf);
        y = cat(1,traj.zf);
    case{'xyz','3d'}
        plot_long_trajectories(traj,1);
        pause
        hl = get(gca,'children');
        keeptraj = zeros(size(hl));
        for i = 1:length(hl)
            keeptraj(i) = get(hl(i),'UserData'); 
        end
        traj = traj(keeptraj);
        return
    otherwise
        error('wrong second input: xy,yz,xz')
end

trajid = cat(1,traj.trajid);

figure, hold on
plot(x,y,'.');
[xv,yv] = ginput;
in = inpolygon(x,y,xv,yv);
id = find(in);
if ~isempty(id)
    plot(x(id),y(id),'ro');
    hold off
    cleanout = [];
    for i = 1:length(traj)
        if ismember(traj(i).trajid(1),trajid(id))
            cleanout = cat(1,cleanout,i);
        end
    end
    traj(cleanout) = [];
end


