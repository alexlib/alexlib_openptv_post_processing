function ptv = update(traj, ptv, columns)

% global traj;
% global ptv;

for i=1:length(traj)
    ptv(traj(i,1),traj(i,2),columns+1)=1;
end