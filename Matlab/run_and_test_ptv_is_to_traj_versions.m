function run_and_test_ptv_is_to_traj_versions
% inside ptv2traj_v2 (should be renamed)
% we changed spline smoothing routines and played with their parameters
% the one chosen is from Beat's work: 21 points at max, otherwise 5 pieces of least
% square spline
traj = ptv_is_to_traj('/Users/alex/Documents/PTV/0104_water/');
trajlen = zeros(length(traj),1);
for i = 1:length(traj), trajlen(i) = length(traj(i).xf); end
ind = find(trajlen > 200);
figure
hold on
for j = 1:length(ind), i = ind(j);  quiver3(traj(i).xf,traj(i).yf,traj(i).zf,traj(i).uf,traj(i).vf,traj(i).wf); end;
ax = cat(1,traj.axf);
figure, plot(ax)
figure, nhist(ax,100)
