uf = cat(2,traj.uf);
vf = cat(2,traj.vf);
wf = cat(2,traj.wf);


nhist(uf(:))
hold on
nhist(vf(:))
nhist(wf(:))

Rs = uf.*vf + uf.*wf + wf.*vf;

nhist(Rs(:),1000)