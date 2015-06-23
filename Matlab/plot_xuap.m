function quiver_xuap(xuap)
% QUIVER_XUAP(XUAP)
% plots as a 3D quiver the XUAP structure (a single frame)
% 

quiver3(xuap.x,xuap.y,xuap.z,xuap.u,xuap.v,xuap.w);

