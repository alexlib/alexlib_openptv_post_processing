function plot_pdf_piv_data(matfile,ylocation,dt)
% PLOT_PDF_PIV_DATA(MATFILE,YLOCATION)
% Usage:
% 
% >> plot_pdf_piv_data('../2011/piv_grid_2011/90.mat',1200);
% 
% Author: Alex Liberzon (c) 2012 Turbulence Structure Laboratory
%


pix2mm = 8.983;

load(matfile)
[~,ind] = min(abs(y(:,1) - ylocation)); % find closest
u= u*(1000/pix2mm/dt);
v= v*(1000/pix2mm/dt);
tmpu = u(ind,:,1:500);
tmpv = v(ind,:,1:500);

figure, 
subplot(121)
nhist(tmpu(:),101);
xlabel('u, mm/s')

subplot(122)
nhist(tmpv(:),101);
xlabel('v, mm/s')

