
figure
plot(xcorr(newtraj(187).uf,'coeff'),'r-')
hold on
plot(xcorr(newtraj(187).vf,'coeff'),'g--')
plot(xcorr(newtraj(187).wf,'coeff'),'b-.')


figure
plot(xcorr(newtraj(187).axf(3:end-2),'coeff'),'r-')
hold on
plot(xcorr(newtraj(187).ayf(3:end-2),'coeff'),'g--')
plot(xcorr(newtraj(187).azf(3:end-2),'coeff'),'b-.')