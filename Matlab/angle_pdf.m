m_vel=(m_u.^2+m_v.^2+m_w.^2).^0.5;
ij=find(m_vel==0);
abscostheta3=abs((m_eig_vec3x.*m_u+m_eig_vec3y.*m_v+m_eig_vec3z.*m_w)...
    ./(m_eig_vec3x.^2+m_eig_vec3y.^2+m_eig_vec3z.^2).^0.5...
    ./(m_u.^2+m_v.^2+m_w.^2).^0.5);
theta_radian3=(acos(abscostheta3));
theta_degree3=rad2deg(theta_radian3);
theta_degree3(ij)=NaN;

ind=find(isnan(theta_degree3(:)));
theta_degree3(ind)=[];
mu3=mean(abscostheta3)
sigma3=std(abscostheta3)
most_frequent3=mode(abscostheta3)
x3=floor(min(abscostheta3)):ceil(max(abscostheta3));
pdfNormal3 = normpdf(x3, mu3, sigma3);

ij=find(m_vel==0);
abscostheta2=abs((m_eig_vec2x.*m_u+m_eig_vec2y.*m_v+m_eig_vec2z.*m_w)...
    ./(m_eig_vec2x.^2+m_eig_vec2y.^2+m_eig_vec2z.^2).^0.5...
    ./(m_u.^2+m_v.^2+m_w.^2).^0.5);
theta_radian2=(acos(abscostheta2));
theta_degree2=rad2deg(theta_radian2);
theta_degree2(ij)=NaN;

ind=find(isnan(theta_degree2(:)));
theta_degree2(ind)=[];
mu2=mean(abscostheta2)
sigma2=std(abscostheta2)
most_frequent2=mode(abscostheta2)
x2=floor(min(abscostheta2)):ceil(max(abscostheta2));
pdfNormal2 = normpdf(x2, mu2, sigma2);

ij=find(m_vel==0);
abscostheta1=abs((m_eig_vec1x.*m_u+m_eig_vec1y.*m_v+m_eig_vec1z.*m_w)...
    ./(m_eig_vec1x.^2+m_eig_vec1y.^2+m_eig_vec1z.^2).^0.5...
    ./(m_u.^2+m_v.^2+m_w.^2).^0.5);
theta_radian1=(acos(abscostheta1));
theta_degree1=rad2deg(theta_radian1);
theta_degree1(ij)=NaN;

ind=find(isnan(theta_degree1(:)));
theta_degree1(ind)=[];
mu1=mean(abscostheta1)
sigma1=std(abscostheta1)
most_frequent1=mode(abscostheta1)
x1=floor(min(abscostheta1)):ceil(max(abscostheta1));
pdfNormal1 = normpdf(x1, mu1, sigma1);


figure;
subplot(411);
plot(x3,pdfNormal3,'b');grid;xlabel('pdf for angle between third eigen vector and velocity vector','FontName','Times New Roman','Fontsize',14)
subplot(412);
plot(x2,pdfNormal2,'g');grid;xlabel('pdf for angle between second eigen vector and velocity vector','FontName','Times New Roman','Fontsize',14)
subplot(413);
plot(x1,pdfNormal1,'r');grid;xlabel('pdf for angle between first eigen vector and velocity vector','FontName','Times New Roman','Fontsize',14)
subplot(414);
plot(x3,pdfNormal3,'b',x2,pdfNormal2,'g',x1,pdfNormal1,'r');grid;xlabel('Combined Pdf','FontName','Times New Roman','Fontsize',14)
legend('third','second','first',2)



plot(x3,pdfNormal3,'b');grid;hold on;
plot(x2,pdfNormal2,'g');grid;hold on;
plot(x1,pdfNormal1,'r');grid;
hold off
gtext('Combined Pdf')

%pdfNormal = normpdf(x, mu, sigma);
%plot(x, pdfNormal)
