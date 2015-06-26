function res=do_script(name,first,frames,col)

clear all;
%frames=250;
%first=5700;

u_vect=zeros(10*17*frames,15);
v_vect=zeros(10*17*frames,15);
w_vect=zeros(10*17*frames,15);
u_vect=NaN;
v_vect=NaN;
w_vect=NaN;

last=first+frames-1;




count=zeros(15,1);
count=count+1;

for timestep=first:last 

    timestep
    filename_in=strcat(name,num2str(timestep));
    D=load(filename_in);

    y = D(:,2);
    u = D(:,4);
    v = D(:,5);
    w = D(:,6);

    for i = 8:3:50%5:3:17%
        ind = find(abs(y-i*0.001)<0.001 & isnan(u)==0);
        le=length(ind);
        i_ind=(i-5)/3;%(i-2)/3;%
        u_vect(count(i_ind):count(i_ind)+le-1,i_ind)=u(ind);
        v_vect(count(i_ind):count(i_ind)+le-1,i_ind)=v(ind);
        w_vect(count(i_ind):count(i_ind)+le-1,i_ind)=w(ind);
        count(i_ind)=count(i_ind)+le;
    end
end

for i = 8:3:50%5:3:17%
    i_ind=(i-5)/3;%(i-2)/3;%
    uu_vect = (u_vect(:,i_ind)-nanmean(u_vect(:,i_ind))).*(u_vect(:,i_ind)-nanmean(u_vect(:,i_ind)));
    vv_vect = (v_vect(:,i_ind)-nanmean(v_vect(:,i_ind))).*(v_vect(:,i_ind)-nanmean(v_vect(:,i_ind)));
    ww_vect = (w_vect(:,i_ind)-nanmean(w_vect(:,i_ind))).*(w_vect(:,i_ind)-nanmean(w_vect(:,i_ind)));
    uv_vect = (u_vect(:,i_ind)-nanmean(u_vect(:,i_ind))).*(v_vect(:,i_ind)-nanmean(v_vect(:,i_ind)));
    uw_vect = (u_vect(:,i_ind)-nanmean(u_vect(:,i_ind))).*(w_vect(:,i_ind)-nanmean(w_vect(:,i_ind)));
    vw_vect = (v_vect(:,i_ind)-nanmean(v_vect(:,i_ind))).*(w_vect(:,i_ind)-nanmean(w_vect(:,i_ind)));
    
%     uu_vect = (u_vect(:,i_ind)).*(u_vect(:,i_ind));
%     vv_vect = (v_vect(:,i_ind)).*(v_vect(:,i_ind));
%     ww_vect = (w_vect(:,i_ind)).*(w_vect(:,i_ind));
%     uv_vect = (u_vect(:,i_ind)).*(v_vect(:,i_ind));
%     uw_vect = (u_vect(:,i_ind)).*(w_vect(:,i_ind));
%     vw_vect = (v_vect(:,i_ind)).*(w_vect(:,i_ind));
    
    uu = nanmean(uu_vect);
    vv = nanmean(vv_vect);
    ww = nanmean(ww_vect);
    uv = nanmean(uv_vect);
    uw = nanmean(uw_vect);
    vw = nanmean(vw_vect);

    Rs = [uu,uv,uw;uv,vv,vw;uw,vw,ww];

    [av_II(i_ind),av_III(i_ind)] = anisotropy(Rs);
end



hf=figure(1);
axis([-.04 .08 0 .35]);
hold on, grid on, box on
plot_Lumley_triangle(hf);
xlabel('III')
ylabel('II')
plot(av_III,-av_II,'r');
scatter(av_III(15),-av_II(15),'r')

figure(2);hold on;
y=8:3:50;%5:3:17%
y=60-y;
plot(y,-av_II,'r')
xlabel('y');
ylabel('II');

figure(3);hold on;
y=8:3:50;%5:3:17%
y=60-y;
plot(y,av_III,'r')
xlabel('y');
ylabel('III');





