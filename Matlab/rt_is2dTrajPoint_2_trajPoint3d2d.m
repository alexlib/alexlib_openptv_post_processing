clear all
tic
first=19700;
last=20100;


Nmax_particles=300;
N_timesteps=200;
N_col=28;

F=repmat(NaN,[Nmax_particles N_col N_timesteps]);

for n=first:N_timesteps
    f1=['E:\PTV\Working_folder\WF_1\3d2d_info_diagnosis\rt_is2d_',num2Str(n),'.mat'];
    load (f1)
    si=size(a1_a,1);

    F(1:si,1,n)=a1_a;
    F(1:si,2,n)=a2_a;
    F(1:si,3,n)=a3_a;
    F(1:si,4,n)=a4_a;

    F(1:si,5,n)=b1_a;
    F(1:si,6,n)=b2_a;
    F(1:si,7,n)=b3_a;
    F(1:si,8,n)=b4_a;

    F(1:si,9,n)=c1_a;
    F(1:si,10,n)=c2_a;
    F(1:si,11,n)=c3_a;
    F(1:si,12,n)=c4_a;

    F(1:si,13,n)=d1_a;
    F(1:si,14,n)=d2_a;
    F(1:si,15,n)=d3_a;
    F(1:si,16,n)=d4_a;

    F(1:si,17,n)=e1_a;
    F(1:si,18,n)=e2_a;
    F(1:si,19,n)=e3_a;
    F(1:si,20,n)=e4_a;

    F(1:si,21,n)=f1_a;
    F(1:si,22,n)=f2_a;
    F(1:si,23,n)=f3_a;
    F(1:si,24,n)=f4_a;

    F(1:si,25,n)=g1_a;
    F(1:si,26,n)=g2_a;
    F(1:si,27,n)=g3_a;
    F(1:si,28,n)=g4_a;

end



for i=first:last-N_timesteps %-------------------'-'
                             % Why i starts from 2 ?

    if mod(i,5)==0
        i
    end

    name=['D:\PTV\Working_folder\WF_1\res\trajPoint.',num2Str(i)];
    f=load(name);
    s=size(f);
    clear x y z;
    if s(1,1)>0
        x=f(:,1);
        y=f(:,2);
        z=f(:,3);
        u=f(:,4);
        v=f(:,5);
        w=f(:,6);
        ax=f(:,7);
        ay=f(:,8);
        az=f(:,9);

        age=f(:,10);
        row=f(:,11);
        time=i+age; %----------------------


        a1_t=repmat(NaN,size(x));
        a2_t=repmat(NaN,size(x));
        a3_t=repmat(NaN,size(x));
        a4_t=repmat(NaN,size(x));

        b1_t=repmat(NaN,size(x));
        b2_t=repmat(NaN,size(x));
        b3_t=repmat(NaN,size(x));
        b4_t=repmat(NaN,size(x));

        c1_t=repmat(NaN,size(x));
        c2_t=repmat(NaN,size(x));
        c3_t=repmat(NaN,size(x));
        c4_t=repmat(NaN,size(x));

        d1_t=repmat(NaN,size(x));
        d2_t=repmat(NaN,size(x));
        d3_t=repmat(NaN,size(x));
        d4_t=repmat(NaN,size(x));

        e1_t=repmat(NaN,size(x));
        e2_t=repmat(NaN,size(x));
        e3_t=repmat(NaN,size(x));
        e4_t=repmat(NaN,size(x));

        f1_t=repmat(NaN,size(x));
        f2_t=repmat(NaN,size(x));
        f3_t=repmat(NaN,size(x));
        f4_t=repmat(NaN,size(x));

        g1_t=repmat(NaN,size(x));
        g2_t=repmat(NaN,size(x));
        g3_t=repmat(NaN,size(x));
        g4_t=repmat(NaN,size(x));

        for k=1:size(x,1)

            a1_t(k)=F(row(k),1,time(k)-i+1);
            a2_t(k)=F(row(k),2,time(k)-i+1);
            a3_t(k)=F(row(k),3,time(k)-i+1);
            a4_t(k)=F(row(k),4,time(k)-i+1);

            b1_t(k)=F(row(k),5,time(k)-i+1);
            b2_t(k)=F(row(k),6,time(k)-i+1);
            b3_t(k)=F(row(k),7,time(k)-i+1);
            b4_t(k)=F(row(k),8,time(k)-i+1);

            c1_t(k)=F(row(k),9,time(k)-i+1);
            c2_t(k)=F(row(k),10,time(k)-i+1);
            c3_t(k)=F(row(k),11,time(k)-i+1);
            c4_t(k)=F(row(k),12,time(k)-i+1);

            d1_t(k)=F(row(k),13,time(k)-i+1);
            d2_t(k)=F(row(k),14,time(k)-i+1);
            d3_t(k)=F(row(k),15,time(k)-i+1);
            d4_t(k)=F(row(k),16,time(k)-i+1);

            e1_t(k)=F(row(k),17,time(k)-i+1);
            e2_t(k)=F(row(k),18,time(k)-i+1);
            e3_t(k)=F(row(k),19,time(k)-i+1);
            e4_t(k)=F(row(k),20,time(k)-i+1);

            f1_t(k)=F(row(k),21,time(k)-i+1);
            f2_t(k)=F(row(k),22,time(k)-i+1);
            f3_t(k)=F(row(k),23,time(k)-i+1);
            f4_t(k)=F(row(k),24,time(k)-i+1);

            g1_t(k)=F(row(k),25,time(k)-i+1);
            g2_t(k)=F(row(k),26,time(k)-i+1);
            g3_t(k)=F(row(k),27,time(k)-i+1);
            g4_t(k)=F(row(k),28,time(k)-i+1);

        end

        f1=['E:\PTV\Working_folder\WF_1\3d2d_info_diagnosis\rt_is2dTrajPoint_2_trajPoint3d2d_',num2Str(i),'.mat'];

        save(f1,'x','y','z','u','v','w','ax','ay','az','age','row', ...
            'a1_t','a2_t','a3_t','a4_t','b1_t','b2_t','b3_t','b4_t', ...
            'c1_t','c2_t','c3_t','c4_t','d1_t','d2_t','d3_t','d4_t', ...
            'e1_t','e2_t','e3_t','e4_t','f1_t','f2_t','f3_t','f4_t', ...
            'g1_t','g2_t','g3_t','g4_t');



    else
        f1=['E:\PTV\Working_folder\WF_1\3d2d_info_diagnosis\rt_is2dTrajPoint_2_trajPoint3d2d_',num2Str(i),'.mat'];

        save(f1,'');
    end


    f1=['E:\PTV\Working_folder\WF_1\3d2d_info_diagnosis\rt_is2d_',num2Str(N_timesteps+i),'.mat'];
    load (f1)
    si=size(a1_a,1);

    F(1:si,1,N_timesteps+1)=a1_a;
    F(1:si,2,N_timesteps+1)=a2_a;
    F(1:si,3,N_timesteps+1)=a3_a;
    F(1:si,4,N_timesteps+1)=a4_a;

    F(1:si,5,N_timesteps+1)=b1_a;
    F(1:si,6,N_timesteps+1)=b2_a;
    F(1:si,7,N_timesteps+1)=b3_a;
    F(1:si,8,N_timesteps+1)=b4_a;

    F(1:si,9,N_timesteps+1)=c1_a;
    F(1:si,10,N_timesteps+1)=c2_a;
    F(1:si,11,N_timesteps+1)=c3_a;
    F(1:si,12,N_timesteps+1)=c4_a;

    F(1:si,13,N_timesteps+1)=d1_a;
    F(1:si,14,N_timesteps+1)=d2_a;
    F(1:si,15,N_timesteps+1)=d3_a;
    F(1:si,16,N_timesteps+1)=d4_a;

    F(1:si,17,N_timesteps+1)=e1_a;
    F(1:si,18,N_timesteps+1)=e2_a;
    F(1:si,19,N_timesteps+1)=e3_a;
    F(1:si,20,N_timesteps+1)=e4_a;

    F(1:si,21,N_timesteps+1)=f1_a;
    F(1:si,22,N_timesteps+1)=f2_a;
    F(1:si,23,N_timesteps+1)=f3_a;
    F(1:si,24,N_timesteps+1)=f4_a;

    F(1:si,25,N_timesteps+1)=g1_a;
    F(1:si,26,N_timesteps+1)=g2_a;
    F(1:si,27,N_timesteps+1)=g3_a;
    F(1:si,28,N_timesteps+1)=g4_a;

    F(:,:,1)=[];%-----------------------------------------


end
toc

