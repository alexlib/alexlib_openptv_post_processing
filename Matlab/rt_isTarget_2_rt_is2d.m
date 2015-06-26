clear all;
tic

first=19700;
last=20100;

for n=first:last
    if mod(n,5)==0
        n
    end

    
    name_rt=['D:\PTV\Working_folder\WF_1\res\rt_is.',num2Str(n)];
    fid_rt = fopen(name_rt);
    number_rt = fscanf(fid_rt, '%i', [1 1]);
    f=fscanf(fid_rt, '%i %f %f %f %g %g %g %g', [8 inf]);
    f=f';
    fclose(fid_rt);
    si_rt=size(f);

    if si_rt(1,1)>0
        id=round(f(:,1));
        x=f(:,2);
        y=f(:,3);
        z=f(:,4);
        c1=round(f(:,5));
        c2=round(f(:,6));
        c3=round(f(:,7));
        c4=round(f(:,8));
    end

    %%%%%%%%%%%% Target files %%%%%%%%%%
    name_t1=['D:\PTV\Working_folder\WF_1\img_for_diagnosis\cam1.',num2Str(n)];
    name_t2=['D:\PTV\Working_folder\WF_1\img_for_diagnosis\cam2.',num2Str(n)];
    name_t3=['D:\PTV\Working_folder\WF_1\img_for_diagnosis\cam3.',num2Str(n)];
    name_t4=['D:\PTV\Working_folder\WF_1\img_for_diagnosis\cam4.',num2Str(n)];

   


    for k=1:4

        switch k
            case 1
                name=name_t1;
            case 2
                name=name_t2;
            case 3
                name=name_t3;
            case 4
                name=name_t4;
        end

        name=[name '_targets'];

        fid_target = fopen(name);
        number_target = fscanf(fid_target, '%i', [1 1]);
        g=fscanf(fid_target, '%i %f %f %g %g %g %g %g', [8 inf]);
        g=g';
        fclose(fid_target);
        si_target=size(g);


        switch k
            case 1
                if si_target(1,1)>0
                    idt1=round(g(:,1));
                    xt1=g(:,2);
                    yt1=g(:,3);
                    n_px1=g(:,4);
                    nx_px1=g(:,5);
                    ny_px1=g(:,6);
                    sumg1=g(:,7);
                    ptr_rt1=g(:,8);
                end
            case 2
                if si_target(1,1)>0
                    idt2=round(g(:,1));
                    xt2=g(:,2);
                    yt2=g(:,3);
                    n_px2=g(:,4);
                    nx_px2=g(:,5);
                    ny_px2=g(:,6);
                    sumg2=g(:,7);
                    ptr_rt2=g(:,8);
                end
            case 3
                if si_target(1,1)>0
                    idt3=round(g(:,1));
                    xt3=g(:,2);
                    yt3=g(:,3);
                    n_px3=g(:,4);
                    nx_px3=g(:,5);
                    ny_px3=g(:,6);
                    sumg3=g(:,7);
                    ptr_rt3=g(:,8);
                end
            case 4
                if si_target(1,1)>0
                    idt4=round(g(:,1));
                    xt4=g(:,2);
                    yt4=g(:,3);
                    n_px4=g(:,4);
                    nx_px4=g(:,5);
                    ny_px4=g(:,6);
                    sumg4=g(:,7);
                    ptr_rt4=g(:,8);
                end
        end
    end


   


    row_t_cam1=c1(1:si_rt(1,1))+1;
    row_t1=row_t_cam1;
    i0_1=find(row_t_cam1==0);
    row_t1(i0_1)=1;%1 place kora hochhe karon 0 gulo k rekhe dile pore error show korbe (real + int or logical error)...
    % pore abar ei 1 gulo kei nan dea replace kore deya
    % hochhe.

    a1_a=n_px1(row_t1);
    a1_a(i0_1)=NaN;% Ekhane 1 er karone assign kora value k nan dea replace kora holo
    b1_a=nx_px1(row_t1);
    b1_a(i0_1)=NaN;
    c1_a=ny_px1(row_t1);
    c1_a(i0_1)=NaN;
    d1_a=sumg1(row_t1);
    d1_a(i0_1)=NaN;
    e1_a=ptr_rt1(row_t1);
    e1_a(i0_1)=NaN;
    f1_a=xt1(row_t1);
    f1_a(i0_1)=NaN;
    g1_a=yt1(row_t1);
    g1_a(i0_1)=NaN;




    row_t_cam2=c2(1:si_rt(1,1))+1;
    row_t2=row_t_cam2;
    i0_2=find(row_t_cam2==0);
    row_t2(i0_2)=1;

    a2_a=n_px2(row_t2);
    a2_a(i0_2)=NaN;
    b2_a=nx_px2(row_t2);
    b2_a(i0_2)=NaN;
    c2_a=ny_px2(row_t2);
    c2_a(i0_2)=NaN;
    d2_a=sumg2(row_t2);
    d2_a(i0_2)=NaN;
    e2_a=ptr_rt2(row_t2);
    e2_a(i0_2)=NaN;
    f2_a=xt2(row_t2);
    f2_a(i0_2)=NaN;
    g2_a=yt2(row_t2);
    g2_a(i0_2)=NaN;

    row_t_cam3=c3(1:si_rt(1,1))+1;
    row_t3=row_t_cam3;
    i0_3=find(row_t_cam3==0);
    row_t3(i0_3)=1;

    a3_a=n_px3(row_t3);
    a3_a(i0_3)=NaN;
    b3_a=nx_px3(row_t3);
    b3_a(i0_3)=NaN;
    c3_a=ny_px3(row_t3);
    c3_a(i0_3)=NaN;
    d3_a=sumg3(row_t3);
    d3_a(i0_3)=NaN;
    e3_a=ptr_rt3(row_t3);
    e3_a(i0_3)=NaN;
    f3_a=xt3(row_t3);
    f3_a(i0_3)=NaN;
    g3_a=yt3(row_t3);
    g3_a(i0_3)=NaN;


    row_t_cam4=c4(1:si_rt(1,1))+1;
    row_t4=row_t_cam4;
    i0_4=find(row_t_cam4==0);
    row_t4(i0_4)=1;

    a4_a=n_px4(row_t4);
    a4_a(i0_4)=NaN;
    b4_a=nx_px4(row_t4);
    b4_a(i0_4)=NaN;
    c4_a=ny_px4(row_t4);
    c4_a(i0_4)=NaN;
    d4_a=sumg4(row_t4);
    d4_a(i0_4)=NaN;
    e4_a=ptr_rt4(row_t4);
    e4_a(i0_4)=NaN;
    f4_a=xt4(row_t4);
    f4_a(i0_4)=NaN;
    g4_a=yt4(row_t4);
    g4_a(i0_4)=NaN;

    f1=['E:\PTV\Working_folder\WF_1/3d2d_info_diagnosis\rt_is2d_',num2Str(n),'.mat'];
    save(f1,'a1_a','a2_a','a3_a','a4_a','b1_a','b2_a','b3_a','b4_a', ...
        'c1_a','c2_a','c3_a','c4_a','d1_a','d2_a','d3_a','d4_a', ...
        'e1_a','e2_a','e3_a','e4_a', ...
        'f1_a','f2_a','f3_a','f4_a', ...
        'g1_a','g2_a','g3_a','g4_a');

end
toc
