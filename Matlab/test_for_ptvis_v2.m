sss=zeros(1,8);
%first=1300;
%last=1300;
write_rt=0;

for m=1:1
    direc2=['F:\Exp_07112010\Exp',num2str(m),'_07112010','\image'];
    direc=['F:\Exp_07112010\Exp',num2str(m),'_07112010','\res'];
    d = [direc,'/ptv_is*'];
    
    
    dd=dir(d)
    filename1=dd(1).name
    %filename2=dd(end).name
    
    
    
    first = eval(filename1(8:end));
    if first==1
        last = size(dd,1);
    else
        filename2=dd(end).name
        last=eval(filename2(8:end))
    end
    
  
    
    for n=first:last
        if mod(n,5)==0
            n
        end
        
        
        name_rt=[direc,'\rt_is.',num2Str(n)];
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
            %c3=round(f(:,7));
            %c4=round(f(:,8));
        end
        
        %%%%%%%%%%%% Target files %%%%%%%%%%
        name_t1=[direc2,'\filt_Cam1.',num2Str(n)];
        name_t2=[direc2,'\filt_Cam2.',num2Str(n)];
        
        
        
        for k=1:2
            
            switch k
                case 1
                    name=name_t1;
                case 2
                    name=name_t2;
                    
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
                    
            end
        end
 
        row_t_cam1=c1(1:si_rt(1,1))+1;
        row_t1=row_t_cam1;
        i0_1=find(row_t_cam1==0);
        row_t1(i0_1)=1;%1 place kora hochhe karon 0 gulo k rekhe dile pore error show korbe (real + int or logical error)...
        % pore abar ei 1 gulo kei nan dea replace kore deya
        % hochhe.
        
        total_pix1=n_px1(row_t1);
        total_pix1(i0_1)=NaN;% Ekhane 1 er karone assign kora value k nan dea replace kora holo
        x_pix1=nx_px1(row_t1);
        x_pix1(i0_1)=NaN;
        y_pix1=ny_px1(row_t1);
        y_pix1(i0_1)=NaN;
        grv1=sumg1(row_t1);
        grv1(i0_1)=NaN;
        rt_id1=ptr_rt1(row_t1);
        rt_id1(i0_1)=NaN;
  
        row_t_cam2=c2(1:si_rt(1,1))+1;
        row_t2=row_t_cam2;
        i0_2=find(row_t_cam2==0);
        row_t2(i0_2)=1;
        
        total_pix2=n_px2(row_t2);
        total_pix2(i0_2)=NaN;
        x_pix2=nx_px2(row_t2);
        x_pix2(i0_2)=NaN;
        y_pix2=ny_px2(row_t2);
        y_pix2(i0_2)=NaN;
        grv2=sumg2(row_t2);
        grv2(i0_2)=NaN;
        rt_id2=ptr_rt2(row_t2);
        rt_id2(i0_2)=NaN;
        
        
        
        total_pix=mean([total_pix1 total_pix2],2);
        xpix=mean([x_pix1 x_pix2],2);
        ypix=mean([y_pix1 y_pix2],2);
        grv=mean([grv1 grv2],2);
 
        if write_rt==1
            fid = fopen(['E:\PTV\Working_folder\Exp32b_09112010\res\rt_is_v2.',num2Str(n)],'w');
            fprintf(fid,'%d\n',length(total_pix));
            fprintf(fid,'%4d %9.4f %9.4f %9.4f %5d %5d %9.2f %9.2f %9.2f %9.4f\n',[(1:length(total_pix))',x,y,z,c1,c2,total_pix,xpix,ypix,grv]');
            fclose(fid);
        end
        
        
        
        name_ptv=['D:\PTV\Working_folder\Exp19b_1211201\res\ptv_is.',num2Str(n)];
        fid=fopen(name_ptv,'r');
        total_entry=fscanf(fid,'%i',[1 1]);
        data=fscanf(fid,'%i %i %f %f %f ',[5 total_entry]);
        data=data';
        fclose(fid);
        if ~isempty(data)
            
            past=data(:,1);
            future=data(:,2);
            xptv=data(:,3);
            yptv=data(:,4);
            zptv=data(:,5);
            
           
            
            %%%%%%%%%%
            for kk=1:length(xptv)
                ss=inter_v2(xptv(kk)/1000,yptv(kk)/1000,zptv(kk)/1000);
                ss=ss';
                sss(kk,:)=ss;
            end
            
            x_m=sss(:,1);
            y_m=sss(:,2);
            z_m=sss(:,3);
            
            strain=sss(:,4);
            lambda1=sss(:,5);
            lambda2=sss(:,6);
            lambda3=sss(:,7);
            velocity=sss(:,8);
            
            
            
        end
        
        fid = fopen(['D:\PTV\Working_folder\Exp19b_1211201\res\ptv_is_v2.',num2Str(n)],'w');
        fprintf(fid,'%d\n',length(past));
        fprintf(fid,'%4d %4d %9.4f %9.4f %9.4f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n',[past,future,x_m,y_m,z_m,total_pix,xpix,ypix,grv,strain,lambda1,lambda2,lambda3,velocity]');
        fclose(fid);
        
        
    end
    
 
end