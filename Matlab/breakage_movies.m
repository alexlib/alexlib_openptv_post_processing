figure(1)
first=19720;
last=19740;
img=2;
% first=21021;
% last=21050;

%name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\WF/'];
name_root=['D:\technical\codes\IfU-IPG-PTV\working_folder_colloid_dup/'];

dummy=0;
for i=first:last
    dummy=dummy+1;
    name=[name_root,'img\cam',num2Str(img),'.',num2Str(i)];
    %name=['D:\technical\codes\IfU-IPG-PTV\working_folder_colloid_dup\img\cam2.',num2Str(i)];
    f1=imread(name);
    f1=cat(3,f1,f1,f1);
    
    f2=[name_root,'rt_is_2d/rt_is2d_',num2Str(i),'.mat'];%+first-1
    load(f2,'a1_a','a2_a','a3_a','a4_a','b1_a','b2_a','b3_a','b4_a', ...
        'c1_a','c2_a','c3_a','c4_a','d1_a','d2_a','d3_a','d4_a', ...
        'e1_a','e2_a','e3_a','e4_a', ...
        'f1_a','f2_a','f3_a','f4_a', ...
        'g1_a','g2_a','g3_a','g4_a');
    f3=[name_root,'res/xuap.',num2Str(i)];
    xuap=load(f3);
    u=xuap(:,9);
    v=xuap(:,10);
    w=xuap(:,11);
    vel=(u.^2+v.^2+w.^2).^0.5;
    for j=1:length(f2)
        switch img
            case 1
                b_a=b1_a(j); c_a=c1_a(j); g_a=g1_a(j); f_a=f1_a(j);
            case 2
                b_a=b2_a(j); c_a=c2_a(j); g_a=g2_a(j); f_a=f2_a(j);
            case 3
                b_a=b3_a(j); c_a=c3_a(j); g_a=g3_a(j); f_a=f3_a(j);
            case 4
                b_a=b4_a(j); c_a=c4_a(j); g_a=g4_a(j); f_a=f4_a(j);
        end
              
        if isnan(g_a)==0 & isnan(f_a)==0
            dx=round(c_a);
            dy=round(b_a);
            x_im=round(g_a);
            y_im=round(f_a);
            xc=0;
            for xx=-0.5*dx-0:0.5*dx+0
                xc=xc+1;
                yc=0;
                for yy=-0.5*dy-0:0.5*dy+0
                    yc=yc+1;
                    if vel(j)>0
                        if abs(xx-0.5*dx)<2 | abs(xx+0.5*dx)<2 | abs(yy-0.5*dy)<2 | abs(yy+0.5*dy)<2 %| abs(xx-0)<2 | abs(yy-0)<2
                            f1(round(x_im+xx),round(y_im+yy),1)=0;
                            f1(round(x_im+xx),round(y_im+yy),2)=255;
                            f1(round(x_im+xx),round(y_im+yy),3)=0;
                        end
                    else
                        if abs(xx-0.5*dx)<2 | abs(xx+0.5*dx)<2 | abs(yy-0.5*dy)<2 | abs(yy+0.5*dy)<2 %| abs(xx-0)<2 | abs(yy-0)<2
                            f1(round(x_im+xx),round(y_im+yy),1)=255;
                            f1(round(x_im+xx),round(y_im+yy),2)=255;
                            f1(round(x_im+xx),round(y_im+yy),3)=255;
                        end
                    end
                end
            end
        end
    end
    
    imshow(f1);
    tit=['frame nuber: ',num2Str(i)];
    title(tit);
    A(dummy)=getframe;
end