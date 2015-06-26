%%
%load xuap
clear all
stack=repmat(nan,[100 6 5457]);
load data_rsl_p7mm;
load velocity_component.mat  m_x m_y m_z m_u m_v m_w m_vel;
first=1;last=5457;
j=0;
for i=first:last
    i
    j=j+1;
    name=['E:\PTV\Working_folder\WF_3\res\xuap.',num2Str(i)];
    %name=['D:\Utku\xuap.',num2Str(ii)];

    f=load(name);
    s=size(f);
    if s(1,1)>0
        x=f(:,6);
        y=f(:,7);
        z=f(:,8);
        u=f(:,9);
        v=f(:,10);
        w=f(:,11);

        ind_v=find(abs(v)>0);

        if length(ind_v)>0
            xuap_v=v(ind_v);

            xx=x(ind_v);
            yy=y(ind_v);
            zz=z(ind_v);
            clear inter_v;
            for k=1:length(ind_v);
                
                %int_v(k)=interp3(m_x,m_y,m_z,m_v,xx(ind_v(k)),yy(ind_v(k)),zz(ind_v(k)));
                inter_v(k)=interp3(m_x,m_y,m_z,m_v,xx(k),yy(k),zz(k));
                %                 if size(inter_v,1)>length(ind_v)
                %                     break;
                %                 end
            end
            if size(inter_v',1)>size(xuap_v,1)
                disp('Wrong')
                break;
            else
                diff=abs(inter_v'- xuap_v);
                stack(1:length(ind_v),1,j)=xx;
                stack(1:length(ind_v),2,j)=yy;
                stack(1:length(ind_v),3,j)=zz;

                stack(1:length(ind_v),4,j)=xuap_v;
                stack(1:length(ind_v),5,j)=inter_v';
                stack(1:length(ind_v),6,j)=diff;
            end
        end

    end


end
figure;
yyy=stack(:,2,:);
all_diff=stack(:,6,:);
scatter(yyy(:),all_diff(:))
 figure;
 xuapp=stack(:,4,:);
 scatter(yyy(:),xuapp(:))


%%
%go thorugh xuap and take those for which av field has non zero value


%%
%get interp for x,y,z (from xuap) for, say v-comp.

%%
%store difference


%%
%scatter plot y-diff