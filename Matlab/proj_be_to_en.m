function vec_be_proj=proj_be_to_en(traj_bin, vec_be,ind_time,ind_be,ind_en,max_jump,proj_frame)

% global traj_bin;

vec_be_proj=[];
order=2;

for i=1:length(ind_time)
    A=zeros(max_jump+1,order+1);
    y=zeros(max_jump+1,order+1);
    count=0;
    for frame=traj_bin(ind_be(ind_time(i)),1):traj_bin(ind_be(ind_time(i)),1)+max_jump
        count=count+1;
        for or=0:order
            A(count,or+1)=frame^or;
        end
        for n=1:3
            y(count,n)=traj_bin(ind_be(ind_time(i))+count-1,3+n);
        end
    end
    X=(A'*A)\A'*y;
    
    %loop through jump size until max jump is reached
    proj=[];
    for n=1:n
        proj(n)=0;
        for or=0:order
            proj(n)=proj(n)+X(or+1,n)*proj_frame^or;
        end
    end
    vec_be_proj=[vec_be_proj;proj];
end