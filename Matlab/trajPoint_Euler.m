function res=trajPoint_Euler(first,last,tail);


value=zeros(200,30,2000);
pointInd=zeros(200,1);

for i=first:last
    clear f;
    i
    name=['D:/data/Riso_micro_PTV/trajPoint.',num2Str(i)];
    f=load(name);
    s=size(f);
    if s(1,1)>0
        u=f(:,4);
        v=f(:,5);
        w=f(:,6);
        absu=(u.^2+v.^2+w.^2).^0.5;
        ax=f(:,7);%%%%%%%%%%%
        ay=f(:,8);%%%%%%%%%%%
        az=f(:,9);%%%%%%%%%%%
        absa=(ax.^2+ay.^2+az.^2).^0.5;%%%%%%%%%%%%%
        curv=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3);
        
        
        
        reldiv=f(:,31);
        age=f(:,32);

        ende=0;
        for iii=1:s(1,1)-1
            start=ende+1;
            ende=start;
            leng=1;
            searching=1;
            while ende<s(1,1)-1 & searching==1
                if age(ende+1)==age(ende)+1
                    ende=ende+1;
                    leng=leng+1;
                else
                    searching=0;
                end
            end
            if ende<s(1,1)
                for ii=start:ende
                    if age(ii)>tail & age(ii)<leng-tail+1 & reldiv(ii)<0.2 & absa(ii)^2/8e-4<100 & curv(ii)<3000
                        pointInd(age(ii))=pointInd(age(ii))+1;
                        for attr=1:30
                            value(age(ii),attr,pointInd(age(ii)))=f(ii,attr);
                        end
                        
                    end
                end
            end
        end    
    end
    %write out
    name=['D:/data/Riso_micro_PTV/Euler.',num2Str(i)];
    fid=fopen(name,'w');
    pointInd(1)
    for ii=1:pointInd(1)
       for attr=1:29
           dummy=value(1,attr,ii);
           %[attr dummy]
           fprintf(fid, '%s\t', num2Str(value(1,attr,ii)));
       end
       fprintf(fid, '%s\n', num2Str(value(1,30,ii)));
    end
    fclose(fid);
    %shift upwards
    value(1:199,:,:)=value(2:200,:,:);
    pointInd(1:199,1)=pointInd(2:200,1);
    pointInd(200)=0;
end


