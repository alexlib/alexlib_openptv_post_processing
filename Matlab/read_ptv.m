function [future x_c,y_c,z_c]=read_ptv(name)
fid=fopen(name,'r');
    total_entry=fscanf(fid,'%i',[1 1]);
    data=fscanf(fid,'%i %i %f %f %f ',[5 total_entry]);
    data=data';
    fclose(fid);
    if ~isempty(data)
        f=find(data(:,2)>0);
        if ~isempty(f)
        future=data(:,2);
        x=data(:,3);
        y=data(:,4);
        z=data(:,5);
        future=future(f);
        x_c=x(f);
        y_c=y(f);
        z_c=z(f);
        end
    end