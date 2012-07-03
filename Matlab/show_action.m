for i = 1:100
    fprintf(1,'\\');
    pause(.05)
    fprintf(1,'\b')
    fprintf(1,'|');
    pause(.05)
    fprintf(1,'\b')
    fprintf(1,'/');
    pause(.05)
    fprintf(1,'\b')
    fprintf(1,'--');
    pause(.05)
    fprintf(1,'\b\b');
    if ~mod(i,10),fprintf(1,'.'),end
end
fprintf('\n');
    