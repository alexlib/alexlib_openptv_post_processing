

h = awaitbar(0,'deba','ptv_is_v2');
steps = 100;
for step = 1:steps
    for i=1:10;
        for j=1:10
            
            % computations take place here
            awaitbar(step / steps,h,sprintf('steps=%d i=%d j=%d ',step,i,j));
        end
    end
end
%delete(h)

