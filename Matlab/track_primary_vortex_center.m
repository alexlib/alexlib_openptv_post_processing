% track the center of the primary vortex


% filter, get the norm
tmp = vec2scal(filterf(v20,2,'gauss'),'norm');
[xc,yc] = deal(zeros(length(tmp),11));

for i = 1:length(tmp)
    % take the min of the norm map:
    [minw,minI] = min(tmp(i).w(:));
    [k,l] = ind2sub(size(tmp(1).w),minI);

    % visualize if necessary
    contourf(tmp(i).x,tmp(i).y,tmp(i).w');

    xc(i) = tmp(i).x(k);
    yc(i) = tmp(i).y(l);
end
