function x = smooth2d(x)
% x = filter2(fspecial('average',2),x);
% x = filter2(fspecial('unsharp',.7),x);
x = filter2(fspecial('gaussian',[3 3], .7),x);
