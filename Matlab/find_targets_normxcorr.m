function [i1,j1,i2,j2] = find_targets_normxcorr(im, ball1, ball2)
%FIND_TARGETS_NORMXCORR (IMAGE, BALL1, BALL2)

verbose = 0; % 0 or 1


% image registration using normalized cross-correlation
% C = conv2(single(db2_1a),single(ball1(end:-1:1,end:-1:1)),'same');
% C = normxcorr2_mex(double(ball1),double(im),'same');
C = normxcorr2(uint8(ball1),uint8(im));

shiftx = ceil(size(ball1,1)/2);
shifty = ceil(size(ball1,2)/2);
C = C(shiftx:end-shiftx,shifty:end-shifty);

% find the peak - the "best" location
% [i1,j1] = find(C == max(C(:)));
[max_cc, imax] = max(abs(C(:)));
[i1, j1] = ind2sub(size(C),imax(1));

if verbose
    imshow(im); hold on;
end

% black out the ball 1 position in order not to get double peak: bug found
% on July 12, 2010 using the db2_2.10766

ia = max(1,i1-shiftx);
ib  = min(i1+shiftx,size(C,1));
ja = max(1,j1-shifty);
jb = min(j1+shifty, size(C,2));


im(ia:ib,ja:jb) = 0;

% option 1
% C = conv2(single(db2_1a),single(ball2(end:-1:1,end:-1:1)),'same');

% option2
% C = normxcorr2(ball2,db2_1a);
% shiftx = ceil(size(ball2,1)/2);
% shifty = ceil(size(ball2,2)/2);
% C = C(shiftx:end-shiftx,shifty:end-shifty);


% option 3
% C = normxcorr2_mex(double(ball2),double(im),'same');

C = normxcorr2(uint8(ball2),uint8(im));

shiftx = ceil(size(ball2,1)/2);
shifty = ceil(size(ball2,2)/2);
C = C(shiftx:end-shiftx,shifty:end-shifty);

% find the peak
[max_cc, imax] = max(abs(C(:)));
[i2, j2] = ind2sub(size(C),imax(1));

if verbose
    scatter(j1,i1,50,'ro','fill');
    scatter(j2,i2,50,'go','fill');
    drawnow
end
