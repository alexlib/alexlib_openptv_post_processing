%Read in test images and display them Figure 1
cir = imread('circles.tif');
ic = imread('ic.tif');
icbw = im2bw(ic);
icbwa = im2bw(double(ic)/255, 0.5);

figure('Name', 'Test cases:');
subplot(2,2,1), imshow(cir), title('Circles');
subplot(2,2,2), imshow(ic), title('IC - original');
subplot(2,2,3), imshow(icbw), title('ICBW -  optimal threshold');
subplot(2,2,4), imshow(icbwa), title('ICBW - threshold 0.5');

%BWMORPH bothat Figure 2
cirbothat1 = bwmorph(cir, 'bothat', 1);
cirbothat2 = bwmorph(cir, 'bothat', 5);
cirbothatinf = bwmorph(cir, 'bothat', Inf);
icbwbothat1 = bwmorph(icbw, 'bothat', 1);
icbwbothat2 = bwmorph(icbw, 'bothat', 5);
icbwbothatinf = bwmorph(icbw, 'bothat', Inf);
figure('Name', 'Testing the bothat operation');
subplot(2,3,1), imshow(cirbothat1), title('Circle, N=1');
subplot(2,3,2), imshow(cirbothat2), title('Circle, N=5');
subplot(2,3,3), imshow(cirbothatinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwbothat1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwbothat2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwbothatinf), title('ICBW, N=Inf');

%BWMORPH bridge Figure 3
cirbridge1 = bwmorph(cir, 'bridge', 1);
cirbridge2 = bwmorph(cir, 'bridge', 5);
cirbridgeinf = bwmorph(cir, 'bridge', Inf);
icbwbridge1 = bwmorph(icbw, 'bridge', 1);
icbwbridge2 = bwmorph(icbw, 'bridge', 5);
icbwbridgeinf = bwmorph(icbw, 'bridge', Inf);
figure('Name', 'Testing the bridge operation');
subplot(2,3,1), imshow(cirbridge1), title('Circle, N=1');
subplot(2,3,2), imshow(cirbridge2), title('Circle, N=5');
subplot(2,3,3), imshow(cirbridgeinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwbridge1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwbridge2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwbridgeinf), title('ICBW, N=Inf');

%BWMORPH clean Figure 4
circlean1 = bwmorph(cir, 'clean', 1);
circlean2 = bwmorph(cir, 'clean', 5);
circleaninf = bwmorph(cir, 'clean', Inf);
icbwclean1 = bwmorph(icbw, 'clean', 1);
icbwclean2 = bwmorph(icbw, 'clean', 5);
icbwcleaninf = bwmorph(icbw, 'clean', Inf);
figure('Name', 'Testing the clean operation');
subplot(2,3,1), imshow(circlean1), title('Circle, N=1');
subplot(2,3,2), imshow(circlean2), title('Circle, N=5');
subplot(2,3,3), imshow(circleaninf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwclean1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwclean2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwcleaninf), title('ICBW, N=Inf');

%BWMORPH close Figure 5
circlose1 = bwmorph(cir, 'close', 1);
circlose2 = bwmorph(cir, 'close', 5);
circloseinf = bwmorph(cir, 'close', Inf);
icbwclose1 = bwmorph(icbw, 'close', 1);
icbwclose2 = bwmorph(icbw, 'close', 5);
icbwcloseinf = bwmorph(icbw, 'close', Inf);
figure('Name', 'Testing the close operation');
subplot(2,3,1), imshow(circlose1), title('Circle, N=1');
subplot(2,3,2), imshow(circlose2), title('Circle, N=5');
subplot(2,3,3), imshow(circloseinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwclose1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwclose2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwcloseinf), title('ICBW, N=Inf');

%BWMORPH diag Figure 6
cirdiag1 = bwmorph(cir, 'diag', 1);
cirdiag2 = bwmorph(cir, 'diag', 5);
cirdiaginf = bwmorph(cir, 'diag', Inf);
icbwdiag1 = bwmorph(icbw, 'diag', 1);
icbwdiag2 = bwmorph(icbw, 'diag', 5);
icbwdiaginf = bwmorph(icbw, 'diag', Inf);
figure('Name', 'Testing the diag operation');
subplot(2,3,1), imshow(cirdiag1), title('Circle, N=1');
subplot(2,3,2), imshow(cirdiag2), title('Circle, N=5');
subplot(2,3,3), imshow(cirdiaginf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwdiag1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwdiag2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwdiaginf), title('ICBW, N=Inf');

%BWMORPH dilate Figure 7
cirdilate1 = bwmorph(cir, 'dilate', 1);
cirdilate2 = bwmorph(cir, 'dilate', 5);
cirdilateinf = bwmorph(cir, 'dilate', Inf);
icbwdilate1 = bwmorph(icbw, 'dilate', 1);
icbwdilate2 = bwmorph(icbw, 'dilate', 5);
icbwdilateinf = bwmorph(icbw, 'dilate', Inf);
figure('Name', 'Testing the dilate operation');
subplot(2,3,1), imshow(cirdilate1), title('Circle, N=1');
subplot(2,3,2), imshow(cirdilate2), title('Circle, N=5');
subplot(2,3,3), imshow(cirdilateinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwdilate1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwdilate2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwdilateinf), title('ICBW, N=Inf');

%BWMORPH erode Figure 8
cirerode1 = bwmorph(cir, 'erode', 1);
cirerode2 = bwmorph(cir, 'erode', 5);
cirerodeinf = bwmorph(cir, 'erode', Inf);
icbwerode1 = bwmorph(icbw, 'erode', 1);
icbwerode2 = bwmorph(icbw, 'erode', 5);
icbwerodeinf = bwmorph(icbw, 'erode', Inf);
figure('Name', 'Testing the erode operation');
subplot(2,3,1), imshow(cirerode1), title('Circle, N=1');
subplot(2,3,2), imshow(cirerode2), title('Circle, N=5');
subplot(2,3,3), imshow(cirerodeinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwerode1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwerode2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwerodeinf), title('ICBW, N=Inf');

%BWMORPH fill Figure 9
cirfill1 = bwmorph(cir, 'fill', 1);
cirfill2 = bwmorph(cir, 'fill', 5);
cirfillinf = bwmorph(cir, 'fill', Inf);
icbwfill1 = bwmorph(icbw, 'fill', 1);
icbwfill2 = bwmorph(icbw, 'fill', 5);
icbwfillinf = bwmorph(icbw, 'fill', Inf);
figure('Name', 'Testing the fill operation');
subplot(2,3,1), imshow(cirfill1), title('Circle, N=1');
subplot(2,3,2), imshow(cirfill2), title('Circle, N=5');
subplot(2,3,3), imshow(cirfillinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwfill1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwfill2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwfillinf), title('ICBW, N=Inf');

%BWMORPH hbreak Figure 10
cirhbreak1 = bwmorph(cir, 'hbreak', 1);
cirhbreak2 = bwmorph(cir, 'hbreak', 5);
cirhbreakinf = bwmorph(cir, 'hbreak', Inf);
icbwhbreak1 = bwmorph(icbw, 'hbreak', 1);
icbwhbreak2 = bwmorph(icbw, 'hbreak', 5);
icbwhbreakinf = bwmorph(icbw, 'hbreak', Inf);
figure('Name', 'Testing the hbreak operation');
subplot(2,3,1), imshow(cirhbreak1), title('Circle, N=1');
subplot(2,3,2), imshow(cirhbreak2), title('Circle, N=5');
subplot(2,3,3), imshow(cirhbreakinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwhbreak1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwhbreak2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwhbreakinf), title('ICBW, N=Inf');

%BWMORPH majority Figure 11
cirmajority1 = bwmorph(cir, 'majority', 1);
cirmajority2 = bwmorph(cir, 'majority', 5);
cirmajorityinf = bwmorph(cir, 'majority', Inf);
icbwmajority1 = bwmorph(icbw, 'majority', 1);
icbwmajority2 = bwmorph(icbw, 'majority', 5);
icbwmajorityinf = bwmorph(icbw, 'majority', Inf);
figure('Name', 'Testing the majority operation');
subplot(2,3,1), imshow(cirmajority1), title('Circle, N=1');
subplot(2,3,2), imshow(cirmajority2), title('Circle, N=5');
subplot(2,3,3), imshow(cirmajorityinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwmajority1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwmajority2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwmajorityinf), title('ICBW, N=Inf');

%BWMORPH open Figure 12
ciropen1 = bwmorph(cir, 'open', 1);
ciropen2 = bwmorph(cir, 'open', 5);
ciropeninf = bwmorph(cir, 'open', Inf);
icbwopen1 = bwmorph(icbw, 'open', 1);
icbwopen2 = bwmorph(icbw, 'open', 5);
icbwopeninf = bwmorph(icbw, 'open', Inf);
figure('Name', 'Testing the open operation');
subplot(2,3,1), imshow(ciropen1), title('Circle, N=1');
subplot(2,3,2), imshow(ciropen2), title('Circle, N=5');
subplot(2,3,3), imshow(ciropeninf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwopen1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwopen2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwopeninf), title('ICBW, N=Inf');

%BWMORPH remove Figure 13
cirremove1 = bwmorph(cir, 'remove', 1);
cirremove2 = bwmorph(cir, 'remove', 5);
cirremoveinf = bwmorph(cir, 'remove', Inf);
icbwremove1 = bwmorph(icbw, 'remove', 1);
icbwremove2 = bwmorph(icbw, 'remove', 5);
icbwremoveinf = bwmorph(icbw, 'remove', Inf);
figure('Name', 'Testing the remove operation');
subplot(2,3,1), imshow(cirremove1), title('Circle, N=1');
subplot(2,3,2), imshow(cirremove2), title('Circle, N=5');
subplot(2,3,3), imshow(cirremoveinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwremove1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwremove2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwremoveinf), title('ICBW, N=Inf');

%BWMORPH shrink Figure 14
cirshrink1 = bwmorph(cir, 'shrink', 1);
cirshrink2 = bwmorph(cir, 'shrink', 5);
cirshrinkinf = bwmorph(cir, 'shrink', Inf);
icbwshrink1 = bwmorph(icbw, 'shrink', 1);
icbwshrink2 = bwmorph(icbw, 'shrink', 5);
icbwshrinkinf = bwmorph(icbw, 'shrink', Inf);
figure('Name', 'Testing the shrink operation');
subplot(2,3,1), imshow(cirshrink1), title('Circle, N=1');
subplot(2,3,2), imshow(cirshrink2), title('Circle, N=5');
subplot(2,3,3), imshow(cirshrinkinf), title('Circle, N=Inf');
subplot(2,3,4), imshow(icbwshrink1), title('ICBW, N=1');
subplot(2,3,5), imshow(icbwshrink2), title('ICBW, N=5');
subplot(2,3,6), imshow(icbwshrinkinf), title('ICBW, N=Inf');