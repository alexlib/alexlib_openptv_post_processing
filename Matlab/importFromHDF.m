A2 = hdfread('E:\Rotating\Matlab_codes\script_graveyard\test.hdf', '/1', 'Fields', 'ax,ay,az', 'FirstRecord',1 ,'NumRecords',60550);
figure, nhist(A2{1},1000,'r-')
hold on
nhist(A2{2},1000,'g--')
nhist(A2{3},1000,'b-.')


A2 = hdfread('E:\Rotating\Matlab_codes\script_graveyard\test.hdf', '/1', 'Fields', 'x,y,z,u,v,w', 'FirstRecord',1 ,'NumRecords',60550);
figure, nhist(A2{1},1000,'r-')
hold on
nhist(A2{2},1000,'g--')
nhist(A2{3},1000,'b-.')