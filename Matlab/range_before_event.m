function [xs,zs] = range_before_event (x_of_event,z_of_event,u_of_event,w_of_event,R) 

% A is a rotation  matrix with  an angle of 180 degrees
 %A = [-1,0,0;0,1,0;0,0,-1];
 %instead of A we can use the following equation:
 xs = (x_of_event-u_of_event*R/sqrt((u_of_event)^2+(w_of_event)^2));
 zs = (z_of_event-w_of_event*R/sqrt((u_of_event)^2+(w_of_event)^2));