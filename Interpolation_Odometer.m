close all;
clear all;
clc;

odometer=load('D:\Study\ENGO 623\Project 3 and 4\2011-04-22 Kingston East Long\2011-04-22 Kingston East Long\ProPakRover\CarChip\CarChip_Speed_interpolated.mat');
 odo_Vel_1Hz = odometer.CarChip_Speed_1HZ;
x = odometer.CarChip_second_1HZ; 
y = odo_Vel_1Hz; 
xi = x(1):1/100:x(length(x)) ;
yi = interp1(x,y,xi); 
plot(x,y,'o',xi,yi);