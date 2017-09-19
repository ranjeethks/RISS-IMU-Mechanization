% SYnchrozization: From 490977th second between Xbow data, INSPVAS, Carship
% Reference: 1st values of INS_PVA
% g, Rn and Rm are updated in every loop, 
% bias values from first project, no bias or SF removed for Xbow gyro Wz
% in 1 Hz rate
% nova_wz(te)*cos(p)*cos(r)   0.0205345245148081   357.582036317179 delta phi, delta azi
% nova_wz(te)   0.0208595034299677   357.535182879638, delta phi, delta azi
% bias not removed for odometer       
% done for both xbow , novatel data

clear all;
close all;
clc;
format long g;

gravity	 = 9.805209209982110;


we = 7.2921151467e-5; % rad/sec
a = 6378137; % semi major in m, axis WGS84
e2 = 0.00669437999014; % eccentricity square,  Department of Defense World Geodetic System 1984, Its Definition and Relationships with Local Geodetic Systems http://home.online.no/~sigurdhu/WGS84_Eng.html

a1=9.7803267714; 
a4=-0.0000030876910891;
a2=0.0052790414; 
a5=0.0000000043977311;
a3=0.0000232718; 
a6=0.0000000000007211;

g_p1 = 9.80520920998301;
% Initial g
xbow_fx_bias =  4.6*1e-3*g_p1; %converted mg to g 
xbow_fx_sf = -0.02/100; %converted % to actual SF

xbow_fy_bias =  7.25*1e-3*g_p1; %converted mg to g 
xbow_fy_sf = -0.066/100; %converted % to actual SF

xbow_fz_bias =  -11.65*1e-3*g_p1; %converted mg to g 
xbow_fz_sf = -0.025/100; %converted % to actual SF

xbow_wx_bias = 0.079*pi/180; % No Conversion is needed, already in deg/ sec;
xbow_wx_sf = -0.05/100; %converted % to actual SF

xbow_wy_bias = 0.088*pi/180; % No Conversion is needed, already in deg/ sec;
xbow_wy_sf = 0.17/100; %converted % to actual SF 

xbow_wz_bias = 0.198*pi/180; % No Conversion is needed, already in deg/ sec;
xbow_wz_sf = 0.17/100; %converted % to actual SF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nova_fx_bias =  -0.48*1e-3*g_p1; %converted mg to m/s^2 
nova_fx_sf = 65/1000000; %converted ppm to actual SF


nova_fy_bias =  0.165*1e-3*g_p1; %converted mg to m/s^2 
nova_fy_sf = 45/1000000; %converted ppm to actual SF


nova_wz_bias = -2.3/(180*3600/pi); % converted deg/hr to rad/sec;
nova_wz_sf = -20000/1000000; %converted ppm to actual SF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

novatel_gps_pos = load('BESTGPSPOS.mat');
novatel_gps_vel = load('BESTGPSVEL.mat');
ins_PVA = load('INSPVAS.mat');
raw_nova_imu = load('RAWIMUS.mat');
raw_xbow_imu=load('Xbow_Apr_22_11.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_t=1;
%car_vel_100hz = load('charchip_1Hz_intrap_vel.mat');
%car_time_100hz = load('charchip_1Hz_intrap_time.mat');

car_chip_1hz = load('CarChip_Speed_interpolated.mat');
car_vel_1hz = car_chip_1hz.CarChip_Speed_1HZ;
car_time_1hz = car_chip_1hz.CarChip_second_1HZ;

car_acc_1hz(1) = car_vel_1hz(1);
m=2:length(car_vel_1hz);
n=1:length(car_vel_1hz)-1;
%o=2:length(car_vel_1hz);
acc_in = 2:length(car_vel_1hz);

car_acc_1hz(acc_in)=(car_vel_1hz(m)-car_vel_1hz(n))/delta_t;
%car_at_1hz(acc_in)=(car_time_1hz(m)-car_time_1hz(n))/delta_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nova_fx_raw = raw_nova_imu.f.x;
nova_fy_raw = -(raw_nova_imu.f.y);
nova_wz_raw = raw_nova_imu.w.z;



nova_fx_100Hz = (nova_fx_raw - nova_fx_bias)/(1+nova_fx_sf);  
nova_fy_100Hz = (nova_fy_raw - nova_fy_bias)/(1+nova_fy_sf);  
nova_wz_100Hz = (nova_wz_raw - nova_wz_bias)/(1+nova_wz_sf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xbow_fx_raw = raw_xbow_imu.f.x;
xbow_fy_raw = (raw_xbow_imu.f.y);
xbow_wz_raw = raw_xbow_imu.w.z;



xbow_fx_100Hz = (xbow_fx_raw - xbow_fx_bias)/(1+xbow_fx_sf);  
xbow_fy_100Hz = (xbow_fy_raw - xbow_fy_bias)/(1+xbow_fy_sf);  
xbow_wz_100Hz = (xbow_wz_raw);% - xbow_wz_bias)/(1+xbow_wz_sf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind=1:100:length(nova_fx_100Hz);
nova_fx=nova_fx_100Hz(ind);

ind=1:100:length(nova_fx_100Hz);
nova_fy=nova_fy_100Hz(ind);

ind=1:100:length(nova_fx_100Hz);
nova_wz=nova_wz_100Hz(ind);


% sync_from_xbow=92116;
% ind=sync_from_xbow:100:length(xbow_fx_100Hz);
% nova_fx=xbow_fx_100Hz(ind);
% 
% ind=sync_from_xbow:100:length(xbow_fx_100Hz);
% nova_fy=xbow_fy_100Hz(ind);
% 
% ind=sync_from_xbow:100:length(xbow_fx_100Hz);
% nova_wz=xbow_wz_100Hz(ind);

%ins_PVA.INS_second(1) - raw_xbow_imu.IMU_second(sync_from_xbow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi(1) = ins_PVA.INS_Lat(1)*pi/180; % at time 490977.001
lambda(1) = ins_PVA.INS_Long(1)*pi/180; % at time 490977.001
height(1) = ins_PVA.INS_Alt(1); % at time 490977.001
g = a1*(1+a2*sin(phi(1))*sin(phi(1))+a3*sin(phi(1))*sin(phi(1))*sin(phi(1))*sin(phi(1)))+(a4+a5*sin(phi(1))*sin(phi(1)))*height(1)+a6*height(1)*height(1);  % at time 490977.001
v_e(1)=ins_PVA.INS_ve(1);
v_n(1)=ins_PVA.INS_vn(1);
v_u(1)=ins_PVA.INS_vu(1);
R_N = a./sqrt(1-e2*sin(phi(1)).*sin(phi(1)));% Normal radius
R_M = (a*(1-e2))./((1-e2*sin(phi(1)).*sin(phi(1))).^(1.5));% Meridian radius

p=(asin((nova_fy(1)-car_acc_1hz(1))/g));
pitch(1)=p*180/pi;

r=-(asin((nova_fx(1)+(car_vel_1hz(1)*nova_wz(1)))/(g*cos(p))));
roll(1)=r*180/pi;

Azi(1)=ins_PVA.INS_Azi(1);
az = Azi(1)*pi/180;

plot_len=1800;
 for te=2:plot_len;
     p=(asin((nova_fy(te)-car_acc_1hz(te))/g)); % prev g, because prev phi, h is known
     pitch(te)=p*180/pi;
     
     r=-(asin((nova_fx(te)+(car_vel_1hz(te)*nova_wz(te)))/(g*cos(p)))); % prev g, because prev phi, h is known
     roll(te)=r*180/pi;
     
     az = -(nova_wz(te)*(cos(p)*cos(r)) - we*sin(phi(te-1)) - ((v_e(te-1)*tan(phi(te-1)))/(R_N + height(te-1))))*delta_t + Azi(te-1)*pi/180; % first term present, next terms previous

     Azi(te)=  az*180/pi;
     if (Azi(te) > 360)
         Azi(te)= Azi(te)-360;
     end
     if (Azi(te) < 0)
         Azi(te)= Azi(te)+360;
     end
     %Azi(te)=ins_PVA.INS_Azi(te); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     v_e(te)= car_vel_1hz(te)*sin(az)*cos(p);
     v_n(te)= car_vel_1hz(te)*cos(az)*cos(p);
     v_u(te)= car_vel_1hz(te)*sin(p);
     
     height(te)=height(te-1)+v_u(te)*delta_t;
     phi(te)=phi(te-1)+ v_n(te)*delta_t/(R_M + height(te));
     lambda(te)=lambda(te-1)+ v_e(te)*delta_t/((R_N + height(te))*cos(phi(te)));
     
     g = a1*(1+a2*sin(phi(te))*sin(phi(te))+a3*sin(phi(te))*sin(phi(te))*sin(phi(te))*sin(phi(te)))+(a4+a5*sin(phi(1))*sin(phi(1)))*height(1)+a6*height(1)*height(1);  % at time 490977.001
     R_N = a./sqrt(1-e2*sin(phi(te)).*sin(phi(te)));% Normal radius
     R_M = (a*(1-e2))./((1-e2*sin(phi(te)).*sin(phi(te))).^(1.5));% Meridian radius

     

     
 end
 n_phi=phi;
 n_lambda=lambda;
 n_height=height;
 
 
 
 sync_from_xbow=92116;
ind=sync_from_xbow:100:length(xbow_fx_100Hz);
nova_fx=xbow_fx_100Hz(ind);

ind=sync_from_xbow:100:length(xbow_fx_100Hz);
nova_fy=xbow_fy_100Hz(ind);

ind=sync_from_xbow:100:length(xbow_fx_100Hz);
nova_wz=xbow_wz_100Hz(ind);

ins_PVA.INS_second(1) - raw_xbow_imu.IMU_second(sync_from_xbow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi(1) = ins_PVA.INS_Lat(1)*pi/180; % at time 490977.001
lambda(1) = ins_PVA.INS_Long(1)*pi/180; % at time 490977.001
height(1) = ins_PVA.INS_Alt(1); % at time 490977.001
g = a1*(1+a2*sin(phi(1))*sin(phi(1))+a3*sin(phi(1))*sin(phi(1))*sin(phi(1))*sin(phi(1)))+(a4+a5*sin(phi(1))*sin(phi(1)))*height(1)+a6*height(1)*height(1);  % at time 490977.001
v_e(1)=ins_PVA.INS_ve(1);
v_n(1)=ins_PVA.INS_vn(1);
v_u(1)=ins_PVA.INS_vu(1);
R_N = a./sqrt(1-e2*sin(phi(1)).*sin(phi(1)));% Normal radius
R_M = (a*(1-e2))./((1-e2*sin(phi(1)).*sin(phi(1))).^(1.5));% Meridian radius

p=(asin((nova_fy(1)-car_acc_1hz(1))/g));
pitch(1)=p*180/pi;

r=-(asin((nova_fx(1)+(car_vel_1hz(1)*nova_wz(1)))/(g*cos(p))));
roll(1)=r*180/pi;

Azi(1)=ins_PVA.INS_Azi(1);
az = Azi(1)*pi/180;

plot_len=1800;
 for te=2:plot_len;
     p=(asin((nova_fy(te)-car_acc_1hz(te))/g)); % prev g, because prev phi, h is known
     pitch(te)=p*180/pi;
     
     r=-(asin((nova_fx(te)+(car_vel_1hz(te)*nova_wz(te)))/(g*cos(p)))); % prev g, because prev phi, h is known
     roll(te)=r*180/pi;
     
     az = -(nova_wz(te)*(cos(p)*cos(r)) - we*sin(phi(te-1)) - ((v_e(te-1)*tan(phi(te-1)))/(R_N + height(te-1))))*delta_t + Azi(te-1)*pi/180; % first term present, next terms previous

     Azi(te)=  az*180/pi;
     if (Azi(te) > 360)
         Azi(te)= Azi(te)-360;
     end
     if (Azi(te) < 0)
         Azi(te)= Azi(te)+360;
     end
     %Azi(te)=ins_PVA.INS_Azi(te); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     v_e(te)= car_vel_1hz(te)*sin(az)*cos(p);
     v_n(te)= car_vel_1hz(te)*cos(az)*cos(p);
     v_u(te)= car_vel_1hz(te)*sin(p);
     
     height(te)=height(te-1)+v_u(te)*delta_t;
     phi(te)=phi(te-1)+ v_n(te)*delta_t/(R_M + height(te));
     lambda(te)=lambda(te-1)+ v_e(te)*delta_t/((R_N + height(te))*cos(phi(te)));
     
     g = a1*(1+a2*sin(phi(te))*sin(phi(te))+a3*sin(phi(te))*sin(phi(te))*sin(phi(te))*sin(phi(te)))+(a4+a5*sin(phi(1))*sin(phi(1)))*height(1)+a6*height(1)*height(1);  % at time 490977.001
     R_N = a./sqrt(1-e2*sin(phi(te)).*sin(phi(te)));% Normal radius
     R_M = (a*(1-e2))./((1-e2*sin(phi(te)).*sin(phi(te))).^(1.5));% Meridian radius

     

     
 end
 x_phi=phi;
 x_lambda=lambda;
 x_height=height;
  %plot(pitch); hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'r');
  %plot(roll); hold on; plot(ins_PVA.INS_Roll(1:plot_len),'r');
  %stem(Azi); hold on; plot(ins_PVA.INS_Azi(1:plot_len),'r');
 %plot(phi*180/pi); hold on; plot(ins_PVA.INS_Lat(1:plot_len),'r');
 %plot(lambda*180/pi); hold on; plot(ins_PVA.INS_Long(1:plot_len),'r');
  
 h=figure;
 
 plot(x_lambda*180/pi, x_phi*180/pi, 'r','LineWidth',2); hold on; plot(n_lambda*180/pi, n_phi*180/pi, 'b','LineWidth',2); hold on; plot(ins_PVA.INS_Long(1:plot_len),ins_PVA.INS_Lat(1:plot_len),'g','LineWidth',3); 
 %hold on;  plot3(novatel_gps_pos.GP_Long,novatel_gps_pos.GP_Lat, novatel_gps_pos.GP_Alt,'b--','LineWidth',3);
 grid on;
 lg=legend('Un-aided 3D RISS Mechanization with MEMS IMU','Un-aided 3D RISS Mechanization with Tactical IMU','Reference');
 gt1=findobj(lg,'type','text');
 set(gt1,'fontname','--','fontweight','bold');
    
 xlabel('Longitude (Degrees)','fontweight','bold','fontsize',12);
 ylabel('Latitude (Degrees)','fontweight','bold','fontsize',12);
 %zlabel('Height (meters)','fontweight','bold','fontsize',12);

 dummy=(ins_PVA.INS_Lat(1:plot_len))';
diff(1)=max(abs(((phi(1,1:plot_len))*180/pi)-dummy));
 plot(((phi(1,1:plot_len)*180/pi)-dummy));
 dummy=(ins_PVA.INS_Azi(1:plot_len))';
diff(2)=max(abs((Azi(1,1:plot_len))-dummy));
diff'
 
 

 
%  hh=figure;
%  
%  plot3(lambda*180/pi, phi*180/pi,  'r','LineWidth',2); hold on;  plot(ins_PVA.INS_Long(1:plot_len),ins_PVA.INS_Lat(1:plot_len),'g','LineWidth',3); 
%  %hold on;  plot(novatel_gps_pos.GP_Long,novatel_gps_pos.GP_Lat,'b--','LineWidth',3);
%  grid on;
%  lg=legend('Un-aided 3D RISS Mechanization','Reference','GPS Position');
%  gt1=findobj(lg,'type','text');
%  set(gt1,'fontname','--','fontweight','bold');
%     
%  xlabel('Longitude (Degrees)','fontweight','bold','fontsize',12);
%  ylabel('Latitude (Degrees)','fontweight','bold','fontsize',12);
%  zlabel('Height (meters)','fontweight','bold','fontsize',12);
% 
%  dummy=(ins_PVA.INS_Lat(1:plot_len))';
% max((phi(1,1:plot_len)*180/pi)-dummy)
%  
%  dummy=(ins_PVA.INS_Azi(1:plot_len))';
% max((Azi(1,1:plot_len))-dummy)
%  
 

 
 
 