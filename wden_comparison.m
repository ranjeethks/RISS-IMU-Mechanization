%Ranjeeth KS, University of Calgary

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

xbow_wz_bias =0.03*pi/180; % 0.0198/ previously no bias/0.198/ 0.03
xbow_wz_sf = 0.17/100; %3.17/ previously no SF/0.17/ 0.17


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nova_fx_bias =  -0.48*1e-3*g_p1; %converted mg to m/s^2 
nova_fx_sf = 65/1000000; %converted ppm to actual SF


nova_fy_bias =  0.165*1e-3*g_p1; %converted mg to m/s^2 
nova_fy_sf = 45/1000000; %converted ppm to actual SF


nova_wz_bias = 2.5/(180*3600/pi); % converted deg/hr to rad/sec; 20.3/-2.3
nova_wz_sf = -15000/1000000; %converted ppm to actual SF; -1500/-20000

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bias_odo=mean(car_vel_1hz(1:length(ins_PVA.INS_vn))-(sqrt((ins_PVA.INS_vn).*(ins_PVA.INS_vn)+(ins_PVA.INS_ve).*(ins_PVA.INS_ve))));
car_vel_1hz=car_vel_1hz-(-0.2); %previously no bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


car_acc_1hz(1) = car_vel_1hz(1);
m=2:length(car_vel_1hz);
n=1:length(car_vel_1hz)-1; 
%o=2:length(car_vel_1hz);
acc_in = 2:length(car_vel_1hz);

car_acc_1hz(acc_in)=(car_vel_1hz(m)-car_vel_1hz(n))/delta_t;
%car_at_1hz(acc_in)=(car_time_1hz(m)-car_time_1hz(n))/delta_t;

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
xbow_wz_100Hz = (xbow_wz_raw - xbow_wz_bias)/(1+xbow_wz_sf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shift=00; % samples
 x_shift = round(shift/100); % in seconds
%  car_acc_1hz=[(car_acc_1hz((1+x_shift):length(car_acc_1hz)))  (car_acc_1hz(1:x_shift))];
%  car_vel_1hz=[(car_vel_1hz((1+x_shift):length(car_vel_1hz)))'  (car_vel_1hz(1:x_shift))'];

% nova_fx_100Hz = wden(nova_fx_100Hz, 'rigrsure', 's', 'one', 6, 'db5');
% nova_fy_100Hz = wden(nova_fy_100Hz, 'rigrsure', 's', 'one', 6, 'db5');
% nova_wz_100Hz = wden(nova_wz_100Hz, 'rigrsure', 's', 'one',6, 'db5');

% xbow_fx_100Hz = wden(xbow_fx_100Hz, 'rigrsure', 's', 'one', 6, 'db5');
% xbow_fy_100Hz = wden(xbow_fy_100Hz, 'rigrsure', 's', 'one', 6, 'db5');
% xbow_wz_100Hz = wden(xbow_wz_100Hz, 'rigrsure', 's', 'one',6, 'db5');


sync_from_xbow=1+shift;
ind=sync_from_xbow:100:length(nova_fx_100Hz);
nova_fx=nova_fx_100Hz(ind);

ind=sync_from_xbow:100:length(nova_fx_100Hz);
nova_fy=nova_fy_100Hz(ind);

ind=sync_from_xbow:100:length(nova_fx_100Hz);
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



phi(1) = ins_PVA.INS_Lat(1+x_shift)*pi/180; % at time 490977.001
lambda(1) = ins_PVA.INS_Long(1+x_shift)*pi/180; % at time 490977.001
height(1) = ins_PVA.INS_Alt(1+x_shift); % at time 490977.001
g = a1*(1+a2*sin(phi(1))*sin(phi(1))+a3*sin(phi(1))*sin(phi(1))*sin(phi(1))*sin(phi(1)))+(a4+a5*sin(phi(1))*sin(phi(1)))*height(1)+a6*height(1)*height(1);  % at time 490977.001
v_e(1)=ins_PVA.INS_ve(1+x_shift);
v_n(1)=ins_PVA.INS_vn(1+x_shift);
v_u(1)=ins_PVA.INS_vu(1+x_shift);
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
%      if (Azi(te) > 270)
%          Azi(te)= Azi(te)-360;
%      end
%      if (Azi(te) < 0)
%          Azi(te)= Azi(te)+360;
%      end
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
 n_phi=phi*180/pi;
 n_lambda=lambda*180/pi;
 n_height=height;
 n_v_e = v_e;
  n_v_n = v_n;
   n_v_u = v_u;
   n_pitch = pitch;
   n_roll = roll;
   n_Azi = Azi;


  nova_fx_100Hz = wden(nova_fx_100Hz, 'rigrsure', 's', 'one', 6, 'db5');
nova_fy_100Hz = wden(nova_fy_100Hz, 'rigrsure', 's', 'one', 6, 'db5');
nova_wz_100Hz = wden(nova_wz_100Hz, 'rigrsure', 's', 'one',6, 'db5'); 
   
 sync_from_xbow=1+shift;
ind=sync_from_xbow:100:length(nova_fx_100Hz);
nova_fx=nova_fx_100Hz(ind);

ind=sync_from_xbow:100:length(nova_fx_100Hz);
nova_fy=nova_fy_100Hz(ind);

ind=sync_from_xbow:100:length(nova_fx_100Hz);
nova_wz=nova_wz_100Hz(ind);


 shift=0000; % samples
 x_shift = shift/100; % in seconds
%  car_acc_1hz=[(car_acc_1hz((1+x_shift):length(car_acc_1hz)))  (car_acc_1hz(1:x_shift))];
%  car_vel_1hz=[(car_vel_1hz((1+x_shift):length(car_vel_1hz)))'  (car_vel_1hz(1:x_shift))'];
  %stem(Azi); hold on; plot(ins_PVA.INS_Azi(1:plot_len),'r');

 
%  sync_from_xbow=92116+shift;
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

phi(1) = ins_PVA.INS_Lat(1+ x_shift)*pi/180; % at time 490977.001
lambda(1) = ins_PVA.INS_Long(1+ x_shift)*pi/180; % at time 490977.001
height(1) = ins_PVA.INS_Alt(1+ x_shift); % at time 490977.001
g = a1*(1+a2*sin(phi(1))*sin(phi(1))+a3*sin(phi(1))*sin(phi(1))*sin(phi(1))*sin(phi(1)))+(a4+a5*sin(phi(1))*sin(phi(1)))*height(1)+a6*height(1)*height(1);  % at time 490977.001
v_e(1)=ins_PVA.INS_ve(1+ x_shift);
v_n(1)=ins_PVA.INS_vn(1+ x_shift);
v_u(1)=ins_PVA.INS_vu(1+ x_shift);
R_N = a./sqrt(1-e2*sin(phi(1)).*sin(phi(1)));% Normal radius
R_M = (a*(1-e2))./((1-e2*sin(phi(1)).*sin(phi(1))).^(1.5));% Meridian radius

p=(asin((nova_fy(1)-car_acc_1hz(1))/g));
pitch(1)=p*180/pi;

r=-(asin((nova_fx(1)+(car_vel_1hz(1)*nova_wz(1)))/(g*cos(p))));
roll(1)=r*180/pi;

Azi(1)=ins_PVA.INS_Azi(1)-10;
az = Azi(1)*pi/180;

plot_len=1800;
 for te=2:plot_len;
     p=(asin((nova_fy(te)-car_acc_1hz(te))/g)); % prev g, because prev phi, h is known
     pitch(te)=p*180/pi;
     
     r=-(asin((nova_fx(te)+(car_vel_1hz(te)*nova_wz(te)))/(g*cos(p)))); % prev g, because prev phi, h is known
     roll(te)=r*180/pi;
     
     az = -(nova_wz(te)*(cos(p)*cos(r)) - we*sin(phi(te-1)) - ((v_e(te-1)*tan(phi(te-1)))/(R_N + height(te-1))))*delta_t + Azi(te-1)*pi/180; % first term present, next terms previous

     Azi(te)=  az*180/pi;
%      if (Azi(te) > 270)
%          Azi(te)= Azi(te)-360;
%      end
%      if (Azi(te) < 0)
%          Azi(te)= Azi(te)+360;
%      end
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
%      R_M=6190000;
%      R_N=6190000;

     
 end
 x_phi=phi*180/pi;
 x_lambda=lambda*180/pi;
 x_height=height;
 x_v_e = v_e;
  x_v_n = v_n;
   x_v_u = v_u;
   x_pitch = pitch;
   x_roll = roll;
   x_Azi = Azi;

   x_azi_err=x_Azi-(ins_PVA.INS_Azi(1:plot_len))';
   n_azi_err=n_Azi-(ins_PVA.INS_Azi(1:plot_len))';
   
   for te=1:plot_len;
       if x_azi_err(te)>270
           x_azi_err(te)=x_azi_err(te)-360;
       elseif x_azi_err(te)<-270
           x_azi_err(te)=x_azi_err(te)+360;
       end
       
       if n_azi_err(te)>270
           n_azi_err(te)=n_azi_err(te)-360;
       elseif n_azi_err(te)<-270
           n_azi_err(te)=n_azi_err(te)+360;
       end
       
       if (x_Azi(te) > 360)
          x_Azi(te)= x_Azi(te)-360;
      end
      if (n_Azi(te) > 360)
           n_Azi(te)= n_Azi(te)-360;
      end
       
   end
       
   
%   figure;
%    subplot(3,1,1);
%    %title('Attitude error plot ','fontweight','bold','fontsize',12);
%   plot(x_phi-(ins_PVA.INS_Lat(1:plot_len))','r','LineWidth',2); hold on; plot(n_phi-(ins_PVA.INS_Lat(1:plot_len))','b','LineWidth',2); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   grid on;
%   lg=legend('Latitude error MEMS IMU','Latitude error Tactical IMU');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Latitude error (degrees)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,2);
%   plot(x_lambda-(ins_PVA.INS_Long(1:plot_len))','r','LineWidth',2); hold on; plot(n_lambda-(ins_PVA.INS_Long(1:plot_len))','b','LineWidth',2); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   grid on;
%   lg=legend('Longitude error MEMS IMU','Longitude error Tactical IMU');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Longitude error (degrees)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,3);
%   plot(x_height-(ins_PVA.INS_Alt(1:plot_len))','r','LineWidth',2); hold on; plot(n_height-(ins_PVA.INS_Alt(1:plot_len))','b','LineWidth',2);
%   grid on;
%   lg=legend('Height error MEMS IMU','Height error Tactical IMU');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Height error (m)','fontweight','bold','fontsize',10);
 
  
%   figure;
%    subplot(3,1,1);
%    %title('Attitude error plot ','fontweight','bold','fontsize',12);
%   plot(x_v_e-(ins_PVA.INS_ve(1:plot_len))','r','LineWidth',2); hold on; plot(n_v_e-(ins_PVA.INS_ve(1:plot_len))','b','LineWidth',2); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   grid on;
%   lg=legend('Velocity East error MEMS IMU','Velocity East error Tactical IMU');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Velocity East error (m/s)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,2);
%   plot(x_v_n-(ins_PVA.INS_vn(1:plot_len))','r','LineWidth',2); hold on; plot(n_v_n-(ins_PVA.INS_vn(1:plot_len))','b','LineWidth',2); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   grid on;
%   lg=legend('Velocity North error MEMS IMU','Velocity North error Tactical IMU');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Velocity North error (m/s)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,3);
%   plot(x_v_u-(ins_PVA.INS_vu(1:plot_len))','r','LineWidth',2); hold on; plot(n_v_u-(ins_PVA.INS_vu(1:plot_len))','b','LineWidth',2);
%   grid on;
%   lg=legend('Velocity Up error MEMS IMU','Velocity Up error Tactical IMU');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Velocity Up error (m/s)','fontweight','bold','fontsize',10);
  
% figure;
%    subplot(3,1,1);
   plot(n_pitch-(ins_PVA.INS_Pitch(1:plot_len))','b','LineWidth',2); hold on; plot(x_pitch-(ins_PVA.INS_Pitch(1:plot_len))','r','LineWidth',2); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
  grid on;
  lg=legend('Pitch error before de-noising','Pitch error after de-noising');
  gt1=findobj(lg,'type','text');
  set(gt1,'fontname','--','fontweight','bold');
  
  xlabel('time (seconds)','fontweight','bold','fontsize',10);
  ylabel('Pitch error (degrees)','fontweight','bold','fontsize',10);
     title('Pitch error plot before and after de-noising ','fontweight','bold','fontsize',12);

%   
%   subplot(3,1,2);
%   plot(x_roll-(ins_PVA.INS_Roll(1:plot_len))','r','LineWidth',2); hold on; plot(n_roll-(ins_PVA.INS_Roll(1:plot_len))','b','LineWidth',2); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   grid on;
%   lg=legend('Roll error MEMS IMU','Roll error Tactical IMU');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Roll error (degrees)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,3);
%   plot(x_azi_err,'r','LineWidth',2); hold on; plot(n_azi_err,'b','LineWidth',2); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   grid on;
%   lg=legend('Azimuth error MEMS IMU','Azimuth error Tactical IMU');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Azimuth error (degrees)','fontweight','bold','fontsize',10);  



% figure;
%    subplot(3,1,1);
%   plot(x_pitch,'r','LineWidth',2); hold on; plot(n_pitch,'b','LineWidth',2); hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g','LineWidth',2);
%   grid on;
%   lg=legend('Pitch MEMS IMU','Pitch  Tactical IMU','Pitch Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Pitch  (degrees)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,2);
%   plot(x_roll,'r','LineWidth',2); hold on; plot(n_roll,'b','LineWidth',2); hold on; plot(ins_PVA.INS_Roll(1:plot_len),'g','LineWidth',2);
%   grid on;
%   lg=legend('Roll MEMS IMU','Roll  Tactical IMU','Roll Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Roll  (degrees)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,3);
%   plot(x_Azi,'r','LineWidth',2); hold on; plot(n_Azi,'b','LineWidth',2); hold on; plot(ins_PVA.INS_Azi(1:plot_len),'g','LineWidth',2);
%   grid on;
%   lg=legend('Azimuth MEMS IMU','Azimuth  Tactical IMU','Azimuth Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Azimuth  (degrees)','fontweight','bold','fontsize',10);


% figure;
%    subplot(3,1,1);
%   plot(x_v_e,'r','LineWidth',2); hold on; plot(n_v_e,'b','LineWidth',2); hold on; plot(ins_PVA.INS_ve(1:plot_len),'g','LineWidth',2);
%   grid on;
%   lg=legend('East Velocity MEMS IMU','East Velocity  Tactical IMU','East Velocity Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('East Velocity (m/s)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,2);
%   plot(x_v_n,'r','LineWidth',2); hold on; plot(n_v_n,'b','LineWidth',2); hold on; plot(ins_PVA.INS_vn(1:plot_len),'g','LineWidth',2);
%   grid on;
%   lg=legend('North Velocity MEMS IMU','North Velocity  Tactical IMU','North Velocity Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('North Velocity  (m/s)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,3);
%   plot(x_v_u,'r','LineWidth',2); hold on; plot(n_v_u,'b','LineWidth',2); hold on; plot(ins_PVA.INS_vu(1:plot_len),'g','LineWidth',2);
%   grid on;
%   lg=legend('Up Velocity MEMS IMU','Up Velocity  Tactical IMU','Up Velocity Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Up Velocity  (m/s)','fontweight','bold','fontsize',10);
%   



% figure;
%    subplot(3,1,1);
%   plot(x_phi,'r','LineWidth',2); hold on; plot(n_phi,'b','LineWidth',2); hold on; plot(ins_PVA.INS_Lat(1:plot_len),'g','LineWidth',2);
%   grid on;
%   lg=legend('Lattitude MEMS IMU','Lattitude  Tactical IMU','Lattitude Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Lattitude (degrees)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,2);
%   plot(x_lambda,'r','LineWidth',2); hold on; plot(n_lambda,'b','LineWidth',2); hold on; plot(ins_PVA.INS_Long(1:plot_len),'g','LineWidth',2);
%   grid on;
%   lg=legend('Longitude MEMS IMU','Longitude  Tactical IMU','Longitude Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Longitude  (degrees)','fontweight','bold','fontsize',10);
%   
%   subplot(3,1,3);
%   plot(x_height,'r','LineWidth',2); hold on; plot(n_height,'b','LineWidth',2); hold on; plot(ins_PVA.INS_Alt(1:plot_len),'g','LineWidth',2);
%   grid on;
%   lg=legend('Altitude MEMS IMU','Altitude  Tactical IMU','Up Velocity Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
%   
%   xlabel('time (seconds)','fontweight','bold','fontsize',10);
%   ylabel('Altitude  (m)','fontweight','bold','fontsize',10);
  
  
  
%   figure;
%   plot(x_v_e-(ins_PVA.INS_ve(1:plot_len))','r'); hold on; plot(n_v_e-(ins_PVA.INS_ve(1:plot_len))','b'); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   figure;  
%   plot(x_v_n-(ins_PVA.INS_vn(1:plot_len))','r'); hold on; plot(n_v_n-(ins_PVA.INS_vn(1:plot_len))','b'); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   figure;  
%   plot(x_v_u-(ins_PVA.INS_vu(1:plot_len))','r'); hold on; plot(n_v_u-(ins_PVA.INS_vu(1:plot_len))','b');
%   
%   
%   figure;
%   plot(x_phi-(ins_PVA.INS_Lat(1:plot_len))','r'); hold on; plot(n_phi-(ins_PVA.INS_Lat(1:plot_len))','b'); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   figure;  
%   plot(x_lambda-(ins_PVA.INS_Long(1:plot_len))','r'); hold on; plot(n_lambda-(ins_PVA.INS_Long(1:plot_len))','b'); %hold on; plot(ins_PVA.INS_Pitch(1:plot_len),'g');
%   figure;  
%   plot(x_height-(ins_PVA.INS_Alt(1:plot_len))','r'); hold on; plot(n_height-(ins_PVA.INS_Alt(1:plot_len))','b');
%   
% %   plot(roll); hold on; plot(ins_PVA.INS_Roll(1:plot_len),'r');
% %   plot(Azi); hold on; plot(ins_PVA.INS_Azi(1:plot_len),'r');
%  %plot(phi*180/pi); hold on; plot(ins_PVA.INS_Lat(1:plot_len),'r');
%  %plot(lambda*180/pi); hold on; plot(ins_PVA.INS_Long(1:plot_len),'r');
%   
%  
%  hh=figure;
%  plot(x_lambda, x_phi, 'r','LineWidth',2); hold on; plot(n_lambda, n_phi, 'b','LineWidth',2); hold on; plot(ins_PVA.INS_Long(1:plot_len),ins_PVA.INS_Lat(1:plot_len),'g','LineWidth',3); 
%  %hold on;  plot3(novatel_gps_pos.GP_Long,novatel_gps_pos.GP_Lat, novatel_gps_pos.GP_Alt,'b--','LineWidth',3);
%   grid on;
%   lg=legend('Un-aided 3D-RISS Mechanization with MEMS IMU','Un-aided 3D RISS Mechanization with Tactical IMU','Reference');
%   gt1=findobj(lg,'type','text');
%   set(gt1,'fontname','--','fontweight','bold');
% %     
%   xlabel('Longitude (Degrees)','fontweight','bold','fontsize',12);
%   ylabel('Latitude (Degrees)','fontweight','bold','fontsize',12);
% %  %zlabel('Height (meters)','fontweight','bold','fontsize',12);
%   title('Position plot: 3D RISS Mechanizations ','fontweight','bold','fontsize',12);
% %print(hh,'-djpeg','-r500','pos_plot_3d_riss');
%  dummy=(ins_PVA.INS_Lat(1:plot_len))';
% diff(1)=max(abs(((phi(1,1:plot_len))*180/pi)-dummy));
%  %figure(3);
%  %plot(((phi(1,1:plot_len)*180/pi)-dummy));
%  dummy=(ins_PVA.INS_Azi(1:plot_len))';
% diff(2)=max(abs((Azi(1,1:plot_len))-dummy));
%  dummy=(ins_PVA.INS_Long(1:plot_len))';
% diff(3)=max(abs((lambda(1,1:plot_len))*180/pi-dummy));
% diff'
 
 

 
%  hh=figure;
%  
%  plot3(x_lambda*180/pi, x_phi*180/pi, x_height,  'r','LineWidth',2); hold on; plot3(n_lambda*180/pi, n_phi*180/pi, n_height,  'b','LineWidth',2); hold on; plot3(ins_PVA.INS_Long(1:plot_len),ins_PVA.INS_Lat(1:plot_len),ins_PVA.INS_Alt(1:plot_len),'g','LineWidth',3); 
%  %hold on;  plot(novatel_gps_pos.GP_Long,novatel_gps_pos.GP_Lat,'b--','LineWidth',3);
%  grid on;
%  %lg=legend('Un-aided 3D RISS Mechanization','Reference','GPS Position');
%  %gt1=findobj(lg,'type','text');
%  %set(gt1,'fontname','--','fontweight','bold');
%     
%  xlabel('Longitude (Degrees)','fontweight','bold','fontsize',12);
%  ylabel('Latitude (Degrees)','fontweight','bold','fontsize',12);
%  zlabel('Height (meters)','fontweight','bold','fontsize',12);
% title('3D plot','fontweight','bold','fontsize',12);
% %print(hh,'-djpeg','-r500','3D_plot_3d_riss');
%  dummy=(ins_PVA.INS_Lat(1:plot_len))';
% max((phi(1,1:plot_len)*180/pi)-dummy)
%  
%  dummy=(ins_PVA.INS_Azi(1:plot_len))';
% max((Azi(1,1:plot_len))-dummy)
%  
 

 
 
 