figure;
   subplot(3,1,1);
  plot(x_phi,'r','LineWidth',2); hold on; plot(n_phi,'b','LineWidth',2); hold on; plot(ins_PVA.INS_Lat(1:plot_len),'g','LineWidth',2);
  grid on;
  lg=legend('Lattitude MEMS IMU','Lattitude  Tactical IMU','Lattitude Reference');
  gt1=findobj(lg,'type','text');
  set(gt1,'fontname','--','fontweight','bold');
  
  xlabel('time (seconds)','fontweight','bold','fontsize',10);
  ylabel('Lattitude (degrees)','fontweight','bold','fontsize',10);
  
  subplot(3,1,2);
  plot(x_lambda,'r','LineWidth',2); hold on; plot(n_lambda,'b','LineWidth',2); hold on; plot(ins_PVA.INS_Long(1:plot_len),'g','LineWidth',2);
  grid on;
  lg=legend('Longitude MEMS IMU','Longitude  Tactical IMU','Longitude Reference');
  gt1=findobj(lg,'type','text');
  set(gt1,'fontname','--','fontweight','bold');
  
  xlabel('time (seconds)','fontweight','bold','fontsize',10);
  ylabel('Longitude  (degrees)','fontweight','bold','fontsize',10);
  
  subplot(3,1,3);
  plot(x_height,'r','LineWidth',2); hold on; plot(n_height,'b','LineWidth',2); hold on; plot(ins_PVA.INS_Alt(1:plot_len),'g','LineWidth',2);
  grid on;
  lg=legend('Up Velocity MEMS IMU','Up Velocity  Tactical IMU','Up Velocity Reference');
  gt1=findobj(lg,'type','text');
  set(gt1,'fontname','--','fontweight','bold');
  
  xlabel('time (seconds)','fontweight','bold','fontsize',10);
  ylabel('Height  (m)','fontweight','bold','fontsize',10);