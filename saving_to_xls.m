 xbow_result = [x_phi x_lambda]*180/pi;
novatel_result = [n_phi n_lambda]*180/pi;
ref_result = [ins_PVA.INS_Lat(1:plot_len) ins_PVA.INS_Long(1:plot_len)];
result_xls (:,1) = [(x_phi) (n_phi) (ins_PVA.INS_Lat(1:plot_len))'];
result_xls (:,2) = [(x_lambda) (n_lambda) (ins_PVA.INS_Long(1:plot_len))'];