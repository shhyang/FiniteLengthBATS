clear

filename = 'K128M1.mat';

load(filename)

flopt = BATSFLopt(M,K,q,h);

% Soliton Distribution

dd_sol = solitonDist(K);

flopt.setDegreeDist(dd_sol);
stop_sol = flopt.FixedRec_acc(400,'BP');
[err_sol,eco_sol] = flopt.errorProb(stop_sol);

save(filename,'dd_sol','stop_sol','err_sol','eco_sol','-append');

inac_sol = flopt.FixedRec_acc(400,'inac');
ein_sol = flopt.expinac(inac_sol);

save(filename,'inac_sol','ein_sol','-append');

% Robust Solition Distribution

[dd_ros,c_ros,delta_ros] = flopt.robustSolitonOpt(150,100,'BP');
flopt.setDegreeDist(dd_ros);
stop_ros = flopt.FixedRec_acc(400,'BP');
[err_ros,eco_ros] = flopt.errorProb(stop_ros);

save(filename,'dd_ros', 'c_ros', 'delta_ros', 'stop_ros','err_ros','eco_ros','-append');

% Optimize BP decoding

[dd_150, dd_150_diff] = flopt.ddOpt(150,dd_sol,100,'BP');

flopt.setDegreeDist(dd_150);
stop_150 = flopt.FixedRec_acc(400,'BP');
[err_150,eco_150] = flopt.errorProb(stop_150);

save(filename,'dd_150','dd_150_diff','stop_150','err_150','eco_150','-append');

% [dd_160, dd_160_diff] = flopt.ddOpt(160,dd_sol,100,'BP');
% 
% flopt.setDegreeDist(dd_160);
% stop_160 = flopt.FixedRec_acc(400,'BP');
% [err_160,eco_160] = flopt.errorProb(stop_160);
% 
% save(filename,'dd_160','dd_160_diff','stop_160','err_160','eco_160','-append');


% Optimize Inactivation decoding

[dd_130, dd_130_diff] = flopt.ddOpt(130,dd_sol,100,'inac');

flopt.setDegreeDist(dd_130);
inac_130 = flopt.FixedRec_acc(400,'inac');
ein_130 = flopt.expinac(inac_130);

save(filename,'dd_130','dd_130_diff','inac_130','ein_130','-append');

% [dd_135, dd_135_diff] = flopt.ddOpt(135,dd_sol,100,'inac');
% 
% flopt.setDegreeDist(dd_135);
% inac_135 = flopt.FixedRec_acc(400,'inac');
% ein_135 = flopt.expinac(inac_135);
% 
% save(filename,'dd_135','dd_135_diff','inac_135','ein_135','-append');
