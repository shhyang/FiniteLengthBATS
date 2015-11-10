addpath expmv

clear

filename = 'K256M16.mat';

load(filename)

asy = BATSAsymp(M,K,q,h);
flopt = BATSFLopt(M,K,q,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Asymptotic optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Run Asymptotic Optimization!');
dd_asy = asy.opt(0.01);

flopt.setDegreeDist(dd_asy);
stop_asy = flopt.FixedRec_acc(400,'BP');
[err_asy,eco_asy] = flopt.errorProb(stop_asy);

save(filename,'dd_asy','stop_asy','err_asy','eco_asy','-append');

inac_asy = flopt.FixedRec_acc(400,'inac');
ein_asy = flopt.expInac(inac_asy);

save(filename,'inac_asy','ein_asy','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Run Finite-Length BP Optimization!');
[dd_40, dd_40_diff] = flopt.ddOpt(40,dd_asy,100,'BP');

flopt.setDegreeDist(dd_40);
stop_40 = flopt.FixedRec_acc(400,'BP');
[err_40,eco_40] = flopt.errorProb(stop_40);

save(filename,'dd_40','dd_40_diff','stop_40','err_40','eco_40','-append');

% Poisson number of batches

flopt.setDegreeDist(dd_40);
stop_40_p = flopt.PoissonRec_acc(400,0.5);
[err_40_p,eco_40_p] = flopt.errorProb(stop_40_p);
save(filename,'stop_40_p','err_40_p','eco_40_p','-append');

% [dd_45, dd_45_diff] = flopt.ddOpt(45,dd_asy,100,'BP');
% 
% flopt.setDegreeDist(dd_45);
% stop_45 = flopt.FixedRec_acc(400);
% [err_45,eco_45] = flopt.errorProb(stop_45);
% 
% save(filename,'dd_45','dd_45_diff','stop_45','err_45','eco_45','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inactivation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Run Finite-Length Inactivation Optimization!');

[dd_25, dd_25_diff] = flopt.ddOpt(25,dd_asy,100,'inac');

flopt.setDegreeDist(dd_25);
inac_25 = flopt.FixedRec_acc(400,'inac');
ein_25 = flopt.expInac(inac_25);

save(filename,'dd_25','dd_25_diff','inac_25','ein_25','-append');

% Poisson number of batches

flopt.setDegreeDist(dd_25);
inac_25_p = flopt.PoissonRec_acc(400,0.5,'inac');
ein_25_p = flopt.expInac(inac_25_p);

save(filename,'inac_25_p','ein_25_p','-append');

% [dd_30, dd_30_diff] = flopt.ddOpt(30,dd_asy,100,'inac');
% 
% flopt.setDegreeDist(dd_30);
% inac_30 = flopt.FixedRec_acc(400,'inac');
% ein_30 = flopt.expInac(inac_30);
% 
% save(filename,'dd_30','dd_30_diff','inac_30','ein_30','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximize error exponent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dd_mee,~] = flopt.maxErrExponent();

flopt.setDegreeDist(dd_mee);
stop_mee = flopt.FixedRec_acc(400,'BP');
[err_mee,eco_mee] = flopt.errorProb(stop_mee);
inac_mee = flopt.FixedRec_acc(400,'inac');
ein_mee = flopt.expInac(inac_mee);

save(filename,'dd_mee','stop_mee','err_mee','eco_mee','inac_mee','ein_mee','-append');

% end of the script

