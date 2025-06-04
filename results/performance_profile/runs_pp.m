clear;

load('BIQ_all_NAL.mat');
perf('table_BIQ_pp',BIQ_all_NAL, 1);

load('QKP_all.mat');
perf('table_QKP_pp',QKP_all,1);
