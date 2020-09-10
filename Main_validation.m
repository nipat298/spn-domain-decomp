simulation ='E';
SettingsFile = 'sample_el';

%% basic SPN
disp('Running for basic SPN')
func = @FwdSolver;
fcs_flag=0;
[mesh_spn,Jd,Fluence_spn]=Datagen_validation(simulation,SettingsFile,func,fcs_flag);
Jd_spn = Jd{1};
%%
% fname=[SettingsFile,'45_spn_hg_1207.mat'];
% save(fname)