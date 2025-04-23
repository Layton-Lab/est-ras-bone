function name2 = postprocess_name(name)
% use to postprocess parameter names for figures
if strcmp(name, 'eRen')
    name2 = 'e_{renin}';
elseif strcmp(name, 'Bone2Plas_Ca0')
    name2 = '\Phi_{Bone2PlasCa0}';
elseif strcmp(name, 'T15')
    name2 = '\Phi_{Plas2BoneCa0}';
elseif strcmp(name, 'alpha_ct')
    name2 = '\alpha_{Ctriol, renin}';
elseif strcmp(name, 'EmaxMeffOC')
    name2 = 'E_{max, M-OC}';
elseif strcmp(name, 'EmaxL')
    name2 = 'E_{max, L}';
elseif strcmp(name, 'IPhos')
    name2 = 'I_{Phos}';
elseif strcmp(name, 'delta_AT1RL')
    name2 = '\delta_{AT1R, L}';
elseif strcmp(name, 'alpha_CaSR')
    name2 = '\alpha_{CaSR, renin}';
elseif strcmp(name, 'h_AT2R')
    name2 = 'h_{AT2R}';
elseif strcmp(name, 'h_AngII')
    name2 = 'h_{AngII}';
elseif strcmp(name, 'KPTH_renin')
    name2 = '\delta_{PTH, renin}';
elseif strcmp(name, 'KCaSR')
    name2 = '\delta_{CaSR, renin}';
elseif strcmp(name, 'h_AT1R')
    name2 = 'h_{AT1R}';
elseif strcmp(name, 'eAGT')
    name2 = 'e_{AGT}';
elseif strcmp(name, 'eACE')
    name2 = 'e_{ACE}';
elseif strcmp(name, 'eAT1R')
    name2 = 'e_{AT1R}';
elseif strcmp(name, 'eAT2R')
    name2 = 'e_{AT2R}';
elseif strcmp(name, 'c_ACE')
    name2 = 'c_{ACE}';
elseif strcmp(name, 'c_AT1R')
    name2 = 'c_{AT1R}';
elseif strcmp(name, 'c_AT2R')
    name2 = 'c_{AT2R}';
elseif strcmp(name, 'koutL')
    name2 = 'k_{out,L}';
elseif strcmp(name, 'k3')
    name2 = 'k_3';
elseif strcmp(name, 'kinRNK')
    name2 = 'k_{in, RNK}';
elseif strcmp(name, 'k4')
    name2 = 'k_4';
elseif strcmp(name,'koutRNK')
    name2 = 'k_{out,RNK}';
elseif strcmp(name, 'EmaxLpth')
    name2 = 'E_{max,L-PTH}';
elseif strcmp(name, 'KCaPTH')
    name2 = 'K_{Ca,PTH}';
elseif strcmp(name, 'kinLbase')
    name2 = 'k_{in, L}';
elseif strcmp(name, 'kdeg_AOH')
    name2 = 'k_{deg, AOH}';
elseif strcmp(name, 'Vmax_PTH_AOH')
    name2 = 'V_{max, PTH-AOH}';
elseif strcmp(name, 'kbase_AOH')
    name2 = 'k_{base, AOH}';
elseif strcmp(name, 'PhosMax')
    name2 = 'V_{max, Phos}';
elseif strcmp(name, 'KAOH_Phos')
    name2 = '\delta_{AOH, Phos}';
elseif strcmp(name, 'X_PRCPRA')
    name2 = 'X_{PRCPRA}';
elseif strcmp(name, 'Nrs')
    name2 = 'N_{rs}';
elseif strcmp(name, 'h_renin')
    name2 = 'h_{renin}';
elseif strcmp(name, 'rho_AT1RL')
    name2 = '\rho_{AT1R,L}';
elseif strcmp(name, 'EC50MeffOC')
    name2 = '\delta_{RNKL-RNK, OC}';
elseif strcmp(name, 'KCt_renin')
    name2 = '\delta_{Ctriol, renin}';
elseif strcmp(name, 'VmaxCaEff')
    name2 = 'V_{max, Ca-PTH}';
elseif strcmp(name, 'T16')
    name2 = 'V_{max, PTH-reabs}';
elseif strcmp(name, 'T70')
    name2 = 'k_{prod, PTG}';
elseif strcmp(name, 'gamCaPTH')
    name2 = '\gamma_{Ca,PTH}';
elseif strcmp(name, 'deltaAOH_PTH')
    name2 = '\delta_{AOH,PTH}';
else
    name2 = regexprep(name, '_', '\\_');
end % if
end % function