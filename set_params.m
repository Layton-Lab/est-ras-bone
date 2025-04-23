function pars = set_params()
% Set parameters
pars.OC0 = 0.00115398; 
pars.OB0 = 0.00501324; 

%% Volumes
pars.Vp = 14; % plasma volume
pars.Vint = 24; % intracellular volume

%% PTH model
% PTmax
pars.PTout = 0.0001604; % OpenBoneMin.cpp
pars.CtriolMax = 2;  % OpenBoneMin.cpp
pars.CtriolMin = 0.9; % OpenBoneMin.cpp
pars.gamCtriol = 12.5033; % OpenBoneMin.cpp CtriolPTgam
pars.KPTCtriol = 68.3805; % this value was solved using OpenBoneMin.cpp equations for Ctriol50%exp(log(INparenCtriol) / pars.gamCtriol); % OpenBoneMin.cpp Ctriol50
pars.alpha_511 = 1.5178; % Gaweda et al 2021
pars.gam511 =  2.2864; % Gaweda et al 2021
pars.delta511 = 0.7; %0.8; %0.9; % Gaweda et al 2021


% PTG
pars.kdegPTG = 0.02; % 2 * T70 * 0.85 + 2 * 0.15 * T70; solved from OpenBoneMin.cpp 
pars.T70 = 0.01; % parameter from OpenBoneMin.cpp
pars.T71 = 0.03; % parameter from OpenBoneMin.cpp
pars.gamCaEff = 0.9; % parameter from OpenBoneMin.cpp (ScaEffGam)
pars.CaConc0 = 2.3; % baseline CaConc
% this par depends on Ctriol level
pars.VmaxCaEff = 90; % hardcoded into T72 in OpenBoneMin.cpp (set to 90)

% PTH
pars.kdegPTH = 100/14; % kout from OpenBoneMin.cpp
pars.alphaPTH = 6249.09; %6249.09; % T58 from OpenBoneMin.cpp
pars.rhoPTH = 96.25; % T61 from OpenBoneMin.cpp
pars.gamCaPTH = 11.7387; % T59 from OpenBoneMin.cpp
pars.KCaPTH = 1.7796; % solved from OpenBoneMin.cpp (T60)

%% Calcitriol model
pars.kdeg_AOH = 0.05; %
pars.kbase_AOH = 6.45;  % 
pars.kdeg_Ctriol = 0.1; % 

% PTH impact on AOH
pars.deltaAOH_PTH = 2; 
pars.gamAOH_PTH = 0.75; 
pars.Vmax_PTH_AOH = 1.9037;

% Phos impact on AOH (Phos inhibits AOH)
pars.gamAOH_Phos = 2; 
pars.KAOH_Phos = 1.3; 
pars.PhosMax = 1.8155;

%% Calcium homeostasis
% Renal handling (OpenBoneMin)
pars.GFR0 = 100.0/16.667;
pars.maxTmESTkid = 0.923737;
pars.Reabs50 = 1.57322;
pars.T16 = 1.3; 
pars.PTHconc0 = 2.8;


% Gut2Plas_Ca
CaDay = 40;  % from OpenBoneMin
pars.ICa = CaDay/24; % calcium intake, CaDay/24 from OpenBoneMin

pars.KCtriol = 90; 
pars.Vmax_CtriolGut_Ca = 0.35; 
pars.GutAbs_Ca0 = 0.1;


pars.gamCtriolGut = 2;

%% Phosphate homeostasis
% Renal handling (OpenBoneMin)
pars.T46 = 1.142;

% Bone transport
pars.eta_Bone_Phos = 0.464; % OpenBoneMin -- phosphate scale factor

% Gut2Plas
PhosDay = 22.6; % mmol/day 
pars.IPhos = PhosDay/24; % Phosphate intake
pars.Vmax_CtriolGut_Phos = 0.3; % 
pars.GutAbs_Phos0 = 0.5; % basal gut fractional Phos absorption

% ECC 2 Intra Phosphate
eta = 51.8;
ECCPhosConcBase = 1.22;
INTPhosConcBase = 100;
pars.eta_Ecc2Int = eta; % OpenBoneMin (T49)
pars.eta_Int2Ecc = eta *(ECCPhosConcBase/INTPhosConcBase); % OpenBoneMin (T55)

% TERI
pars.TERICL = 62.2; % from OpenBoneMin
pars.TERIVC = 94.4; % from OpenBoneMin

% Bone compartment (Peterson & Riggs update)
pars.Da = 0.0292; % OpenBoneMin.cpp
pars.E0RANKL = 3.80338; % OpenBoneMin.cpp
pars.EmaxL = 0.469779; 
pars.LsurvOCgam =3.09023;

RANKL0 = 0.4; 
pars.k2 = 0.112013;
pars.k3 = 0.00000624;
pars.k4 = 0.112013;
PiL0 = (pars.k3/pars.k4)*RANKL0;
EC50survInPar = (pars.E0RANKL - pars.EmaxL)*(PiL0^pars.LsurvOCgam/(pars.E0RANKL - 1)) - PiL0^pars.LsurvOCgam;
pars.EC50surv = exp(log(EC50survInPar)/pars.LsurvOCgam);

pars.PicOCgam = 1.0168;

FracPicOC = 0.878215;
pars.Pic0 = 0.22814;
pars.E0PicOC = FracPicOC*pars.Pic0;
pars.EmaxPicOC = 1.9746;

TGFBact0 = pars.Pic0;
EC50PicOCparen = (pars.EmaxPicOC*TGFBact0^pars.PicOCgam/(pars.Pic0 - pars.E0PicOC)) - TGFBact0^pars.PicOCgam;
pars.EC50PicOC = exp(log(EC50PicOCparen)/pars.PicOCgam);



pars.E0Meff = 0.388267;
pars.EmaxMeffOC = 3.15667;
pars.kinOCgam = 6; 
L_init = 0.4;
RNK_init = 10.0;
M0 = pars.k3*RNK_init*L_init/pars.k4;
pars.EC50MeffOC = exp(log(M0^pars.kinOCgam*pars.EmaxMeffOC/(1-pars.E0Meff) - M0^pars.kinOCgam)/pars.kinOCgam);

pars.FracOBfast = 0.797629;

ROB1_init = 0.00104122;
pars.ROB0 = ROB1_init;


pars.OBfast0 = pars.OB0*pars.FracOBfast;

pars.kb = 0.000605516;
pars.bigDb = pars.kb*pars.OB0*pars.Pic0/pars.ROB0;
pars.PicOBgam = 0.122313;
pars.EmaxPicOB = 0.251636;
FracPicOB = 0.000244818;
pars.E0PicOB = FracPicOB*pars.Pic0;
EC50PicOBparen = (pars.EmaxPicOB*TGFBact0^pars.PicOBgam/(pars.Pic0 - pars.E0PicOB)) - TGFBact0^pars.PicOBgam;
pars.EC50PicOB = exp(log(EC50PicOBparen)/pars.PicOBgam);

E0RUNX2kbEffFACT = 1.01;
RUNkbMaxFact = 0.638114;
pars.E0RUNX2kbEff= E0RUNX2kbEffFACT*pars.kb;
pars.RUNkbMax = pars.E0RUNX2kbEff*RUNkbMaxFact;
pars.RUNkbGAM = 3.67798;
RUNX20 = 10;
INparen = (pars.RUNkbMax * RUNX20^pars.RUNkbGAM) / (pars.E0RUNX2kbEff - pars.kb) - RUNX20^pars.RUNkbGAM;
pars.RUNkb50 = exp(log(INparen)/pars.RUNkbGAM);
MultPicOBkb = 3.11842;
pars.E0PicOBkb = MultPicOBkb*pars.Pic0;
FracPic0kb = 0.764028;
pars.EmaxPicOBkb = FracPic0kb*pars.Pic0;
pars.PicOBgamkb = 2.92375;
EC50PicOBparenKb = ((pars.E0PicOBkb - pars.EmaxPicOBkb)*TGFBact0^pars.PicOBgamkb) / (pars.E0PicOBkb - pars.Pic0)  - TGFBact0^pars.PicOBgamkb;
pars.EC50PicOBkb = exp(log(EC50PicOBparenKb)/pars.PicOBgamkb);
pars.Frackb = 0.313186;

RNK0 = RNK_init;
pars.koutRNK = 0.00323667;
pars.kinRNKgam = 0.151825;
pars.kinRNK = (pars.koutRNK*RNK0 + pars.k3*RNK0*RANKL0 - pars.k4*M0) / TGFBact0^pars.kinRNKgam;

pars.koutL = 0.00293273;
pars.kinLbase = pars.koutL*RANKL0;
pars.k1 = 0.00000624;

pars.OsteoEffectGam = 0.173833;
pars.EmaxLpth = 1.30721;
pars.PTH50 = 0.5; %pars.EmaxLpth*3.85 - 3.85;

pars.koutTGF0 = 0.0000298449;
pars.OCtgfGAM = 0.593891;
pars.OBtgfGAM = 0.0111319;

TGFB0 = pars.Pic0*1000;
pars.kinTGF = pars.koutTGF0*TGFB0;

pars.Dr = pars.kb*pars.OB0/pars.Pic0;
FracPicROB = 0.883824;
pars.E0PicROB = FracPicROB*pars.Pic0;
pars.EmaxPicROB = 3.9745;
pars.PicROBgam = 1.80968;
EC50PicROBparen= (pars.EmaxPicROB*TGFBact0^pars.PicROBgam / (pars.Pic0 - pars.E0PicROB)) - TGFBact0^pars.PicROBgam;
pars.EC50PicROB = exp(log(EC50PicROBparen)/pars.PicROBgam);
pars.bigDb = pars.kb*pars.OB0*pars.Pic0/pars.ROB0;

pars.bcl2Kout = 0.693;

OPG0 = 4; % O_init
pars.kO = 15.8885;
pars.pObase = pars.kO*OPG0;
pars.opgPTH50 = 2.225183665235689;

RX2Kout0 = 0.693;
RX2_init = 10.0;
RX20 = RX2_init;
pars.RX2Kin = RX2Kout0*RX20;
pars.E0rx2Kout = 0.125;
pars.EmaxPTHRX2x = 5;
pars.EC50PTHRX2x = 17.362700711839036;

pars.crebKout = 0.00279513;
CREB_init = 10.0;
CREB0 = CREB_init;
pars.crebKin0 = pars.crebKout*CREB0;
pars.E0crebKin = 0.5;
pars.EmaxPTHcreb = 3.39745;
pars.EC50PTHcreb = 12.894716821674301;

HApMRT = 3.60609;
pars.kLShap = 1/HApMRT;

Q_0 = 100.0;
T13 = (CaDay/24)/Q_0;
FracJ14 = 0.107763;
pars.Bone2Plas_Ca0 = T13*Q_0*(1-FracJ14);

pars.Vmax_Ca_OC = 0.543488 * Q_0 * FracJ14; 
pars.OCgam = 1.6971;

pars.OC50 = exp(log((pars.Vmax_Ca_OC*(pars.OC0^pars.OCgam)/T13) - (pars.OC0^pars.OCgam))/pars.OCgam); 

pars.MOCratio0 = M0/pars.OC0;
pars.MOCratioGam = 0.603754;

pars.T15 = CaDay/(pars.CaConc0*pars.Vp*24);
pars.FracJ15 = 0.114376;


%% RAS parameters
pars.h_renin = 12/60; % hours %12; % mins
pars.h_AGT = 10; % hours (Lo et al 2011)
pars.k_AGT = 610.39 * 60; % fmol/ml/hour % from female human get_pars in BP-regulation

pars.X_PRCPRA = 61/60 * 60; % fmol/hour/pg %61/60.0; %fmol/min/pg
pars.KAGT0 = 513525.02;%520385; % AGT impact, new parameter
pars.Nrs = 0.8*60; %60; % Leete & Layton  % Hallow et al 2014
pars.A = 0.0102;  % A_AT1-renin, Hallow et al 2014
pars.B = 0.95;  % B_AT1-renin in BP-regulation code, Hallow et al 2014
pars.AT1R0 = 3.75; 

pars.c_ACE = 1.4079 * 60; % (female BP-regulation) 1/hr 0.88492; % 1/min
pars.c_Chym = 0.1482 * 60; % (female BP-regulation) 1/hr % 0.09315; % 1/min
pars.c_NEP = 0.060759 * 60; %(female BP-regulation) 1/hr %0.038189; % 1/min

pars.c_ACE2 =  0.0037603 * 60; % (female, BP-regulation) 1/hr 0.0078009; % 1/min
pars.c_IIIV =  0.038644 * 60; % (female, BP-regulation) 1/hr 0.25056; % 1/min
pars.c_AT1R = 0.027089 * 60; % (female, BP-regulation) 1/hr 0.17008; % 1/min
pars.c_AT2R = 0.038699 * 60; % (female, BP-regulation) 1/hr %0.065667; % 1/min

pars.h_AngI = 0.5/60; % hours %0.5; % min
pars.h_AngII = 0.66/60; % hours %0.66; % min

pars.h_Ang17 = 30/60; % hours %30; % min
pars.h_AngIV = 0.5/60; % hours %0.5; % min
pars.h_AT1R  = 12/60; % hours 12; % min
pars.h_AT2R  = 12/60; % hours %12; % min

%% RAS impact on bone
% RANKL effect
pars.rho_AT1RL = 0.95;  % min value
pars.alpha_AT1RL = 5;  % max value (8-fold from Shimizu in vitro)
pars.gamma_AT1RL = 3;  % steepness
pars.delta_AT1RL = 16; % half max

% OPG effect
pars.rho_AT1RO = 0.90; 
pars.alpha_AT1RO = 2.5;
pars.gamma_AT1RO = 2; 
pars.delta_AT1RO = 14; 

%% RAS impact on PTH
pars.rho_AT1RPTH = 0.7;  % fit to Grant 1992 data
pars.alpha_AT1RPTH = 3.2; % fit to Grant 1992 data
pars.gamAT1R_PTH = 1.5;  % fit to Grant 1992 data
pars.KAT1R_PTH = 10; % fit to Grant 1992 data

%% CaSR impact on renin
pars.rho_CaSR = 0.25;
pars.alpha_CaSR = 1.2;
pars.gamCaSR = 15;
pars.KCaSR = 2.6;

%% PTH impact on renin
pars.rho_PTHrenin =  0.5; 
pars.alpha_PTHrenin = 2.5;  
pars.gamPTH_renin = 3;
pars.KPTH_renin = 4; 

%% Calcitriol impact on renin
pars.alpha_ct = 1.11; 
pars.rho_ct= 0.8; 
pars.gamCt_renin = 1; 
pars.KCt_renin = 180; 

%% Estrogen effects on bone
pars.E2scalePicB1 = 0.0000116832;
pars.tgfbGAM = 0.0374;
pars.tgfbactGAM = 0.045273;
pars.robGAM_EST = 0.16;

%% BMD
pars.gamOB = 0.0793;
koutBMDfnBAS = 0.000005651;
pars.koutBMDfn = koutBMDfnBAS;
gamOCfnBAS = 0.3101;
pars.gamOCfn = gamOCfnBAS;
pars.OCbase = 0.001714042543982; 

%% Estrogen impact on RAS
pars.eRen = 0.15;  
pars.eAGT = 0.15; 
pars.eACE = 5; 
pars.eAT1R = 0.5; 
pars.eAT2R = 0.05; 
end