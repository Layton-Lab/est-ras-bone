function dydt = mod_eqns(t,y,params,varargin)
% Model equations for ca-phosphate model

%% Get variable inputs
do_PTH = true; % do PTH compartment
do_Ctriol = true; % do calcitriol compartment
%do_FGF = false; % do FGF compartment
do_Ca = true; % do Ca compartment
do_Phos = true; % do Phos compartment
do_RAS = true; % do RAS compartment
do_bone = true; % bone compartment

% effects
do_LAT1REff = true;
do_OAT1REff = true;
do_CtriolRsec = true;
do_PTHRsec = true;
do_AT1RPTH = true;
do_CaSRRsec = true;

% BMD
do_BMD = true;

% disease states
%do_CKD = false; % do lower GFR

% drugs
do_ACEi = false;
do_ARB = false;

% estrogen
do_EST = false; % change estrogen with time
do_EST_RAS = do_EST; % estrogen impact on RAS

% infusion experiments
%do_TERIiv = false;
do_ANGII_inf = false;

for ii = 1:2:length(varargin)
    if strcmp(varargin{ii}, 'do_PTH')
        do_PTH = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_Ctriol')
        do_Ctriol = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_Ca')
        do_Ca  = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_Phos')
        do_Phos = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_RAS')
        do_RAS = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_bone')
        do_bone = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_ANGII_inf')
        do_ANGII_inf = varargin{ii+1};
    %elseif strcmp(varargin{ii}, 'do_TERIiv')
        %do_TERIiv = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_EST')
        temp = varargin{ii+1};
        do_EST = temp(1);
        if length(temp) == 2
            do_EST_RAS = temp(2);
        else
            do_EST_RAS = do_EST;
        end
    elseif strcmp(varargin{ii}, 'do_AT1R_bone')
        temp = varargin{ii+1};
        do_LAT1REff = temp(1);
        do_OAT1REff = temp(2);
    elseif strcmp(varargin{ii}, 'do_CtriolPTH_Rsec')
        temp = varargin{ii+1};
        do_CtriolRsec = temp(1);
        do_PTHRsec = temp(2);
    elseif strcmp(varargin{ii}, 'do_AT1RPTH')
        do_AT1RPTH = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_BMD')
        do_BMD = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_ACEi')
        temp = varargin{ii+1};
        do_ACEi = temp(1);
        age_ACEi = temp(2);
        pct_inhib_ACEi = temp(3);
    elseif strcmp(varargin{ii}, 'do_ARB')
        temp = varargin{ii+1};
        do_ARB = temp(1);
        age_ARB = temp(2);
        pct_inhib_ARB = temp(3);
    else
        fprintf('What is this varargin input? %s \n', varargin{ii})
        error('varagin input not available')
    end
end

%% Set variable names
% PTH model
PTmax = y(1); % PT gland max capacity, unitless
PTG   = y(2); % PTH gland pool, unitless (OpenBoneMin.cpp = S)
PTH   = y(3); % Circulating PTH, pmol
% Calcitriol
AOH    = y(4); % 1alpha-hydroxylase, mmol/h 
Ctriol = y(5); % calcitriol, pmol
% FGF23
% FGF23  = y(6); % FGF23 in blood
% Calcium
Ca     = y(7); % extracellular calcium, mmol
% Phosphate
PhosECC = y(8); % extracelllular phosphate, mmol
PhosInt = y(9); % intracellular phosphate, mmol

% TERI compartment
TERISC  = y(10); % exogenous PTH (subcutaneous compartment)

% Bone compartment
OC      = y(11); % osteoclast
OBfast  = y(12); % fast osteoblast
OBslow  = y(13); % slow osteoblast
M       = y(14); % RANK-RANKL
RNK     = y(15); % RANK
L       = y(16); % RANKL
TGFBact = y(17); % active TGFB
TGFB    = y(18); % latent TGFB
ROB1    = y(19); % responding osteoblast
BCL2    = y(20); % BCL2
OPG     = y(21); % OPG
N       = y(22); % RANK-OPG
RX2     = y(23); % RUNX2
CREB    = y(24); % CREB
HAp     = y(25); % hydroxyapatite

% RAS compartment
PRC       = y(26); % 
AGT       = y(27); % pmol/ml
AngI      = y(28); % fmol/ml = pmol/L
AngII     = y(29); % fmol/ml= pmol/L
Ang17     = y(30); % fmol/ml= pmol/L
AngIV     = y(31); % fmol/ml= pmol/L
AT1R      = y(32); % fmol/ml= pmol/L
AT2R      = y(33); % fmol/ml= pmol/L

% femoral neck BMD
BMDfn     = y(34); % normalized

% TERI central compartment
TERI_CENT = y(35);

%% Set parameter names
% NOTE: update updatepars.m to get parameter list from set_params()
OC0 = params(1);
OB0 = params(2);
Vp = params(3);
Vint = params(4);
PTout = params(5);
CtriolMax = params(6);
CtriolMin = params(7);
gamCtriol = params(8);
KPTCtriol = params(9);
alpha_511 = params(10);
gam511 = params(11);
delta511 = params(12);
kdegPTG = params(13);
T70 = params(14);
T71 = params(15);
gamCaEff = params(16);
CaConc0 = params(17);
VmaxCaEff = params(18);
kdegPTH = params(19);
alphaPTH = params(20);
rhoPTH = params(21);
gamCaPTH = params(22);
KCaPTH = params(23);
kdeg_AOH = params(24);
kbase_AOH = params(25);
kdeg_Ctriol = params(26);
deltaAOH_PTH = params(27);
gamAOH_PTH = params(28);
Vmax_PTH_AOH = params(29);
gamAOH_Phos = params(30);
KAOH_Phos = params(31);
PhosMax = params(32);
GFR0 = params(33);
maxTmESTkid = params(34);
Reabs50 = params(35);
T16 = params(36);
PTHconc0 = params(37);
ICa = params(38);
KCtriol = params(39);
Vmax_CtriolGut_Ca = params(40);
GutAbs_Ca0 = params(41);
gamCtriolGut = params(42);
T46 = params(43);
eta_Bone_Phos = params(44);
IPhos = params(45);
Vmax_CtriolGut_Phos = params(46);
GutAbs_Phos0 = params(47);
eta_Ecc2Int = params(48);
eta_Int2Ecc = params(49);
TERICL = params(50);
TERIVC = params(51);
Da = params(52);
E0RANKL = params(53);
EmaxL = params(54);
LsurvOCgam = params(55);
k2 = params(56);
k3 = params(57);
k4 = params(58);
EC50surv = params(59);
PicOCgam = params(60);
Pic0 = params(61);
E0PicOC = params(62);
EmaxPicOC = params(63);
EC50PicOC = params(64);
E0Meff = params(65);
EmaxMeffOC = params(66);
kinOCgam = params(67);
EC50MeffOC = params(68);
FracOBfast = params(69);
ROB0 = params(70);
OBfast0 = params(71);
kb = params(72);
bigDb = params(73);
PicOBgam = params(74);
EmaxPicOB = params(75);
E0PicOB = params(76);
EC50PicOB = params(77);
E0RUNX2kbEff = params(78);
RUNkbMax = params(79);
RUNkbGAM = params(80);
RUNkb50 = params(81);
E0PicOBkb = params(82);
EmaxPicOBkb = params(83);
PicOBgamkb = params(84);
EC50PicOBkb = params(85);
Frackb = params(86);
koutRNK = params(87);
kinRNKgam = params(88);
kinRNK = params(89);
koutL = params(90);
kinLbase = params(91);
k1 = params(92);
OsteoEffectGam = params(93);
EmaxLpth = params(94);
PTH50 = params(95);
koutTGF0 = params(96);
OCtgfGAM = params(97);
OBtgfGAM = params(98);
kinTGF = params(99);
Dr = params(100);
E0PicROB = params(101);
EmaxPicROB = params(102);
PicROBgam = params(103);
EC50PicROB = params(104);
bcl2Kout = params(105);
kO = params(106);
pObase = params(107);
opgPTH50 = params(108);
RX2Kin = params(109);
E0rx2Kout = params(110);
EmaxPTHRX2x = params(111);
EC50PTHRX2x = params(112);
crebKout = params(113);
crebKin0 = params(114);
E0crebKin = params(115);
EmaxPTHcreb = params(116);
EC50PTHcreb = params(117);
kLShap = params(118);
Bone2Plas_Ca0 = params(119);
Vmax_Ca_OC = params(120);
OCgam = params(121);
OC50 = params(122);
MOCratio0 = params(123);
MOCratioGam = params(124);
T15 = params(125);
FracJ15 = params(126);
h_renin = params(127);
h_AGT = params(128);
k_AGT = params(129);
X_PRCPRA = params(130);
KAGT0 = params(131);
Nrs = params(132);
A = params(133);
B = params(134);
AT1R0 = params(135);
c_ACE = params(136);
c_Chym = params(137);
c_NEP = params(138);
c_ACE2 = params(139);
c_IIIV = params(140);
c_AT1R = params(141);
c_AT2R = params(142);
h_AngI = params(143);
h_AngII = params(144);
h_Ang17 = params(145);
h_AngIV = params(146);
h_AT1R = params(147);
h_AT2R = params(148);
rho_AT1RL = params(149);
alpha_AT1RL = params(150);
gamma_AT1RL = params(151);
delta_AT1RL = params(152);
rho_AT1RO = params(153);
alpha_AT1RO = params(154);
gamma_AT1RO = params(155);
delta_AT1RO = params(156);
rho_AT1RPTH = params(157);
alpha_AT1RPTH = params(158);
gamAT1R_PTH = params(159);
KAT1R_PTH = params(160);
rho_CaSR = params(161);
alpha_CaSR = params(162);
gamCaSR = params(163);
KCaSR = params(164);
rho_PTHrenin = params(165);
alpha_PTHrenin = params(166);
gamPTH_renin = params(167);
KPTH_renin = params(168);
alpha_ct = params(169);
rho_ct = params(170);
gamCt_renin = params(171);
KCt_renin = params(172);
E2scalePicB1 = params(173);
tgfbGAM = params(174);
tgfbactGAM = params(175);
robGAM_EST = params(176);
gamOB = params(177);
koutBMDfn = params(178);
gamOCfn = params(179);
OCbase = params(180);
eRen = params(181);
eAGT = params(182);
eACE = params(183);
eAT1R = params(184);
eAT2R = params(185);


%% Model equations
dydt = zeros(length(y),1);

% Concentrations
CtriolConc  = Ctriol/Vp;
PTHConc     = (PTH + TERI_CENT)/Vp; %PTH/Vp;
CaConc      = Ca/Vp;
ECCPhosConc = PhosECC/Vp;
IntPhosConc = PhosInt/Vint;

% Estrogen
if do_EST
    EST = get_estrogen(t); 
else
    EST = 1;
end

% renal insufficiency
GFR = GFR0;
% if do_CKD
%     GFR = get_GFR(t, GFR0);
%     rensuff = GFR/GFR0;
%     alpha = 0.5;
%     m = (1 - alpha)/(0.9 - 0.16);
%     temp =  m * (rensuff - 0.16) + alpha;
%     ScaEff_AOH = min(1, temp);
% else
%     GFR = GFR0;
%     ScaEff_AOH=1;
% end

%% PTH model

    % PT max (PT gland max capacity) 
    CtriolPTeff = CtriolMax - (CtriolMax - CtriolMin) * HillFunc(CtriolConc, gamCtriol, KPTCtriol);
    % Phosphate effect from Gaweda et al 2021
    PhosPTeff = alpha_511 * HillFunc(ECCPhosConc, gam511, delta511);

    % PTG (PTH gland pool)
    CaEff_PTG = VmaxCaEff *(CaConc0/CaConc).^gamCaEff; % from OpenBoneMin.cpp T72 
    T73 = T71 * (CtriolConc - CaEff_PTG);
    T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73)); % tanh

    PTGprod = T70 * (0.85 * (1 - T74) + 0.15); % From OpenBoneMin.cpp (T76 equation)

    % PTH
    CaEff_PTH = alphaPTH - (alphaPTH - rhoPTH) * HillFunc(CaConc, gamCaPTH, KCaPTH); % T63 from OpenBoneMin.cpp
    FCTD = (PTG / 0.5) * PTmax;
    
    if do_AT1RPTH
        AT1REff_PTH = rho_AT1RPTH + (alpha_AT1RPTH - rho_AT1RPTH)*HillFunc(AT1R,gamAT1R_PTH,KAT1R_PTH);
    else
        AT1REff_PTH = rho_AT1RPTH + (alpha_AT1RPTH - rho_AT1RPTH)*HillFunc(3.5475,gamAT1R_PTH,KAT1R_PTH);
    end
    
    %FCTD = 0.8859;
    %CaEff_PTH = 285.6317;
    %AT1REff_PTH = 1.1392;
    PTH_secretion = AT1REff_PTH*CaEff_PTH*FCTD;


    % TERI PK
    k_TERI = TERICL / TERIVC;
    TERIPKin = TERISC * k_TERI; % first order rate subq doing into plasma
    % d(TERISC)/dt
    dydt(10) = - TERIPKin;

    % d(TERI_CENT)/dt -- from OpenBoneMin
    %dydt(35) = TERIPKin - TERI_CENT * TERIKA;
    dydt(35) = TERIPKin - TERI_CENT * kdegPTH;

if do_PTH
    % d(PTmax)/dt
    dydt(1) = PTout * CtriolPTeff * PhosPTeff - PTout * PTmax;

   
    % d(PTG)/dt
    dydt(2) = PTGprod - kdegPTG * PTG; 

    % d(PTH)/dt
    %dydt(3) = PTH_secretion - kdegPTH*PTH + TERIPKin; 
    dydt(3) = PTH_secretion - kdegPTH*PTH;

end % do_PTH

%% Calcitriol model
    % PTH stimulates AOH
    PTHEffAOH = HillFunc(PTHConc, gamAOH_PTH, deltaAOH_PTH);
    etaAOH_PTH = Vmax_PTH_AOH * PTHEffAOH;

    % Phosphate inhibits AOH
    %PhosEff = PhosEff0 - (PhosEff0 * PhosEffMax)*HillFunc(ECCPhosConc, gamAOH_Phos, KAOH_Phos);
    %etaAOH_Phos = min(1, PhosEff);
    %etaAOH_Phos = min(1, alpha_AOHPhos - (alpha_AOHPhos - rho_AOHPhos) * PhosEffAOH);
    etaAOH_Phos = min(1, PhosMax * (KAOH_Phos^gamAOH_Phos / (KAOH_Phos^gamAOH_Phos + ECCPhosConc^gamAOH_Phos)));
if do_Ctriol
    % d(AOH)/dt
    dydt(4) =  kbase_AOH * etaAOH_PTH * etaAOH_Phos - kdeg_AOH * AOH;

    % d(Ctriol)/dt 
    dydt(5) = AOH - kdeg_Ctriol * Ctriol;
end

%% Calcium regulation
    % Urine Calcium (from OpenBoneMin)
    CaFilt = 0.6 * 0.5 * GFR * CaConc;

    % estrogen impact on kidney
    mtmEST = (1-maxTmESTkid)/(1-0.1);
    tmEST = 1 - mtmEST + mtmEST*EST;

    ReabsMax = tmEST * (0.3*GFR*CaConc0 - 0.149997)*(Reabs50 + CaConc0) / CaConc0;
 
    % Effect of PTH on calcium reabsorption
    T17 = PTHconc0*T16 - PTHconc0;

    ReabsPTHeff = (T16*PTHConc)/(PTHConc + T17);

    % PTH-sensitivie calcium reabsorption in kidney
    CaReabsActive =  (ReabsMax*CaConc/(Reabs50 + CaConc))*ReabsPTHeff;

    Urine_Ca = CaFilt - CaReabsActive;


    % Bone2Plas_Ca
    OCEff = HillFunc(OC, OCgam, OC50);
    MOCratioEff = ((M/OC)/MOCratio0)^MOCratioGam;
    Bone2Plas_Ca_OCdepend = Vmax_Ca_OC * OCEff*MOCratioEff;

    Bone2Plas_Ca = Bone2Plas_Ca0 + Bone2Plas_Ca_OCdepend;

    % Plas2Bone_Ca
    OB = OBslow + OBfast;
    %J15 = (T15*P*(1-FracJ15) + T15*P*FracJ15*HAp);
    Plas2Bone_Ca = T15 * Ca * (1 - FracJ15) + T15 * Ca * FracJ15 * HAp;

    % Gut2Plas_Ca
    Gut_abs_CtriolEff = (CtriolConc^gamCtriolGut)./(CtriolConc^gamCtriolGut + KCtriol^gamCtriolGut);
    Gut_frac_abs_Ca = GutAbs_Ca0 + Vmax_CtriolGut_Ca * Gut_abs_CtriolEff;
    Gut2Plas_Ca = ICa*Gut_frac_abs_Ca;

if do_Ca
    % d(Ca)/dt
    dydt(7) = Bone2Plas_Ca - Plas2Bone_Ca - Urine_Ca + Gut2Plas_Ca; 
end % do_Ca

%% Phosphate regulation
    % Urine phosphate (from OpenBoneMin)
    T47 = T46*0.88*GFR;
    J48 = max(0, 0.88*GFR*ECCPhosConc - T47);
    Urine_Phos = J48;   

    % Bone2Plas_Phos
    Bone2Plas_Phos = eta_Bone_Phos * Bone2Plas_Ca;
    % Plas2Bone_Phos
    Plas2Bone_Phos = eta_Bone_Phos * Plas2Bone_Ca;

    % Gut2Plas_Phos
    Gut_frac_abs_Phos = GutAbs_Phos0 + Vmax_CtriolGut_Phos * Gut_abs_CtriolEff;
    Gut2Plas_Phos = IPhos * Gut_frac_abs_Phos;

    % Intra2Extra_Phos % intra 2 extracellular and vise versa
    Ecc2Int_Phos = eta_Ecc2Int * ECCPhosConc;
    Int2Ecc_Phos = eta_Int2Ecc * IntPhosConc;

if do_Phos   
    % d(PhosECC)/dt
    dydt(8) = Bone2Plas_Phos - Plas2Bone_Phos ...
                        - Urine_Phos + Gut2Plas_Phos ...
                        - Ecc2Int_Phos + Int2Ecc_Phos;
    % d(PhosInt)/dt
    dydt(9) = Ecc2Int_Phos - Int2Ecc_Phos;
end % do_Phos

% Bone compartment (from peterson & riggs)
if do_bone
% TGFB act induces OC apoptosis
PicOC = E0PicOC + ((EmaxPicOC*TGFBact^PicOCgam)/(TGFBact^PicOCgam + EC50PicOC^PicOCgam));

% RANKL-RANK inhibits OC decay
PiL = M/10;
LsurvOC = E0RANKL - (E0RANKL - EmaxL)*(PiL^LsurvOCgam/(PiL^LsurvOCgam + EC50surv^LsurvOCgam));

KLSoc = Da*PicOC*LsurvOC;

% RANKL-RANKL increases OC differentiation
MeffOC = E0Meff + EmaxMeffOC * HillFunc(M, kinOCgam, EC50MeffOC);

kinOC2 = Da*Pic0*MeffOC*OC0; 

% d(OC)/dt
dydt(11) = kinOC2 - KLSoc*OC;

PicOB = E0PicOB + EmaxPicOB*TGFBact^PicOBgam / (TGFBact^PicOBgam + EC50PicOB^PicOBgam);
if BCL2 > 105
    RUNX2 = BCL2 - 90;
else
    RUNX2 = 10;
end
RUNX2kbPrimeEff = RUNkbMax*RUNX2^RUNkbGAM / (RUNX2^RUNkbGAM + RUNkb50^RUNkbGAM);
PicOBkb = E0PicOBkb - (E0PicOBkb  - EmaxPicOBkb)*TGFBact^PicOBgamkb / (TGFBact^PicOBgamkb + EC50PicOBkb^PicOBgamkb);
PicOBkbEff = (PicOBkb/Pic0)*(1 / (EST^E2scalePicB1));
kbprime = E0RUNX2kbEff*PicOBkbEff - RUNX2kbPrimeEff;
kbslow = kbprime*Frackb;
kbfast = (kb*OB0 + kbslow*OBfast0 - kbslow*OB0) / OBfast0;
Frackb2 = kbfast/kbprime;
% d(OBfast)/dt
dydt(12) = (bigDb/PicOB)*ROB1*FracOBfast*Frackb2  - kbfast*OBfast;

% d(OBslow)/dt
dydt(13) = (bigDb/PicOB)*ROB1*(1-FracOBfast)*Frackb - kbslow*OBslow;

% d(M)/dt
dydt(14) = k3*RNK*L - k4*M;

% d(RNK)/dt
dydt(15) = kinRNK*TGFBact^kinRNKgam - koutRNK*RNK - k3*RNK*L + k4*M;

OsteoEffect = (OB/OB0)^OsteoEffectGam;
LpthEff = EmaxLpth*(PTHConc) / (PTH50 + (PTHConc));

if do_LAT1REff
    LAT1REff = rho_AT1RL + (alpha_AT1RL - rho_AT1RL) * HillFunc(AT1R, gamma_AT1RL, delta_AT1RL);
else
    LAT1REff = rho_AT1RL + (alpha_AT1RL - rho_AT1RL) * HillFunc(3.372731659146763, gamma_AT1RL, delta_AT1RL);
end

kinL = kinLbase*(OsteoEffect)*LpthEff*LAT1REff;

% d(L)/dt
dydt(16) = kinL - koutL*L - k1*OPG*L + k2*N - k3*RNK*L + k4*M; % - kdenosl*DENMOL*L;

koutTGFeqn = koutTGF0*TGFB*((OC/OC0)^OCtgfGAM);
koutTGFact = koutTGF0 * 1000;
TGFact_EST = EST^tgfbactGAM;
% d(TGFBact)/dt
dydt(17) = koutTGFeqn * TGFact_EST - koutTGFact*TGFBact;

% d(TGFB)/dt
TGF_EST = (1/EST)^tgfbGAM;

dydt(18) = kinTGF*(((OB/OB0)^OBtgfGAM))*TGF_EST - koutTGFeqn * TGFact_EST;

PicROB = E0PicROB + EmaxPicROB*TGFBact^PicROBgam/(TGFBact^PicROBgam + EC50PicROB^PicROBgam);
ROBin = Dr*PicROB;

PicOB = E0PicOB + EmaxPicOB*TGFBact^PicOBgam / (TGFBact^PicOBgam + EC50PicOB^PicOBgam);
KPT =1*(bigDb/PicOB);
% d(ROB1)/dt
dydt(19) = ROBin * (1/EST)^robGAM_EST - KPT*ROB1;

bcl2Kin = RX2*CREB*0.693;
% d(BCL2)/dt
dydt(20) = bcl2Kin - bcl2Kout*BCL2;

pO = pObase*(ROB1/ROB0)*((PTHConc+(opgPTH50*(ROB1/ROB0)))/(2*PTHConc));

% Impact of AT1R on OPG
if do_OAT1REff
    OAT1REff = rho_AT1RO + (alpha_AT1RO - rho_AT1RO)*HillFunc(AT1R, gamma_AT1RO, delta_AT1RO);
else
    OAT1REff = rho_AT1RO + (alpha_AT1RO - rho_AT1RO)*HillFunc(3.372731659146763, gamma_AT1RO, delta_AT1RO);
end

% d(OPG)/dt
dydt(21) = pO * OAT1REff - k1*OPG*L + k2*N - kO*OPG;

% d(N)/dt
dydt(22) = k1*OPG*L - k2*N;

RX2Kout = E0rx2Kout + EmaxPTHRX2x*PTHConc/(PTHConc+EC50PTHRX2x);
% d(RX2)/dt
dydt(23) = RX2Kin - RX2Kout*RX2;

crebKin =crebKin0* (E0crebKin + EmaxPTHcreb*PTHConc/(PTHConc+EC50PTHcreb));
% d(CREB)/dt
dydt(24) = crebKin - crebKout*CREB;

% Hydroxyapatite
dydt(25) = (kLShap/OB0)*OB - kLShap * HAp;
end % do_bone

%% RAS model
if do_RAS
    
    % Do estrogen impacts on RAS
    if do_EST_RAS
        % Estrogen inhibits renin secretion (inhibition)
        etaEST_Ren_max = (eRen + 1)/eRen; 
        etaEST_Ren =  etaEST_Ren_max*g_minus(EST/eRen);
        % Estrogen impact on AGT (activation)
        etaEST_AGTmax = eAGT + 1;
        etaEST_AGT = etaEST_AGTmax * g_plus(EST/eAGT);
        % Estrogen impact on ACE (inhibition)
        etaEST_ACEmax = (eACE + 1)/eACE;
        etaEST_ACE = etaEST_ACEmax * g_minus(EST/eACE);
        % Estrogen impact on AT1R receptor binding (inhibition)
        etaEST_AT1Rmax = (eAT1R + 1)/eAT1R;
        etaEST_AT1R = etaEST_AT1Rmax * g_minus(EST/eAT1R);
        % Estrogen impact on AT2R receptor binding (activation)
        etaEST_AT2Rmax = (eAT2R + 1);
        etaEST_AT2R = etaEST_AT2Rmax * g_plus(EST/eAT2R);
        %etaEST_AT2R = 1;
    else
        etaEST_AGT = 1;
        etaEST_Ren = 1;
        etaEST_AT1R = 1;
        etaEST_ACE = 1;
        etaEST_AT2R = 1;
    end

    eta_exp = A - B * log10(AT1R/AT1R0);
    eta_AT1 = 10.^eta_exp;
    % CaSR impact on renin secretion (inhibition)
    if do_CaSRRsec
        eta_CaSR = alpha_CaSR - (alpha_CaSR - rho_CaSR)*HillFunc(CaConc, gamCaSR, KCaSR);
    else
        eta_CaSR = alpha_CaSR - (alpha_CaSR - rho_CaSR)*HillFunc(CaConc0, gamCaSR, KCaSR);
    end

    % PTH impact on renin secretion (activation)
    if do_PTHRsec
        eta_PTH = rho_PTHrenin + (alpha_PTHrenin - rho_PTHrenin)*HillFunc(PTHConc, gamPTH_renin, KPTH_renin);
    else
        eta_PTH = rho_PTHrenin + (alpha_PTHrenin - rho_PTHrenin)*HillFunc(2.8889, gamPTH_renin, KPTH_renin);
    end
    % Calcitriol impact on renin secretion (inhibition)
    if do_CtriolRsec
        eta_Ctriol = alpha_ct - (alpha_ct - rho_ct)*HillFunc(CtriolConc, gamCt_renin, KCt_renin);
    else
        CtriolConc0 = 94.8497;
        eta_Ctriol = alpha_ct - (alpha_ct - rho_ct)*HillFunc(CtriolConc0, gamCt_renin, KCt_renin);
    end
    R_sec = Nrs * eta_AT1 * eta_CaSR * eta_PTH * eta_Ctriol * etaEST_Ren;
    % d(PRC)/dt
    dydt(26) = R_sec - (log(2)/h_renin) * PRC;
    
    % d(AGT)/dt
    %PRA = PRC * X_PRCPRA; % From BP code (Hallow 2014)
    PRA = PRC * X_PRCPRA * 2 * (AGT/(AGT + KAGT0));
   

    dydt(27) = etaEST_AGT * k_AGT - PRA - (log(2)/h_AGT)*AGT;
    

    % d(AngI)/dt
    if do_ACEi
        t_yrs = t/(365*24);
        if t_yrs > age_ACEi
            if t_yrs < age_ACEi + 1
                % gradual decrease for numerics
                pct_ACEi = (t_yrs - age_ACEi)*pct_inhib_ACEi;
            else
                pct_ACEi = pct_inhib_ACEi;
            end
        else
            pct_ACEi = 0;
        end
    else
        pct_ACEi = 0;
    end
    AngItoAngII = (etaEST_ACE*(1 - pct_ACEi)*c_ACE + c_Chym)*AngI;
    AngItoAng17 = c_NEP * AngI;
    dydt(28) = PRA - AngItoAngII - AngItoAng17 - (log(2)/h_AngI)*AngI;
    
    if do_ANGII_inf
        kinf0 = 130;
        % Graded infusion (Grant 1992)
        if t < 0
            k_inf2 = 0;
        elseif t < 1
            if t < (20/60)
                k_inf2 = kinf0;
            elseif t < (40/60)
                k_inf2 = kinf0*3;
            else
                k_inf2 = kinf0 * 10;
            end
        else
            k_inf2 = 0;
        end
    else
        k_inf2 = 0;
    end

    if do_ARB
        t_yrs = t/(365*24);
        if t_yrs > age_ARB
            if t_yrs < age_ARB + 1
                % gradual decrease for numerics
                pct_ARB = (t_yrs - age_ARB)*pct_inhib_ARB; 
            else
                pct_ARB = pct_inhib_ARB;
            end
        else
            pct_ARB = 0;
        end
    else
        pct_ARB = 0;
    end
    % d(AngII)/dt
    AngIItoAT1R = etaEST_AT1R*(1 - pct_ARB)*c_AT1R * AngII;
    dydt(29) = AngItoAngII ...
               - (c_ACE2 + c_IIIV +...
               + etaEST_AT2R*c_AT2R)*AngII ...
               - AngIItoAT1R...
               - (log(2)/h_AngII)*AngII ...
           + k_inf2;

    % d(AngI7)/dt
    dydt(30) = AngItoAng17 + c_ACE2 * AngII...
        - log(2)/h_Ang17 * Ang17;
    
    % d(AngIV)/dt
    dydt(31) = c_IIIV * AngII - log(2)/h_AngIV * AngIV;
    
    % d(AT1R)/dt
    dydt(32) = AngIItoAT1R - log(2)/h_AT1R * AT1R;
    
    % d(AT2R)/dt
    dydt(33) = etaEST_AT2R*c_AT2R * AngII - log(2)/h_AT2R * AT2R;
end


if do_BMD
    % femoral neck BMD
    kinBMDfn =  koutBMDfn;
    if t < 0
        dydt(34) = 0;
    else
        dydt(34) = kinBMDfn * (OB/OB0)^gamOB - koutBMDfn * (OC/OCbase)^gamOCfn * BMDfn;
    end
end
end % mod_eqns