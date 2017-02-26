function [t,final,finalc]=ABN_mathematical_model_new()

%Original description of model from Kwong et al., PNAS (2015), 112(41):12627. 
%Constants taken from Kwong et al., unless noted

close all;
clc;

% Calls the ODE function below  over 200 minutes with the initial 
% concentration in the blood of 2.5 uM and 0 uM everywhere else
options = odeset('AbsTol',1e-20,'RelTol',1e-20);
[t,y] = ode15s(@(t,y) model(t,y),0:200, [2.5e-6, 0, 0,0,0],options);

%Tumor bearing values from ODE
Np=y(:,1); %Nanoparticle in plasma
Nt=y(:,2); %Nanoparticle in tumor
Rt=y(:,3); %Reporter in tumor
Rp=y(:,4); %Reporter in plasma
Ru=y(:,5); %Reporter in urine 

% Plot urine reporter
plot(t,Ru,'b');
hold on;

%Control (-tumor)
[t2,y2] = ode15s(@(t,y) model2(t,y),0:200, [2.5e-6, 0, 0,0,0],options);
  
finalc = y2(:,5); %Reporter in urine for control sample
final = Ru; %Reporter in urine for tumor sample

plot(t2,finalc,'m');
legend('Model Tumor','Model Control','Location', 'Best')


%defining ODE: calls a +tumor case
function dydt = model(t,y)

%clearing constants
mps_np = 4.62e-02; %m^-1 Reflect clearance rate from Kwon/Dudani et al. (0.25 hr blood half life.)  
mps_r = 0.0064; %m^-1 Reporter clearance
k_np_filter = 1.66e-05; %m^-1 Urinary filtration of NP
k_r_filter = 0.032; %m^-1 Urinary filtration of reporter
k_absorb = 2.89; % [M]/min Fit to account for reporter loss.

%permeability constants
p_np_tissue = 1.4e-04*1.18; %m^-1 Transport of NP into tumor. 1.18 term from improved accumulation (Fig. 2, Kwon/Dudani et al.). 
p_r_tissue = 0.09; %m^-1 Transport of reporter

%%%PLASMA
%plasma enzyme kinetic constants
k_cat_bckg = 0.0659/29.15; % [1/min] Plasma cleavage Kcat. 29.5 term from decreased thrombin cleavage (Fig. 2, Kwon/Dudani et al.).
k_m_bckg = 1.0063e-05; % [M] Plasma cleavage Km

%TUMOR
%tumor enzyme kinetic constants 
k_cat = 0.5*3.99; % [1/min] MMP9 cleavage kcat. 3.99 term from increased MMP9 cleavage (Fig. 2, Kwon/Dudani et al.).
k_m = 2.13e-06; % [M] MMP9 cleavage Km

%enzymes
e_nonspecific = 4e-6; % [M] 
e_tissue = 7.1624e-7; % [M] 
e_plasma = 7.1623e-8;% [M] Plasma = 1/10th [Tissue], see Hori et al., STM (2011), 3(109):109.

%%%%State Conditions
C_np_plasma = y(1);
C_np_tissue = y(2);
C_r_tissue = y(3);
C_r_plasma = y(4);
C_r_urine = y(5);

%Nanoparticle concentration in plasma: function of circulation time, clearance into
%tumor, cleaving due to free, nonspecific proteases, and cleaving due to secreted
%enzymes from tumor
dC_np_plasmadt = -(mps_np+k_np_filter)*C_np_plasma - (p_np_tissue*(C_np_plasma-C_np_tissue))- (k_cat_bckg*e_nonspecific*C_np_plasma)/(k_m_bckg+C_np_plasma) - (k_cat*e_plasma*C_np_plasma)/(k_m+C_np_plasma);

%Nanoparticle concentration in tissue: function of clearance into
%tumor and cleaving due to secreted enzymes
dC_np_tissuedt = (p_np_tissue*(C_np_plasma-C_np_tissue))- (k_cat*e_tissue*C_np_tissue)/(k_m+C_np_tissue);

%Reporter concentration in tissue: function of cleaving of nps and clearance into
%blood
dC_r_tissuedt =(k_cat*e_tissue*C_np_tissue)/(k_m+C_np_tissue) - (p_r_tissue*(C_r_tissue));

%Reporter concentration in plasma: function of clearance into blood from
%tissue, natural decay, and cleavage in blood of np, and clearance into
%urine
dC_r_plasmadt = (p_r_tissue*(C_r_tissue)) - (k_r_filter+mps_r+k_absorb)*C_r_plasma + (k_cat_bckg*e_nonspecific*C_np_plasma)/(k_m_bckg+C_np_plasma) + (k_cat*e_plasma*C_np_plasma)/(k_m+C_np_plasma);

%reporter concentration in urine: function of clearance into urine from
%blood
dC_r_urinedt = k_r_filter*(C_r_plasma);

dydt = [dC_np_plasmadt; dC_np_tissuedt; dC_r_tissuedt; dC_r_plasmadt; dC_r_urinedt];
return 

%Control (-tumor case)
function dydt = model2(t,y)

%clearing constants
mps_np = 4.62e-02; %m^-1 Reflect clearance rate from Kwon/Dudani et al. 0.25 hr blood half life. 
mps_r = 0.0064; %m^-1 Reporter clearance
k_np_filter = 1.7e-05; %m^-1 Urinary filtration of NP
k_r_filter = 0.031757007798849; %m^-1 Urinary filtration of reporter
k_absorb = 2.8915; % [M]/min Fit to account for reporter loss.

%Multiply tumor terms by 0
%permeability constants
p_np_tissue = 0*1.4e-04; %m^-1 Transport of NP into tumor
p_r_tissue = 0*0.09; %m^-1 Transport of reporter

%%%PLASMA
%plasma enzyme kinetic constants
k_cat_bckg = 0.0659/29.15; % [1/min] Plasma cleavage Kcat. 29.5 term from decreased thrombin cleavage (Fig. 2).
k_m_bckg = 1.0063e-05; % [M] Plasma cleavage Km

%TUMOR
%tumor enzyme kinetic constants 
k_cat = 0.5*3.99; % [1/min] MMP9 cleavage kcat. 
k_m = 2.13e-06; % [M] MMP9 cleavage Km

%Multiply tumor terms by 0
%enzymes
e_nonspecific = 4e-6; % [M] 
e_tissue = 0*7.1624e-7; % [M] 
e_plasma = 0*7.1623e-8;% [M] Plasma = 1/10th [Tissue], see Hori et al., STM (2011), 3(109):109.

%%%%State Conditions
C_np_plasma = y(1);
C_np_tissue = y(2);
C_r_tissue = y(3);
C_r_plasma = y(4);
C_r_urine = y(5);

%Nanoparticle concentration in plasma: function of circulation time, clearance into
%tumor, cleaving due to free, nonspecific proteases
dC_np_plasmadt = -(mps_np+k_np_filter)*C_np_plasma - (p_np_tissue*(C_np_plasma-C_np_tissue))- (k_cat_bckg*e_nonspecific*C_np_plasma)/(k_m_bckg+C_np_plasma) - (k_cat*e_plasma*C_np_plasma)/(k_m+C_np_plasma);

%Nanoparticle concentration in tissue: function of clearance into
%tumor and cleaving due to secreted enzymes... in control remains at 0.
dC_np_tissuedt = (p_np_tissue*(C_np_plasma-C_np_tissue))- (k_cat*e_tissue*C_np_tissue)/(k_m+C_np_tissue);

%reporter concentration in tissue: function of cleaving of nps and clearance into
%blood... in control remains at 0.
dC_r_tissuedt =(k_cat*e_tissue*C_np_tissue)/(k_m+C_np_tissue) - (p_r_tissue*(C_r_tissue));

%reporter concentration in plasma: natural decay, cleavage in blood of np, and clearance into
%urine
dC_r_plasmadt = (p_r_tissue*(C_r_tissue)) - (k_r_filter+mps_r+k_absorb)*C_r_plasma + (k_cat_bckg*e_nonspecific*C_np_plasma)/(k_m_bckg+C_np_plasma) + (k_cat*e_plasma*C_np_plasma)/(k_m+C_np_plasma);

%reporter concentration in urine: function of clearance into urine from
%blood
dC_r_urinedt = k_r_filter*(C_r_plasma);

dydt = [dC_np_plasmadt; dC_np_tissuedt; dC_r_tissuedt; dC_r_plasmadt; dC_r_urinedt];
return 
