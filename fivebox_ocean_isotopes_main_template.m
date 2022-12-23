%==============================================
%INITIAL CONDITIONS and time span of simulation
%To set initial conditions
%===============================================
% See Figure 3a in Fripiat et al.(2023) for the architecture of the 5-box model
% N for nitrate concentration (in µmol l-1); 
%deep for deep ocean; AZ for surface water in the Antarctic Zone; SAZ for surface water in the Subantarctic Zone; 
%LLML for surface water in the low-latitude area; pycn for pycnocline
N_deep_ini          = 31; 
N_AZ_ini            = 31;
N_SAZ_ini           = 31;
N_pycn_ini          = 31;
N_LLML_ini          = 31;

R15ref = 0.0036782; % isotopic ratio of Air N2 (i.e., the international reference for N isotopes)

% d15N for delta values
d15N_deep_ini        = 5.4; 
d15N_AZ_ini          = 5.4;
d15N_SAZ_ini         = 5.4;
d15N_pycn_ini        = 5.4;
d15N_LLML_ini        = 5.4;

%15N for 15N concentration (in µmol l-1)
N15_deep_ini         =((((d15N_deep_ini)/1000)+1)*R15ref)*N_deep_ini; 
N15_AZ_ini           =((((d15N_AZ_ini)/1000)+1)*R15ref)*N_AZ_ini;
N15_SAZ_ini          =((((d15N_SAZ_ini)/1000)+1)*R15ref)*N_SAZ_ini;
N15_pycn_ini         =((((d15N_pycn_ini)/1000)+1)*R15ref)*N_pycn_ini;
N15_LLML_ini         =((((d15N_LLML_ini)/1000)+1)*R15ref)*N_LLML_ini;

%Box volume in l
V_deep               = 9.0000e20; 
V_AZ                 = 6.1e18;
V_SAZ                = 9.76e18;
V_pycn               = 2.7700e+20;
V_LLML               = 5.9000e+19;
V_tot                = V_deep+V_AZ+V_SAZ+V_pycn+V_LLML;

%Array with initial conditions
x0 = [N_deep_ini, N_AZ_ini, N_SAZ_ini, N_pycn_ini, N_LLML_ini, N15_deep_ini, N15_AZ_ini, N15_SAZ_ini, N15_pycn_ini,N15_LLML_ini];
tspan = (0:50000); %Model Spin-Up to reach the equilibrium

%=================
%SOLVING THE ODE: to run the model 
%=================

options = odeset('RelTol',1e-7,'AbsTol',1e-7);
[t,x] = ode15s(@fiveboxes_ocean_isotopes_ode,tspan,x0,options);

%=================
%To create the final matrix
%=================
% to convert 15N concentration into delta value
d15N=zeros(length(tspan),5);
d15N(:,1)     = (((x(:,6)./x(:,1))./R15ref)-1)*1000;
d15N(:,2)     = (((x(:,7)./x(:,2))./R15ref)-1)*1000;
d15N(:,3)     = (((x(:,8)./x(:,3))./R15ref)-1)*1000;
d15N(:,4)     = (((x(:,9)./x(:,4))./R15ref)-1)*1000;
d15N(:,5)     = (((x(:,10)./x(:,5))./R15ref)-1)*1000;
% to check the mass and isotopic balance 
conservation=zeros(length(tspan),1);
conservation(:,1)=(x(:,1).*V_deep+x(:,2).*V_AZ+x(:,3).*V_SAZ+x(:,4).*V_pycn+x(:,5).*V_LLML)/V_tot;
conservation(:,2)=(x(:,1).*d15N(:,1).*V_deep+x(:,2).*d15N(:,2).*V_AZ+x(:,3).*d15N(:,3).*V_SAZ+x(:,4).*d15N(:,4).*V_pycn...
    +x(:,5).*d15N(:,5).*V_LLML)./(V_deep.*x(:,1)+V_AZ.*x(:,2)+V_SAZ.*x(:,3)+V_pycn.*x(:,4)+V_LLML.*x(:,5));

% To create a matrix with all interesting model results
% 1 = N_deep, 2 = N_AZ, 3 = N_SAZ, 4 = N_pycn, 5 = N_LLML, 6 = d15N_deep,
% 7 = d15N_AZ, 8 = d15N_SAZ, 9 = d15N_pycn, 10 = d15N_LLML,
% 11 = mass balance, 12 = isotopic balance
Final = [x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),d15N(:,1),d15N(:,2),d15N(:,3),d15N(:,4),d15N(:,5),conservation(:,1),conservation(:,2)];
clearvars -except keep Final 

%%%%%%%%%%%%%% plot results
figure(1);
hold on; plot(Final(:,1),'k','LineWidth',2);
hold on; plot(Final(:,4),'r','LineWidth',2);

figure(2);
hold on;plot(Final(:,6),'k','LineWidth',2);
hold on;plot(Final(:,9),'r','LineWidth',2);

