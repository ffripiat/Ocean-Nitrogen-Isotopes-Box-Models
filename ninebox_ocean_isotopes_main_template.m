%==============================================
%INITIAL CONDITIONS and time span of simulation
%To set the initial conditions
%===============================================
% See Figure 3b in Fripiat et al.(2023) for the architecture of the 9-box model
% N for nitrate concentration (in µmol l-1); 
N_1_ini          = 31; % deep ocean Indo-Pac-South Atlantic; 
N_2_ini          = 31; % deep ocean North Atlantic;
N_3_ini          = 31; % Pyncocline Indo-Pac-South Atlantic
N_4_ini          = 31; % Pyncocline North Atlantic
N_5_ini          = 31; % AZ surface
N_6_ini          = 31; % SAZ surface
N_7_ini          = 31; % LLML Indo-Pac-South Atlantic
N_8_ini          = 31; % LLML North Atlantic
N_9_ini          = 31; % North Atl surface

R15ref = 0.0036782; % isotopic ratio of Air N2 (i.e., the international reference for N isotopes)

% d15N for delta values
d15N_1_ini        = 5.4; 
d15N_2_ini        = 5.4;
d15N_3_ini        = 5.4;
d15N_4_ini        = 5.4;
d15N_5_ini        = 5.4;
d15N_6_ini        = 5.4;
d15N_7_ini        = 5.4;
d15N_8_ini        = 5.4;
d15N_9_ini        = 5.4;

%15N for 15N concentration
N15_1_ini         =((((d15N_1_ini)/1000)+1)*R15ref)*N_1_ini; 
N15_2_ini         =((((d15N_2_ini)/1000)+1)*R15ref)*N_2_ini;
N15_3_ini         =((((d15N_3_ini)/1000)+1)*R15ref)*N_3_ini;
N15_4_ini         =((((d15N_4_ini)/1000)+1)*R15ref)*N_4_ini;
N15_5_ini         =((((d15N_5_ini)/1000)+1)*R15ref)*N_5_ini;
N15_6_ini         =((((d15N_6_ini)/1000)+1)*R15ref)*N_6_ini;
N15_7_ini         =((((d15N_7_ini)/1000)+1)*R15ref)*N_7_ini;
N15_8_ini         =((((d15N_8_ini)/1000)+1)*R15ref)*N_8_ini;
N15_9_ini         =((((d15N_9_ini)/1000)+1)*R15ref)*N_9_ini;

%Box volume in l
V_1               = 7.9200e+20;
V_2               = 1.0800e+20;
V_3               = 2.4376e+20;
V_4               = 3.3240e+19;
V_5               = 6.1e18;
V_6               = 9.76e18;
V_7               = 5.1920e+19;
V_8               = 7.0800e+18;
V_9               = 3e18;
V_tot             = V_1+V_2+V_3+V_4+V_5+V_6+V_7+V_8+V_9;

%Array with initial conditions
x0 = [N_1_ini, N_2_ini, N_3_ini, N_4_ini, N_5_ini, N_6_ini, N_7_ini, N_8_ini, N_9_ini, N15_1_ini, N15_2_ini, N15_3_ini, N15_4_ini, N15_5_ini, N15_6_ini, N15_7_ini, N15_8_ini, N15_9_ini];  
tspan = (0:50000); %Model Spin-Up to reach the equilibrium

%=================
%SOLVING THE ODE: to run the model 
%=================

options = odeset('RelTol',1e-7,'AbsTol',1e-7);
[t,x] = ode15s(@ninebox_ocean_isotopes_ode,tspan,x0,options);

%=================
%To create the final matrix
%=================
% to convert 15N concentration into delta value
d15N=zeros(length(tspan),9);
d15N(:,1)     = (((x(:,10)./x(:,1))./R15ref)-1)*1000;
d15N(:,2)     = (((x(:,11)./x(:,2))./R15ref)-1)*1000;
d15N(:,3)     = (((x(:,12)./x(:,3))./R15ref)-1)*1000;
d15N(:,4)     = (((x(:,13)./x(:,4))./R15ref)-1)*1000;
d15N(:,5)     = (((x(:,14)./x(:,5))./R15ref)-1)*1000;
d15N(:,6)     = (((x(:,15)./x(:,6))./R15ref)-1)*1000;
d15N(:,7)     = (((x(:,16)./x(:,7))./R15ref)-1)*1000;
d15N(:,8)     = (((x(:,17)./x(:,8))./R15ref)-1)*1000;
d15N(:,9)     = (((x(:,18)./x(:,9))./R15ref)-1)*1000;

% to calculate the value for the global deep ocean, pycnocline and LLML
lump=zeros(length(tspan),6);
lump(:,1)      = (x(:,1).*V_1+x(:,2).*V_2)./(V_1+V_2);
lump(:,2)      = (x(:,3).*V_3+x(:,4).*V_4)./(V_3+V_4);
lump(:,3)      = (x(:,7).*V_7+x(:,8).*V_8)./(V_7+V_8);
lump(:,4)      = (x(:,1).*d15N(:,1).*V_1+x(:,2).*d15N(:,2).*V_2)./(x(:,1).*V_1+x(:,2).*V_2);
lump(:,5)      = (x(:,3).*d15N(:,3).*V_3+x(:,4).*d15N(:,4).*V_4)./(x(:,3).*V_3+x(:,4).*V_4);
lump(:,6)      = (x(:,7).*d15N(:,7).*V_7+x(:,8).*d15N(:,8).*V_8)./(x(:,7).*V_7+x(:,8).*V_8);

% to check the mass and isotopic balance 
conservation=zeros(length(tspan),2);
conservation(:,1)=(x(:,1).*V_1+x(:,2).*V_2+x(:,3).*V_3+x(:,4).*V_4+x(:,5).*V_5+x(:,6).*V_6+x(:,7).*V_7+x(:,8).*V_8+x(:,9).*V_9)/V_tot;
conservation(:,2)=(x(:,1).*d15N(:,1).*V_1+x(:,2).*d15N(:,2).*V_2+x(:,3).*d15N(:,3).*V_3+x(:,4).*d15N(:,4).*V_4+x(:,5).*d15N(:,5).*V_5+x(:,6).*d15N(:,6).*V_6+x(:,7).*d15N(:,7).*V_7+x(:,8).*d15N(:,8).*V_8+x(:,9).*d15N(:,9).*V_9)./(x(:,1).*V_1+x(:,2).*V_2+x(:,3).*V_3+x(:,4).*V_4+x(:,5).*V_5+x(:,6).*V_6+x(:,7).*V_7+x(:,8).*V_8+x(:,9).*V_9);

% To create a matrix with all interesting model results
% 1 = conc Deep Indo-Pacific-South Atlantic; 2 = conc deep North Atlantic; 
% 3 = conc pycnocline Indo-Pacific-South Atlantic; 4 = conc pycnocline North Atlantic; 
% 5 = conc AZ surface; 6 = conc SAZ surface; 7 = conc LLML Indo-Pacific-South Atlantic; 
% 8 = conc LLML North Atlantic; 9 = conc North Atl surface; 
% 10 = d15N Deep Indo-Pacific-South Atlantic; 11 = d15N deep North Atlantic; 
% 12 = d15N pycnocline Indo-Pacific-South Atlantic; 13 = d15N pycnocline North Atlantic; 
% 14 = d15N AZ surface; 15 = d15N SAZ surface; 16 = d15N LLML Indo-Pacific-South Atlantic; 
% 17 = d15N LLML North Atlantic; 18 = d15N North Atl surface; 
% 19 = conc Deep; 20 = conc pycnocline; 21 = conc LLML
% 22 = d15N Deep; 23 = d15N pycnocline; 24 = d15N LLML
%25 = mass balance; 26 = isotopic balance
Final = [x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6),x(:,7),x(:,8),x(:,9),d15N(:,1),d15N(:,2),d15N(:,3),d15N(:,4),d15N(:,5),d15N(:,6),d15N(:,7),d15N(:,8),d15N(:,9),lump(:,1),lump(:,2),lump(:,3),lump(:,4),lump(:,5),lump(:,6),conservation(:,1),conservation(:,2)];
clearvars -except keep Final 

%%%%%%%%%%%%%% plot results
figure(1);
hold on; plot(Final(:,19),'k-','LineWidth',2)
hold on; plot(Final(:,20),'r-','LineWidth',2)
figure(2);
hold on; plot(Final(:,22),'k-','LineWidth',2)
hold on; plot(Final(:,23),'r-','LineWidth',2)
