%Sensitivity experiments for the five-box model: 
%To test the sensitivity of the model to parameters, we vary parameters over a range well beyond literature estimates
%To run first the "sensitivity_generation.m" to create a matrix "sensivity" in which each row corresponds to a combination of parameters 

% A = final matrix in which each row corresponds to steady-state model simulation for a given combination of parameters 
% The following number indicate the corresponding column:
%1 = conc deep; 2 = conc AZ; 3 = conc SAZ; 4 = conc pycn; 5 = conc LLML; 
%6 = 15N deep; 7 = 15N AZ; 8 = 15N SAZ; 9 = 15N pycn; 10 = 15N LLML
%11 = d15N deep; 12 = d15N AZ; 13 = d15N SAZ; 14 = d15N pycn; 15 = d15N LLML
%16 = Ws; 17 = Ms; 18 = Wd; 19 = Md; 20 = Mll; 21 = x_Wn; 22 = K_AZ; 23 =
%K_SAZ; 24 = R_SAZ; 25 = R_LL_P; 26 = x_deni; 27 = deni; 28 = overturning ratio
%29 = SO total supply; 30 = SO Adv supply; 31 = SO diff supply; 32 = deep
%total supply; 33 = deep adv supply; 34 = deep diff supply; 35 = total
%supply; 36 = pycnocline recipe; 37 = diff/total ratio; 38 = EP_LL; 39 =
%R_LL; 40 = EP SAZ; 41 = R_SAZ; 42 = EP; 43 = R; 44 = pycn - deep d15N diff;
%45 = mass balance; 46 = isotope balance; 47 = denitrification; 

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

%Array with initial conditions
x0 = [N_deep_ini, N_AZ_ini, N_SAZ_ini, N_pycn_ini, N_LLML_ini, N15_deep_ini, N15_AZ_ini, N15_SAZ_ini, N15_pycn_ini,N15_LLML_ini];
tspan = (0:20000); %Here what is important is to reach the equilibrium (i.e.,no variation with time)

%=======================
%Parameters declaration
%=======================
%see figure 3a in Fripiat et al. (2023) for a description of the parameters
%box volume in l
V_deep               = 9.0000e20; 
V_AZ                 = 6.1e18;
V_SAZ                = 9.76e18;
V_pycn               = 2.7700e+20;
V_LLML               = 5.9000e+19;
V_tot                = V_deep+V_AZ+V_SAZ+V_pycn+V_LLML;

%Constant parameters
K_LLML  = 0.99999;%degree of consumption in the low-latitude surface water
R_AZ    = 1;
eps     = 5.5; %isotope effect of nitrate assimilation
eps_deni_P= 7.25;%net isotope effect of denitrification in the pycnocline (both water column and sedimentary denitrification) 
eps_deni_D= 0;% net isotope effect of denitrification in the deep ocean (only sedimentary denitrification) 
R15ref = 0.0036782;% isotopic ratio of Air N2 (i.e., the international reference for N isotopes)
d15N_N2 = -1;% d15N of N2 fixation

%create a matrix of zero value with a number of rows being equal to the sensitivity matrix (i.e., number of model runs) and a number of column being equal to the model variables
xend = zeros(length(sensitivity(:,1)),length(x0));

%the following loop run the model for each of the scenario (i.e.,combination of parameters) corresponding to a raw in the "sensitivity" matrix
% Column in the "sensitivity": 1 = Ws, Ms=2, Wd=3, Md = 4, Mll = 5, x_Wn = 6, K_AZ = 7, K_SAZ = 8, R_SAZ = 9, R_LL_P = 10, x_deni = 11, deni = 12, MAZ = 13
for i=1:length(sensitivity(:,1))
    Ws      = sensitivity(i,1); 
    Wd      = sensitivity(i,3); 
    Wn      = Ws+Wd; 
    x_Wn    = sensitivity(i,6); % fraction of the advective output passing through the low-latitude mixed layer
    Wnll    = x_Wn.*Wn; 
    Wnp     = (1-x_Wn).*Wn; 
    MAZ     = sensitivity(i,13); 
    Ms      = sensitivity(i,2); 
    Md      = sensitivity(i,4); 
    Mll     = sensitivity(i,5); 
    K_AZ    = sensitivity(i,7); %degree of consumption in AZ
    K_SAZ   = sensitivity(i,8); %degree of consumption in PFZ-SAZ
    R_SAZ_P = sensitivity(i,9); 
    R_SAZ_D = 1-R_SAZ_P; 
    R_LL_P  = sensitivity(i,10); 
    R_LL_D = 1-R_LL_P; 
    x_deni  = sensitivity(i,11); %fraction of denitrification occuring in the pycnocline (i.e., the remaining fraction occurs in the deep ocean)
    deni    = sensitivity(i,12); % total denitrification in µmol yr-1
    deni_P  = x_deni.*deni; % denitrification in the pycnocline (both water column and sedimentary denitrification) = N2 fixation in the pycnocline
    deni_D  = (1-x_deni).*deni; % denitrification in the deep ocean (only sedimentary denitrification) = N2 fixation in the deep
    N15_N2_P  = ((((d15N_N2)/1000)+1)*R15ref).*deni_P; %15N input by N2 fixation in µmol yr-1 (in the pycnocline)
    N15_N2_D  = ((((d15N_N2)/1000)+1)*R15ref).*deni_D; %15N input by N2 fixation in µmol yr-1 (in the deep ocean)
    %=================
    %SOLVING THE ODE for each scenario (i.e., row) of the "sensitivity" matrix
    %=================
    options = odeset('RelTol',1e-7,'AbsTol',1e-7);
    [t,x] = ode15s(@(t,x)fivebox_ocean_isotopes_sensitivity_ode(t,x,V_deep,V_AZ, V_SAZ, V_pycn, V_LLML, Ws, Wd, Ms, Md, Mll, Wnll, Wnp, MAZ, K_AZ, K_SAZ, K_LLML, R_AZ, R_SAZ_P, R_SAZ_D, R_LL_P, R_LL_D, deni_P, deni_D, eps, eps_deni_D, eps_deni_P, N15_N2_P, N15_N2_D),tspan,x0,options);
    % to convert 15N concentration in delta values
    d15N=zeros(length(tspan),5);
    d15N(:,1)     = (((x(:,6)./x(:,1))./R15ref)-1)*1000;
    d15N(:,2)     = (((x(:,7)./x(:,2))./R15ref)-1)*1000;
    d15N(:,3)     = (((x(:,8)./x(:,3))./R15ref)-1)*1000;
    d15N(:,4)     = (((x(:,9)./x(:,4))./R15ref)-1)*1000;
    d15N(:,5)     = (((x(:,10)./x(:,5))./R15ref)-1)*1000;
    % For each scenario, to take steady-state solution (i.e., last row in x) and insert it in the xend matrix 
    xend(i,1:10)=x(length(tspan),1:10);
    xend(i,11:15)=d15N(length(tspan),1:5);
end
   
% to create the final matrix which includes model variables/parameters/offline calcualtion
A = [xend,sensitivity]; % merging the model output (only equilibrium solution) with sensitivity matrix (i.e., model parameters)
% the following corresponds to offline calculation 
A(:,28) = sensitivity(:,1)./(sensitivity(:,1)+sensitivity(:,3));%overturning ratio
A(:,29) = (sensitivity(:,1)+sensitivity(:,2)).*A(:,3).*14./(10^6*10^12); % SO total nitrate supply into the pycnocline (Tg N yr-1)
A(:,30) = (sensitivity(:,1)).*A(:,3).*14./(10^6*10^12); % SO advective nitrate supply into the pycnocline (Tg N yr-1)
A(:,31) = (sensitivity(:,2)).*A(:,3).*14./(10^6*10^12); % SO diffusive nitrate supply into the pycnocline (Tg N yr-1)
A(:,32) = (sensitivity(:,3)+sensitivity(:,4)).*A(:,1).*14./(10^6*10^12); % deep total nitrate supply into the pycnocline (Tg N yr-1)
A(:,33) = (sensitivity(:,3)).*A(:,1).*14./(10^6*10^12); % deep advective nitrate supply into the pycnocline (Tg N yr-1)
A(:,34) = (sensitivity(:,4)).*A(:,1).*14./(10^6*10^12); % deep diffusive nitrate supply into the pycnocline (Tg N yr-1)
A(:,35) = A(:,29)+A(:,32); % total nitrate supply into the pycnocline (Tg)
A(:,36) = A(:,29)./A(:,35); %pycnocline recipe 
A(:,37) = (A(:,31)+A(:,34))./A(:,35); %diffusive/total nitrate supply ratio into the pycnocline
A(:,38) = ((sensitivity(:,6).*(sensitivity(:,1)+sensitivity(:,3))+sensitivity(:,5)).*A(:,4)*14/(10^6*10^12)).*K_LLML; % export production in low-latitude surface water(Tg N yr-1)
A(:,39) = ((sensitivity(:,6).*(sensitivity(:,1)+sensitivity(:,3))+sensitivity(:,5)).*A(:,4)*14/(10^6*10^12)).*K_LLML.*sensitivity(:,10); % remineralization in the pycnocline from export production in the low-latitude surface water(Tg N yr-1)
A(:,40) = (sensitivity(:,1).*A(:,2)+sensitivity(:,2).*A(:,4))*14/(10^6*10^12).*sensitivity(:,8); % export production in PFZ-SAZ (Tg N yr-1)
A(:,41) = (sensitivity(:,1).*A(:,2)+sensitivity(:,2).*A(:,4))*14/(10^6*10^12).*sensitivity(:,8).*sensitivity(:,9); % remineralization in the pycnocline from export production in the PFZ-SAZ (Tg N yr-1)
A(:,42) = A(:,38)+A(:,40);% Export production from both low-latitude surface water and PFZ-SAZ (Tg N yr-1)
A(:,43) = A(:,39)+A(:,41);% Remineralization in the pycnocline from export production from both low-latitude surface water and PFZ-SAZ (Tg N yr-1)
A(:,44) = A(:,14)-A(:,11);% difference between deep ocean and pycnocline d15N
A(:,45) = (A(:,1).*V_deep+A(:,2).*V_AZ+A(:,3).*V_SAZ+A(:,4).*V_pycn+A(:,5).*V_LLML)/V_tot; % mass balance
A(:,46) = (A(:,1).*A(:,11).*V_deep+A(:,2).*A(:,12).*V_AZ+A(:,3).*A(:,13).*V_SAZ+A(:,4).*A(:,14).*V_pycn...
    +A(:,5).*A(:,15).*V_LLML)./(V_deep.*A(:,1)+V_AZ.*A(:,2)+V_SAZ.*A(:,3)+V_pycn.*A(:,4)+V_LLML.*A(:,5)); % isotopic balance
A(:,47) = sensitivity(:,12).*14./(10^6*10^12);%denitrifcation in Tg N yr-1

%to clean up the workspace
clearvars -except keep A sensitivity

%to plot results...can be adapted 
figure(1);
hold on; plot(A(:,36),A(:,14),'ko','MarkerFaceColor','k')
hold on; plot(A(:,36),A(:,11),'ro','MarkerFaceColor','r')
figure(2);
hold on; plot(A(:,36),A(:,4),'ko','MarkerFaceColor','k')
hold on; plot(A(:,36),A(:,1),'ro','MarkerFaceColor','r')
