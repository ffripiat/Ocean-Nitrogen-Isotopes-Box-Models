%Sensitivity experiments for the nine-box model: 
%To test the sensitivity of the model to parameters, we vary parameters over a range well beyond literature estimates
%To run first the "ninebox_sensitivity_generation.m" to create a matrix "sensivity" in which each row corresponds to a combination of parameters 

% A = matrix with all model input (model variables, parameters and offline calculation)
% The following number indicate the corresponding column
% 1 = conc Deep Indo-Pacific-South Atlantic; 2 = conc deep North Atlantic; 
% 3 = conc pycnocline Indo-Pacific-South Atlantic; 4 = conc pycnocline North Atlantic; 
% 5 = conc AZ surface; 6 = conc SAZ surface; 7 = conc LLML Indo-Pacific-South Atlantic; 
% 8 = conc LLML North Atlantic; 9 = conc North Atl surface; 
% 10 = 15Nconc Deep Indo-Pacific-South Atlantic; 11 = 15Nconc deep North Atlantic; 
% 12 = 15Nconc pycnocline Indo-Pacific-South Atlantic; 13 = 15Nconc pycnocline North Atlantic; 
% 14 = 15Nconc AZ surface; 15 = 15Nconc SAZ surface; 16 = 15Nconc LLML Indo-Pacific-South Atlantic; 
% 17 = 15Nconc LLML North Atlantic; 18 = 15Nconc North Atl surface; 
% 19 = 15Nconc Deep Indo-Pacific-South Atlantic; 20 = 15Nconc deep North Atlantic; 
% 21 = 15Nconc pycnocline Indo-Pacific-South Atlantic; 22 = 15Nconc pycnocline North Atlantic; 
% 23 = 15Nconc AZ surface; 24 = 15Nconc SAZ surface; 25 = 15Nconc LLML Indo-Pacific-South Atlantic; 
% 26 = 15Nconc LLML North Atlantic; 27 = 15Nconc North Atl surface; 
% 28 = conc Deep; 29 = conc pycnocline; 30 = conc LLML
% 31 = d15N Deep; 32 = d15N pycnocline; 33 = d15N LLML
%34 = Ws; 35 = Ms; 36 = Wd; 37 = Md; 38 = Mll; 39 = X3; 40 = K_AZ; 41 =
%K_SAZ; 42 = R_SAZ; 43 = R_LL_P; 44 = x_deni; 45 = deni; 46 = MAZ; 47 = X4
%48 = M34; 49 = M29; 50 = M12; 51 = K9
%52 = overturning ratio; %53 = SO total supply; 54 = SO Adv supply; 55 = SO diff supply; 56 = deep
%total supply; 57 = deep adv supply; 58 = deep diff supply; 59 = total
%supply; 60 = pycnocline recipe; 61 = diff/total ratio; 62 = isotope difference;
%63 = mass balance; 64 = isotopic balance;

%==============================================
%INITIAL CONDITIONS and time span of simulation
%To set the initial conditions:
%===============================================
% See Figure 3b in Fripiat et al.(2023) for the architecture of the 9-box model
% N for nitrate concentration (in µmol l-1); 
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

%Array with initial conditions
x0 = [N_1_ini, N_2_ini, N_3_ini, N_4_ini, N_5_ini, N_6_ini, N_7_ini, N_8_ini, N_9_ini, N15_1_ini, N15_2_ini, N15_3_ini, N15_4_ini, N15_5_ini, N15_6_ini, N15_7_ini, N15_8_ini, N15_9_ini]; 
tspan = (0:20000); %Model Spin-Up to reach the equilibrium

%=======================
%Parameters declariation
%=======================
%fixed model parameters (i.e., not changing)
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

%Constant parameters
fracA             = 0.12;%fraction for North Atlantic
fracP             = 0.88;%fraction for Ind-Pac-South Atlantic
fracdeni          = 0.9;%fraction of denitrification occuring in the Ind-Pac-South Atlantic
fracN2fix         = 0.83;%fraction of N2 fixation occuring in the Ind-Pac-South Atlantic
K_LL              = 0.99999;%degree of consumption in the low-latitude surface water
R_5               = 1; 
R_9               = 1; 
eps               = 5.5; %isotope effect of nitrate assimilation
eps_deni_3        = 8.0; % net isotope effect of denitrification in the pycnocline (both water column and sedimentary denitrification)
eps_deni_4        = 0; % net isotope effect of denitrification in the North Atlantic pycnocline (only sedimentary denitrification)  
eps_deni_1        = 0; % net isotope effect of denitrification in the Ind-Pac-South Atlantic deep ocean (only sedimentary denitrification) 
eps_deni_2        = 0; % net isotope effect of denitrification in the North Atlantic deep ocean (only sedimentary denitrification)  
R15ref            = 0.0036782;% isotopic ratio of Air N2 (i.e., the international reference for N isotopes)
d15N_N2           = -1;% d15N of N2 fixation

%create a matrix of zero value with a number of rows being equal to the sensitivity matrix (i.e., number of model runs) and a number of column being equal to the model variables
xend = zeros(length(sensitivity(:,1)),length(x0));

%the following loop run the model for each of the scenario (i.e.,combination of parameters) corresponding to a raw in the "sensitivity" matrix
for i=1:length(sensitivity(:,1))
    Ws             = sensitivity(i,1); 
    Wd             = sensitivity(i,3); 
    W15            = Ws;
    WN             = Ws+Wd;
    X3             = sensitivity(i,6);
    X4             = sensitivity(1,14);
    W13            = fracP.*Wd;
    W24            = fracA.*Wd;  
    W34            = (1-X3)*(W15+W13);
    W37            = X3*(W15+W13);
    W42            = (1-X4)*(W24+W34);
    W48            = X4*(W24+W34);
    W21            = WN-W24;
    M36            = sensitivity(i,2); 
    Md             = sensitivity(i,4); 
    M13            = fracP.*Md;
    M24            = fracA.*Md;
    Mll            = sensitivity(i,5); 
    M37            = fracP.*Mll;
    M48            = fracA.*Mll;
    M15            = sensitivity(i,13);
    M34            = sensitivity(i,15);
    M29            = sensitivity(i,16);
    M12            = sensitivity(i,17);
    K_5            = sensitivity(i,7); %degree of consumption in AZ
    K_6            = sensitivity(i,8); %degree of consumption in PFZ-SAZ
    K_9            = sensitivity(i,18); %degree of consumption in North Atlantic
    R_6_P          = sensitivity(i,9); 
    R_6_D          = 1-R_6_P; 
    R_LL_P         = sensitivity(i,10); 
    R_LL_D         = 1-R_LL_P; 
    x_deni         = sensitivity(i,11); %fraction of denitrification occuring in the pycnocline (i.e., the remaining fraction occurs in the deep ocean)
    deni           = sensitivity(i,12); % total denitrification in µmol yr-1
    deni_3         = fracdeni.*x_deni.*deni; % denitrification in the Ind-Pac-South Atlantic pycnocline
    deni_4         = (1-fracdeni).*x_deni.*deni; % denitrification in the North Atlantic pycnocline
    deni_1         = fracP.*(1-x_deni).*deni; % denitrification in the Ind-Pac-South Atlantic deep ocean 
    deni_2         = fracA.*(1-x_deni).*deni; % denitrification in the North Atlantic deep ocean 
    N2fix_3        = fracN2fix.*x_deni.*deni; % N2 fixation in the Ind-Pac-South Atlantic pycnocline
    N2fix_4        = (1-fracN2fix).*x_deni.*deni; % N2 fixation in the North Atlantic pycnocline
    N2fix_1        = deni_1;% N2 fixation in the Ind-Pac-South Atlantic deep ocean
    N2fix_2        = deni_2;% N2 fixation in the North Atlantic deep ocean
    N15_N2_3       = ((((d15N_N2)/1000)+1)*R15ref).*N2fix_3; %15N input by N2 fixation in µmol yr-1 (in the Ind-Pac-South Atlantic pycnocline)
    N15_N2_4       = ((((d15N_N2)/1000)+1)*R15ref).*N2fix_4; %15N input by N2 fixation in µmol yr-1 (in the North Atlantic pycnocline)
    N15_N2_1       = ((((d15N_N2)/1000)+1)*R15ref).*N2fix_1; %15N input by N2 fixation in µmol yr-1 (in the Ind-Pac-South Atlantic deep ocean)
    N15_N2_2       = ((((d15N_N2)/1000)+1)*R15ref).*N2fix_2; %15N input by N2 fixation in µmol yr-1 (in the North Atlantic deep ocean )
    %=================
    %SOLVING THE ODE for each scenario (i.e., row) of the "sensitivity" matrix
    %=================
    options = odeset('RelTol',1e-7,'AbsTol',1e-7);
    [t,x] = ode15s(@(t,x)ninebox_ocean_isotopes_sensitivity_ode(t,x,V_1,V_2, V_3, V_4,V_5, V_6,V_7, V_8,V_9,fracA,fracP,X3,X4,fracdeni,fracN2fix, Ws,Wd, W15,W24,W21,W13,W34,W37,W42,W48,M36,Md,M13,M24, Mll,M37,M48, M15,M34, M29,M12,K_5, K_6,K_9, K_LL, R_5, R_6_P, R_6_D, R_LL_P, R_LL_D,R_9, x_deni,deni_3,deni_4,deni_1,deni_2,N2fix_3,N2fix_4,N2fix_1,N2fix_2, eps, eps_deni_3, eps_deni_4,eps_deni_1, eps_deni_2, R15ref, d15N_N2, N15_N2_1, N15_N2_2,N15_N2_3, N15_N2_4),tspan,x0,options);
    % to convert 15N concentration in delta values
    d15N=zeros(length(tspan),9);
    d15N(:,1)      = (((x(:,10)./x(:,1))./R15ref)-1)*1000;
    d15N(:,2)      = (((x(:,11)./x(:,2))./R15ref)-1)*1000;
    d15N(:,3)      = (((x(:,12)./x(:,3))./R15ref)-1)*1000;
    d15N(:,4)      = (((x(:,13)./x(:,4))./R15ref)-1)*1000;
    d15N(:,5)      = (((x(:,14)./x(:,5))./R15ref)-1)*1000;
    d15N(:,6)      = (((x(:,15)./x(:,6))./R15ref)-1)*1000;
    d15N(:,7)      = (((x(:,16)./x(:,7))./R15ref)-1)*1000;
    d15N(:,8)      = (((x(:,17)./x(:,8))./R15ref)-1)*1000;
    d15N(:,9)      = (((x(:,18)./x(:,9))./R15ref)-1)*1000;
    % to calculate the value for the global deep ocean, pycnocline and LLML
    lump=zeros(length(tspan),6);
    lump(:,1)      = (x(:,1).*V_1+x(:,2).*V_2)./(V_1+V_2);
    lump(:,2)      = (x(:,3).*V_3+x(:,4).*V_4)./(V_3+V_4);
    lump(:,3)      = (x(:,7).*V_7+x(:,8).*V_8)./(V_7+V_8);
    lump(:,4)      = (x(:,1).*d15N(:,1).*V_1+x(:,2).*d15N(:,2).*V_2)./(x(:,1).*V_1+x(:,2).*V_2);
    lump(:,5)      = (x(:,3).*d15N(:,3).*V_3+x(:,4).*d15N(:,4).*V_4)./(x(:,3).*V_3+x(:,4).*V_4);
    lump(:,6)      = (x(:,7).*d15N(:,7).*V_7+x(:,8).*d15N(:,8).*V_8)./(x(:,7).*V_7+x(:,8).*V_8);
    % For each scenario, to take steady-state solution (i.e., last row in x) and insert it in the xend matrix 
    xend(i,1:18)=x(length(tspan),1:18);
    xend(i,19:27)=d15N(length(tspan),1:9);
    xend(i,28:33)=lump(length(tspan),1:6);
end

% to create the final matrix which includes model variables/parameters/offline calcualtion
A = [xend,sensitivity]; % merging the model output (only equilibrium solution) with sensitivity matrix (i.e., model parameters)
% the following corresponds to offline calculation 
A(:,52) = sensitivity(:,1)./(sensitivity(:,1)+sensitivity(:,3)); %overturning ratio
A(:,53) = (sensitivity(:,1)+sensitivity(:,2)).*A(:,6).*14./(10^6*10^12); % SO total nitrate supply into the pycnocline (Tg N yr-1)
A(:,54) = (sensitivity(:,1)).*A(:,6).*14./(10^6*10^12); % SO advective nitrate supply into the pycnocline (Tg N yr-1)
A(:,55) = (sensitivity(:,2)).*A(:,6).*14./(10^6*10^12); % SO diffusive nitrate supply into the pycnocline (Tg N yr-1)
A(:,56) = (sensitivity(:,3)+sensitivity(:,4)).*fracP.*A(:,1).*14./(10^6*10^12)+(sensitivity(:,3)+sensitivity(:,4)).*fracA.*A(:,2).*14./(10^6*10^12); % deep total nitrate supply into the pycnocline (Tg N yr-1)
A(:,57) = sensitivity(:,3).*fracP.*A(:,1).*14./(10^6*10^12)+sensitivity(:,3).*fracA.*A(:,2).*14./(10^6*10^12); % deep advective nitrate supply into the pycnocline (Tg N yr-1)
A(:,58) = sensitivity(:,4).*fracP.*A(:,1).*14./(10^6*10^12)+sensitivity(:,4).*fracA.*A(:,2).*14./(10^6*10^12); % deep diffusive nitrate supply into the pycnocline (Tg N yr-1)
A(:,59) = A(:,53)+A(:,56); % total nitrate supply into the pycnocline (Tg N yr-1)
A(:,60) = A(:,53)./A(:,59); %pycnocline recipe 
A(:,61) = (A(:,55)+A(:,58))./A(:,59); %diffusive/total nitrate supply ratio into the pycnocline
A(:,62) = A(:,32)-A(:,31);% difference between deep ocean and pycnocline d15N
A(:,63) = (A(:,1).*V_1+A(:,2).*V_2+A(:,3).*V_3+A(:,4).*V_4+A(:,5).*5+A(:,6).*6+A(:,7).*V_7+A(:,8).*V_8+A(:,9).*V_9)/V_tot;% mass balance
A(:,64) = (A(:,1).*A(:,19).*V_1+A(:,2).*A(:,20).*V_2+A(:,3).*A(:,21).*V_3+A(:,4).*A(:,22).*V_4...
    +A(:,5).*A(:,23).*V_5+A(:,6).*A(:,24).*V_6+A(:,7).*A(:,25).*V_7+A(:,8).*A(:,26).*V_8+A(:,9).*A(:,27).*V_9)./(A(:,1).*V_1+A(:,2).*V_2+A(:,3).*V_3+A(:,4).*V_4+A(:,5).*V_5+A(:,6).*V_6+A(:,7).*V_7+A(:,8).*V_8+A(:,9).*V_9); % isotopic balance

%to clean up the workspace
clearvars -except keep A sensitivity

%%%%%%%%%%%%%% plot results
figure(1);
hold on; plot(A(:,60),A(:,28),'ko','MarkerFaceColor','k')
hold on; plot(A(:,60),A(:,29),'ro','MarkerFaceColor','r')
figure(2);
hold on; plot(A(:,60),A(:,31),'ko','MarkerFaceColor','k')
hold on; plot(A(:,60),A(:,32),'ro','MarkerFaceColor','r')
