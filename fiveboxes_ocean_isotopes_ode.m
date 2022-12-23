function dx = fiveboxes_ocean_isotopes_imbalance_ode(t,x,V_deep,V_AZ, V_SAZ, V_pycn, V_LLML, Ws, Wd, Ms, Md, Mll, Wnll, Wnp, MAZ, K_AZ, K_SAZ, K_LLML, R_AZ, R_SAZ_P, R_SAZ_D, R_LL_P, R_LL_D, x_deni, deni_P, deni_D, eps, eps_deni_D, eps_deni_P, R15ref, d15N_N2, N15_N2_P, N15_N2_D)

%=======================
%Parameters declaration
%=======================
%see figure 3a in Fripiat et al. (2023) for a descriptin of the parameters
%box volume in l
V_deep               = 9.0000e20; 
V_AZ                 = 6.1e18;
V_SAZ                = 9.76e18;
V_pycn               = 2.7700e+20;
V_LLML               = 5.9000e+19;

% flow in Sv (10^6 m^3 s^-1) converted in l yr-1 using *10^6*10^3*365*60*60*24
Ws                   = 31*10^6*10^3*365*60*60*24; %31 % Southern ocean upwelling (i.e., upper cell)
Wd                   = 6*10^6*10^3*365*60*60*24;%6 % Upwelling from the deep ocean to the pycnocline
Ms                   = 52*10^6*10^3*365*60*60*24;% 52 % isopycnal mixing between PFZ-SAZ and pycnocline
Md                   = 24*10^6*10^3*365*60*60*24; %24 % diapycnal mixing between deep ocean and pycnocline
Mll                  = 24*10^6*10^3*365*60*60*24; % diapycnal mixing between pycnocline and low-latitude surface water
Wnll                 = 30*10^6*10^3*365*60*60*24; % Advective output from the pycnocline passing through the low-latitude surface water
Wnp                  = 7*10^6*10^3*365*60*60*24; % Advective output from the pycnocline directly from the pycnocline
MAZ                  = 20*10^6*10^3*365*60*60*24; % Vertical mixing in the AZ (i.e., mimicking the lower overturning cell + mixing)

%degree of consumption in AZ, SAZ and Low-Latitude mixed layer
K_AZ    = 0.12;
K_SAZ   = 0.32;
K_LLML  = 0.9999;

%degree of remineralization
R_AZ    = 1; 
R_SAZ_P = 0.22;
R_SAZ_D = 0.78;
R_LL_P = 0.85;
R_LL_D = 0.15;

eps     = 5.5; %Isotope effect of nitrate assimilation
x_deni  = 0.99999; %fraction of denitrification occuring in the pycnocline (i.e., the remaining fraction occurs in the deep ocean)
deni_P    = (150/14*10^6*10^12)*x_deni; %µmol yr-1 % denitrification (= N2 fixation) in the pycnocline %% the first number is in Tg N yr-1 (being converted in µmol yr-1 with /14*10^6*10^12)
deni_D    = (150/14*10^6*10^12)*(1-x_deni); %µmol yr-1 % denitrification (= N2 fixation) in the deep ocean %% the first number is in Tg N yr-1 (being converted in µmol yr-1 with /14*10^6*10^12)
eps_deni_P= 7.25; % net isotope effect of denitrification in the pycnocline (both water column and sedimentary denitrification 
eps_deni_D= 0; % net isotope effect of denitrification in the deep ocean (only sedimentary denitrification) 
R15ref = 0.0036782; % isotopic ratio of Air N2 (i.e., the international reference for N isotopes)
d15N_N2 = -1; % d15N of N2 fixation
N15_N2_P  = ((((d15N_N2)/1000)+1)*R15ref)*deni_P; %15N input by N2 fixation in µmol yr-1 (in the pycnocline)
N15_N2_D  = ((((d15N_N2)/1000)+1)*R15ref)*deni_D; %15N input by N2 fixation in µmol yr-1 (in the deep ocean)

%================
% Equation being used by the ODE
%================
dx=zeros(10,1);
% nitrate concentration with 1 for the deep ocean, 2 for the AZ, 3 for the
% SAZ, 4 for the pycnocline, and 5 for the low-latitude mixed layer
dx(1)= MAZ/V_deep*(x(2)-x(1))+Wnll/V_deep*x(5)+Wnp/V_deep*x(4)-Ws/V_deep*x(1)-Wd/V_deep*x(1)+Md/V_deep*(x(4)-x(1))+(Ws+MAZ)/V_deep*x(1)*K_AZ*R_AZ+Ws/V_deep*x(2)*K_SAZ*R_SAZ_D+Ms/V_deep*x(4)*K_SAZ*R_SAZ_D+(Wnll+Mll)/V_deep*x(4)*K_LLML*R_LL_D+deni_D/V_deep-deni_D/V_deep;
dx(2)= MAZ/V_AZ*(x(1)-x(2))+Ws/V_AZ*(x(1)-x(2))-(Ws+MAZ)/V_AZ*x(1)*K_AZ;
dx(3)= Ws/V_SAZ*(x(2)-x(3))+Ms/V_SAZ*(x(4)-x(3))-Ws/V_SAZ*x(2)*K_SAZ-Ms/V_SAZ*x(4)*K_SAZ;
dx(4)= Ws/V_pycn*x(3)+Ms/V_pycn*(x(3)-x(4))+Wd/V_pycn*x(1)+Md/V_pycn*(x(1)-x(4))-Wnp/V_pycn*x(4)-Wnll/V_pycn*x(4)+Mll/V_pycn*(x(5)-x(4))+(Wnll+Mll)/V_pycn*x(4)*K_LLML*R_LL_P+Ws/V_pycn*x(2)*K_SAZ*R_SAZ_P+Ms/V_pycn*x(4)*K_SAZ*R_SAZ_P+deni_P/V_pycn-deni_P/V_pycn;
dx(5)= Wnll/V_LLML*(x(4)-x(5))+Mll/V_LLML*(x(4)-x(5))-(Wnll+Mll)/V_LLML*x(4)*K_LLML;
% 15N concentration in nitrate with 6 for the deep ocean, 7 for the AZ, 8 for the
% SAZ, 9 for the pycnocline, and 10 for the low-latitude mixed layer
dx(6)= MAZ/V_deep*(x(2)*x(7)/x(2)-x(1)*x(6)/x(1))+Wnll/V_deep*x(5)*x(10)/x(5)+Wnp/V_deep*x(4)*x(9)/x(4)-Ws/V_deep*x(1)*x(6)/x(1)-Wd/V_deep*x(1)*x(6)/x(1)+Md/V_deep*(x(4)*x(9)/x(4)-x(1)*x(6)/x(1))+(Ws+MAZ)/V_deep*x(1)*x(6)/x(1)*K_AZ*R_AZ*(1-(1-K_AZ)^(1-eps/1000))/K_AZ+Ws/V_deep*x(2)*x(7)/x(2)*K_SAZ*R_SAZ_D*(1-(1-K_SAZ)^(1-(eps/4)/1000))/K_SAZ+Ms/V_deep*x(4)*x(9)/x(4)*K_SAZ*R_SAZ_D*(1-(1-K_SAZ)^(1-(eps/4)/1000))/K_SAZ+(Wnll+Mll)/V_deep*x(4)*x(9)/x(4)*K_LLML*R_LL_D*(1-(1-K_LLML)^(1-eps/1000))/K_LLML+deni_D/V_deep*N15_N2_D/deni_D-deni_D/V_deep*x(6)./x(1)*(1-eps_deni_D/1000);
dx(7)= MAZ/V_AZ*(x(1)*x(6)/x(1)-x(2)*x(7)/x(2))+Ws/V_AZ*(x(1)*x(6)/x(1)-x(2)*x(7)/x(2))-(Ws+MAZ)/V_AZ*x(1)*x(6)/x(1)*K_AZ*(1-(1-K_AZ)^(1-eps/1000))/K_AZ;
dx(8)= Ws/V_SAZ*(x(2)*x(7)/x(2)-x(3)*x(8)/x(3))+Ms/V_SAZ*(x(4)*x(9)/x(4)-x(3)*x(8)/x(3))-Ws/V_SAZ*x(2)*x(7)/x(2)*K_SAZ*(1-(1-K_SAZ)^(1-(eps/4)/1000))/K_SAZ-Ms/V_SAZ*x(4)*x(9)/x(4)*K_SAZ*(1-(1-K_SAZ)^(1-(eps/4)/1000))/K_SAZ;
dx(9)= Ws/V_pycn*x(3)*x(8)/x(3)+Ms/V_pycn*(x(3)*x(8)/x(3)-x(4)*x(9)/x(4))+Wd/V_pycn*x(1)*x(6)/x(1)+Md/V_pycn*(x(1)*x(6)/x(1)-x(4)*x(9)/x(4))-Wnp/V_pycn*x(4)*x(9)/x(4)-Wnll/V_pycn*x(4)*x(9)/x(4)+Mll/V_pycn*(x(5)*x(10)/x(5)-x(4)*x(9)/x(4))+(Wnll+Mll)/V_pycn*x(4)*x(9)/x(4)*K_LLML*R_LL_P*(1-(1-K_LLML)^(1-eps/1000))/K_LLML+Ws/V_pycn*x(2)*x(7)/x(2)*K_SAZ*R_SAZ_P*(1-(1-K_SAZ)^(1-(eps/4)/1000))/K_SAZ+Ms/V_pycn*x(4)*x(9)/x(4)*K_SAZ*R_SAZ_P*(1-(1-K_SAZ)^(1-(eps/4)/1000))/K_SAZ+deni_P/V_pycn*N15_N2_P/deni_P-deni_P/V_pycn*x(9)./x(4)*(1-eps_deni_P/1000);
dx(10)= Wnll/V_LLML*(x(4)*x(9)/x(4)-x(5)*x(10)/x(5))+Mll/V_LLML*(x(4)*x(9)/x(4)-x(5)*x(10)/x(5))-(Wnll+Mll)/V_LLML*x(4)*x(9)/x(4)*K_LLML*(1-(1-K_LLML)^(1-eps/1000))/K_LLML;

