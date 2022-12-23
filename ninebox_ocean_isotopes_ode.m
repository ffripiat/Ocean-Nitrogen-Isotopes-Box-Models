function dx = ninebox_ocean_isotopes_ode(t,x,V_1,V_2, V_3, V_4,V_5, V_6,V_7, V_8,V_9,fracA,fracP,X3,X4,fracdeni3,fracN2fix3, Ws,Wd, W15,W24,W21,W13,W34,W37,W42,W48,M36,Md,M13,M24, Mll,M37,M48, M15,M34, M29,M12,K_5, K_6,K_9, K_LL, R_5, R_6_P, R_6_D, R_LL_P, R_LL_D,R_9, deni_3,deni_4,deni_1,deni_2,N2fix_3,N2fix_4, N2fix_1,N2fix_2,eps, eps_deni_3, eps_deni_4, eps_deni_1, eps_deni_2, R15ref, d15N_N2, N15_N2_3, N15_N2_4, N15_N2_1, N15_N2_2)

%=======================
%Parameters declaration
%=======================
%see figure 3b in Fripiat et al. (2023) for a description of the parameters
%box volume in l
V_1               = 7.9200e+20;
V_2               = 1.0800e+20;
V_3               = 2.4376e+20;
V_4               = 3.3240e+19;
V_5               = 6.1e18;
V_6               = 9.76e18;
V_7               = 5.1920e+19;
V_8               = 7.0800e+18;
V_9               = 3e18;

%parameters which are related to the nine-box model
fracA             = 0.12;%fraction for North Atlantic
fracP             = 0.88;%fraction for Ind-Pac-South Atlantic
X3                = 0.82;%fraction of outflow through LLML in the Ind-Pac-South Atlantic
X4                = 0.05;%fraction of outflow through LLML in the North Atlantic 
fracdeni          = 0.9; %fraction of denitrification occuring in the Ind-Pac-South Atlantic
fracN2fix         = 0.83; %fraction of N2 fixation occuring in the Ind-Pac-South Atlantic


% flow in Sv (10^6 m^3 s^-1) converted in l yr-1 using *10^6*10^3*365*60*60*24
Ws                   = 31*10^6*10^3*365*60*60*24; 
Wd                   = 6*10^6*10^3*365*60*60*24;
W15                  = Ws;
WN                   = Ws+Wd;
W13                  = fracP.*Wd;
W24                  = fracA.*Wd;
W34                  = (1-X3)*(W15+W13);
W37                  = X3*(W15+W13);
W42                  = (1-X4)*(W24+W34);
W48                  = X4*(W24+W34);
W21                  = WN-W24;
M36                  = 52*10^6*10^3*365*60*60*24;
Md                   = 24*10^6*10^3*365*60*60*24; 
M13                  = fracP.*Md;
M24                  = fracA.*Md;
Mll                  = 24*10^6*10^3*365*60*60*24; 
M37                  = fracP.*Mll;
M48                  = fracA.*Mll;
M15                  = 20*10^6*10^3*365*60*60*24; 
M34                  = 20*10^6*10^3*365*60*60*24;
M29                  = 20*10^6*10^3*365*60*60*24;
M12                  = 20*10^6*10^3*365*60*60*24;

%degree of consumption in AZ, SAZ Low-Latitude mixed layer, and North Atlantic
K_5    = 0.16;
K_6    = 0.32;
K_LL   = 0.9999;
K_9    = 0.6;

%degree of remineralization
R_5    = 1; 
R_9    = 1; 
R_6_P  = 0.22;
R_6_D  = 0.78;
R_LL_P = 0.85;
R_LL_D = 0.15;

eps            = 5; %isotope effect of nitrate assimilation
x_deni         = 0.9999; %fraction of denitrification occuring in the pycnocline (i.e., the remaining fraction occurs in the deep ocean)
deni           = (150/14*10^6*10^12);%µmol yr-1 % denitrification (= N2 fixation) in the pycnocline 
deni_3         = fracdeni.*x_deni.*deni; % denitrification in the Ind-Pac-South Atlantic pycnocline
deni_4         = (1-fracdeni).*x_deni.*deni; % denitrification in the North Atlantic pycnocline
deni_1         = fracP.*(1-x_deni).*deni; % denitrification in the Ind-Pac-South Atlantic deep ocean
deni_2         = fracA.*(1-x_deni).*deni; % denitrification in the North Atlantic deep ocean
N2fix_3        = fracN2fix.*x_deni.*deni; % N2 fixation in the Ind-Pac-South Atlantic pycnocline
N2fix_4        = (1-fracN2fix).*x_deni.*deni; % N2 fixation in the North Atlantic pycnocline
N2fix_1        = deni_1; % N2 fixation in the Ind-Pac-South Atlantic deep ocean
N2fix_2        = deni_2; % N2 fixation in the North Atlantic deep ocean
eps_deni_3     = 8; % net isotope effect of denitrification in the Ind-Pac-South Atlantic pycnocline (both water column and sedimentary denitrification 
eps_deni_4     = 0; % net isotope effect of denitrification in the North Atlantic pycnocline (only sedimentary denitrification) 
eps_deni_1     = 0; % net isotope effect of denitrification in the Ind-Pac-South Atlantic deep ocean (only sedimentary denitrification) 
eps_deni_2     = 0; % net isotope effect of denitrification in the North Atlantic deep ocean (only sedimentary denitrification) 
R15ref         = 0.0036782; % isotopic ratio of Air N2 (i.e., the international reference for N isotopes)
d15N_N2        = -1; % d15N of N2 fixation
N15_N2_3       = ((((d15N_N2)/1000)+1)*R15ref)*N2fix_3; %15N input by N2 fixation in µmol yr-1 (in the Ind-Pac-South Atlantic pycnocline)
N15_N2_4       = ((((d15N_N2)/1000)+1)*R15ref)*N2fix_4; %15N input by N2 fixation in µmol yr-1 (in the North Atlantic pycnocline)
N15_N2_1       = ((((d15N_N2)/1000)+1)*R15ref).*N2fix_1; %15N input by N2 fixation in µmol yr-1 (in the Ind-Pac-South Atlantic deep ocean) 
N15_N2_2       = ((((d15N_N2)/1000)+1)*R15ref).*N2fix_2; %15N input by N2 fixation in µmol yr-1 (in the North Atlantic deep ocean )

%================
% Equation being used by the ODE
%================
dx=zeros(18,1);
% 1 = Deep Indo-Pacific-South Atlantic; 2 = Deep North Atlantic; 3 = pycnocline Indo-Pacific-South Atlantic; 
% 4 = pycnocline North Atlantic; 5 = AZ surface; 6 = SAZ surface; 7 = LLML Indo-Pacific-South Atlantic; 
% 8 = LLML North Atlantic; 9 = North Atl surface

%N concentration
dx(1) = W21/V_1*x(2)-W13/V_1*x(1)-W15/V_1*x(1)+M13/V_1*(x(3)-x(1))+M15/V_1*(x(5)-x(1))+M12/V_1*(x(2)-x(1))+(W15+M15)/V_1*x(1)*K_5*R_5+W15/V_1*x(5)*K_6*R_6_D+M36/V_1*x(3)*K_6*R_6_D+(W37+M37)/V_1*x(3)*K_LL*R_LL_D+N2fix_1/V_1-deni_1/V_1;
dx(2) = W42/V_2*x(4)+(W48+W37)/V_2*x(9)-W21/V_2*x(2)-W24/V_2*x(2)+M24/V_2*(x(4)-x(2))+M29/V_2*(x(9)-x(2))+M12/V_2*(x(1)-x(2))+(M48+W48)/V_2*x(4)*K_LL*R_LL_D+W37/V_2*x(7)*K_LL*R_LL_D+M29/V_2*x(2)*K_9*R_9+(W48+W37)/V_2*x(8)*K_9*R_9+N2fix_2/V_2-deni_2/V_2;
dx(3) = W15/V_3*x(6)+W13/V_3*x(1)-W34/V_3*x(3)-W37/V_3*x(3)+M13/V_3*(x(1)-x(3))+M36/V_3*(x(6)-x(3))+M34/V_3*(x(4)-x(3))+M37/V_3*(x(7)-x(3))+W15/V_3*x(5)*K_6*R_6_P+M36/V_3*x(3)*K_6*R_6_P+(W37+M37)/V_3*x(3)*K_LL*R_LL_P+N2fix_3/V_3-deni_3/V_3;
dx(4) = W24/V_4*x(2)+W34/V_4*x(3)-W42/V_4*x(4)-W48/V_4*x(4)+M34/V_4*(x(3)-x(4))+M24/V_4*(x(2)-x(4))+M48/V_4*(x(8)-x(4))+(W48+M48)/V_4*x(4)*K_LL*R_LL_P+W37/V_4*x(7)*K_LL*R_LL_P+N2fix_4/V_4-deni_4/V_4;
dx(5) = W15/V_5*(x(1)-x(5))+M15/V_5*(x(1)-x(5))-(W15+M15)/V_5*x(1)*K_5;
dx(6) = W15/V_6*(x(5)-x(6))+M36/V_6*(x(3)-x(6))-W15/V_6*x(5)*K_6-M36/V_6*x(3)*K_6;
dx(7) = W37/V_7*(x(3)-x(7))+M37/V_7*(x(3)-x(7))-(W37+M37)/V_7*x(3)*K_LL;
dx(8) = W48/V_8*x(4)+W37/V_8*x(7)-W48/V_8*x(8)-W37/V_8*x(8)+M48/V_8*(x(4)-x(8))-(W48+M48)/V_8*x(4)*K_LL-W37/V_8*x(7)*K_LL;
dx(9) = (W48+W37)/V_9*(x(8)-x(9))+M29/V_9*(x(2)-x(9))-(W48+W37)/V_9*x(8)*K_9-M29/V_9*x(2)*K_9;

% 15N concentration 
dx(10) = W21/V_1*x(2)*x(11)/x(2)-W13/V_1*x(1)*x(10)/x(1)-W15/V_1*x(1)*x(10)/x(1)+M13/V_1*(x(3)*x(12)/x(3)-x(1)*x(10)/x(1))+M15/V_1*(x(5)*x(14)/x(5)-x(1)*x(10)/x(1))+M12/V_1*(x(2)*x(11)/x(2)-x(1)*x(10)/x(1))+(W15+M15)/V_1*x(1)*x(10)/x(1)*K_5*R_5*(1-(1-K_5)^(1-eps/1000))/K_5+W15/V_1*x(5)*x(14)/x(5)*K_6*R_6_D*(1-(1-K_6)^(1-(eps/4)/1000))/K_6+M36/V_1*x(3)*x(12)/x(3)*K_6*R_6_D*(1-(1-K_6)^(1-(eps/4)/1000))/K_6+(W37+M37)/V_1*x(3)*x(12)/x(3)*K_LL*R_LL_D*(1-(1-K_LL)^(1-eps/1000))/K_LL+N2fix_1/V_1*N15_N2_1/N2fix_1-deni_1/V_1*x(10)/x(1)*(1-eps_deni_1/1000);
dx(11) = W42/V_2*x(4)*x(13)/x(4)+(W48+W37)/V_2*x(9)*x(18)/x(9)-W21/V_2*x(2)*x(11)/x(2)-W24/V_2*x(2)*x(11)/x(2)+M24/V_2*(x(4)*x(13)/x(4)-x(2)*x(11)/x(2))+M29/V_2*(x(9)*x(18)/x(9)-x(2)*x(11)/x(2))+M12/V_2*(x(1)*x(10)/x(1)-x(2)*x(11)/x(2))+(M48+W48)/V_2*x(4)*x(13)/x(4)*K_LL*R_LL_D*(1-(1-K_LL)^(1-eps/1000))/K_LL+W37/V_2*x(7)*x(16)/x(7)*K_LL*R_LL_D*(1-(1-K_LL)^(1-eps/1000))/K_LL+M29/V_2*x(2)*x(11)/x(2)*K_9*R_9*(1-(1-K_9)^(1-eps/1000))/K_9+(W48+W37)/V_2*x(8)*x(17)/x(8)*K_9*R_9*(1-(1-K_9)^(1-eps/1000))/K_9+N2fix_2/V_2*N15_N2_2/N2fix_2-deni_2/V_2*x(11)/x(2)*(1-eps_deni_2/1000);
dx(12) = W15/V_3*x(6)*x(15)/x(6)+W13/V_3*x(1)*x(10)/x(1)-W34/V_3*x(3)*x(12)/x(3)-W37/V_3*x(3)*x(12)/x(3)+M13/V_3*(x(1)*x(10)/x(1)-x(3)*x(12)/x(3))+M36/V_3*(x(6)*x(15)/x(6)-x(3)*x(12)/x(3))+M34/V_3*(x(4)*x(13)/x(4)-x(3)*x(12)/x(3))+M37/V_3*(x(7)*x(16)/x(7)-x(3)*x(12)/x(3))+W15/V_3*x(5)*x(14)/x(5)*K_6*R_6_P*(1-(1-K_6)^(1-(eps/4)/1000))/K_6+M36/V_3*x(3)*x(12)/x(3)*K_6*R_6_P*(1-(1-K_6)^(1-(eps/4)/1000))/K_6+(W37+M37)/V_3*x(3)*x(12)/x(3)*K_LL*R_LL_P*(1-(1-K_LL)^(1-eps/1000))/K_LL+N2fix_3/V_3*N15_N2_3/N2fix_3-deni_3/V_3*x(12)/x(3)*(1-eps_deni_3/1000);
dx(13) = W24/V_4*x(2)*x(11)/x(2)+W34/V_4*x(3)*x(12)/x(3)-W42/V_4*x(4)*x(13)/x(4)-W48/V_4*x(4)*x(13)/x(4)+M34/V_4*(x(3)*x(12)/x(3)-x(4)*x(13)/x(4))+M24/V_4*(x(2)*x(11)/x(2)-x(4)*x(13)/x(4))+M48/V_4*(x(8)*x(17)/x(8)-x(4)*x(13)/x(4))+(W48+M48)/V_4*x(4)*x(13)/x(4)*K_LL*R_LL_P*(1-(1-K_LL)^(1-eps/1000))/K_LL+W37/V_4*x(7)*x(16)/x(7)*K_LL*R_LL_P*(1-(1-K_LL)^(1-eps/1000))/K_LL+N2fix_4/V_4*N15_N2_4/N2fix_4-deni_4/V_4*x(13)/x(4)*(1-eps_deni_4/1000);
dx(14) = W15/V_5*(x(1)*x(10)/x(1)-x(5)*x(14)/x(5))+M15/V_5*(x(1)*x(10)/x(1)-x(5)*x(14)/x(5))-(W15+M15)/V_5*x(1)*x(10)/x(1)*K_5*(1-(1-K_5)^(1-eps/1000))/K_5;
dx(15) = W15/V_6*(x(5)*x(14)/x(5)-x(6)*x(15)/x(6))+M36/V_6*(x(3)*x(12)/x(3)-x(6)*x(15)/x(6))-W15/V_6*x(5)*x(14)/x(5)*K_6*(1-(1-K_6)^(1-(eps/4)/1000))/K_6-M36/V_6*x(3)*x(12)/x(3)*K_6*(1-(1-K_6)^(1-(eps/4)/1000))/K_6;
dx(16) = W37/V_7*(x(3)*x(12)/x(3)-x(7)*x(16)/x(7))+M37/V_7*(x(3)*x(12)/x(3)-x(7)*x(16)/x(7))-(W37+M37)/V_7*x(3)*x(12)/x(3)*K_LL*(1-(1-K_LL)^(1-eps/1000))/K_LL;
dx(17) = W48/V_8*x(4)*x(13)/x(4)+W37/V_8*x(7)*x(16)/x(7)-W48/V_8*x(8)*x(17)/x(8)-W37/V_8*x(8)*x(17)/x(8)+M48/V_8*(x(4)*x(13)/x(4)-x(8)*x(17)/x(8))-(W48+M48)/V_8*x(4)*x(13)/x(4)*K_LL*(1-(1-K_LL)^(1-eps/1000))/K_LL-W37/V_8*x(7)*x(16)/x(7)*K_LL*(1-(1-K_LL)^(1-eps/1000))/K_LL;
dx(18) = (W48+W37)/V_9*(x(8)*x(17)/x(8)-x(9)*x(18)/x(9))+M29/V_9*(x(2)*x(11)/x(2)-x(9)*x(18)/x(9))-(W48+W37)/V_9*x(8)*x(17)/x(8)*K_9*(1-(1-K_9)^(1-eps/1000))/K_9-M29/V_9*x(2)*x(11)/x(2)*K_9*(1-(1-K_9)^(1-eps/1000))/K_9;

