size=100;%number of scenario (10000 should take 1 hour to run)
%see figure 3a in Fripiat et al. (2023) for a description of the parameters
%To create a "sensivity" matrix in which each row corresponds to a combination of parameters 

a = 0.1*10^3*10^6*3600*24*365;
b = 100*10^3*10^6*3600*24*365;
r = (b-a).*rand(size,1) + a;
Ws = r; 

a = 0.1*10^3*10^6*3600*24*365;
b = 100*10^3*10^6*3600*24*365;
r = (b-a).*rand(size,1) + a;
Ms = r; 

a = 0.1*10^3*10^6*3600*24*365;
b = 100*10^3*10^6*3600*24*365;
r = (b-a).*rand(size,1) + a;
Wd = r; 

a = 0.1*10^3*10^6*3600*24*365;
b = 100*10^3*10^6*3600*24*365;
r = (b-a).*rand(size,1) + a;
Md = r; 

a = 0.1*10^3*10^6*3600*24*365;
b = 50*10^3*10^6*3600*24*365;
r = (b-a).*rand(size,1) + a;
Mll = r; 

a = 0.5;
b = 1;
r = (b-a).*rand(size,1) + a;
x_Wn = r; 

a = 0.23;
b = 0.23;
r = (b-a).*rand(size,1) + a;
K_AZ = r; %degree of consumption in AZ

a = 0.35;
b = 0.35;
r = (b-a).*rand(size,1) + a;
K_SAZ = r; %degree of consumption in PFZ-SAZ

a = 0.0;
b = 0.5;
r = (b-a).*rand(size,1) + a;
R_SAZ = r; 

a = 0.80;
b = 0.90;
r = (b-a).*rand(size,1) + a;
R_LL_P = r; 

a = 0.9999;
b = 0.9999;
r = (b-a).*rand(size,1) + a;
x_deni = r; %fraction of denitrification occuring in the pycnocline (i.e., the remaining fraction occurs in the deep ocean)

a = 100/14*10^6*10^12;
b = 250/14*10^6*10^12;
r = (b-a).*rand(size,1) + a;
deni = r; % total denitrification in µmol yr-1

a = 20*10^3*10^6*3600*24*365;
b = 20*10^3*10^6*3600*24*365;
r = (b-a).*rand(size,1) + a;
MAZ = r; 

%to generate the sensitivy matrix (each row = one model scenario)
sensitivity = zeros(size,13);
sensitivity(:,1)   = Ws;
sensitivity(:,2)   = Ms;%Ms
sensitivity(:,3)   = Wd;
sensitivity(:,4)   = Md;%Md
sensitivity(:,5)   = Mll;
sensitivity(:,6)   = x_Wn;
sensitivity(:,7)   = K_AZ;
sensitivity(:,8)   = K_SAZ;%K_SAZ
sensitivity(:,9)   = R_SAZ;
sensitivity(:,10)  = R_LL_P;
sensitivity(:,11)  = x_deni;
sensitivity(:,12)  = deni;
sensitivity(:,13)  = MAZ;
clearvars -except keep sensitivity
