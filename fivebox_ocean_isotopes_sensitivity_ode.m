function dx = fivebox_ocean_isotopes_sensitivity_ode(t,x,V_deep,V_AZ, V_SAZ, V_pycn, V_LLML, Ws, Wd, Ms, Md, Mll, Wnll, Wnp, MAZ, K_AZ, K_SAZ, K_LLML, R_AZ, R_SAZ_P, R_SAZ_D, R_LL_P, R_LL_D, deni_P, deni_D, eps, eps_deni_D, eps_deni_P, N15_N2_P, N15_N2_D)

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
end
