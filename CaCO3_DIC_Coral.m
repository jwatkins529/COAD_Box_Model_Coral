function [alpha_c_spec,alpha_o_spec,alpha63_spec,alpha64_spec,alpha65_spec,F]=CaCO3_DIC(TC,Ksp,pH, A, ~, ~, B1)
%%
%=========================Input conditions=================================
TK = TC+273.15;
omega = A*B1/Ksp;
%%
%==========================Rp-omega relationship===========================
krate=exp(11.54-8690/TK);
F=krate*(omega-1)^(1.7);
Salinity=35;
A = A;
logOmega = 0.0005:0.05:8;
Omega = exp(logOmega);
%%
%==========================Initializing====================================
S=zeros(length(Omega),1);
B_r=zeros(length(Omega),1);
B0_r=zeros(length(Omega),1);
B1_r=zeros(length(Omega),1);
B2_r=zeros(length(Omega),1);
DIC=zeros(length(Omega),1);
PA=zeros(length(Omega),1);
PB1=zeros(length(Omega),1);
PB2=zeros(length(Omega),1);
uA=zeros(length(Omega),1);
uB=zeros(length(Omega),1);
uC=zeros(length(Omega),1);
iA=zeros(length(Omega),1);
iB=zeros(length(Omega),1);
iC=zeros(length(Omega),1);
rho_c=zeros(length(Omega),1);
y_o=zeros(length(Omega),1);
R_c=zeros(length(Omega),1);
r_c=zeros(length(Omega),1);
delta_B0 = zeros(length(Omega),1);
delta_B1 = zeros(length(Omega),1);
delta_B2 = zeros(length(Omega),1);
Delta_bw = zeros(length(Omega),1);
Delta_cw = zeros(length(Omega),1);
Delta_cb = zeros(length(Omega),1);
Delta_bb = zeros(length(Omega),1);
r_B1 = zeros(length(Omega),1);
r_B2 = zeros(length(Omega),1);
r_o = zeros(length(Omega),1);
AttRate = zeros(length(Omega),1);
DetRate = zeros(length(Omega),1);
term1 = zeros(length(Omega),1);
term2 = zeros(length(Omega),1);
RR63 = zeros(length(Omega),1);
D63_EIC = zeros(length(Omega),1);
alpha63 = zeros(length(Omega),1);
Delta63 = zeros(length(Omega),1);
Delta47 = zeros(length(Omega),1);
RR64 = zeros(length(Omega),1);
D64_EIC = zeros(length(Omega),1);
Delta64 = zeros(length(Omega),1);
Delta48 = zeros(length(Omega),1);
alpha64 = zeros(length(Omega),1);
RR65 = zeros(length(Omega),1);
D65_EIC = zeros(length(Omega),1);
Delta65 = zeros(length(Omega),1);
Delta49 = zeros(length(Omega),1);
alpha65 = zeros(length(Omega),1);
D=zeros(length(Omega),1);
%===========================New Arrays=====================================
d2 = zeros(length(Omega),1);
B3 = zeros(length(Omega),1);
B4 = zeros(length(Omega),1);
B5 = zeros(length(Omega),1);
B6 = zeros(length(Omega),1);
alpha_o = zeros(length(Omega),1);
alpha_c = zeros(length(Omega),1);
%%
for i = 1: length(Omega)
%%
%================================pK's======================================
H = 10^-pH;
%-----------Millero at al. (2006) seawater---------------------------------
pK1o = 6320.813/TK+19.568224*log(TK)-126.34048;       
pK1 = 13.4191*Salinity^0.5+0.0331*Salinity-5.33*10^(-5)*Salinity^2-(530.123*Salinity^0.5+6.103*Salinity)/TK-2.0695*Salinity^0.5*log(TK)+pK1o;
pK2o = 5143.692/TK+14.613358*log(TK)-90.18333;        
pK2 = 21.0894*Salinity^0.5+0.1248*Salinity-3.687*10^-4*Salinity^2-(772.483*Salinity^0.5+20.051*Salinity)/TK-3.3336*Salinity^0.5*log(TK)+pK2o;
K1 = 10^(-pK1);
K2 = 10^(-pK2);
Ks = Ksp;%2*10^(log10(10^-8.486)+(-0.77712+0.0028426*TK_r+178.34/TK_r)*(Salinity_r^0.5)+(-0.07711)*Salinity_r+0.0041249*(Salinity_r^1.5));%2.0 prefactor is for aragonite
phi = (1+H/K2+H*H/K1/K2)/(1+H/K1+K2/H);                                     %Ratio of HCO3/CO3 in bulk solution
theta = 10^(8.6-pH);  

%-------------Millero et al. (2007) - NaCl solutions-----------------------
% AA1 = 35.2911*mNaCl^0.5+0.8491*mNaCl-0.32*mNaCl^1.5+0.055*mNaCl^2;
% BB1 = -1583.09*mNaCl^0.5;
% CC1 = -5.4366*(mNaCl^0.5);
% pK1o = -114.3106+5773.67/TK+17.779524*log(TK);
% pK1 = AA1+BB1/TK+CC1*log(TK)+pK1o;                                           
% AA2 = 38.2746*mNaCl^0.5+1.6057*mNaCl-0.647*mNaCl^1.5+0.113*mNaCl^2;
% BB2 = -1738.16*mNaCl^0.5;
% CC2 = -6.0346*mNaCl^0.5;
% pK2o = -83.2997+4821.38/TK+13.5962*log(TK);
% pK2 = AA2+BB2/TK+CC2*log(TK)+pK2o;                                         
% Ks = 10^-8.48;                                                              %Jacobsen (1974) calculated by Ellen using PHREEQC for freshwater
% K1 = 10^(-pK1);
% K2 = 10^(-pK2);
% theta = 10^(8.6-pH);                                                        %Ratio of HCO3/CO3 on calcite surface
% phi = (1+H/K2+H*H/K1/K2)/(1+H/K1+K2/H);                                     %Ratio of HCO3/CO3 in bulk solution

%%
%===========================Speciation=====================================
B1_r(i) = Omega(i)*Ks/A;% [CO32-] (M)
DIC(i) = B1_r(i)*(1+H/K2+(H^2)/(K1*K2));
B2_r(i) = DIC(i)/(1+H/K1+K2/H);% [HCO3-] (M)
B0_r(i) = DIC(i)/(1+K1/H+K1*K2/(H^2)); %[CO2 + H2CO3] (M)
B_r(i) = B1_r(i)+B2_r(i);
%%
%===========================Input parameters===============================
boltz=1.38065E-23;                                                          %Boltzmann constant
epsilon = 0.67E-20;                                                         %kink formation energy 0.62 to 12 is acceptable
gamma = 1.2E-19;                                                            %edge work 0.1 to 1.3 is acceptable
kA1 = 3E6; 
% epsilon = 1e-20;%0.67E-20;                                                    %kink formation energy 0.62 to 12 is acceptable
% gamma = 2.62e-20;%1.2E-19;                                                    %edge work 0.1 to 1.3 is acceptable
% kA1 = 1.7e5;%3E6;                                                      %Qicui
kA2 = kA1;                                                                  %Wolthers et al. (2012)
kB1 = 2*kA1*(1+theta)/(1+phi);                                              %Wolthers et al. (2012)
kB2 = kB1;                                                                  %Wolthers et al. (2012)
vA1 = 2E3;                                                                  %Wolthers et al. (2012)
vA2 = vA1;                                                                  %Wolthers et al. (2012)
%%
%========================Dependent parameters==============================
a = 3.199E-10;                                                              %Closest spacing between A and B sites
d = 27100;                                                                  %Molar density of calcite
S(i) = (A*B1_r(i)/Ks)^0.5;                                                    %Saturation ratio for calcite
k_bar_A = kA1+theta*kA2;                                                    %Rate coefficient for A attachment
k_bar_B = kB1+phi*kB2;                                                      %Rate coefficient for B1 and B2 attachment
v_bar_A = vA1+vA2;                                                          %Rate coefficient for A detachment
vB1 = Ks*k_bar_A*k_bar_B/(v_bar_A*(1+theta));                               %B1 detachment frequency
vB2 = vB1;                                                                  %Wolthers et al. (2012)
v_bar_B = vB1+theta*vB2;                                                    %Wolthers et al. (2012)
PB1(i) = (k_bar_B*B1_r(i)+v_bar_A)/(k_bar_A*A+v_bar_B+(1+theta)*(k_bar_B*B1_r(i)+v_bar_A));% Probability that a given site is a B1 site
PA(i) = 1-(1+theta)*PB1(i);                                                 %Probability that a given site is an A site
PB2(i) = 1-PA(i)-PB1(i);                                                    %Probability that a given site is a B2 site
uA(i) = k_bar_A*A*PB1(i)-v_bar_A*PA(i);                                     %Kink propagation rate
uB(i) = k_bar_B*B1_r(i)*PA(i)-v_bar_B*PB1(i);                                 %Kink propagation rate
uC(i) = uA(i)+uB(i);                                                        %Kink propagation rate
iA(i) = 2*exp(-2*epsilon/(boltz*TK))*(S(i)^2-1)*(v_bar_B*k_bar_A*A/(k_bar_A*A+v_bar_B));%Rate of kink formation on B sites
iB(i) = 2*exp(-2*epsilon/(boltz*TK))*(S(i)^2-1)*(v_bar_A*k_bar_B*B_r(i)/(k_bar_B*B_r(i)+v_bar_A));%Rate of kink formation on A sites
iC(i) = (iA(i)+iB(i))/2;                                                    %Net rate of kink formation
rho_c(i) = (2*iC(i)/(uA(i)+uB(i)))^0.5;                                     %Steady state kink density
y_o(i) = 19*a*gamma/(boltz*TK*log(S(i)));                                   %Step spacing
R_c(i) = rho_c(i)*uC(i)*a^2*d/y_o(i);                                       %Calcite growth rate
%%
%======Oxygen Isotopes============
d18Ow = 0.3;
rVSMOW = 0.0020052;
rw = (d18Ow/1000+1)*rVSMOW;
B2_r(i) = phi*B1_r(i);                                                          %OLD NOTE: don't mix B2 from K's with B2 from phi because they are slightly different
alpha_bw = exp(2590/TK^2 + 0.00189);                                        %Beck et al. (2005)
alpha_cw = exp(2390/TK^2 - 0.00270);                                        %Beck et al. (2005)
alpha_xw = exp((22.5*1000/TK-46.1)/1000);                                   %Wang et al. (2013) - aragonite
%alpha_xw = 0.001+exp((17747/TK-29.777)/1000);                              %Watkins et al. (2013) - calcite
    alpha_O_eq_1 = alpha_xw/alpha_cw;                                       %EFF
    alpha_O_eq_2 = alpha_xw/alpha_bw;                                       %EFF
    alpha_O_f_1 = 0.9995;                                                   %KFF
    alpha_O_f_2 = 0.9966;                                                   %KFF
    alpha_O_b_1 = alpha_O_f_1/alpha_O_eq_1;                                 
    alpha_O_b_2 = alpha_O_f_2/alpha_O_eq_2; 
    alpha_O_b2_b1 = alpha_bw/alpha_cw;
alpha_eq_c_w = alpha_xw;
Delta_cw(i) = 1000*log(alpha_cw);
Delta_bw(i) = 1000*log(alpha_bw);
%%
%=========Carbon isotopes ===============
d13C_DIC = 0.6;
rDIC = 0.01118;                                                             %This is PDB, but the value is arbitrary
%Fractionation factors given by Zhang et al. (1995) as in Zeebe & WG (2001)
epsilon_bg = -0.1141*TC + 10.78;                                            %Zhang et al. (1995) HCO3 - CO2(g)
%epsilon_cg = -0.052*TC + 7.22;                                              %Not used.  Zhang et al. (1995) CO3 - CO2(g)
%epsilon_cb = (epsilon_cg - epsilon_bg)/(1+epsilon_bg/1000);                %Zhang et al. (1995) CO3 - HCO3
%epsilon_db = (epsilon_dg-epsilon_bg)/(1+epsilon_bg/1000);                  %Zhang et al. (1995) CO2 - HCO3
epsilon_cb = -867/TK+2.52;                                                  %Mook (1986) - 0.4 permil difference between CO3 and HCO3
epsilon_db = -9866/TK+24.12;                                                %Mook (1986) - CO2(aq)-HCO3 (is nearly identical to that in Zhang et al. 1995)
delta_B2(i) = (d13C_DIC*(B0_r(i)+B1_r(i)+B2_r(i))-(epsilon_db*B0_r(i)+epsilon_cb*B1_r(i)))/((1+epsilon_db/1000)*B0_r(i)+B2_r(i)+(1+epsilon_cb/1000)*B1_r(i));
delta_B1(i) = delta_B2(i)*(1+epsilon_cb/1000)+epsilon_cb;
delta_B0(i) = delta_B2(i)*(1+epsilon_db/1000)+epsilon_db;%CO2(aq)
alpha_gb = 1/(epsilon_bg/1000+1);
alpha_gx = exp((-2.4612+(7.6663*1000/TK)-(2.988*1000000/(TK^2)))/1000);
alpha_xb = 0.0005+alpha_gb/alpha_gx;
alpha_cb = epsilon_cb/1000+1;
alpha_bc = 1/alpha_cb;
alpha_xc = 0.0005+alpha_xb/alpha_cb;
Delta_cb(i) = 1000*log(alpha_cb);
Delta_bb(i) = 0;
    alpha_C_eq_1 = alpha_xc;                                                %EFF
    alpha_C_eq_2 = alpha_xb;                                                %EFF
    alpha_C_f_1 = 1.000+0.0005;                                                     %KFF
    alpha_C_f_2 = 1.000+0.0005;                                                        %KFF
    alpha_C_b_1 = alpha_C_f_1/alpha_C_eq_1;
    alpha_C_b_2 = alpha_C_f_2/alpha_C_eq_2;
%%
%=====================Isotopologue k's and nu's ===========================
r_B2(i) = (delta_B2(i)/1000+1)*rDIC;
r_B1(i) = (delta_B1(i)/1000+1)*rDIC;
d2(i) = alpha_C_eq_2*r_B2(1);                                               %Old note: =r_c^eq; replaces alpha_eq_c_w for O isotopes
k1 = kB1;                                                                   %Mass 60    
k2 = kB2;                                                                   %Mass 60
k3 = kB1*alpha_C_f_1;                                                       %Mass 61 
k4 = kB2*alpha_C_f_2;                                                       %Mass 61 
k5 = kB1*alpha_O_f_1;                                                       %Mass 62 
k6 = kB2*alpha_O_f_2;                                                       %Mass 62 
v1 = vB1;                                                                   %Mass 60 
v2 = vB2;                                                                   %Mass 60 
v3 = vB1*alpha_C_b_1;                                                       %Mass 61 
v4 = vB2*alpha_C_b_2;                                                       %Mass 61 
v5 = vB1*alpha_O_b_1;                                                       %Mass 62 
v6 = vB2*alpha_O_b_2;                                                       %Mass 62 
%%
%===================Isotopologue concentrations ===========================
B3(i) = B1_r(i)*(delta_B1(i)/1000+1)*rDIC;                                    %Mass 61
B4(i) = B2_r(i)*(delta_B2(i)/1000+1)*rDIC;                                    %Mass 61
B5(i) = 3*alpha_cw*rw*B1_r(i);                                                %Mass 62
B6(i) = 3*alpha_bw*rw*B2_r(i);                                                %Mass 62
%===============Carbon isotope composition of calcite======================
D(i) = ((k3*B3(i)+k4*B4(i))/(k1*B1_r(i)+k2*B2_r(i))) * (v1*PB1(i)+v2*PB2(i)) / (d2(i)*(v3-v4))-v4*(PB1(i)+PB2(i))/(v3-v4);
uB(i) = k1*B1_r(i)*PA(i)+k2*B2_r(i)*PA(i)-v1*PB1(i)-v2*PB2(i);
r_c(i) = (k3*B3(i)*PA(i)+k4*B4(i)*PA(i))/(uB(i)+v3*D(i)+v4*PB1(i)+v4*PB2(i)-v4*D(i));
alpha_c(i) = r_c(i)/r_B1(i);
%===============Oxygen isotope composition of calcite======================
D(i) = ((k5*B5(i)+k6*B6(i))/(k1*B1_r(i)+k2*B2_r(i))) * (v1*PB1(i)+v2*PB2(i))/(3*alpha_eq_c_w*rw*(v5-v6))-v6*(PB1(i)+PB2(i))/(v5-v6);
uB(i) = k1*B1_r(i)*PA(i)+k2*B2_r(i)*PA(i)-v1*PB1(i)-v2*PB2(i);
r_o(i) = (k5*B5(i)*PA(i)+k6*B6(i)*PA(i))/(uB(i)+v5*D(i)+v6*PB1(i)+v6*PB2(i)-v6*D(i));
alpha_o(i) = r_o(i)/(B5(i)/B1_r(i));                                      %Calcite relative to CO32-
%=====================Clumped isotopes=====================================
AFF = 0;%0.280;                                                                %25C Acid fractionation factor from Tripati et al. (2015)
AFF2 = 0;                                                                   %AFF for double clumped
AFF3 = 0;                                                                   %AFF for triple clumped
%D63_B1 = 9.03e-6*TC^2-(3.13e-3)*TC+7.23e-1 - 0.28;                         %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
%D63_B2 = 9.39e-6*TC^2-(3.31e-3)*TC+7.91e-1 - 0.28;                         %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
D63_B1 = 43187/(TK^2)-34.833/TK+0.0007;                                     %Hill et al. 2020
D63_B2 = 43655/(TK^2)-23.643/TK-0.0088;                                     %Hill et al. 2020
D64_B1 = 23492/(TK^2)-52.842/TK+0.0304;                                     %Hill et al. 2020
D64_B2 = 21842/(TK^2)-50.457/TK+0.0291;                                     %Hill et al. 2020
D65_B1 = 112667/(TK^2)-123.11/TK+0.0304;                                    %Hill et al. 2020
D65_B2 = 112026/(TK^2)-97.208/TK+0.009;                                     %Hill et al. 2020
D63eq = 43159/(TK^2)-25.095/TK-0.0078;                                      %Calcite; Hill et al. 2020
D64eq = 23566/(TK^2)-52.319/TK+0.0297;                                      %Calcite; Hill et al. 2020
D65eq = 112667/(TK^2)-102.28/TK+0.012;                                      %Calcite; Hill et al. 2020
%===============Clumped isotope composition of calcite=====================
factor = 0;
alpha_63_f_1 = alpha_O_f_1*alpha_C_f_1+factor;
alpha_63_f_2 = alpha_O_f_2*alpha_C_f_2+factor;
uB(i) = k1*B1_r(i)*PA(i)+k2*B2_r(i)*PA(i)-v1*PB1(i)-v2*PB2(i);
uBc = alpha_O_eq_1*alpha_C_eq_1*(D63eq/1000+1)*(1+phi)/(alpha_63_f_1*(D63_B1/1000+1)+phi*alpha_63_f_2*alpha_O_b2_b1*alpha_bc*(D63_B2/1000+1));%Coefficient for uB
AttRate(i) = k1*B1_r(i)*PA(i)+k2*B2_r(i)*PA(i);
DetRate(i) = v1*PB1(i)+v2*PB2(i);
term1(i) = (alpha_xc*uB(i)*(1+phi)/(alpha_C_f_1 + phi*alpha_C_f_2*alpha_bc) + DetRate(i));
term2(i) = (alpha_O_eq_1*uB(i)*(1+phi)/(alpha_O_f_1 + phi*alpha_O_f_2*alpha_O_b2_b1) + DetRate(i));
RR63(i) = term1(i)*term2(i)*(D63eq/1000+1)/(AttRate(i)*(DetRate(i)+uBc*uB(i)));
Delta63(i) = (RR63(i)-1)*1000;
Delta47(i) = Delta63(i) + AFF;

D63_EIC(i) = B2_r(i)/(B1_r(i)+B2_r(i))*D63_B2 + B1_r(i)/(B1_r(i)+B2_r(i))*D63_B1;
alpha63(i) = RR63(i)/(D63_EIC(i)/1000+1);

%===============Double clumped composition of calcite=====================
factor = 0;
alpha_64_f_1 = alpha_O_f_1*alpha_O_f_1+factor;
alpha_64_f_2 = alpha_O_f_2*alpha_O_f_2+factor;
uB(i) = k1*B1_r(i)*PA(i)+k2*B2_r(i)*PA(i)-v1*PB1(i)-v2*PB2(i);
uBc = alpha_O_eq_1*alpha_O_eq_1*(D64eq/1000+1)*(1+phi)/(alpha_64_f_1*(D64_B1/1000+1)+phi*alpha_64_f_2*alpha_O_b2_b1*alpha_O_b2_b1*(D64_B2/1000+1));%Coefficient for uB
AttRate(i) = k1*B1_r(i)*PA(i)+k2*B2_r(i)*PA(i);
DetRate(i) = v1*PB1(i)+v2*PB2(i);
term1(i) = (alpha_O_eq_1*uB(i)*(1+phi)/(alpha_O_f_1 + phi*alpha_O_f_2*alpha_O_b2_b1) + DetRate(i));
term2(i) = (alpha_O_eq_1*uB(i)*(1+phi)/(alpha_O_f_1 + phi*alpha_O_f_2*alpha_O_b2_b1) + DetRate(i));
RR64(i) = term1(i)*term2(i)*(D64eq/1000+1)/(AttRate(i)*(DetRate(i)+uBc*uB(i)));
Delta64(i) = (RR64(i)-1)*1000;
Delta48(i) = Delta64(i) + AFF2;

D64_EIC(i) = B2_r(i)/(B1_r(i)+B2_r(i))*D64_B2 + B1_r(i)/(B1_r(i)+B2_r(i))*D64_B1;
alpha64(i) = RR64(i)/(D64_EIC(i)/1000+1);
%===============Triple clumped composition of calcite=====================
factor = 0;
alpha_65_f_1 = alpha_O_f_1*alpha_O_f_1*alpha_C_f_1+factor;
alpha_65_f_2 = alpha_O_f_2*alpha_O_f_2*alpha_C_f_2+factor;
uB(i) = k1*B1_r(i)*PA(i)+k2*B2_r(i)*PA(i)-v1*PB1(i)-v2*PB2(i);
uBc = alpha_O_eq_1*alpha_O_eq_1*alpha_C_eq_1*(D65eq/1000+1)*(1+phi)/(alpha_65_f_1*(D65_B1/1000+1)+phi*alpha_65_f_2*alpha_O_b2_b1*alpha_O_b2_b1*alpha_bc*(D65_B2/1000+1));%Coefficient for uB
AttRate(i) = k1*B1_r(i)*PA(i)+k2*B2_r(i)*PA(i);
DetRate(i) = v1*PB1(i)+v2*PB2(i);
term1(i) = (alpha_xc*uB(i)*(1+phi)/(alpha_C_f_1 + phi*alpha_C_f_2*alpha_bc) + DetRate(i));
term2(i) = (alpha_O_eq_1*uB(i)*(1+phi)/(alpha_O_f_1 + phi*alpha_O_f_2*alpha_O_b2_b1) + DetRate(i));
RR65(i) = term1(i)*term2(i)^2*(D65eq/1000+1)/(AttRate(i)^2*(DetRate(i)+uBc*uB(i)));
Delta65(i) = (RR65(i)-1)*1000;
Delta49(i) = Delta65(i) + AFF3;

D65_EIC(i)= B2_r(i)/(B1_r(i)+B2_r(i))*D65_B2+ B1_r(i)/(B1_r(i)+B2_r(i))*D65_B1;
alpha65(i) = RR65(i)/(D65_EIC(i)/1000+1);%now refers to r/r* of CaCO3 relative to r/r* of EIC
% %=============================================
end
%%
alpha_c_spec = interp1(R_c, alpha_c,F);
alpha_o_spec = interp1(R_c, alpha_o,F);
alpha63_spec = interp1(R_c, alpha63,F);
alpha64_spec = interp1(R_c, alpha64,F);
alpha65_spec = interp1(R_c, alpha65,F);
end