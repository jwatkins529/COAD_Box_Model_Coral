clear all
close all
%%
FAlk = 1e-6*[0:0.05:1.2]';

xalpha_o = zeros(length(FAlk),1);
xalpha_c = zeros(length(FAlk),1);
xalpha63 = zeros(length(FAlk),1);
xalpha64 = zeros(length(FAlk),1);
xR_c = zeros(length(FAlk),1);
xSAcarb = zeros(length(FAlk),1);
xDIC = zeros(length(FAlk),1);
xH = zeros(length(FAlk),1);
xchi = zeros(length(FAlk),1);
xomega = zeros(length(FAlk),1);
xCa = zeros(length(FAlk),1);
xJCaCO3 = zeros(length(FAlk),1);
xalpha_cEIC = zeros(length(FAlk),1);
xalpha_oEIC = zeros(length(FAlk),1);
xR_oCaCO3 = zeros(length(FAlk),1);
xR_cCaCO3 = zeros(length(FAlk),1);
xC266 = zeros(length(FAlk),1);
xC268 = zeros(length(FAlk),1);
xC366 = zeros(length(FAlk),1);
xC386 = zeros(length(FAlk),1);
xC288 = zeros(length(FAlk),1);
xB2666 = zeros(length(FAlk),1);
xB2668 = zeros(length(FAlk),1);
xB3666 = zeros(length(FAlk),1);
xB3866 = zeros(length(FAlk),1);
xB2886 = zeros(length(FAlk),1);
xC2666 = zeros(length(FAlk),1);
xC2668 = zeros(length(FAlk),1);
xC3666 = zeros(length(FAlk),1);
xC3866 = zeros(length(FAlk),1);
xC2886 = zeros(length(FAlk),1);
xE2666 = zeros(length(FAlk),1);
xE2668 = zeros(length(FAlk),1);
xE3666 = zeros(length(FAlk),1);
xE3866 = zeros(length(FAlk),1);
xE2886 = zeros(length(FAlk),1);
xd18O_CO2 = zeros(length(FAlk),1);
xd18O_HCO3 = zeros(length(FAlk),1);
xd18O_CO3 = zeros(length(FAlk),1);
xd18O_EIC = zeros(length(FAlk),1);
xd18O_CaCO3 = zeros(length(FAlk),1);
xD18O_HCO3 = zeros(length(FAlk),1);
xD18O_CO3 = zeros(length(FAlk),1);
xD18O_EIC = zeros(length(FAlk),1);
xD18OCaCO3 = zeros(length(FAlk),1);
xd13C_CO2 = zeros(length(FAlk),1);
xd13C_HCO3 = zeros(length(FAlk),1);
xd13C_CO3 = zeros(length(FAlk),1);
xd13C_EIC = zeros(length(FAlk),1);
xd13C_CaCO3 = zeros(length(FAlk),1);
xD63_HCO3 = zeros(length(FAlk),1);
xD63_CO3 = zeros(length(FAlk),1);
xD63_EIC = zeros(length(FAlk),1);
xD63_CaCO3 = zeros(length(FAlk),1);
xD64_HCO3 = zeros(length(FAlk),1);
xD64_CO3 = zeros(length(FAlk),1);
xD64_EIC = zeros(length(FAlk),1);
xD64_CaCO3 = zeros(length(FAlk),1);
xF = zeros(length(FAlk),1);
xFAlk = zeros(length(FAlk),1);
xpH = zeros(length(FAlk),1);
xCA = zeros(length(FAlk),1);
xTA = zeros(length(FAlk),1);
xRatio = zeros(length(FAlk),1);
xRatiohydration = zeros(length(FAlk),1);
xRatiohydroxylation = zeros(length(FAlk),1);
xd18O_offset = zeros(length(FAlk),1);
xD47_offset = zeros(length(FAlk),1);
xD48_offset = zeros(length(FAlk),1);
xfCO2m = zeros(length(FAlk),1);

for ii = 1:length(FAlk)
    disp([FAlk(ii)])
[o,time,out]=COAD_Box_Model_Coral('doit',FAlk(ii));

%%
dt = zeros(length(out(:,1)),1);
alpha_o = zeros(length(out(:,1)),1);
alpha_c = zeros(length(out(:,1)),1);
alpha63 = zeros(length(out(:,1)),1);
alpha64 = zeros(length(out(:,1)),1);
SAcarb = zeros(length(out(:,1)),1);
DIC = zeros(length(out(:,1)),1);
H = zeros(length(out(:,1)),1);
chi = zeros(length(out(:,1)),1);
chi13 = zeros(length(out(:,1)),1);
chi18 = zeros(length(out(:,1)),1);
chi63 = zeros(length(out(:,1)),1);
omega = zeros(length(out(:,1)),1);
Ca = zeros(length(out(:,1)),1);
JCaCO3 = zeros(length(out(:,1)),1);
alpha_cEIC = zeros(length(out(:,1)),1);
alpha_oEIC = zeros(length(out(:,1)),1);
R_oCaCO3 = zeros(length(out(:,1)),1);
R_cCaCO3 = zeros(length(out(:,1)),1);
C266 = zeros(length(out(:,1)),1);
C268 = zeros(length(out(:,1)),1);
C366 = zeros(length(out(:,1)),1);
C386 = zeros(length(out(:,1)),1);
C288 = zeros(length(out(:,1)),1);
B2666 = zeros(length(out(:,1)),1); 
B2668 = zeros(length(out(:,1)),1);
B3666 = zeros(length(out(:,1)),1);
B3866 = zeros(length(out(:,1)),1);
B2886 = zeros(length(out(:,1)),1);
C2666 = zeros(length(out(:,1)),1); 
C2668 = zeros(length(out(:,1)),1);
C3666 = zeros(length(out(:,1)),1);
C3866 = zeros(length(out(:,1)),1);
C2886 = zeros(length(out(:,1)),1);
E2666 = zeros(length(out(:,1)),1); 
E2668 = zeros(length(out(:,1)),1);
E3666 = zeros(length(out(:,1)),1);
E3866 = zeros(length(out(:,1)),1);
E2886 = zeros(length(out(:,1)),1);
d18O_CO2 = zeros(length(out(:,1)),1);
d18O_HCO3 = zeros(length(out(:,1)),1);
d18O_CO3 = zeros(length(out(:,1)),1);
d18O_EIC = zeros(length(out(:,1)),1);
d18O_CaCO3 = zeros(length(out(:,1)),1);
D18O_HCO3 = zeros(length(out(:,1)),1);
D18O_CO3 = zeros(length(out(:,1)),1);
D18O_EIC = zeros(length(out(:,1)),1);
D18OCaCO3 = zeros(length(out(:,1)),1);
d13C_CO2 = zeros(length(out(:,1)),1);
d13C_HCO3 = zeros(length(out(:,1)),1);
d13C_CO3 = zeros(length(out(:,1)),1);
d13C_EIC = zeros(length(out(:,1)),1);
d13C_CaCO3 = zeros(length(out(:,1)),1);
D63_HCO3 = zeros(length(out(:,1)),1);
D63_CO3 = zeros(length(out(:,1)),1);
D63_EIC = zeros(length(out(:,1)),1);
D63_CaCO3 = zeros(length(out(:,1)),1);
D64_HCO3 = zeros(length(out(:,1)),1);
D64_CO3 = zeros(length(out(:,1)),1);
D64_EIC = zeros(length(out(:,1)),1);
D64_CaCO3 = zeros(length(out(:,1)),1);
F = zeros(length(out(:,1)),1);
TA = zeros(length(out(:,1)),1);
pH = zeros(length(out(:,1)),1);


for i = 1:length(out(:,1))
    DIC(i) = out(i,1)+out(i,3);
    a = 1;
    KN = 0;
    NT = 0;
    K1 = o.K1;
    K2 = o.K2;
    Kw = o.Kw;
    Ksp = o.Ksp;
    TA(i) = out(i,6);
    b = KN+TA(i)+K1;
    c = TA(i)*KN-KN*NT-Kw+TA(i)*K1+KN*K1+K1*K2...
        -DIC(i)*K1;
    d = TA(i)*KN*K1+TA(i)*K1*K2-Kw*KN-KN*NT*K1-Kw*K1+KN*K1*K2...
        -DIC(i)*KN*K1-2*DIC(i)*K1*K2;
    e = TA(i)*KN*K1*K2-Kw*KN*K1-KN*NT*K1*K2-Kw*K1*K2...
        -DIC(i)*KN*2*K1*K2;
    f = -K1*K2*Kw*KN;
    p = [a b c d e f];
    r = roots(p);
    H(i) = max(real(r));
    pH(i) = -log10(H(i));
    chi(i) = 1/(1+K2/H(i));
    chi13(i) = 1/(1+K2*o.alpha_cb/H(i));
    chi18(i) = 1/(1+K2*o.alpha_cw/o.alpha_bw/H(i));
    chi63(i) = 1/(1+o.K2_63/H(i));
    Ca(i) = out(i,5);
    C266(i) = out(i,1);
    C268(i) = out(i,2);
    C366(i) = out(i,7);
    C386(i) = out(i,9); 
    C288(i) = out(i,11);    
    E2666(i) = out(i,3); 
    E2668(i) = out(i,4);
    E3666(i) = out(i,8);
    E3866(i) = out(i,10);
    E2886(i) = out(i,12);
    B2666(i) = out(i,3)*chi(i);
    B2668(i) = out(i,4)*chi18(i);
    B3666(i) = out(i,8)*chi13(i);
    B3866(i) = out(i,10)/(1+o.K2_63/H(i));
    B2886(i) = out(i,12)/(1+o.K2_64/H(i));
    C2666(i) = out(i,3)*(1-chi(i)); 
    C2668(i) = out(i,4)*(1-chi18(i));
    C3666(i) = out(i,8)*(1-chi13(i));
    C3866(i) = out(i,10)-B3866(i);
    C2886(i) = out(i,12)-B2886(i);
    d13C_HCO3(i) =(B3666(i)/B2666(i)/0.01118-1)*1000; 
    d13C_CO3(i) = (C3666(i)/C2666(i)/0.01118-1)*1000; 
    d13C_EIC(i) = (E3666(i)/E2666(i)/0.01118-1)*1000;
    d18O_CO2(i) = ((0.5*out(i,2)/out(i,1))/o.rVSMOW-1)*1000;                  
    d18O_HCO3(i) = ((1/3*B2668(i)/B2666(i))/o.rVSMOW-1)*1000;               
    d18O_CO3(i) =  ((1/3*C2668(i)/C2666(i))/o.rVSMOW-1)*1000;               
    d18O_EIC(i) = ((1/3*E2668(i)/E2666(i))/o.rVSMOW-1)*1000;                
    D18O_HCO3(i) = 1000*log((d18O_HCO3(i)+1000)/(o.d18Ow+1000));
    D18O_CO3(i) = 1000*log((d18O_CO3(i)+1000)/(o.d18Ow+1000));
    D63_HCO3(i) =(B3866(i)*B2666(i)/(B3666(i)*B2668(i))-1)*1000;            %63K = R/R*
    D63_CO3(i) = (C3866(i)*C2666(i)/(C3666(i)*C2668(i))-1)*1000;            %63K = R/R*
    D63_EIC(i) = (E3866(i)*E2666(i)/(E3666(i)*E2668(i))-1)*1000;            %63K = R/R*
    D64_HCO3(i) =(3*B2886(i)*B2666(i)/(B2668(i)*B2668(i))-1)*1000;          %64K = 3*R/R*
    D64_CO3(i) = (3*C2886(i)*C2666(i)/(C2668(i)*C2668(i))-1)*1000;          %64K = 3*R/R*
    D64_EIC(i) = (3*E2886(i)*E2666(i)/(E2668(i)*E2668(i))-1)*1000;          %64K = 3*R/R*
    omega(i) = Ca(i)*C2666(i)/o.Ksp; 
  
    [alpha_c(i),alpha_o(i),alpha63(i),alpha64(i),~,F(i)]=CaCO3_DIC_Coral(o.TC,o.Ksp,pH(i), out(i,5), out(i,1), out(i,3)*chi(i), out(i,3)*(1-chi(i)));

    if omega(i) < 1
        JCaCO3(i) = 0;
        else
        JCaCO3(i) = 1e-3*F(i);                %F moles/m2/s
    end
% 
    alpha_oEIC(i) = alpha_o(i)*((1-chi18(i))/(1-chi(i)));
    R_oCaCO3(i) = alpha_oEIC(i)*(1/3*E2668(i)/E2666(i));
    alpha_cEIC(i) = alpha_c(i)*((1-chi13(i))/(1-chi(i)));
%     R_oCaCO3(i) = alpha_o(i)*(1/3*C2668(i)/C2666(i));
    R_cCaCO3(i) = alpha_c(i)*C3666(i)/C2666(i);
    d18O_CaCO3(i) = ((R_oCaCO3(i)/o.rVSMOW-1)*1000);  
    d13C_CaCO3(i) = (R_cCaCO3(i)/0.01118-1)*1000;   
    D63_CaCO3(i) = (alpha63(i)*(D63_EIC(i)/1000+1)-1)*1000;
    D64_CaCO3(i) = (alpha64(i)*(D64_EIC(i)/1000+1)-1)*1000;
end
for i = 2:length(out(:,1))
     dt(i) = time(i)-time(i-1);
end


%%
xalpha_oEIC(ii) = alpha_oEIC(end);
xalpha_cEIC(ii) = alpha_cEIC(end);
xalpha_o(ii) = alpha_o(end);                                                    %CaCO3 relative to CO32-
xalpha_c(ii) = alpha_c(end);                                                    %CaCO3 relative to CO32-
xalpha63(ii) = alpha63(end);
xalpha64(ii) = alpha64(end);
xSAcarb(ii) = SAcarb(end);
xDIC(ii) = DIC(end);
xH(ii) = H(end) ;
xpH(ii) = -log10(xH(ii)) ;
xTA(ii) = TA(end);
xchi(ii) = chi(end) ;
xomega(ii) = omega(end);
xCa(ii) = Ca(end);
xJCaCO3(ii) = JCaCO3(end);
xalpha_cEIC(ii) =alpha_cEIC(end);
xalpha_oEIC(ii) =alpha_oEIC(end);
xR_oCaCO3(ii) = R_oCaCO3(end);
xR_cCaCO3(ii) = R_cCaCO3(end);
xC266(ii) = C266(end);
xC268(ii) =C268(end);
xC366(ii) = C366(end) ;
xC386(ii) = C386(end);
xC288(ii) = C288(end);
xB2666(ii) = B2666(end);
xB2668(ii) = B2668(end);
xB3666(ii) = B3666(end);
xB3866(ii) = B3866(end);
xB2886(ii) = B2886(end);
xC2666(ii) = C2666(end);
xC2668(ii) = C2668(end);
xC3666(ii) = C3666(end);
xC3866(ii) = C3866(end);
xC2886(ii) = C2886(end);
xE2666(ii) = E2666(end);
xE2668(ii) = E2668(end);
xE3666(ii) = E3666(end);
xE3866(ii) = E3866(end);
xE2886(ii) = E2886(end);
xd18O_CO2(ii) = d18O_CO2(end);
xd18O_HCO3(ii) = d18O_HCO3(end) ;
xd18O_CO3(ii) = d18O_CO3(end);
xd18O_EIC(ii) = d18O_EIC(end);
xd18O_CaCO3(ii) = d18O_CaCO3(end);
xD18O_HCO3(ii) = D18O_HCO3(end);
xD18O_CO3(ii) = D18O_CO3(end);
xD18O_EIC(ii) = D18O_EIC(end);
xD18OCaCO3(ii) = D18OCaCO3(end);
xd13C_CO2(ii) = d13C_CO2(end) ;
xd13C_HCO3(ii) = d13C_HCO3(end) ;
xd13C_CO3(ii) = d13C_CO3(end) ;
xd13C_EIC(ii) = d13C_EIC(end) ;
xd13C_CaCO3(ii) =d13C_CaCO3(end);
xD63_HCO3(ii) = D63_HCO3(end) ;
xD63_CO3(ii) = D63_CO3(end);
xD63_EIC(ii) = D63_EIC(end) ;
xD63_CaCO3(ii) = D63_CaCO3(end);
xD64_HCO3(ii) = D64_HCO3(end);
xD64_CO3(ii) = D64_CO3(end);
xD64_EIC(ii) = D64_EIC(end) ;
xD64_CaCO3(ii) = D64_CaCO3(end);
xFAlk(ii) = FAlk(ii);
xF(ii) = F(end);
xCA(ii) = o.prefactor;
xRatiohydration(ii) = (o.kf1*C266(end)/(o.kb1*B2666(end)*H(end)));
xRatiohydroxylation(ii) = (o.kf4*C266(end)*o.Kw/H(end))/(o.kb4*B2666(end));
xRatio(ii) = (o.kf1*C266(end)+o.kf4*C266(end)*o.Kw/H(end))/(o.kb1*B2666(end)*H(end)+o.kb4*B2666(end));
AFF = 0.1962;                                                               %from Lucarelli et al. (2022)
AFF2 = 0.131;   
d18Oeq = o.alpha_xw*(o.d18Ow+1000)-1000;
xd18O_offset(ii) = xd18O_CaCO3(ii)-d18Oeq;
D63eq1 = (43159/(o.TK^2)-25.095/o.TK-0.0078);                               %Calcite; Hill et al. 2020
%D47eq2 = 0.0391*10^6/(o.TK^2)+0.154;                                       %Anderson et al. (2021)
D64eq1 = 23566/(o.TK^2)-52.319/o.TK+0.0297;                                 %Calcite; Hill et al. 2020


xD47_offset(ii) = xD63_CaCO3(ii)-D63eq1;
xD48_offset(ii) = xD64_CaCO3(ii)-D64eq1;
end
%%
% Load data
data = load('Davies_cold.txt'); 
      d13C_Davies_cold = data(:,1);
      d18O_Davies_cold = data(:,2)*0.97001-29.99;  
      D47_Davies_cold = data(:,3);                                                 
      D48_Davies_cold = data(:,4); 
      d18O_offset_Davies_cold = data(:,5);                                                 
      D47_offset_Davies_cold = data(:,6);
      D48_offset_Davies_cold = data(:,7);

data2 = load('Davies_warm.txt'); 
      d13C_Davies_warm = data2(:,1);
      d18O_Davies_warm = data2(:,2)*0.97001-29.99;  
      D47_Davies_warm = data2(:,3);                                                 
      D48_Davies_warm = data2(:,4); 
      d18O_offset_Davies_warm = data2(:,5);                                                 
      D47_offset_Davies_warm = data2(:,6);
      D48_offset_Davies_warm = data2(:,7);                                        

%%
figure 
set(gcf,'DefaultAxesLineWidth',2)
set(gcf,'DefaultLineLineWidth',2)
% set(gcf,'Position',[500, 500, 500, 1000])                                   %[a, b, L, W] (a,b) is the lower left corner. L is width. W is height
set(gcf,'Position',[1000, 0, 1000, 1000])

subplot(4,3,1)
plot(xFAlk,xpH)
set(gca,'yscale','log')
xlabel('Alkalinity Pump')
ylabel('pH')
hold on

subplot(4,3,2)
plot(xFAlk,xomega)
xlabel('Alkalinity Pump')
ylabel('Omega')
%axis([8 12.5 5 35])

subplot(4,3,3)
plot(xFAlk,xC2666*1000)
xlabel('Alkalinity Pump')
ylabel('CO_3^{2-} (mmol/L)')
%axis([8 12.5 5 35])

subplot(4,3,4)
plot(xFAlk,xCa*1000)
xlabel('Alkalinity Pump')
ylabel('Ca^{2+} (mmol/L)')

subplot(4,3,5)
plot(xFAlk,(xE2666+xC266)*1000)
xlabel('Alkalinity Pump')
ylabel('DIC (mmol/L)')

subplot(4,3,6)
plot(xFAlk,xF)
% set(gca,'yscale','log')
xlabel('Alkalinity Pump')
ylabel('R_p (mol/m^2/s)')

subplot(4,3,7)
plot(xFAlk,xd18O_CaCO3.*0.97001-29.99)
hold on
plot(xFAlk,xd18O_CO3.*0.97001-29.99)
hold on
plot(xFAlk,xd18O_HCO3.*0.97001-29.99)
hold on
plot(xFAlk,xd18O_EIC.*0.97001-29.99, '--')
hold on
plot(xFAlk,xd18O_CO2.*0.97001-29.99, '-.')
% set(gca,'yscale','log')
xlabel('Alkalinity Pump')
ylabel('\delta 18O')
legend('CaCO_3','CO_3^{2-}','HCO_3^-','EIC','CO_2')


subplot(4,3,8)
plot(xFAlk,xd13C_CaCO3)
hold on
plot(xFAlk,xd13C_CO3, xFAlk,xd13C_HCO3)
% set(gca,'yscale','log')
xlabel('Alkalinity Pump')
ylabel('\delta 13C') 
legend('CaCO_3','CO_3^{2-}','HCO_3^-')


subplot(4,3,9)
plot(xd18O_CaCO3.*0.97001-29.99,xd13C_CaCO3)
hold on
plot(xd18O_CO3.*0.97001-29.99,xd13C_CO3)
hold on
plot(xd18O_HCO3.*0.97001-29.99,xd13C_HCO3)
xlabel('\delta 18O')
ylabel('\delta 13C')
legend('CaCO_3','CO_3^{2-}','HCO_3^-')


subplot(4,3,10)
plot(xFAlk,xD63_CaCO3+AFF, xFAlk,xD63_CO3+AFF)
ax = gca;
xlabel('Alkalinity Pump')
ylabel('\Delta 47')
legend('CaCO_3','CO_3^{2-}')

subplot(4,3,11)
plot(xFAlk,xD64_CaCO3+AFF2, xFAlk,xD64_CO3+AFF2)
xlabel('Alkalinity Pump')
ylabel('\Delta 48')
legend('CaCO_3','CO_3^{2-}')

subplot(4,3,12)
plot(log10(xF),xalpha_o, log10(xF),xalpha_c,'--')
xlabel('log10R_p (mol/m^2/s)')
ylabel('\alpha calcite-CO_3')
legend('\delta^{18}O','\delta^{13}C')

%%
figure
subplot(2,2,1)
plot(xd18O_CaCO3.*0.97001-29.99,xd13C_CaCO3, xd18O_EIC.*0.97001-29.99,xd13C_EIC, '--')  
hold on 
plot(xd18O_CaCO3.*0.97001-29.99,xd13C_CaCO3)  
hold on 
scatter(d18O_Davies_cold,d13C_Davies_cold, 'bo')
hold on
scatter(d18O_Davies_warm,d13C_Davies_warm, 'ro')
pbaspect([1 1 1])

xlabel('\delta 18O')
ylabel('\delta 13C')
%xlim([-8 7])
%ylim([-10 6])


%%
subplot(2,2,2)
plot(xD64_CaCO3+AFF2,xD63_CaCO3+AFF)
yD64_CaCO3=xD64_CaCO3+AFF2;
yD63_CaCO3=xD63_CaCO3+AFF;
hold on 
scatter(D48_Davies_cold,D47_Davies_cold,'bo')
hold on 
scatter(D48_Davies_warm,D47_Davies_warm,'ro')
xlabel('\Delta 48')
ylabel('\Delta 47')
TempC = (-10:5:40)';
TempC = [TempC];
    for i = 1:length(TempC)
        TempK(i) = TempC(i)+273.15;
        D47eq(i) = (43159/(TempK(i)^2)-25.095/TempK(i)-0.0078)+AFF;             %Calcite; Hill et al. 2020
        D48eq(i) = (23566/(TempK(i)^2)-52.319/TempK(i)+0.0297)+AFF2;            %Calcite; Hill et al. 2020
    end
    
plot(D48eq,D47eq,'ko-')
xlim([0.1 0.33])
ylim([0.55 0.80])
pbaspect([2 2 2])

subplot(2,2,3)
plot(xD48_offset,xD47_offset, D48_offset_Davies_warm, D47_offset_Davies_warm,'ro',D48_offset_Davies_cold, D47_offset_Davies_cold,'bo')
xlabel('\Delta_{64} disequilibrium')
ylabel('\Delta_{63} disequilibrium')
 xlim([-0.2 0.01])
 ylim([-.09 0.12])
grid on
hold on