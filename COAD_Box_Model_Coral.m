%% %File: COAD_Box_Model_Corals.m
%Created: 10/14/2021
%Last modified: 6/3/2025
%Author: James M. Watkins
%Department of Earth Sciences, University of Oregon
%watkins4@uoregon.edu --> jamwatkins@ucdavis.edu as of July 1, 2025
%% Description
%This script is builds on the COAD Box Model of Watkins and Devriendt (2022). 
% 
% If you use it, please cite: 
% Watkins, J., Jia, Q, Zhang, S, Devriendt, L., and Chen, S. 2025, Model for 
% dual-clumped isotopes in corals and the conditions of
% biomineralization, Geochemistry, Geophysics, Geosystems, v. xx
% https://doi.org/tbd 

%%
function varargout = COAD_Box_Model_Coral(varargin)
   [varargout{1:nargout}] = feval(varargin{:});
end

function o=setopts(FAlk)
    %----------Set parameters----------------------------------------------
      o.FAlk = FAlk;                                                              %moles/m2/s
   o.TC = 9.4;                                                                 %Cold - Davies et al. (2022)
   o.pH = 8.3;                                                                 %Cold - Davies et al. (2022) - Chen used 8.1
   o.d18Ow = 0.6457;                                                           %Cold - Davies et al. (2022)   
   o.d13C_DIC = 0.69;                                                          %Cold - Davies et al. (2022)
   o.offset = -1;                                                              %Cold - Davies et al. (2022)
%      o.TC = 28;                                                                %Warm - Davies et al. (2022)
%      o.pH = 8.3;                                                               %Warm - Davies et al. (2022)
%      o.d18Ow = 1.148;                                                          %Warm - Davies et al. (2022)                                                                
%      o.d13C_DIC = 1.72;                                                        %Warm - Davies et al. (2022)
%      o.offset = 5;                                                             %Warm - Davies et al. (2022)
%     o.TC = 5.5;                                                               %Cold - Chen et al. (2018)
%     o.pH = 8.3;                                                               %Cold - Chen et al. (2018)
%     o.d18Ow = 0.3;                                                            %Cold - Chen et al. (2018)                                                                
%     o.d13C_DIC = 0.6;                                                         %Cold - Chen et al. (2018)
%     o.offset = 0;                                                             %Cold - Chen et al. (2018)
%     
    o.fCa = 1;
    o.tau = 150;                                                                %inverse of seawater residence time
    o.Dcell = 30e-6;                                                            %cell permeability m/s
    o.z = 10.3e-6;                                                                %thickness of calcifying fluid zone m
    o.TK = o.TC+273.15;
    o.S = 34.2;
    o.Ca = 10e-3; 
    o.DIC = 2.0e-3;                                                             %moles/kg-soln
    o.rVSMOW = 2005.2/1e6;
    o.krate = exp(11.54-8690/o.TK);
    o.prefactor=2000;
    %-----------Millero at al. (2006) - seawater---------------------------
    o.pK1o = 6320.813/o.TK+19.568224*log(o.TK)-126.34048;       
    o.pK1 = 13.4191*o.S^0.5+0.0331*o.S-5.33*10^-5*o.S^2-(530.123*o.S^0.5+6.103*o.S)/o.TK-2.0695*o.S^0.5*log(o.TK)+o.pK1o;
    o.pK2o = 5143.692/o.TK+14.613358*log(o.TK)-90.18333;        
    o.pK2 = 21.0894*o.S^0.5+0.1248*o.S-3.687*10^-4*o.S^2-(772.483*o.S^0.5+20.051*o.S)/o.TK-3.3336*o.S^0.5*log(o.TK)+o.pK2o;
    o.lnKw = 148.96502-13847.26/o.TK-23.6521*log(o.TK)+(118.67/o.TK-5.977+1.0495*log(o.TK))*o.S^0.5-0.01615*o.S;%p258
    o.lnK0 = 9345.17/o.TK-60.2409+23.3585*log(o.TK/100)+o.S*(0.023517-0.00023656*o.TK+0.0047036*(o.TK/100)^2);%p.257
    o.K1 = 10^(-o.pK1);
    o.K2 = 10^(-o.pK2);
    o.Kw = exp(o.lnKw); 
    o.Ksp = 1.8*10^(log10(10^-8.486)+(-0.77712+0.0028426*o.TK+178.34/o.TK)*(o.S^0.5)+(-0.07711)*o.S+0.0041249*(o.S^1.5));%2.0 prefactor is for aragonite
    o.K0 = exp(o.lnK0);
%   o.Ksp = 7.4*1e-9;                                                           %aragonite
    %----------Speciate----------------------------------------------------
    o.H = 10^-o.pH;
    o.OH = o.Kw/o.H;
    o.CO2 = o.DIC/(1+o.K1/o.H+o.K1*o.K2/(o.H^2));                               %Eq. 1.1.9 of Zeebe and Wolf-Gladrow (2001)
    o.HCO3 = o.DIC/(1+o.H/o.K1+o.K2/o.H);                                       %Eq. 1.1.10 of Zeebe and Wolf-Gladrow (2001)
    o.CO3 = o.DIC/(1+o.H/o.K2+o.H^2/(o.K1*o.K2));                               %Eq. 1.1.11 of Zeebe and Wolf-Gladrow (2001)	
    o.EIC = o.HCO3+o.CO3;
    o.chi = 1/(1+o.K2/o.H);
    o.omega = o.Ca*o.CO3/o.Ksp;
    o.TA = o.CO2*(o.K1/o.H+2*o.K1*o.K2/o.H^2)+o.Kw/o.H-o.H;
    %----------Oxygen isotopes---------------------------------------------
    o.rw = (o.d18Ow/1000+1)*o.rVSMOW;                                           %Definition
    o.alpha_dw = exp(2520/o.TK^2 + 0.01212);                                    %Beck et al. 2005
    o.alpha_bw = exp(2590/o.TK^2 + 0.00189);                                    %Beck et al. 2005
    o.alpha_cw = exp(2390/o.TK^2 - 0.00270);                                    %Beck et al. 2005
    o.alpha_xw = exp((17747/o.TK-29.777)/1000);                                 %Watkins et al. (2014); very similar to Coplen (2007)
    o.eOH = -4.4573+10.3255e3/o.TK -0.5976e6/(o.TK^2);                          %Zeebe (2020)
    o.alpha_OH = 1/(o.eOH/1000+1);                                              %Zeebe (2020)
     %o.alpha_OH = exp(0.00048*o.TK-0.1823);
    o.rOH = o.alpha_OH*o.rw;                                                    %Definition
    o.chi18 = 1/(1+o.K2*o.alpha_cw/o.alpha_bw/o.H);                             %Chen et al. (2018)
    %----------Carbon isotopes---------------------------------------------
    o.rDIC = (o.d13C_DIC/1000+1)*0.01118;                                       %Definition
    o.epsilon_gb = -9483/o.TK+23.89;                                            %Mook (1986)
    o.epsilon_cb = -867/o.TK+2.52;                                              %Mook (1986)
    o.epsilon_db = -9866/o.TK+24.12;                                            %Mook (1986)
    o.delta_HCO3 = (o.d13C_DIC*o.DIC-(o.epsilon_db*o.CO2+o.epsilon_cb*o.CO3))/((1+o.epsilon_db/1000)*o.CO2+o.HCO3+(1+o.epsilon_cb/1000)*o.CO3);%Zeebe and Wolf-Gladrow (2001)
    o.delta_CO3 = o.delta_HCO3*(1+o.epsilon_cb/1000)+o.epsilon_cb;              %Zeebe and Wolf-Gladrow (2001)
    o.delta_CO2 = o.delta_HCO3*(1+o.epsilon_db/1000)+o.epsilon_db;              %Zeebe and Wolf-Gladrow (2001)
    o.alpha_gb = o.epsilon_gb/1000+1;
    o.alpha_cb = o.epsilon_cb/1000+1;
    o.alpha_db = o.epsilon_db/1000+1;
    o.alpha_gx = exp((-2.4612+(7.6663*1000/o.TK)-(2.988*1000000/(o.TK^2))-2)/1000);
    o.alpha_xb = o.alpha_gb/o.alpha_gx;
    o.alpha_xc = o.alpha_xb/o.alpha_cb;
    %----------Clumped isotopes--------------------------------------------
    o.D47_CO2eq = 26447/(o.TK^2)+285.51/o.TK-0.3004;                            %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
    %o.D63_HCO3eq = 9.39e-6*o.TC^2-(3.31e-3)*o.TC+7.91e-1 - 0.28;               %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
    %o.D63_CO3eq = 9.03e-6*o.TC^2-(3.13e-3)*o.TC+7.23e-1 - 0.28;                %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
    o.D63_HCO3eq = 43655/(o.TK^2)-23.643/o.TK-0.0088;                           %Hill et al. 2020
    o.D63_CO3eq = 43187/(o.TK^2)-34.833/o.TK+0.0007;                            %Hill et al. 2020
    o.RR47_CO2 = 1+o.D47_CO2eq/1000;                                            %Table 4 - Watkins and Devriendt (2022)
    o.RR63_HCO3 = 1+o.D63_HCO3eq/1000;                                          %Table 4 - Watkins and Devriendt (2022)
    o.RR63_CO3 = 1+o.D63_CO3eq/1000;                                            %Table 4 - Watkins and Devriendt (2022)
    o.alpha63_cb = (1+o.D63_CO3eq/1000)/(1+o.D63_HCO3eq/1000);                  %Table 4 - Watkins and Devriendt (2022)
    o.K2_63 = o.K2*o.alpha63_cb*o.alpha_cb*o.alpha_cw/o.alpha_bw;               %Table 4 - Watkins and Devriendt (2022)
    %----------Double clumped isotopes-------------------------------------
    o.D48_CO2eq = 29306/(o.TK^2)+93.885/o.TK-0.2914;                            %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
    o.D64_HCO3eq = 21842/(o.TK^2)-50.457/o.TK+0.0291;                           %Hill et al. 2020
    o.D64_CO3eq = 23492/(o.TK^2)-52.842/o.TK+0.0304;                            %Hill et al. 2020
    o.RR48_CO2 = 1+o.D48_CO2eq/1000;                                            %Table 5 - Watkins and Devriendt (2022)
    o.RR64_HCO3 = 1+o.D64_HCO3eq/1000;                                          %Table 5 - Watkins and Devriendt (2022)
    o.RR64_CO3 = 1+o.D64_CO3eq/1000;                                            %Table 5 - Watkins and Devriendt (2022)
    o.alpha64_cb = (1+o.D64_CO3eq/1000)/(1+o.D64_HCO3eq/1000);                  %Table 5 - Watkins and Devriendt (2022)
    o.K2_64 = o.K2*o.alpha64_cb*(o.alpha_cw/o.alpha_bw)^2;                      %Table 5 - Watkins and Devriendt (2022)
    %----------Forward k's-------------------------------------------------
    o.kf1o = 10^(329.85-110.541*log10(o.TK)-(17265.4/o.TK));                    %Uchikawa and Zeebe (2012)
    o.kf1 = o.prefactor*o.kf1o;                                                 %Michaelis-Menten kinetics
    o.af1 = o.kf1*0.991;%1.0000;                                                      %Yumol et al. (2020)
    o.bf1 = o.kf1*0.997;%0.993;                                                       %Yumol et al. (2020) - corrected on 3/6/24
    o.cf1 = o.kf1*0.9824;                                                       %Yumol et al. (2020)
    o.kf4 = 10^(13.635-2895/o.TK);                                              %Uchikawa and Zeebe (2012)
    o.af4 = o.kf4*0.970;                                                        %Christensen et al. (2021) corrected and extrapolated to 5C 
    o.bf4 = o.kf4*1.000;                                                        %Christensen et al. (2021) & Clark et al. (1992)
    o.cf4 = o.kf4/1.019;                                                        %Christensen et al. (2021)
 
  %o.alpha13eq = o.alpha_xc;
  %o.alpha18eq = o.alpha_xw/o.alpha_cw;
   %o.af1 = o.kf1/1.007;            %H2O --> HCO3               %Chen
%   o.bf1 = o.kf1/1.010;            %CO2 --> HCO3               %Chen
%   o.cf1 = o.kf1/1.013;            %CO2 --> HCO3 (carbon)      %Chen TYPO in Chen's paper?
%   o.kf4 = 10^(13.635-2895/o.TK);                              %Uchikawa and Zeebe (2012)
%   o.af4 = o.kf4/1.007;            %H2O --> HCO3               %Chen
%   o.bf4 = o.kf4/1.010;            %CO2 --> HCO3               %Chen
%   o.cf4 = o.kf4/1.011;            %CO2 --> HCO3 (carbon)      %Chen TYPO? in Chen's paper?
    
    %----------Backward k's------------------------------------------------
    o.kb1 = o.kf1/o.K1;
    o.ab1 = o.af1/(o.alpha_bw*o.K1);
    o.bb1 = o.bf1/(o.alpha_bw/o.alpha_dw*o.K1);
    o.cb1 = o.cf1*o.alpha_db/o.K1;
    o.kb4 = o.kf4*o.Kw/o.K1;
    o.ab4 = o.af4*o.Kw/(o.alpha_bw/o.alpha_OH*o.K1);
    o.bb4 = o.bf4*o.Kw/(o.alpha_bw/o.alpha_dw*o.K1);
    o.cb4 = o.cf4*o.Kw*o.alpha_db/o.K1;
    %-------Single clumped-------------------------------------------------
    o.KIE_p1 = (4613.8393/o.TK^2-4.0389/o.TK-0.185)/1000+1;                     %Guo (2020) 
    o.KIE_p4 = (-902.7635/o.TK^2+157.1718/o.TK-0.533)/1000+1;                   %Guo (2020) o.KIE_p4 = 1+0.020/1000
    o.KIE_s1 = (-5705.688/o.TK^2-41.5925/o.TK-0.015)/1000+1;                    %Guo (2020)
    o.KIE_s4 = (-11771.2832/o.TK^2-62.7060/o.TK+.168)/1000+1;                   %Guo (2020)
    o.pf1 = o.KIE_p1*o.cf1*o.af1/o.kf1;                                         %Table 4 - Watkins and Devriendt (2022) 
    o.sf1 = o.KIE_s1*o.cf1*o.bf1/o.kf1;                                         %Table 4 - Watkins and Devriendt (2022) 
    o.pf4 = o.KIE_p4*o.cf4*o.af4/o.kf4;                                         %Table 4 - Watkins and Devriendt (2022) 
    o.sf4 = o.KIE_s4*o.cf4*o.bf4/o.kf4;                                         %Table 4 - Watkins and Devriendt (2022)    
    o.pb1 = o.pf1*o.alpha_db/(o.K1*o.RR63_HCO3*o.alpha_bw);                     %Table 4 - Watkins and Devriendt (2022)        
    o.sb1 = o.sf1*o.RR47_CO2*o.alpha_db/(o.K1*o.RR63_HCO3*o.alpha_bw/o.alpha_dw);%Table 4 - Watkins and Devriendt (2022)       
    o.pb4 = o.pf4*o.Kw*o.alpha_db*o.alpha_OH/(o.K1*o.RR63_HCO3*o.alpha_bw);     %Table 4 - Watkins and Devriendt (2022)     
    o.sb4 = o.sf4*o.Kw*o.RR47_CO2*o.alpha_db/(o.K1*o.RR63_HCO3*o.alpha_bw/o.alpha_dw);%Table 4 - Watkins and Devriendt (2022)
    %======Double clumped======================================================
    o.KIE_p1p = (13249.5324/o.TK^2-37.8964/o.TK+0.027)/1000+1;                  %Guo (2020) 
    o.KIE_p4p = (5859.1625/o.TK^2-3.7964/o.TK-0.197)/1000+1;                    %Guo (2020)
    o.KIE_s1p = (-18411.4121/o.TK^2-3.7575/o.TK+0.074)/1000+1;                  %Guo (2020)
    o.KIE_s4p = (-12333.8137/o.TK^2+8.6005/o.TK+0.024)/1000+1;                  %Guo (2020)
    o.pf1p = o.KIE_p1p*o.bf1*o.af1/o.kf1;                                       %Table 5 - Watkins and Devriendt (2022) 
    o.sf1p = o.KIE_s1p*o.bf1*o.bf1/o.kf1;                                       %Table 5 - Watkins and Devriendt (2022) 
    o.pf4p = o.KIE_p4p*o.bf4*o.af4/o.kf4;                                       %Table 5 - Watkins and Devriendt (2022) 
    o.sf4p = o.KIE_s4p*o.bf4*o.bf4/o.kf4;                                       %Table 5 - Watkins and Devriendt (2022) 
    o.pb1p = o.pf1p/(o.K1*o.RR64_HCO3*(o.alpha_bw/o.alpha_dw)*o.alpha_bw);      %Table 5 - Watkins and Devriendt (2022)                    
    o.sb1p = o.sf1p*o.RR48_CO2/(o.K1*o.RR64_HCO3*(o.alpha_bw/o.alpha_dw)^2);    %Table 5 - Watkins and Devriendt (2022)              
    o.pb4p = o.pf4p/(o.K1/o.Kw*o.RR64_HCO3*(o.alpha_bw/o.alpha_dw)*o.alpha_bw/o.alpha_OH);%Table 5 - Watkins and Devriendt (2022)  
    o.sb4p = o.sf4p*o.RR48_CO2/(o.K1/o.Kw*o.RR64_HCO3*(o.alpha_bw/o.alpha_dw)^2);%Table 5 - Watkins and Devriendt (2022)               
    %----------chi's-------------------------------------------------------
    o.chi18 = 1/(1+o.K2*o.alpha_cw/o.alpha_bw/o.H);                             %Table 3 - Watkins and Devriendt (2022)
    o.chi13 = 1/(1+o.K2*o.alpha_cb/o.H);                                        %Table 3 - Watkins and Devriendt (2022)
    o.chi63 = 1/(1+o.K2_63/o.H);                                                %Table 4 - Watkins and Devriendt (2022)
    o.chi64 = 1/(1+o.K2_64/o.H);                                                %Table 5 - Watkins and Devriendt (2022)
    %----------Isotope ratios----------------------------------------------
    o.oCO2 = o.alpha_dw*o.rw;
    o.oHCO3 = o.alpha_bw*o.rw;
    o.oCO3 = o.alpha_cw*o.rw;
    o.cCO2 = (o.delta_CO2/1000+1)*0.01118;
    o.cHCO3 = (o.delta_HCO3/1000+1)*0.01118;
    o.cCO3 = (o.delta_CO3/1000+1)*0.01118;
    %----------Isotopologue concentrations---------------------------------
    o.C6OH = o.OH;
    o.C8OH = o.rOH*o.C6OH;
    o.C266 = o.CO2;
    o.C268 = 2*o.oCO2*o.C266;
    o.C366 = o.cCO2*o.C266;
    o.B2666 = o.HCO3;
    o.B2668 = 3*o.oHCO3*o.B2666;
    o.B3666 = o.cHCO3*o.B2666;
    o.C2666 = o.CO3;
    o.C2668 = 3*o.oCO3*o.C2666;
    o.C3666 = o.cCO3*o.C2666;  
    o.E2666 = o.B2666+o.C2666;
    o.E2668 = o.B2668+o.C2668;
    o.E3666 = o.B3666+o.C3666;
    o.C386 = o.RR47_CO2*o.C366*o.C268/o.C266;                
    o.B3866 = o.RR63_HCO3*o.B3666*o.B2668/o.B2666;  
    o.C3866 = o.RR63_CO3*o.C3666*o.C2668/o.C2666;          
    o.E3866 = o.B3866/o.chi63;     
    o.C288 = 1/4*o.RR48_CO2*o.C268*o.C268/o.C266;
    o.B2886 = 1/3*o.RR64_HCO3*o.B2668*o.B2668/o.B2666;
    o.C2886 = 1/3*o.RR64_CO3*o.C2668*o.C2668/o.C2666;
    o.E2886 = o.B2886/o.chi64;
    %----------Seawater leak-----------------------------------------------
    o.DICsw = 2.0e-3;
    o.Casw = 10.3e-3;
    o.d13Csw = o.d13C_DIC;%
    o.d18Osw = o.d18Ow;
    o.pHsw = o.pH;
    o.R13sw = (o.d13Csw/1000+1)*0.01118;
    o.R18sw = (o.d18Osw/1000+1)*0.0020052;
    o.Hsw = 10^-o.pHsw;
    o.OHsw = o.Kw/o.Hsw;
    o.CO2sw = o.DICsw/(1+o.K1/o.Hsw+o.K1*o.K2/(o.Hsw^2));
    o.HCO3sw = o.DICsw/(1+o.Hsw/o.K1+o.K2/o.Hsw);
    o.CO3sw = o.DICsw/(1+o.Hsw/o.K2+o.Hsw^2/(o.K1*o.K2));
    o.Alksw = o.CO2sw*(o.K1/o.Hsw+2*o.K1*o.K2/(o.Hsw^2))+o.Kw/o.Hsw-o.Hsw;
    o.EICsw = o.HCO3sw+o.CO3sw;
    o.delta_HCO3sw = (o.d13Csw*o.DICsw-(o.epsilon_db*o.CO2sw+o.epsilon_cb*o.CO3sw))/((1+o.epsilon_db/1000)*o.CO2sw+o.HCO3sw+(1+o.epsilon_cb/1000)*o.CO3sw);
    o.delta_CO3sw = o.delta_HCO3sw*(1+o.epsilon_cb/1000)+o.epsilon_cb;
    o.delta_CO2sw = o.delta_HCO3sw*(1+o.epsilon_db/1000)+o.epsilon_db;
    o.cHCO3sw = (o.delta_HCO3sw/1000+1)*0.01118;
    o.cCO2sw = (o.delta_CO2sw/1000+1)*0.01118;
    o.cCO3sw = (o.delta_CO3sw/1000+1)*0.01118;
    o.oCO2sw = o.alpha_dw*o.R18sw;
    o.oHCO3sw = o.alpha_bw*o.R18sw;
    o.oCO3sw = o.alpha_cw*o.R18sw;
    o.C266sw = o.CO2sw;
    o.C268sw = 2*o.oCO2sw*o.C266sw;
    o.B2666sw = o.HCO3sw;
    o.B2668sw = 3*o.oHCO3sw*o.B2666sw;
    o.C2666sw = o.CO3sw;
    o.C2668sw = 3*o.oCO3sw*o.C2666sw;
    o.E2666sw = o.B2666sw+o.C2666sw;
    o.E2668sw = o.B2668sw+o.C2668sw;
    o.C366sw = o.cCO2sw*o.C266sw;
    o.B3666sw = o.cHCO3sw*o.B2666sw;
    o.C3666sw = o.cCO3sw*o.C2666sw;
    o.E3666sw = o.B3666sw+o.C3666sw;
    o.C288sw = 1/4*o.RR48_CO2*o.C268sw*o.C268sw/o.C266sw;
    o.C386sw = o.RR47_CO2*o.C366sw*o.C268sw/o.C266sw;                
    %----------chi's-------------------------------------------------------
    o.chi63sw = 1/(1+o.K2_63/o.Hsw);                                                %Table 4 - Watkins and Devriendt (2022)
    o.chi64sw = 1/(1+o.K2_64/o.Hsw);                                                %Table 5 - Watkins and Devriendt (2022)    
    o.B3866sw = o.RR63_HCO3*o.B3666sw*o.B2668sw/o.B2666sw;  
    o.B2886sw = 1/3*o.RR64_HCO3*o.B2668sw*o.B2668sw/o.B2666sw;
    o.E3866sw = o.B3866sw/o.chi63sw;     
    o.E2886sw = o.B2886sw/o.chi64sw;
    
    %----------Cellular CO2------------------------------------------------
    o.d13Ccell = o.delta_CO2sw+o.offset;                                            %Value is ~-11 VPDB. Organic matter is ~-22.                                              
    o.R13cell = (o.d13Ccell/1000+1)*0.01118;
    o.CO2cell = 13e-6;             %Chen et al. (2018);
    o.cCO2cell = (o.d13Ccell/1000+1)*0.01118;
    o.oCO2cell = o.oCO2sw;
    o.C266cell = o.CO2cell;
    o.C268cell = 2*o.oCO2cell*o.C266cell;
    o.C366cell = o.cCO2cell*o.C266cell;
    o.C386cell = o.RR47_CO2*o.C366cell*o.C268cell/o.C266cell;                
    o.C288cell = 1/4*o.RR48_CO2*o.C268cell*o.C268cell/o.C266cell;

end

function [o, time, out]=doit(FAlk)
    close all
    o=setopts(FAlk);
    [time, out]=driver(o);
end

function [time, out]=driver(o)
    t=[1e-8 1e4];     % time (s)
    C2660 = o.C266;         %1 = 12C16O16O
    C2680 = o.C268;         %2 = 12C16O18O
    E26660 = o.E2666;       %3 = E12C16O16O16O
    E26680 = o.E2668;       %4 = E12C16O16O18O
    Ca0 = o.Ca;             %5 = Ca2+
    Alk0 = o.TA;            %6 = total alkalinity    
    C3660 = o.C366;         %7 = 13C16O16O
    E36660 = o.E3666;       %8 = E13C16O16O16O
    C3860 = o.C386;         %9 = clumped CO2
    E38660 = o.E3866;       %10 = clumped EIC
    C2880 = o.C288;         %11 = double clumped CO2
    E28860 = o.E2886;       %12 = double clumped EIC
    
    y0=[C2660,C2680,E26660,E26680,Ca0,Alk0,C3660,E36660,C3860,E38660,C2880,E28860]; %Initial conditions vector
    options = odeset('MaxOrder', 4, 'MaxStep', 1e4, 'RelTol', 10^-4, 'AbsTol',10^-4);
     [time, out]=ode23s(@(t,y0) carbonate(t, y0,o.TC,o.FAlk,o.Dcell,o.z,o.tau,o.fCa,o.kf1,o.af1,o.bf1,o.cf1,o.pf1,o.sf1,o.pf1p,o.sf1p,o.kf4,o.af4,o.bf4,o.cf4,o.pf4,o.sf4,o.pf4p,o.sf4p,o.kb1,o.ab1,o.bb1,o.cb1,o.pb1,o.sb1,o.pb1p,o.sb1p,o.kb4,o.ab4,o.bb4,o.cb4,o.pb4,o.sb4,o.pb4p,o.sb4p,o.rw,o.Ksp,o.K1,o.K2,o.K2_63,o.K2_64,o.Kw,o.rOH,o.alpha_cw,o.alpha_bw,o.alpha_cb,o.Casw,o.Alksw,o.C266sw,o.C268sw,o.C366sw,o.E2666sw,o.E2668sw,o.E3666sw,o.C288sw,o.C386sw,o.E3866sw,o.E2886sw,o.C266cell,o.C268cell,o.C366cell,o.C386cell,o.C288cell), t, y0, options);
end

%   


function dy=carbonate(t,y,TC,FAlk,Dcell,z,tau,fCa,kf1,af1,bf1,cf1,pf1,sf1,pf1p,sf1p,kf4,af4,bf4,cf4,pf4,sf4,pf4p,sf4p,kb1,ab1,bb1,cb1,pb1,sb1,pb1p,sb1p,kb4,ab4,bb4,cb4,pb4,sb4,pb4p,sb4p,rw,Ksp,K1,K2,K2_63,K2_64,Kw,rOH,alpha_cw,alpha_bw,alpha_cb,Casw,Alksw,C266sw,C268sw,C366sw,E2666sw,E2668sw,E3666sw,C288sw,C386sw,E3866sw,E2886sw,C266cell,C268cell,C366cell,C386cell,C288cell);
    dy=zeros(size(y));
    DIC = y(1)+y(3);
    a = 1;
    KN = 0;
    NT = 0;
    TA = y(6);
    b = KN+TA+K1;
    c = TA*KN-KN*NT-Kw+TA*K1+KN*K1+K1*K2...
        -DIC*K1;
    d = TA*KN*K1+TA*K1*K2-Kw*KN-KN*NT*K1-Kw*K1+KN*K1*K2...
        -DIC*KN*K1-2*DIC*K1*K2;
    e = TA*KN*K1*K2-Kw*KN*K1-KN*NT*K1*K2-Kw*K1*K2...
        -DIC*KN*2*K1*K2;
    f = -K1*K2*Kw*KN;
    p = [a b c d e f];
    r = roots(p);
    H = max(real(r));
    C6OH = Kw/H;
    C8OH = rOH*C6OH;
    pH = -log10(H);
    chi = 1/(1+K2/H);
    chi18 = 1/(1+K2*alpha_cw/alpha_bw/H);
    chi13 = 1/(1+K2*alpha_cb/H);
    chi63 = 1/(1+K2_63/H);
    chi64 = 1/(1+K2_64/H);                                      
    omega = y(3)*(1-chi)*y(5)/Ksp;                                              %CO3 = y(3)(1-chi)
    [alpha_c,alpha_o,alpha63,alpha64,alpha65,F]=CaCO3_DIC_Coral(TC, Ksp, pH, y(5), y(1), y(3)*chi, y(3)*(1-chi));
        if omega < 1
            JCaCO3 = 0;
        else
           JCaCO3 = 1e-3*F;                                                     %F in moles/m2/s; 1e3 is for converting to moles/m3
        end
       
        alpha_63EIC = alpha63;
        alpha_64EIC = alpha64;
        alpha_65EIC = alpha65;
        alpha_cEIC = alpha_c*((1-chi13)/(1-chi));
        alpha_oEIC = alpha_o*((1-chi18)/(1-chi));  

    dy(1) = -kf1*y(1)+kb1*y(3)*chi*H - kf4*y(1)*C6OH+kb4*y(3)*chi + Dcell/z*(C266cell-y(1)) + 1/tau*(C266sw-y(1));   
    dy(2) = -bf1*y(2)+2/3*bb1*y(4)*chi18*H - bf4*y(2)*C6OH+2/3*bb4*y(4)*chi18 + Dcell/z*(C268cell-y(2))+ 1/tau*(C268sw-y(2));     
    dy(3) =  kf1*y(1)-kb1*y(3)*chi*H + kf4*y(1)*C6OH-kb4*y(3)*chi + 1/tau*(E2666sw-y(3)) - 1/z*JCaCO3; 
    dy(4) =  af1*y(1)*rw-1/3*ab1*y(4)*chi18*H + bf1*y(2)-2/3*bb1*y(4)*chi18*H + af4*y(1)*C8OH-1/3*ab4*y(4)*chi18 + bf4*y(2)*C6OH-2/3*bb4*y(4)*chi18 + 1/tau*(E2668sw-y(4))-1/z*JCaCO3*y(4)/y(3)*alpha_oEIC;
    dy(5) = 1/tau*(Casw-y(5))+1/z*(0.5*fCa*FAlk/1000-JCaCO3);
    dy(6) = 1/tau*(Alksw-y(6))+1/z*(FAlk/1000-2*JCaCO3);
    dy(7) = -cf1*y(7)+cb1*y(8)*chi13*H - cf4*y(7)*C6OH+cb4*y(8)*chi13 + Dcell/z*(C366cell-y(7)) + 1/tau*(C366sw-y(7));
    dy(8) = cf1*y(7)-cb1*y(8)*chi13*H + cf4*y(7)*C6OH-cb4*y(8)*chi13 + 1/tau*(E3666sw-y(8))- 1/z*JCaCO3*y(8)/y(3)*alpha_cEIC;
    dy(9) =-sf1*y(9)+2/3*sb1*y(10)*chi63*H - sf4*y(9)*C6OH+2/3*sb4*y(10)*chi63 + Dcell/z*(C386cell-y(9)) + 1/tau*(C386sw-y(9));
    dy(10) = sf1*y(9)-2/3*sb1*y(10)*chi63*H + sf4*y(9)*C6OH-2/3*sb4*y(10)*chi63 + pf1*y(7)*rw-1/3*pb1*y(10)*chi63*H + pf4*y(7)*C8OH-1/3*pb4*y(10)*chi63+ 1/tau*(E3866sw-y(10))- 1/z*JCaCO3*y(10)/y(3)*alpha_63EIC*alpha_cEIC*alpha_oEIC; 
    dy(11) = -sf1p*y(11)+1/3*sb1p*y(12)*chi64*H - sf4p*y(11)*C6OH+1/3*sb4p*y(12)*chi64 + Dcell/z*(C288cell-y(11)) + 1/tau*(C288sw-y(11));
    dy(12) = sf1p*y(11)-1/3*sb1p*y(12)*chi64*H + sf4p*y(11)*C6OH-1/3*sb4p*y(12)*chi64 + pf1p*y(2)*rw-2/3*pb1p*y(12)*chi64*H + pf4p*y(2)*C8OH-2/3*pb4p*y(12)*chi64+ 1/tau*(E2886sw-y(12))- 1/z*JCaCO3*y(12)/y(3)*alpha_64EIC*alpha_oEIC*alpha_oEIC;     
 
end

