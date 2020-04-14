function dgdt=Glucose_G3P_NEFA_ODE(time,initial_val,input_data,parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input data required
%
%time - time vector for ODE solver
%
%initial_val - initial values [G-3-P adipose space,NEFA adipose space]
%
%input_data -  struct which must contain following arrays
%              parameters_spill [I_spill,D_spill]
%              parameters_NEFA [pgly,B_ATL,ATL_max,K_ATL]
%              parameters_LPL [K_ad,LPL_delay]
%              I        -  vector of arterial insulin concentrations
%              G_art    - vector of arterial glucose concentrations
%              TG_art   - vector of arterial triglyceride concentrations
%              NEFA_art - vector of arterial NEFA concentrations
%              t        - vector of measured time points.
%
%parameters - vector of parameters being estimated in this ODE system
%[re_ester,pnefa]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define parameters for NEFA dynamics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fractional spill over.
D_spill=input_data.parameters_spill(1);
I_b=input_data.I(1);
%define parameter for LPL lipolysis of plasma TG
K_ad=input_data.parameters_LPL(1);
%define parameter governing lipolysis of TG stored within adipocyte
%pgly=input_data.parameters_ATL(1);
B_ATL=input_data.parameters_ATL(2);
ATL_max=input_data.parameters_ATL(3);
K_ATL=input_data.parameters_ATL(4);

%new parameters - passive diffusion of circulating NEFA and addipose tissue
%NEFA  and Reesterification of adipose tissue NEFA
re_ester=parameters(1);
pnefa=parameters(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define data for dependent inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TG_art=spline(-30:1:300,input_data.TG_art_s,time);
NEFA_art=spline(-30:1:300,input_data.NEFA_art_s,time);
I_art=spline(-30:1:300,input_data.I_art_s,time);

I_LPL=spline(-30:1:300,input_data.I_LPL,time);
I_AT=spline(-30:1:300,input_data.I_AT,time);
G_3_P_production=spline(-30:1:300,input_data.G3P,time);

%initial values
G_3_P=initial_val(1);
NEFA_at=initial_val(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Equation describing the rate of change of adipose tissue G-3-P concentration 
%twice the rate of prodiction - rate of utilisation in reesterification
dG_3_P_dt=2.*G_3_P_production-re_ester.*I_AT.*NEFA_at.*G_3_P;

%Equation describing rate of change of adipose tissue NEFA concentration
%= rate of entrapment of NEFA released by LPL lipolysis + transport of NEFA from plasma along concentration gradient + release of NEFA from lipolysis of stored TG within adipose space - reesterification of adipose NEFA
dNEFA_at_dt= 3.*(1-((D_spill.*(I_b./I_art))./100)).*(K_ad.*TG_art.*I_LPL)-pnefa.*(NEFA_at-(NEFA_art+3.*((D_spill.*(I_b./I_art))./100).*(K_ad.*TG_art.*I_LPL)))+3.*(B_ATL + (ATL_max./(1+(I_AT./K_ATL))))-3.*re_ester.*I_AT.*NEFA_at.*G_3_P;

dgdt=[dG_3_P_dt;dNEFA_at_dt];


