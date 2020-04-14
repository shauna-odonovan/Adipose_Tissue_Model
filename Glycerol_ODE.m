function dgdt=Glycerol_ODE(time,initial_GLY,input_data,parameters_AT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input data required
%
%time - time vector for ODE solver
%
%initial_GLY - initial values [Glycerol adipose space,G-3-P adipose space,NEFA adipose space]
%
%input_data -  struct which must contain following arrays
%              parameters_LPL [K_ad,LPL_delay]
%              I_AT     -  vector of adipose tisue delayed insulin.
%              I_LPL    -  vector of LPL delayed insulin.
%              GLY_art  - vector of arterial glycerol concentrations
%              TG_art   - vector of arterial triglyceride concentrations
%              t        - vector of measured time points.
%
%parameters_AT - vector of parameters being estimated in this ODE system
%[pgly,B_ATL,ATL_max,K_ATL,GLY_delay];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LPL lipolysis
K_ad=input_data.parameters_LPL(1);
%Glycerol transport
pgly=parameters_AT(1);
%ATL lipolysis
B_ATL=parameters_AT(2);
ATL_max=parameters_AT(3);
K_ATL=parameters_AT(4);

%intial values
gly_at=initial_GLY(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate vectos of insulin concentrations and arterial TG and GLY as spline 
%of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_LPL=spline(-30:1:300,input_data.I_LPL,time);
I_AT=spline(-30:1:300,input_data.I_AT,time);
TG_art=spline(-30:1:300,input_data.TG_art_s,time);
GLY_art=spline(-30:1:300,input_data.GLY_art_s,time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d.gly_at=-pgly.*(gly_at-GLY_art+K_ad.*TG_art.*I_LPL) +B_ATL +((ATL_max)./(1+(I_AT)./K_ATL));

%array of results
dgdt=[d.gly_at];