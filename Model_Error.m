function error=Model_Error(parameters,input_data,time)

%define arterial concentrations from input_data 
%(only sampled time points are used in the error calculation)
TG_art=input_data.TG_art;
G_art=input_data.G_art;
GLY_art=input_data.GLY_art;
NEFA_art=input_data.NEFA_art;
I_art=input_data.I;

%set parameter values given vector of parameters from input
%parameters for LPL lipolysis
K_ad=parameters(1);
LPL_delay=parameters(2);
input_data.parameters_LPL=[K_ad,LPL_delay];
%parameters for spillover
D_spill=parameters(3);
input_data.parameters_spill=D_spill;
%parameters for glucose flux
GLUT_1=parameters(4);
GLUT_4=parameters(5);
AT_delay=parameters(6);
input_data.parameters_glu=[GLUT_1,GLUT_4,AT_delay];
%parameters for glycerol flux
pgly=parameters(7);
B_ATL=parameters(8);
ATL_max=parameters(9);
K_ATL=parameters(10);
parameters_GLY=[pgly,B_ATL,ATL_max,K_ATL];
input_data.parameters_ATL=parameters_GLY;
%parameters for NEFA
re_ester=parameters(11);
pnefa=parameters(12);
G3P_delay=parameters(13);
glucose_use=parameters(14);
parameters_nefa=[re_ester,pnefa,G3P_delay,glucose_use];
input_data.parameters_nefa=parameters_nefa;
%generate vector of three fold LPL insulin delay
input_data.delay=LPL_delay;
I_pl=input_data.I;
initial_I=[I_pl(1),I_pl(1),I_pl(1)];
[T,I_delay]=ode15s(@insulin_delay,time,initial_I,[],input_data);
I_LPL=I_delay([1,31,91,151,211,271,331],3)';
input_data.I_LPL=I_delay(:,3)';
%generate vector of three fold  AT insulin delay
input_data.delay=AT_delay;
[T,I_delay]=ode15s(@insulin_delay,time,initial_I,[],input_data);
I_AT=I_delay([1,31,91,151,211,271,331],3)';
input_data.I_AT=I_delay(:,3)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRIGLYCERIDE FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AV-difference in TG due to LPL lipolysis
D_TG=-K_ad.* I_LPL.* TG_art;
%measured TG flux
A_V_TG=-1.*input_data.mean_TG_flux;
%error term SSE/sigma^2
error_LPL=(A_V_TG-D_TG)./input_data.std_TG_flux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPILLOVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set basal insulin
I_b=input_data.I(1);
%setting plasma insulin vector equal to the measured arterial insulin
I_PL=input_data.I;
%computed fractional spill over
spill=(D_spill.*(I_b./I_PL))/100;
%compute error for spill over term
error_spill=(spill(3:7)-input_data.mean_spill(1:5))./input_data.std_spill;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GLUCOSE FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%equation for glucose flux
D_G=-GLUT_1.*G_art-GLUT_4.*I_AT.*G_art;
%arteial venous difference in glucose
A_V_G=-1.*input_data.mean_G_flux;
error_glucose=(A_V_G-D_G)./input_data.std_G_flux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GLYCEROL FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set initial value for adipose glycerol concentration - not measured
initial_GLY=input_data.initial_values(1);
%calculate the time series of arterial glycerol concentration given
%parameters suppllied in the input of function (parameters_GLY)
[T,G]=ode15s(@Glycerol_ODE,time,initial_GLY,[],input_data,parameters_GLY);
gly_at=G([1,31,91,151,211,271,331])';
%equation for glycerol flux
AV_GLY=K_ad.*(TG_art).*(I_LPL)+pgly.*(gly_at-(GLY_art+K_ad.*TG_art.*I_LPL));

%error measure
measured_AV_GLY=input_data.mean_GLY_flux;

error_GLY=(-1*(measured_AV_GLY)-(AV_GLY))./input_data.std_GLY_flux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NEFA FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set initial value for adipose NEFA and G-3-P concentration - not measured
initial_values=input_data.initial_values(2:3);

%generate adipose G-3-P concentration as a delay of gluocse flux
initial_G3P=[initial_values(1),initial_values(1)];
[T,G_g3p]=ode15s(@Glucose_G3P,time,initial_G3P,[],input_data);
input_data.G3P=G_g3p(:,2)';


%Generate arterial NEFA
[T,AT_NEFA]=ode15s(@Glucose_G3P_NEFA_ODE,time,initial_values,[],input_data,parameters_nefa);
NEFA_at=AT_NEFA([1,31,91,151,211,271,331],2)';

%equation for NEFA flux
AV_NEFA= 3.*((D_spill.*(I_b./I_art))./100).*(K_ad.*TG_art.*I_LPL)+pnefa.*(NEFA_at-(NEFA_art+3.*((D_spill.*(I_b./I_art))./100).*(K_ad.*TG_art.*I_LPL)));
%error measure.
measured_AV_NEFA=-1.*input_data.mean_NEFA_flux;
error_nefa=(measured_AV_NEFA-AV_NEFA)./input_data.std_NEFA_flux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%total error = composition of five error measures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error=[error_LPL,error_spill,error_glucose,error_GLY,error_nefa];