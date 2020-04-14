function dddt=Glucose_G3P(time,initial_value,input_data)
%A two compartment delay for computation of concentration of G-3-P in
%adipose tissue. Assuming all glucose is converted to G-6-P eventually this
%implies the concentration of G-3-P for reesterification is a delay of
% a proportion of the gluocse uptake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input data required
%
%time - time vector for ODE solver
%
%initial_value - initial values for adipose G-3-P
%
%
%input_data -  struct which must contain following arrays
%              parameters_glu[GLUT1,GLUT4,AT_delay]
%              parameters_NEFA[re_ester,pnefa,G3P_delay,glucose_use]
%              I_AT     -  vector of adipose tisue delayed insulin.
%              G_art    - vector of arterial glucose concentrations.
%              t        - vector of measured time points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define time series of arterial glucose and insulin from spline data in
%input_data. 
G_art=spline(-30:1:300,input_data.G_art_s,time);
I_AT=spline(-30:1:300,input_data.I_AT,time);
%define necessary parameters
tau=input_data.parameters_nefa(3);
usage=input_data.parameters_nefa(4);
GLUT_1=input_data.parameters_glu(1);
GLUT_4=input_data.parameters_glu(2);
%define intial concentrations of adipose tissue G-3-P in all compartments 
G_11=initial_value(1);
G3P=initial_value(2);
%calculate adipose tissue G-3-P concentraion using two compartmental delay
d_G11=(1/tau).*(usage.*(GLUT_1.*G_art+GLUT_4.*I_AT.*G_art)-G_11);
d_G3P=(1/tau).*(G_11-G3P);

dddt=[d_G11;d_G3P];

