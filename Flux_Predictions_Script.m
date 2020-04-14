%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define dependent inputs necessary for model predictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data structure from file - if available.
load('sample_data.mat');

input_data=sample_data;

%or

%manually input_data strucuture for model - by default set to sample data.
%input_data.TG_art=sample_data.TG_art;
%input_data.G_art=sample_data.G_art;
%input_data.GLY_art=sample_data.GLY_art;
%input_data.NEFA_art=sample_data.NEFA_art;
%input_data.I=sample_data.I;
%input_data.t=sample_data.t;
%input_data.labeling='sample data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define parameter values - by default these are set to the optmial
%parameter set below is for the supplied sample data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters for LPL lipolysis
K_ad=0.01;
LPL_delay=158.13;
%parameters for spillover
D_spill=41.53;
%parameters for glucose flux
GLUT_1=0.0159;
GLUT_4=0.00082;
AT_delay=13.85;
%parameters for glycerol flux
pgly=0.2705;
B_ATL=0.0013;
ATL_max=1.4373;
K_ATL=1.49185;
%parameters for NEFA
re_ester=0.000295;
pnefa=0.0478;
G3P_delay=31.74;
glucose_use=0.8803;

parameters=[K_ad,LPL_delay,D_spill,GLUT_1,GLUT_4,AT_delay,pgly,B_ATL,ATL_max,K_ATL,re_ester,pnefa,G3P_delay,glucose_use];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specify vector of values for the initial adipose metabolite concentrations.
%[adipose glycerol, adipose G-3-P, adipose NEFA]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_values=[0.83,0.02,17.6];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specify time span for model predictions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=-30:1:300;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%predict fluxes for input data and parameters using adipose tissue model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Flux_Predictions(parameters,input_data,initial_values,time,1);
