%load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data must be in form of a structure containing the following
%TG_art          - vector of time series of mean arterial triglyceride 
%                  concentrations.
%mean_TG_flux    - vector of time series of mean triglyceride flux across
%                  adipose tissue.
%std_TG_flux     - vector of time series of standard deviations of triglyceride 
%                  flux across the adipose tissue.
%G_art           - vector of time series of mean arterial glucose concentrations.
%mean_G_flux     - vector of time series of mean glucose flux across
%                  adipose tissue.
%std_G_flux      - vector of time series of standard deviations of glucose
%                  flux across the adipose tissue. 
%GLY_art         - vector of time series of mean arterial glycerol concentrations.
%mean_GLY_flux   - vector of time series of mean glycerol flux across
%                  adipose tissue.
%std_GLY_flux    - vector of time series of standard deviations of glycerol
%                  flux across the adipose tissue.
%NEFA_art        - vector of time series of mean arterial NEFA concentrations.
%mean_NEFA_flux  - vector of time series of mean NEFA flux across adipose
%                  tissue.
%std_NEFA_flux   - vector of time series of standard deviations of NEFA flux
%                  across the adipose tissue.
%mean_spill      - vector of mean estimated fractional spill over.
%std_spill       - vector of standard devitation of estimated fractional 
%                  spill over.
%I               - vector of time series of mean arterial insulin
%                  concentrations.
%t               - vector of sampled time points.
%labeling        - string for labelling plots. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%standard deviations are used to weight the error terms. If fitting for
%individual data it is necessary to modifiy the error terms accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('sample_data.mat')
input_data=sample_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up initial parameter values and upper and lower bounds for search
%Define initial values [adipose tissue glycerol,adipose tissue G-3-P, adipose tissue NEFA] 
%specifiy time span for lsqnonlin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_parameters=[0.00957893298670499,155,43.1522822778729,0.0180935942821531,0.000551839775490694,5.38389460989824,0.267467716134800,0.0278107928952439,1.32863066069390,1.71310505639974,0.001,0.0531891784093285,60,0.9];
lower_bounds=[0.001,100,38,0.001,0,1,0.001,0.001,0.01,0.01,0,0,1,0];
upper_bounds=[0.1,300,50,0.1,0.1,45,1,1,5,5,0.05,1.5,100,1];
%initial values 
initial_values=[0.83,0.02,17.6];
time_span=-30:1:300;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate figure of model fit on supplied data: yes=1 no=0
fig_yn=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation of parameter values from data using lsqnonlin algorithm
%Output of Parameter_Fit_Adipose_Model is a structure containing
%p_opt - Optimal paramter set found using lsqnonlin.
%in addition: resnorm, residual, exitflag, output, lambda, and jacobain.
%which may be of use for subsequent analysis of parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
out=Parameter_Fit_Adipose_Model(input_data,initial_parameters,time_span,lower_bounds,upper_bounds,initial_values,fig_yn);
toc
