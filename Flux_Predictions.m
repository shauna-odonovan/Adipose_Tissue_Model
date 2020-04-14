function out=Flux_Predictions(parameters,input_data,initial_val,time,fig_yn)

%This function will predict adipose tissue fluxes of triglyceride, glucose,
%NEFA,and glycerol for a user supplied parameter set and necessary
%depdendent inputs (time sereies of arterial concentrations of metabolites)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters     - Vector containing values for each of the fourteen model 
%                 parameters 
%input_data     - structure containing the dependent inputs for the adipose
%                 tissue model. Input_data must contain the following.
%                 TG_art    - vector of time series of mean arterial triglyceride 
%                             concentrations.
%                 G_art     - vector of time series of mean arterial glucose 
%                             concentrations. 
%                 GLY_art   - vector of time series of mean arterial glycerol 
%                             concentrations.
%                 NEFA_art  - vector of time series of mean arterial NEFA
%                             concentrations.
%                 I         - vector of time series of mean arterial insulin
%                             concentrations.
%                 t         - vector of sampled time points.
%                 labelling - string for labelling plots. 
%initial_values - Vector of initial adipose tissue metabolite concentrations.
%                 [adipose glycerol, adipose G3P, adipose NEFA];
%fig_yn         - Plot figure of model predictions. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
GLY_delay=parameters(6);
parameters_GLY=[pgly,B_ATL,ATL_max,K_ATL,GLY_delay];
input_data.parameters_ATL=parameters_GLY;
%parameters for NEFA
re_ester=parameters(11);
pnefa=parameters(12);
G3P_delay=parameters(13);
glucose_use=parameters(14);
parameters_nefa=[re_ester,pnefa,G3P_delay,glucose_use];
input_data.parameters_nefa=parameters_nefa;

input_data.initial_values=initial_val;

%generate splines of dependent inputs.
TG_art_ss=spline([-60,-45,-37,-30,0,60,120,180,240,300],[input_data.TG_art(1),input_data.TG_art(1),input_data.TG_art(1),input_data.TG_art],-60:1:300);
TG_art=TG_art_ss(31:end);
input_data.TG_art_s=TG_art;
NEFA_art_ss=spline([-60,-45,-37,-30,0,60,120,180,240,300],[input_data.NEFA_art(1),input_data.NEFA_art(1),input_data.NEFA_art(1),input_data.NEFA_art],-60:1:300);
NEFA_art=NEFA_art_ss(31:end);
input_data.NEFA_art_s=NEFA_art;
G_art_ss=spline([-60,-45,-37,-30,0,60,120,180,240,300],[input_data.G_art(1),input_data.G_art(1),input_data.G_art(1),input_data.G_art],-60:1:300);
G_art=G_art_ss(31:end);
input_data.G_art_s=G_art;
GLY_art_ss=spline([-60,-45,-37,-30,0,60,120,180,240,300],[input_data.GLY_art(1),input_data.GLY_art(1),input_data.GLY_art(1),input_data.GLY_art],-60:1:300);
GLY_art=GLY_art_ss(31:end);
input_data.GLY_art_s=GLY_art;
I_PL_ss=spline([-60,-45,-37,input_data.t],[input_data.I(1),input_data.I(1),input_data.I(1),input_data.I],-60:1:300);
I_art=I_PL_ss(31:361);
input_data.I_art_s=I_art;

%generate vector of three fold LPL insulin delay
input_data.delay=LPL_delay;
I_pl=input_data.I;
initial_I=[I_pl(1),I_pl(1),I_pl(1)];
[T,I_delay]=ode15s(@insulin_delay,time,initial_I,[],input_data);
I_LPL=I_delay(:,3)';
input_data.I_LPL=I_delay(:,3)';

%generate vector of three fold  AT insulin delay
input_data.delay=AT_delay;
[T,I_delay]=ode15s(@insulin_delay,time,initial_I,[],input_data);
I_AT=I_delay(:,3)';
input_data.I_AT=I_delay(:,3)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRIGLYCERIDE FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AV-difference in TG due to LPL lipolysis
D_TG=-K_ad.* I_LPL.* TG_art;

out.TG_flux=D_TG;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPILLOVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set basal insulin
I_b=input_data.I(1);
%computed fractional spill over
spill=(D_spill.*(I_b./I_art))/100;

out.frac_spillover=spill;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GLUCOSE FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%equation for glucose flux
D_G=-GLUT_1.*G_art-GLUT_4.*I_AT.*G_art;

out.G_flux=D_G;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GLYCEROL FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set initial value for adipose glycerol concentration - not measured
initial_GLY=input_data.initial_values(1);
%calculate the time series of arterial glycerol concentration given
%parameters suppllied in the input of function (parameters_GLY)
[T,G]=ode15s(@Glycerol_ODE,time,initial_GLY,[],input_data,parameters_GLY);
gly_at=G';
%equation for glycerol flux
D_GLY=K_ad.*(TG_art).*(I_LPL)+pgly.*(gly_at-(GLY_art+K_ad.*TG_art.*I_LPL));

out.GLY_flux=D_GLY;
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
NEFA_at=AT_NEFA(:,2)';

%equation for NEFA flux
D_NEFA= 3.*((D_spill.*(I_b./I_art))./100).*(K_ad.*TG_art.*I_LPL)+pnefa.*(NEFA_at-(NEFA_art+3.*((D_spill.*(I_b./I_art))./100).*(K_ad.*TG_art.*I_LPL)));

out.NEFA_flux=D_NEFA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if fig_yn==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate figure of model predictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure()
    subplot(2,3,1)
    hold on;
    plot(time,-1.*D_TG,'LineWidth',2);
    legend('model prediction');
    xlabel('Time (mins)')
    ylabel('Triglyceride flux (umol/100 ml tissue/min)')
    t_mess=['Triglyceride Flux : ',input_data.labeling];
    title(t_mess);
    xlim([input_data.t(1)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t;
    
    subplot(2,3,2)
    plot(60:1:300,spill(91:331),'LineWidth',2)
    xlabel('time (mins)')
    ylabel('fractional spillover (%)')
    t_mess=['Fractional Spillover : ',input_data.labeling];
    title(t_mess);
    xlim([input_data.t(1)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t;
    
    subplot(2,3,3)
    plot(time,-1.*D_G,'LineWidth',2);
    xlabel('time (mins)')
    ylabel('Glucose flux (mmol/ 100 ml tissue/min)')
    t_mess=['Glucose flux : ',input_data.labeling];
    title(t_mess);
    xlim([input_data.t(1)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t;

    subplot(2,3,4)
    plot(time,-1.*D_GLY,'LineWidth',2)
    xlabel('time (mins)')
    ylabel('glycerol flux (umol/100 ml tissue/min)')
    t_mess=['Glycerol flux : ',input_data.labeling];
    title(t_mess)
    xlim([input_data.t(1)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t;
    
    subplot(2,3,5)
    plot(time,-1.*D_NEFA,'LineWidth',2);
    xlabel('time (mins)')
    ylabel('NEFA flux (umol/ 100 ml tissue/ min)')
    t_mess=['NEFA flux : ',input_data.labeling];
    title(t_mess)
    xlim([input_data.t(1)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t;
 
end



