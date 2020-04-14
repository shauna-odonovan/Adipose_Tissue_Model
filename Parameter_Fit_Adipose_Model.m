function out=Parameter_Fit_Adipose_Model(input_data,initial_p,time,lb,ub,initial_values,fig_yn)

%initial adipose tissue concentrations of glycerol, G-3-P, and NEFA are
%includeded in the input_data structure for ease of use in other functions.
input_data.initial_values=initial_values;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate spline of measured arterial metabolite concenrations for model
%prediciton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-For effcieny this is done once outside of the lsqnonline step and
% data splines are stored in input_data for use by other functions. 
%-It is necessary to enforce a steady state assumption prior to consumption of the meal
% when constructing the cubic splines to avoid the prdiction of erronous
% curves in the fasting state. This is schived simply by introducing
% additional points before the measured fasting data with the same values at
% time point 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
I_PL=I_PL_ss(31:361);
input_data.I_art_s=I_PL;


%define options for lsqnonline algorithm
lsq_options=optimset('Algorithm','trust-region-reflective','MaxFunEvals',500,'TolX',1e-25);
%run lsqnonlin
[p_opt,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(@Model_Error,initial_p,lb,ub,lsq_options,input_data,time);
%define out put structure
out.p_opt=p_opt;
out.resnorm=resnorm;
out.residual=residual;
out.exitflag=exitflag;
out.output=output;
out.lambda=lambda;
out.jacobian=jacobian;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate figure of model fitting if figure==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig_yn==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %generate data for plotting based on p_opt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_data.initial_GLY=initial_values(1);
    %parameters for LPL lipolysis
    K_ad=p_opt(1);
    LPL_delay=p_opt(2);
    input_data.parameters_LPL=[K_ad,LPL_delay];
    %parameters for spillover
    D_spill=p_opt(3);
    input_data.parameters_spill=D_spill;
    %parameters for glucose flux
    GLUT_1=p_opt(4);
    GLUT_4=p_opt(5);
    AT_delay=p_opt(6);
    input_data.parameters_glu=[GLUT_1,GLUT_4,AT_delay];
    %parameters for glycerol flux
    pgly=p_opt(7);
    B_ATL=p_opt(8);
    ATL_max=p_opt(9);
    K_ATL=p_opt(10);
    parameters_GLY=[pgly,B_ATL,ATL_max,K_ATL];
    input_data.parameters_ATL=parameters_GLY;
    %parameters for NEFA
    re_ester=p_opt(11);
    pnefa=p_opt(12);
    G3P_delay=p_opt(13);
    glucose_use=p_opt(14);
    parameters_nefa=[re_ester,pnefa,G3P_delay,glucose_use];
    input_data.parameters_nefa=parameters_nefa;
    
    
    %initial insulin concentrations in remote compartments
    I_pl=input_data.I;
    initial_I=[I_pl(1),I_pl(1),I_pl(1)];

    %generate vector of three fold delayed LPL delay
    input_data.delay=LPL_delay;
    [T,I_delay]=ode15s(@insulin_delay,time,initial_I,[],input_data);
    I_LPL=I_delay(:,3)';
    input_data.I_LPL=I_LPL;
    %define insulin delay for adipose compartment
    input_data.delay=AT_delay;
    [T,I_delay]=ode15s(@insulin_delay,time,initial_I,[],input_data);
    I_AT=I_delay(:,3)';
    input_data.I_AT=I_AT;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TRIGLYCERIDE FLUX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %AV-difference in TG due to LPL lipolysis
    D_TG=-K_ad.* I_LPL.* TG_art;
    %measured TG flux
    A_V_TG=-1.*input_data.mean_TG_flux;
    %error term SSE/sigma^2
    error_LPL=sum(((A_V_TG-D_TG([1,31,91,151,211,271,331]))./input_data.std_TG_flux).^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SPILLOVER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set basal insulin
    I_b=input_data.I(1);
    %computed fractional spill over
    spill=(D_spill.*(I_b./I_PL))/100;
    %compute error for spill over term
    error_spill=sum(((spill([91,151,211,271,331])-input_data.mean_spill(1:5))./input_data.std_spill).^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %GLUCOSE FLUX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %equation for glucose flux
    D_G=-GLUT_1.*G_art-GLUT_4.*I_AT.*G_art;
    %arteial venous difference in glucose
    A_V_G=-1.*input_data.mean_G_flux;
    %compute error for glucose term
    error_glucose=sum(((A_V_G-D_G([1,31,91,151,211,271,331]))./input_data.std_G_flux).^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %GLYCEROL FLUX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    initial_GLY=initial_values(1);
    [T,G]=ode15s(@Glycerol_ODE,time,initial_GLY,[],input_data,parameters_GLY);
    gly_at=G';
    %generate vectos of insulin concentrations and arterial TG as spline of
    %data
    AV_GLY=K_ad.*(TG_art).*(I_LPL)+pgly.*(gly_at-(GLY_art+K_ad.*TG_art.*I_LPL));
 
    %compute error for glycerol term
    measured_AV_GLY=input_data.mean_GLY_flux;
    error_GLY=sum((((-1*measured_AV_GLY)-AV_GLY([1,31,91,151,211,271,331]))./input_data.std_GLY_flux).^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %NEFA FLUX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    initial_values_NEFA=input_data.initial_values(2:3);
    %generate G-3-P concentration as a delay of gluocse flux
    initial_G3P=[initial_values_NEFA(1),initial_values_NEFA(1)];
    
    [T,G_g3p]=ode15s(@Glucose_G3P,time,initial_G3P,[],input_data);
    input_data.G3P=G_g3p(:,2)';


    %Generate arterial NEFA
    [T,AT_NEFA]=ode15s(@Glucose_G3P_NEFA_ODE,time,initial_values_NEFA,[],input_data,parameters_nefa);
    NEFA_at=AT_NEFA(:,2)';
    AV_NEFA= 3.*((D_spill.*(I_b./I_PL))./100).*(K_ad.*TG_art.*I_LPL)+pnefa.*(NEFA_at-(NEFA_art+3.*((D_spill.*(I_b./I_PL))./100).*(K_ad.*TG_art.*I_LPL)));
    %compute error for NEFA term
    measured_AV_NEFA=-1.*input_data.mean_NEFA_flux;
    error_nefa=sum(((measured_AV_NEFA-AV_NEFA([1,31,91,151,211,271,331]))./input_data.std_NEFA_flux).^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save model predicted fluxes in function output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %out.TG_model=D_TG;
    %out.Spill_model=spill;
    %out.GLU_model=D_G;
    %out.GLY_model=AV_GLY;   
    %out.NEFA_model=AV_NEFA;
    out.error=AV_GLY(1)-input_data.mean_GLY_flux(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure()
    title('lsqnonlin');
    subplot(2,3,1)
    message=['C(p) = ',num2str(error_LPL)];
    errorbar(input_data.t,input_data.mean_TG_flux,input_data.std_TG_flux./4,'rx','LineWidth',1.5);
    hold on;
    plot(time,-1.*D_TG,'LineWidth',2);
    refline(0,0)
    legend('mean data','model prediction');
    text(120,input_data.mean_TG_flux(1),message)
    xlabel('Time (mins)')
    ylabel('Triglyceride flux (umol/100 ml tissue/min)')
    t_mess=['Triglyceride Flux : ',input_data.labeling];
    title(t_mess);
    xlim([input_data.t(1)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t;
    
    subplot(2,3,2)
    errorbar([60,120,180,240,300],input_data.mean_spill,input_data.std_spill,'rx','LineWidth',1.5)
    hold on;
    plot(60:1:300,spill(91:331),'LineWidth',2)
    refline(0,0)
    message=['C(p) =',num2str(error_spill)];
    text(120,input_data.mean_spill(1),message)
    xlabel('time (mins)')
    ylabel('fractional spillover (%)')
    legend('mean data','model prediction');
    t_mess=['Fractional Spillover : ',input_data.labeling];
    title(t_mess);
    xlim([input_data.t(3)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t(3:end);
    
    subplot(2,3,3)
    message=['C(p) = ',num2str(error_glucose)];
    errorbar(input_data.t,input_data.mean_G_flux,input_data.std_G_flux,'rx','LineWidth',1.5);
    hold on;
    plot(time,-1.*D_G,'LineWidth',2);
    refline(0,0)
    text(120,input_data.mean_G_flux(1),message)
    xlabel('time (mins)')
    ylabel('Glucose flux (mmol/ 100 ml tissue/min)')
    t_mess=['Glucose flux : ',input_data.labeling];
    title(t_mess);
    xlim([input_data.t(1)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t;

    subplot(2,3,4)
    errorbar(input_data.t,(input_data.mean_GLY_flux),input_data.std_GLY_flux,'rx','LineWidth',1.5);
    hold on;
    plot(time,-1.*AV_GLY,'LineWidth',2)
    refline(0,0)
    xlabel('time (mins)')
    ylabel('glycerol flux (umol/100 ml tissue/min)')
    t_mess=['Glycerol flux : ',input_data.labeling];
    title(t_mess)
    message=['C(p) = ',num2str(error_GLY)];
    text(120,input_data.mean_GLY_flux(1),message);
    xlim([input_data.t(1)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t;
    
    subplot(2,3,5)
    errorbar(input_data.t,input_data.mean_NEFA_flux,input_data.std_NEFA_flux,'rx','LineWidth',1.5);
    hold on;
    plot(time,-1.*AV_NEFA,'LineWidth',2);
    refline(0,0)
    xlabel('time (mins)')
    ylabel('NEFA flux (umol/ 100 ml tissue/ min)')
    t_mess=['NEFA flux : ',input_data.labeling];
    title(t_mess)
    message=['C(p)= ',num2str(error_nefa)];
    text(120,input_data.mean_NEFA_flux(1),message);
    xlim([input_data.t(1)-20,input_data.t(end)+20]);
    ax=gca;
    ax.XTick=input_data.t;
end
total_error=error_nefa+error_GLY+error_glucose+error_spill+error_LPL;
out.error=total_error;
