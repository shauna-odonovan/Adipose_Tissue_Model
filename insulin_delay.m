function dIdt=insulin_delay(time,initial_values,input_data)
%computes three fold insulin delay for given time delay parameter in input
%data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input data required
%
%time - time vector for ODE solver
%
%initial_values - initial values for each insulin delay compartment.
%                 [compartment_1 value,compartment_2_value,compartment_3 value]
%
%input_data -  struct which must contain following arrays
%              parameters_LPL [K_ad,LPL_delay]
%              I        -  vector of arterial insulin concentrations
%              delay    -  delay parameter to be esimated
%              t        - vector of measured time points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate spline of arterial insulin concentration. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spline_time=[-60,-45,-37,input_data.t];
spline_insulin=[input_data.I(1),input_data.I(1),input_data.I(1),input_data.I];
I_ss=spline(spline_time,spline_insulin,time);
I_pl=I_ss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set parameter and initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau=input_data.delay;
I_11=initial_values(1);
I_21=initial_values(2);
I_delay=initial_values(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d.I_11=(3./tau).*(I_pl-I_11);
d.I_21=(3./tau).*(I_11-I_21);
d.I_delay=(3./tau).*(I_21-I_delay);

dIdt=[d.I_11;d.I_21;d.I_delay];