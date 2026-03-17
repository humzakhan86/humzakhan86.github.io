%% MAE434 - Aircraft Design - Optimization Script
%Humza Khan
% Team: Charlie
% Date: 3/11/2025
%--------------------------------------------------------------------------

%% Project Script Starts Here
clear all; close all; clc;
%--------------------------------------------------------------------------

%% Input Variables
alt_ft = 40000; %ft (cruise altitude)
alt_m = convlength(alt_ft,'ft','m'); %Cruise altitude in meters
[alt_T, alt_a, P, alt_rho] = atmosisa(alt_m);
rho_imperial = alt_rho/515.4;

%velocity calculations
V_cruise = 515; %ft/s
alt_a_fts = alt_a*3.281; %ft/s
M_cruise = V_cruise/alt_a_fts; %mach number
fprintf('<strong>Initial Velocity Calculations :)</strong>\n')
fprintf('Cruise Velocity: <strong>%5.2f</strong> ft/s\n',V_cruise)
fprintf('Speed of Sound at 40,000 ft: <strong>%5.2f</strong> ft/s\n',alt_a_fts)
fprintf('Cruise Mach Number: <strong>%5.2f</strong>\n',M_cruise)
fprintf('--------------------------------------------------------------\n')

%Plane/Aerodynamics/Payload specs
w_pl = 185*18; %lbs (taken from the De Havilland Aircraft - Twin Otter)
w_crew = 185*2; %lbs (avg pssnger mass = 185, we plan to carry 18 passangers)
AR = 8; %Estimated Aspect Ratio
swet_sref = 4.34; % 6.2 and & 5 %CAD Reference surface area ratio
K_ld = 15.5; %Civilian Jets ---> Same Page as table 3.6 Aircraft Design - Raymer
fprintf('<strong>Engine Aircraft/Aerodynamics/Payload Specifications :)</strong>\n')
fprintf('Payload Weight: <strong>%5.2f</strong> lb\n',w_pl)
fprintf('Crew weight: <strong>%5.2f</strong> lb\n',w_crew)
fprintf('Aspect Ratio: <strong>%5.2f</strong>\n',AR)
fprintf('Reference Surface Area Ratio: <strong>%5.2f</strong>\n',swet_sref)
fprintf('--------------------------------------------------------------\n')

%range and loitering
Range = 3000 * 6076; %we're planning to fly around 3000 Nautical Miles
L_time = 45*60; %45 min loiter time, converted to seconds
endurance = Range/V_cruise;
fprintf('<strong>Range & Loiter Parameters :)</strong>\n')
fprintf('Range: <strong>%5.2f</strong> ft\n',Range)
fprintf('Loiter time: <strong>%5.2f</strong> seconds\n',L_time)
fprintf('Endurance: <strong>%5.2f</strong> seconds\n',endurance)
fprintf('--------------------------------------------------------------\n')

%Engine Properties - High-bypass turbofan
C_cruise = 0.5/3600; %1/s %SFC ----> Table 3.3 Aircraft Design - Raymer
C_loiter = 0.4/3600; %1/s %SFC ----> Table 3.3 Aircraft Design - Raymer
fprintf('<strong>Engine Properties - High-bypass turbofan :)</strong>\n')
fprintf('C_loiter: <strong>%e</strong>\n',C_loiter)
fprintf('C_cruise: <strong>%e</strong>\n',C_cruise)
fprintf('--------------------------------------------------------------\n')


%Solving for characteristics
SAM_area = 894+894+59+66+165+373+373+461+297+297+313+489+100; %updated wetted area from CAD
s_wet_body = SAM_area;
s_ref_sam = 856.8075;
swet_sref_updated = s_wet_body/s_ref_sam;
L_D_max = K_ld*sqrt(AR/swet_sref_updated); %loiter l/d ratio
L_D_cruise = 0.866*L_D_max; %cruise l/d ratio for jet
fprintf('<strong> L/D Calculation :)</strong>\n')
fprintf('L/D (MAX): <strong>%5.2f</strong>\n',L_D_max)
fprintf('L/D (CRUISE): <strong>%5.2f</strong>\n',L_D_cruise)
fprintf('--------------------------------------------------------------\n')


%5-stage fuel characteristics
w_10 = 0.97; %Take off
w_21 = 0.985; %Climb
w_32 = exp(-(Range*C_cruise)/(V_cruise*L_D_cruise)); %cruise
w_43 = exp(-(L_time*C_loiter)/(L_D_max)); %Loiter
w_54 = 0.995; %landing (Raymer)
w_50 = w_10*w_21*w_32*w_43*w_54; %total weight decrease ratio
w_f0 = 1.06*(1-w_50); %fuel to weigth ratio -----> 6% residual weight allowance

fprintf('<strong>5-Stage Fuel & Weight Fraction Calculation :)</strong>\n')
fprintf('weight fraction during takeoff: <strong>%5.4f</strong>\n',w_10)
fprintf('weight fraction during climb: <strong>%5.4f</strong>\n',w_21)
fprintf('weight fraction during cruise: <strong>%5.4f</strong>\n',w_32)
fprintf('weight fraction during loiter: <strong>%5.4f</strong>\n',w_43)
fprintf('weight fraction during landing: <strong>%5.4f</strong>\n',w_54)
fprintf('weight fraction (total weight decrease ratio): <strong>%5.4f</strong>\n',w_50)
fprintf('weight fraction (fuel with 6-percent residual weight allowance): <strong>%5.4f</strong>\n',w_f0)
fprintf('--------------------------------------------------------------\n')
%--------------------------------------------------------------------------

%% determining takeoff weight (Sensitivity analysis)
a = 0.32; %table 6.1
b = 0.66; %table 6.1
C1 = -0.13;
C2 = 0.3;
C3 = 0.06;
C4 = -0.05;
C5 = 0.05;
K_vs = 1.0; %constant sweep
guess_to_weight = 12500; %lb %initial weight guess from the Twin Otter
w_drop = 0;
rhs = 1; 
t_w_ratio = 1.18; %can be 1.18 depending on the iteration we use
w_os = 45; %wingloading
M_max = M_cruise; %recalling mach number
tolerance = abs((1 - (guess_to_weight/rhs))*100); %tolerance condition
iteration = 0;

%while loop is going to iterate a solution for the empty weight of the
%aircraft
while tolerance > 0.01
    guess_to_weight = guess_to_weight + 0.1; %increasing weight by 5lb
    w_e0 = (a + ((b*guess_to_weight^C1) * (AR^C2) * (t_w_ratio^C3) * (w_os^C4) * (M_max^C5)))*K_vs; %estimating empty weight ratio based on the initial guess
    rhs = w_pl + w_crew + w_drop + (guess_to_weight*w_f0) + w_e0*guess_to_weight; %updating the right hand of the equation
    tolerance = abs((1 - (guess_to_weight/rhs))*100); %tolerance condition updated
    iteration = iteration + 1;
end

%weight calculations (empty and fuel)
w_fuel = rhs*w_f0;
w_e = rhs - w_fuel;

fprintf('<strong>Weight Calculation Results below :)</strong>\n')
fprintf('The Total Take-off Weight of the aircraft is estimated at: <strong>%5.2f</strong> lbs\n',rhs)
fprintf('The Fuel Weight of the aircraft is estimated at: <strong>%5.2f</strong> lbs\n',w_fuel)
fprintf('The Empty Weight of the aircraft is estimated at: <strong>%5.2f</strong> lbs\n',w_e)
fprintf('--------------------------------------------------------------\n')
%--------------------------------------------------------------------------

%% Wing Design and Parameters
WS = w_os; %lb/ft^2 %table 5.5 transport jet value (we had to reduce it to reduce the take off distance and also that 120 from the book was not reasonable as if was of a full-scale bomber)
taper_ratio = 0.8; %estimated from the CAD Design
s_ref = rhs/WS; %reference area
wing_span = sqrt(AR*s_ref); %wingspan
%root_chord = (2*s_ref)/(wing_span*(1-taper_ratio)); %optimal root chord
%tip_chord = taper_ratio*root_chord; % optimal tip chord
%mean_chord = ((2/3)*root_chord*(1+taper_ratio+taper_ratio^2))/(1+taper_ratio);
mean_chord = s_ref/wing_span;
wing_area = s_ref;
cl_cruise = (2*guess_to_weight)/(rho_imperial*V_cruise^2*wing_area);
s_wet = swet_sref_updated*s_ref; %wing wetted surface area
fprintf('<strong>Wing Area Calculation ASSUMING Rectangular Wings  Results below :)</strong>\n')
fprintf('Estimated Coefficient of Lift (cruise): <strong>%8.4f</strong>\n',cl_cruise)
fprintf('The reference area of the wings is estimated at: <strong>%8.4f</strong> ft^2\n',s_ref)
fprintf('Estimated Wetted Surface Area of the wing: <strong>%8.4f</strong> ft^2\n',s_wet)
fprintf('Estimated Wetted Surface Area of the body: <strong>%8.4f</strong> ft^2\n',s_wet_body)
fprintf('Estimated wing span of the aircraft is: <strong>%8.4f</strong> ft\n',wing_span)
fprintf('Estimated mean aerodymamic chord of the wings is: <strong>%8.4f</strong> ft\n',mean_chord)
fprintf('Estimated Taper Ratio (Not Applied) of the wings is: <strong>%8.4f</strong>\n',taper_ratio)
fprintf('--------------------------------------------------------------\n')
%--------------------------------------------------------------------------


%% Take-off distance calculator
%%SG - initial take off roll, before rotation
alt_landing = 0;
[alt_T_land, alt_a_land, P_land, alt_rho_land] = atmosisa(alt_landing);
rho_imperial_landing = alt_rho_land/515.4;
TW = 1.187; %minimizing takeoff distance
CL_to = 1.5;
CL_max = 2.5; %jet flaps/slats/flaps/
CD0 = 0.03;
e = 0.8;
K = 1/(pi*e*AR); %K-factor from Raymer takeoff distance eqns

vstall = sqrt(2*WS/(rho_imperial_landing*CL_max)); %stall speed calculation


pseudo_speed = sqrt(2*WS/(rho_imperial_landing*CL_to));
Vto = 1.2*pseudo_speed; %take off speed calculation

%from raymer pg 673

mu = 0.05; %rolling resistance
KT = TW - mu;
KA = rho_imperial_landing/((2*WS))*(mu*(CL_to)-(CD0)-(K*CL_to));

SG = (1/(2*32.2*KA))*log((KT+KA*Vto^2)/(KT));
%--------------------------------------------------------------------------
TakeOffDistance = (Vto^2)/(2*32.2*(23000/rhs - mu - 0.01));%SG*1.3... 0.01 is the D/W value
fprintf('<strong>Takeoff distance calculation :)</strong>\n')
fprintf('Estimated stall speed: <strong>%8.4f</strong> ft/s\n',vstall)
fprintf('Estimated takeoff speed: <strong>%8.4f</strong> ft/s\n',Vto)
fprintf('Estimated takeoff distance: <strong>%8.4f</strong> ft\n',TakeOffDistance)
fprintf('--------------------------------------------------------------\n')
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Exam 2 starts here
%--------------------------------------------------------------------------
%% Thrust to Weight Ratio
arbi_time = endurance/2; %arbitrary half time where the plane is guaranteed be in cruise
t_w_cruise = 1/(L_D_cruise); %thrust to weight ratio during cruise
cruise_wt = w_10*guess_to_weight; %weight of aircraft during cruise
t_cruise = t_w_cruise*cruise_wt; %estimated thrust required during cruise
t_w_takeoff = t_w_cruise*w_10*(23000/t_cruise); %math checks out
t_takeoff = t_w_takeoff*rhs; %take off thrust
fprintf('<strong>Thrust-to-weight Ratio (Takeoff) Calculation :)</strong>\n')
fprintf('Thrust to Weight Ratio (takeoff): <strong>%8.4f</strong>\n',t_w_takeoff)
fprintf('Thrust at cruise: <strong>%8.4f</strong> lbs\n',t_cruise)
fprintf('Weight during cruise: <strong>%8.4f</strong> lbs\n',cruise_wt)
fprintf('Thrust at Takeoff: <strong>%8.4f</strong> lbs\n',t_takeoff)
fprintf('--------------------------------------------------------------\n')
fprintf('<strong>Script Ends Here For Now </strong>\n')
%--------------------------------------------------------------------------

%% Wing Loading
%--------------------------------------------------------------------------
%% Script Ends here as of now.