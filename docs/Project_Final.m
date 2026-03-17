%% Initial Orbits Information and Calculations
close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% General Info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu= 398600.4418; % km^3/s^2
R_Earth=6378.137; % km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parking Orbit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_park= 7028; % km, semimajor axis
e_park=0.0; % unitless, eccentricity
i_park=deg2rad(10); % radians, inclination
w_park=deg2rad(0); % radians, argument of perigee
O_park=deg2rad(30); % radians, right acension of acending node
nu_park=deg2rad(0); % radians, true anomaly at time=0

r_park= 7028; % km, constant due to circular
v_park= sqrt(mu*((2/a_park)-(1/a_park))); % km/s, circular orbit
SemiR_park= a_park*(1-e_park^2); % km, semilatis rectum

ApN_park=[ cos(O_park)*cos(w_park)-sin(O_park)*sin(w_park)*cos(i_park), -cos(O_park)*sin(w_park)-sin(O_park)*cos(w_park)*cos(i_park),  sin(O_park)*sin(i_park) ;...
           sin(O_park)*cos(w_park)+cos(O_park)*sin(w_park)*cos(i_park), -sin(O_park)*sin(w_park)+cos(O_park)*cos(w_park)*cos(i_park), -cos(O_park)*sin(i_park) ;...
           sin(w_park)*sin(i_park)                                    ,  cos(w_park)*sin(i_park)                                    ,  cos(i_park)            ];
r_P_park_t0= [r_park*cos(nu_park); r_park*sin(nu_park); 0];
r_N_park_t0= ApN_park*r_P_park_t0;

v_P_park_t0= [-sqrt(mu/SemiR_park)*sin(nu_park); sqrt(mu/SemiR_park)*(e_park+cos(nu_park)); 0];
v_N_park_t0= ApN_park*v_P_park_t0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Target Orbit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_tgt= 42164.0; % km
e_tgt= 0.0;
i_tgt=deg2rad(0); % radians
nu_tgt_t0=deg2rad(180); % radians, true anomaly at time=0
P_tgt   = (2*pi*a_tgt   ^(3/2))/sqrt(mu); % seconds
r_tgt= a_tgt; % km, constant becasue circular
v_tgt= sqrt(mu*((2/a_tgt)-(1/a_tgt))); % km/s, circular orbit
SemiR_tgt= a_tgt*(1-e_tgt^2); % km, semilatis rectum

w_tgt= v_tgt/r_tgt; % Angular Velocity

ApN_tgt=[ 1, 0, 0 ;...
          0, 1, 0 ;...
          0, 0, 1];
r_P_tgt_t0= [r_tgt*cos(nu_tgt_t0); r_tgt*sin(nu_tgt_t0); 0];
r_N_tgt_t0=ApN_tgt*r_P_tgt_t0;

v_P_tgt_t0= [-sqrt(mu/SemiR_tgt)*sin(nu_tgt_t0); sqrt(mu/SemiR_tgt)*(e_tgt+cos(nu_tgt_t0)); 0];
v_N_tgt_t0= ApN_tgt*v_P_tgt_t0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Theoretical Orbit 1 - Compares parking to target orbit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_per_theo1= r_park; % km
r_apo_theo1= r_tgt; % km
a_theo1= .5*(r_park+r_tgt); % km

v_per_theo1= sqrt(mu*((2/r_per_theo1)-(1/a_theo1))); % km/s
v_apo_theo1= sqrt(mu*((2/r_apo_theo1)-(1/a_theo1))); % km/s

delta_v_maneuvor1= v_per_theo1 - v_park;
delta_v_maneuvor3= v_tgt - v_apo_theo1;

%% Maneuvor 1 Orbits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 1 - Maneuvor 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_per_trans1= r_park; % km
r_N_per_trans1= r_N_park_t0; % We do the burn here instantaniously at time=0, or possible k parking orbit periods after t0, same position reguardless
v_P_per_trans1= v_P_park_t0 + [0;.468;0]; % Delta V from single burn to determine new speed at periapsis, km/s
ApN_trans1=ApN_park; % This is absolutely true, I dont have to justify it, I am a genius
v_N_per_trans1= ApN_trans1*v_P_per_trans1; % km/s

v_per_trans1=norm(v_N_per_trans1); % km/s
h_trans1= r_per_trans1*v_per_trans1; % angular momentum
SemiR_trans1= (h_trans1^2)/mu; % km

syms a e
eq1= SemiR_trans1 == a*(1-e^2);
eq2= r_per_trans1 == a*(1-e);
sol_trans1=solve(eq1,eq2,a,e);
a_trans1= double(sol_trans1.a);
e_trans1= double(sol_trans1.e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 2 - Maneuvor 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_per_trans2= r_per_trans1; % km
r_N_per_trans2= r_N_per_trans1; 
v_P_per_trans2= v_P_per_trans1 + [0;.468;0]; % Delta V from single burn to determine new speed at periapsis, km/s
ApN_trans2=ApN_trans1; % This is absolutely true, I dont have to justify it, I am a genius
v_N_per_trans2= ApN_trans2*v_P_per_trans2; % km/s

v_per_trans2=norm(v_N_per_trans2); % km/s
h_trans2= r_per_trans2*v_per_trans2; % angular momentum
SemiR_trans2= (h_trans2^2)/mu; % km

syms a e
eq1= SemiR_trans2 == a*(1-e^2);
eq2= r_per_trans2 == a*(1-e);
sol_trans2=solve(eq1,eq2,a,e);
a_trans2= double(sol_trans2.a);
e_trans2= double(sol_trans2.e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 3 - Maneuvor 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_per_trans3= r_per_trans2; % km
r_N_per_trans3= r_N_per_trans2; 
v_P_per_trans3= v_P_per_trans2 + [0;.468;0]; % Delta V from single burn to determine new speed at periapsis, km/s
ApN_trans3=ApN_trans2; % This is absolutely true, I dont have to justify it, I am a genius
v_N_per_trans3= ApN_trans3*v_P_per_trans3; % km/s

v_per_trans3=norm(v_N_per_trans3); % km/s
h_trans3= r_per_trans3*v_per_trans3; % angular momentum
SemiR_trans3= (h_trans3^2)/mu; % km

syms a e 
eq1= SemiR_trans3 == a*(1-e^2);
eq2= r_per_trans3 == a*(1-e);
sol_trans3=solve(eq1,eq2,a,e);
a_trans3= double(sol_trans3.a);
e_trans3= double(sol_trans3.e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 4 - Maneuvor 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_per_trans4= r_per_trans3; % km
r_N_per_trans4= r_N_per_trans3; 
v_P_per_trans4= v_P_per_trans3 + [0;.468;0]; % Delta V from single burn to determine new speed at periapsis, km/s
ApN_trans4=ApN_trans3; % This is absolutely true, I dont have to justify it, I am a genius
v_N_per_trans4= ApN_trans4*v_P_per_trans4; % km/s

v_per_trans4=norm(v_N_per_trans4); % km/s
h_trans4= r_per_trans4*v_per_trans4; % angular momentum
SemiR_trans4= (h_trans4^2)/mu; % km

syms a e 
eq1= SemiR_trans4 == a*(1-e^2);
eq2= r_per_trans4 == a*(1-e);
sol_trans4=solve(eq1,eq2,a,e);
a_trans4= double(sol_trans4.a);
e_trans4= double(sol_trans4.e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 5 - Maneuvor 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_per_trans5= r_per_trans4; % km
r_N_per_trans5= r_N_per_trans4; 
v_P_per_trans5= v_P_per_trans4 + [0; 0.457328723782171;0]; % Delta V from single burn to determine new speed at periapsis, km/s
ApN_trans5=ApN_trans4; % This is absolutely true, I dont have to justify it, I am a genius
v_N_per_trans5= ApN_trans5*v_P_per_trans5; % km/s

v_per_trans5=norm(v_N_per_trans5); % km/s
h_trans5= r_per_trans5*v_per_trans5; % angular momentum
SemiR_trans5= (h_trans5^2)/mu; % km

syms a e 
eq1= SemiR_trans5 == a*(1-e^2);
eq2= r_per_trans5 == a*(1-e);
sol_trans5=solve(eq1,eq2,a,e);
a_trans5= double(sol_trans5.a);
e_trans5= double(sol_trans5.e);


r_apo_trans5= a_trans5*(1 + e_trans5); % km
v_apo_trans5=  sqrt(mu*((2/r_apo_trans5)-(1/a_trans5))); % km/s

nu_apo_trans5=deg2rad(180); % radians

r_P_apo_trans5= [r_apo_trans5*cos(nu_apo_trans5); r_apo_trans5*sin(nu_apo_trans5); 0];
r_N_apo_trans5= ApN_trans5*r_P_apo_trans5;


v_P_apo_trans5= [-sqrt(mu/SemiR_trans5)*sin(nu_apo_trans5); sqrt(mu/SemiR_trans5)*(e_trans5+cos(nu_apo_trans5)); 0];
v_N_apo_trans5= ApN_trans5*v_P_apo_trans5;

delta_v_maneuvor2= 2*v_apo_trans5*(h_trans5/(r_apo_trans5*v_apo_trans5))*sin(i_park/2); % Required delta v for inclination change

h_N_trans5= cross(r_N_per_trans5,v_N_per_trans5);
h_N_trans5_unit= h_N_trans5/norm(h_N_trans5);

%% Maneuvor 2 Orbits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Theoretical Orbit 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_apo_theo2= r_apo_trans5;
r_N_apo_theo2= r_N_apo_trans5;

ApN_theo2=[ cos(O_park)*cos(w_park)-sin(O_park)*sin(w_park)*cos(0), -cos(O_park)*sin(w_park)-sin(O_park)*cos(w_park)*cos(0),  sin(O_park)*sin(0) ;...
           sin(O_park)*cos(w_park)+cos(O_park)*sin(w_park)*cos(0), -sin(O_park)*sin(w_park)+cos(O_park)*cos(w_park)*cos(0), -cos(O_park)*sin(0) ;...
           sin(w_park)*sin(0)                                    ,  cos(w_park)*sin(0)                                    ,  cos(0)            ];

v_N_apo_theo2= ApN_theo2*v_P_apo_trans5;
a_theo2= a_trans5;

v_apo_theo2=norm(v_N_apo_theo2);


delta_v_maneuvor2_N= v_N_apo_trans5 - v_N_apo_theo2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 6 - Maneuvor 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_apo_trans6= r_apo_trans5; % km
r_N_apo_trans6= r_N_apo_trans5; % km
ApN_trans6= ApN_theo2;

pro_Brun_trans6= .177435223300000;

v_N_apo_trans6= v_N_apo_trans5 + -1*[-0.0124845713703247; 0.0216239119241222; -0.285398607480736] + ApN_trans6*-1*[0; pro_Brun_trans6;0]; % Delta V from single burn, fixes inclination and raises perigee of orbit as much as possible

v_apo_trans6=norm(v_N_apo_trans6); % km/s
h_trans6= r_apo_trans6*v_apo_trans6; % angular momentum
SemiR_trans6= (h_trans6^2)/mu; % km

syms a e 
eq1= SemiR_trans6 == a*(1-e^2);
eq2= r_apo_trans6 == a*(1+e);
sol_trans6=solve(eq1,eq2,a,e);
a_trans6= double(sol_trans6.a);
e_trans6= double(sol_trans6.e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 7 - Maneuvor 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_apo_trans7= r_apo_trans6; % km
r_N_apo_trans7= r_N_apo_trans6; % km
ApN_trans7= ApN_trans6;

v_N_apo_trans7= v_N_apo_trans6 + ApN_trans7*[0; -.468; 0];

v_apo_trans7=norm(v_N_apo_trans7); % km/s
h_trans7= r_apo_trans7*v_apo_trans7; % angular momentum
SemiR_trans7= (h_trans7^2)/mu; % km

syms a e 
eq1= SemiR_trans7 == a*(1-e^2);
eq2= r_apo_trans7 == a*(1+e);
sol_trans7=solve(eq1,eq2,a,e);
a_trans7= double(sol_trans7.a);
e_trans7= double(sol_trans7.e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 8 - Maneuvor 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_apo_trans8= r_apo_trans7; % km
r_N_apo_trans8= r_N_apo_trans7; % km
ApN_trans8= ApN_trans7;

v_N_apo_trans8= v_N_apo_trans7 + ApN_trans8*[0; -.468; 0];

v_apo_trans8=norm(v_N_apo_trans8); % km/s
h_trans8= r_apo_trans8*v_apo_trans8; % angular momentum
SemiR_trans8= (h_trans8^2)/mu; % km

syms a e 
eq1= SemiR_trans8 == a*(1-e^2);
eq2= r_apo_trans8 == a*(1+e);
sol_trans8=solve(eq1,eq2,a,e);
a_trans8= double(sol_trans8.a);
e_trans8= double(sol_trans8.e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 9 - Maneuvor 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_apo_trans9= r_apo_trans8; % km
r_N_apo_trans9= r_N_apo_trans8; % km
ApN_trans9= ApN_trans8;

v_N_apo_trans9= v_N_apo_trans8 + ApN_trans9*[0; -(delta_v_maneuvor3 - (pro_Brun_trans6 + .468 + .468)); 0];

v_apo_trans9=norm(v_N_apo_trans9); % km/s
h_trans9= r_apo_trans9*v_apo_trans9; % angular momentum
SemiR_trans9= (h_trans9^2)/mu; % km

syms a e 
eq1= SemiR_trans9 == a*(1-e^2);
eq2= r_apo_trans9 == a*(1+e);
sol_trans9=solve(eq1,eq2,a,e);
a_trans9= double(sol_trans9.a);
e_trans9= double(sol_trans9.e);

%% Phase Orbit Maneuvors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 10 - Maneuvor 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O_phase=deg2rad(0);
w_phase=deg2rad(30);
i_phase=deg2rad(0);
ApN_phase1=[ cos(O_phase)*cos(w_phase)-sin(O_phase)*sin(w_phase)*cos(i_phase), -cos(O_phase)*sin(w_phase)-sin(O_phase)*cos(w_phase)*cos(i_phase),  sin(O_phase)*sin(i_phase) ;...
            sin(O_phase)*cos(w_phase)+cos(O_phase)*sin(w_phase)*cos(i_phase), -sin(O_phase)*sin(w_phase)+cos(O_phase)*cos(w_phase)*cos(i_phase), -cos(O_phase)*sin(i_phase) ;...
            sin(w_phase)*sin(i_phase)                                    ,  cos(w_phase)*sin(i_phase)                                    ,  cos(i_phase)            ];

r_apo_trans10= r_apo_trans9;
P_trans10= P_tgt + [deg2rad(-.5)]/w_tgt;
a_trans10= (P_trans10*(sqrt(mu)/(2*pi)))^(2/3);
e_trans10= (r_apo_trans10/a_trans10)-1;
SemiR_trans10= a_trans10*(1-e_trans10^2);


v_apo_trans10= sqrt(mu*((2/r_apo_trans10)-(1/a_trans10)));

nu_apo_trans10= deg2rad(180);

r_P_apo_trans10= [r_apo_trans10*cos(nu_apo_trans10); r_apo_trans10*sin(nu_apo_trans10); 0];
r_N_apo_trans10= ApN_phase1*r_P_apo_trans10;

v_P_apo_trans10= [-sqrt(mu/SemiR_trans10)*sin(nu_apo_trans10); sqrt(mu/SemiR_trans10)*(e_trans10+cos(nu_apo_trans10)); 0];
v_N_apo_trans10= ApN_phase1*v_P_apo_trans10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 11 - Maneuvor 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O_phase=deg2rad(0);
w_phase=deg2rad(210);
i_phase=deg2rad(0);
ApN_phase2=[ cos(O_phase)*cos(w_phase)-sin(O_phase)*sin(w_phase)*cos(i_phase), -cos(O_phase)*sin(w_phase)-sin(O_phase)*cos(w_phase)*cos(i_phase),  sin(O_phase)*sin(i_phase) ;...
            sin(O_phase)*cos(w_phase)+cos(O_phase)*sin(w_phase)*cos(i_phase), -sin(O_phase)*sin(w_phase)+cos(O_phase)*cos(w_phase)*cos(i_phase), -cos(O_phase)*sin(i_phase) ;...
            sin(w_phase)*sin(i_phase)                                    ,  cos(w_phase)*sin(i_phase)                                    ,  cos(i_phase)            ];

r_per_trans11= r_apo_trans10;
P_trans11= P_tgt + [deg2rad(.35)]/w_tgt;
a_trans11= (P_trans11*(sqrt(mu)/(2*pi)))^(2/3);
e_trans11= 1-(r_per_trans11/a_trans11);
SemiR_trans11= a_trans11*(1-e_trans11^2);


v_per_trans11= sqrt(mu*((2/r_per_trans11)-(1/a_trans11)));

nu_per_trans11= deg2rad(0);

r_P_per_trans11= [r_per_trans11*cos(nu_per_trans11); r_per_trans11*sin(nu_per_trans11); 0];
r_N_per_trans11= ApN_phase2*r_P_per_trans11;

v_P_per_trans11= [-sqrt(mu/SemiR_trans11)*sin(nu_per_trans11); sqrt(mu/SemiR_trans11)*(e_trans11+cos(nu_per_trans11)); 0];
v_N_per_trans11= ApN_phase2*v_P_per_trans11;

%% Loops for plotting orbits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parking Orbit Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x1
r= r_N_park_t0;                       
rdot= v_N_park_t0;                  
tf = (2*pi*a_park^(3/2))/sqrt(mu);     
dt = 1;                                
nt = ceil(tf/dt);                    
x1 = zeros(6,nt);                    
t = 0;                             
x1(1:3,1) = r;                        
x1(4:6,1) = rdot;                
for k = 1:nt-1
    x_k = x1(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x1(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Target Orbit Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x2
r= r_N_tgt_t0;                  
rdot= v_N_tgt_t0;                    
tf = (2*pi*a_tgt^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x2 = zeros(6,nt);                
t = 0;                             
x2(1:3,1) = r;                       
x2(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x2(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x2(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Theoretical Orbit 1 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x3
r=[r_park;0;0];        
rdot=[0;v_per_theo1;0];                
tf = (2*pi*a_theo1^(3/2))/sqrt(mu);       
dt = 1;                               
nt = ceil(tf/dt);                 
x3 = zeros(6,nt);                
t = 0;                              
x3(1:3,1) = r;                           
x3(4:6,1) = rdot;                         
for k = 1:nt-1
    x_k = x3(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x3(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 1 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x4
r= r_N_per_trans1;                  
rdot= v_N_per_trans1;                    
tf = (2*pi*a_trans1^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x4 = zeros(6,nt);                
t = 0;                             
x4(1:3,1) = r;                       
x4(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x4(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x4(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 2 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x5
r= r_N_per_trans2;                  
rdot= v_N_per_trans2;                    
tf = (2*pi*a_trans2^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x5 = zeros(6,nt);                
t = 0;                             
x5(1:3,1) = r;                       
x5(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x5(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x5(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 3 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x6
r= r_N_per_trans3;                  
rdot= v_N_per_trans3;                    
tf = (2*pi*a_trans3^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x6 = zeros(6,nt);                
t = 0;                             
x6(1:3,1) = r;                       
x6(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x6(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x6(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 4 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x7
r= r_N_per_trans4;                  
rdot= v_N_per_trans4;                    
tf = (2*pi*a_trans4^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x7 = zeros(6,nt);                
t = 0;                             
x7(1:3,1) = r;                       
x7(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x7(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x7(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 5 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x8
r= r_N_per_trans5;                  
rdot= v_N_per_trans5;                    
tf = (2*pi*a_trans5^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x8 = zeros(6,nt);                
t = 0;                             
x8(1:3,1) = r;                       
x8(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x8(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x8(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 6 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x9
r= r_N_apo_trans6;                  
rdot= v_N_apo_trans6;                    
tf = (2*pi*a_trans6^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x9 = zeros(6,nt);                
t = 0;                             
x9(1:3,1) = r;                       
x9(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x9(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x9(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 7 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x11
r= r_N_apo_trans7;                  
rdot= v_N_apo_trans7;                    
tf = (2*pi*a_trans7^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x11 = zeros(6,nt);                
t = 0;                             
x11(1:3,1) = r;                       
x11(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x11(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x11(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 8 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x12
r= r_N_apo_trans8;                  
rdot= v_N_apo_trans8;                    
tf = (2*pi*a_trans8^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x12 = zeros(6,nt);                
t = 0;                             
x12(1:3,1) = r;                       
x12(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x12(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x12(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 9 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x13
r= r_N_apo_trans9;                  
rdot= v_N_apo_trans9;                    
tf = (2*pi*a_trans9^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x13 = zeros(6,nt);                
t = 0;                             
x13(1:3,1) = r;                       
x13(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x13(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x13(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 10 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x14
r= r_N_apo_trans10;                  
rdot= v_N_apo_trans10;                    
tf = (2*pi*a_trans10^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x14 = zeros(6,nt);                
t = 0;                             
x14(1:3,1) = r;                       
x14(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x14(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x14(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer Orbit 11 Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x15
r= r_N_per_trans11;                  
rdot= v_N_per_trans11;                    
tf = (2*pi*a_trans11^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x15 = zeros(6,nt);                
t = 0;                             
x15(1:3,1) = r;                       
x15(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x15(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x15(:,k+1) = x_kp1; 
end 

%% Plot Earth, Orbits, and time points
%%%%%%%%%% Code to plot the Earth %%%%%%%%%%
imData = imread('2_no_clouds_4k.jpg');
[xS,yS,zS] = sphere(50);
earth_radius = 6378.1370;
xSE = earth_radius*xS;
ySE = earth_radius*yS;
zSE = earth_radius*zS;
surface(xSE,ySE,zSE);
axis equal
grid on
xlabel('Inertial x (km)')
ylabel('Inertial y (km)')
zlabel('Inertial z (km)')
ch = get(gca,'children');
set(ch,'facecolor','texturemap','cdata',flipud(imData),'edgecolor','none');

%%%%%%%%%% Plot Orbits %%%%%%%%%%
hold on 

% Theoretical Orbits
%plot3(x3(1,:), x3(2,:), x3(3,:), '-','Color','[0 1 0]', 'Linewidth', 2)
%plot3(x10(1,:), x10(2,:), x10(3,:), '-','Color','[0 1 0]', 'Linewidth', 2)

% Maneuvor 1 Orbits
plot3(x4(1,:), x4(2,:), x4(3,:), '-','Color','[1 0 0]', 'Linewidth', 1)
plot3(x5(1,:), x5(2,:), x5(3,:), '-','Color','[1 0 0]', 'Linewidth', 1)
plot3(x6(1,:), x6(2,:), x6(3,:), '-','Color','[1 0 0]', 'Linewidth', 1)
plot3(x7(1,:), x7(2,:), x7(3,:), '-','Color','[1 0 0]', 'Linewidth', 1)
plot3(x8(1,:), x8(2,:), x8(3,:), '-','Color','[1 0 0]', 'Linewidth', 1)

% Maneuvor 2 Orbits
plot3(x9(1,:), x9(2,:), x9(3,:), '-','Color','[0 1 0]', 'Linewidth', 1)
plot3(x11(1,:), x11(2,:), x11(3,:), '-','Color','[0 1 0]', 'Linewidth', 1)
plot3(x12(1,:), x12(2,:), x12(3,:), '-','Color','[0 1 0]', 'Linewidth', 1)
plot3(x13(1,:), x13(2,:), x13(3,:), '-','Color','[0 1 0]', 'Linewidth', 1)

% Initial Orbits
plot3(x1(1,:), x1(2,:), x1(3,:), '-','Color','[0 0 0]', 'Linewidth', 3)
plot3(x2(1,:), x2(2,:), x2(3,:), '-','Color','[1 1 0]', 'Linewidth', 3)

% Interception Orbits
plot3(x14(1,:), x14(2,:), x14(3,:), '-','Color','[0 0 1]', 'Linewidth', 1)
plot3(x15(1,:), x15(2,:), x15(3,:), '-','Color','[0 0 1]', 'Linewidth', 1)

% Points
plot3(r_N_park_t0(1), r_N_park_t0(2), r_N_park_t0(3), 'o','Color','[1 0 1]', 'Linewidth', 2)                  % Chase position at t0
plot3(r_N_tgt_t0(1), r_N_tgt_t0(2), r_N_tgt_t0(3), 'o','Color','[0 1 1]', 'Linewidth', 2)                     % Target position at t0
plot3(r_N_apo_trans9(1), r_N_apo_trans9(2), r_N_apo_trans9(3), 'o','Color','[.8 0 .8]', 'Linewidth', 2)       % Chase position at t1
plot3(r_N_apo_trans10(1), r_N_apo_trans10(2), r_N_apo_trans10(3), 'o','Color','[.6 0 .6]', 'Linewidth', 2)    % Chase position at t2
plot3(r_N_per_trans11(1), r_N_per_trans11(2), r_N_per_trans11(3), 'o','Color','[.4 0 .4]', 'Linewidth', 2)    % Chase position at t3

%% Keep Track of Time
P_park  = (2*pi*a_park  ^(3/2))/sqrt(mu);
P_tgt   = (2*pi*a_tgt   ^(3/2))/sqrt(mu);
P_trans1= (2*pi*a_trans1^(3/2))/sqrt(mu);
P_trans2= (2*pi*a_trans2^(3/2))/sqrt(mu);
P_trans3= (2*pi*a_trans3^(3/2))/sqrt(mu);
P_trans4= (2*pi*a_trans4^(3/2))/sqrt(mu);
P_trans5= (2*pi*a_trans5^(3/2))/sqrt(mu);
P_trans6= (2*pi*a_trans6^(3/2))/sqrt(mu);
P_trans7= (2*pi*a_trans7^(3/2))/sqrt(mu);
P_trans8= (2*pi*a_trans8^(3/2))/sqrt(mu);
P_trans9= (2*pi*a_trans9^(3/2))/sqrt(mu);
P_trans10= (2*pi*a_trans10^(3/2))/sqrt(mu);
P_trans11= (2*pi*a_trans11^(3/2))/sqrt(mu);

%%%%%%%%%%%%%%%%%%%% t1 math %%%%%%%%%%%%%%%%%%%%
hold on
for k=7:7
    ManeuvorTime_t1= k*P_park + P_trans1 + P_trans2 + P_trans3 + P_trans4 + .5*P_trans5 + P_trans6 + P_trans7 + P_trans8;

    theta_tgt_t1= (w_tgt*ManeuvorTime_t1 + nu_tgt_t0);
    nu_tgt_t1= 2*pi*(theta_tgt_t1/(2*pi) - floor(theta_tgt_t1/(2*pi)));
    nu_tgt_t1_deg= rad2deg(nu_tgt_t1);

    r_P_tgt_t1= [r_tgt*cos(nu_tgt_t1); r_tgt*sin(nu_tgt_t1); 0];
    r_N_tgt_t1=ApN_tgt*r_P_tgt_t1;

    v_P_tgt_t1= [-sqrt(mu/SemiR_tgt)*sin(nu_tgt_t1); sqrt(mu/SemiR_tgt)*(e_tgt+cos(nu_tgt_t1)); 0];
    v_N_tgt_t1= ApN_tgt*v_P_tgt_t1;

%%%%%%%%%% Plotting %%%%%%%%%%
plot3(r_N_tgt_t1(1), r_N_tgt_t1(2), r_N_tgt_t1(3), 'o','Color','[0 .8 .8]', 'Linewidth', 2) % Target Position at t1
end
hold off

%%%%%%%%%%%%%%%%%%%% t2 math %%%%%%%%%%%%%%%%%%%%
ManeuvorTime_t2= ManeuvorTime_t1 + P_trans9 + P_trans10;

    theta_tgt_t2= (w_tgt*ManeuvorTime_t2 + nu_tgt_t0);
    nu_tgt_t2= 2*pi*(theta_tgt_t2/(2*pi) - floor(theta_tgt_t2/(2*pi)));
    nu_tgt_t2_deg= rad2deg(nu_tgt_t2);

    r_P_tgt_t2= [r_tgt*cos(nu_tgt_t2); r_tgt*sin(nu_tgt_t2); 0];
    r_N_tgt_t2=ApN_tgt*r_P_tgt_t2;

%%%%%%%%%% Plotting %%%%%%%%%%
hold on
plot3(r_N_tgt_t2(1), r_N_tgt_t2(2), r_N_tgt_t2(3), 'o','Color','[0 .6 .6]', 'Linewidth', 2) % Target Position at t2
hold off

%%%%%%%%%%%%%%%%%%%% t3 math %%%%%%%%%%%%%%%%%%%%
ManeuvorTime_t3= ManeuvorTime_t2 + P_trans11;

    theta_tgt_t3= (w_tgt*ManeuvorTime_t3 + nu_tgt_t0);
    nu_tgt_t3= 2*pi*(theta_tgt_t3/(2*pi) - floor(theta_tgt_t3/(2*pi)));
    nu_tgt_t3_deg= rad2deg(nu_tgt_t3);

    r_P_tgt_t3= [r_tgt*cos(nu_tgt_t3); r_tgt*sin(nu_tgt_t3); 0];
    r_N_tgt_t3=ApN_tgt*r_P_tgt_t3;

%%%%%%%%%% Plotting %%%%%%%%%%
hold on
plot3(r_N_tgt_t3(1), r_N_tgt_t3(2), r_N_tgt_t3(3), 'o','Color','[0 .4 .4]', 'Linewidth', 2) % Target Position at t3
hold off

%% Delta V and Elapsed Time Calcualtions

% Burn 1
deltaV_Burn1= abs(v_per_trans1 - v_park);
elapst_Burn1= 7*P_park;

% Burn 2
deltaV_Burn2= abs(v_per_trans2 - v_per_trans1);
elapst_Burn2= elapst_Burn1 + P_trans1;

% Burn 3
deltaV_Burn3= abs(v_per_trans3 - v_per_trans2);
elapst_Burn3= elapst_Burn2 + P_trans2;

% Burn 4
deltaV_Burn4= abs(v_per_trans4 - v_per_trans3);
elapst_Burn4= elapst_Burn3 + P_trans3;

% Burn 5
deltaV_Burn5= abs(v_per_trans5 - v_per_trans4);
elapst_Burn5= elapst_Burn4 + P_trans4;

% Burn 6
deltaV_Burn6= abs(norm(v_N_apo_trans6 - v_N_apo_trans5));
elapst_Burn6= elapst_Burn5 + .5*P_trans5;

% Burn 7
deltaV_Burn7= abs(v_apo_trans7 - v_apo_trans6);
elapst_Burn7= elapst_Burn6 + P_trans6;

% Burn 8
deltaV_Burn8= abs(v_apo_trans8 - v_apo_trans7);
elapst_Burn8= elapst_Burn7 + P_trans7;

% Burn 9
deltaV_Burn9= abs(norm(v_N_apo_trans9 - v_N_apo_trans8));
elapst_Burn9= elapst_Burn8 + P_trans8;

% Burn 10
deltaV_Burn10= abs(v_apo_trans10 - v_apo_trans9);
elapst_Burn10= elapst_Burn9 + P_trans9;

% Burn 11
deltaV_Burn11= abs(v_per_trans11 - v_apo_trans10);
elapst_Burn11= elapst_Burn10 + P_trans10;

% Burn 12
deltaV_Burn12= abs(v_per_trans11 - v_apo_trans9);
elapst_Burn12= elapst_Burn11 + P_trans11;

DELTA_V_TOTAL= deltaV_Burn1 + deltaV_Burn2 + deltaV_Burn3 + deltaV_Burn4 + deltaV_Burn5 + deltaV_Burn6 + deltaV_Burn7 + deltaV_Burn8 + deltaV_Burn9 + deltaV_Burn10 + deltaV_Burn11 + deltaV_Burn12;


deltaV_vector= 1000*[deltaV_Burn1  deltaV_Burn2  deltaV_Burn3  deltaV_Burn4  deltaV_Burn5  deltaV_Burn6  deltaV_Burn7  deltaV_Burn8  deltaV_Burn9  deltaV_Burn10  deltaV_Burn11  deltaV_Burn12]';

%% Ground Trace
hold off

EarthRoationalSpeed= 7.292123516990375e-05; % rad/s - calculated from 23hr, 56min, 4s day rather than 24hr day

lat_x1= asin(x1(3,:) ./ a_park ); % Parking Orbit
length_x1= length(lat_x1);
lon_x1=zeros([1,length_x1]);
for i= 1:length_x1
    lon_x1(i)= atan2(x1(2,i),x1(1,i)) - EarthRoationalSpeed*(i + 0);
    if lon_x1(i) < -pi
        lon_x1(i)= lon_x1(i) + 2*pi;
    end
end

lat_x4= asin(x4(3,:) ./ a_trans1 ); % Transfer Orbit 1
length_x4= length(lat_x4);
lon_x4=zeros([1,length_x4]);
for i= 1:length_x4
    lon_x4(i)= atan2(x4(2,i),x4(1,i)) - EarthRoationalSpeed*(i + 7*P_park);
    if lon_x4(i) < -pi
        lon_x4(i)= lon_x4(i) + 2*pi;
    end
end

lat_x5= asin(x5(3,:) ./ a_trans2 ); % Transfer Orbit 2
length_x5= length(lat_x5);
lon_x5=zeros([1,length_x5]);
for i= 1:length_x5
    lon_x5(i)= atan2(x5(2,i),x5(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1);
    if lon_x5(i) < -pi
        lon_x5(i)= lon_x5(i) + 2*pi;
    end
end

lat_x6= asin(x6(3,:) ./ a_trans3 ); % Transfer Orbit 3
length_x6= length(lat_x6);
lon_x6=zeros([1,length_x6]);
for i= 1:length_x6
    lon_x6(i)= atan2(x6(2,i),x6(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1 + P_trans2);
    if lon_x6(i) < -pi
        lon_x6(i)= lon_x6(i) + 2*pi;
    end
end

lat_x7= asin(x7(3,:) ./ a_trans4 ); % Transfer Orbit 4
length_x7= length(lat_x7);
lon_x7=zeros([1,length_x7]);
for i= 1:length_x7
    lon_x7(i)= atan2(x7(2,i),x7(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1 + P_trans2 + P_trans3);
    if lon_x7(i) < -pi
        lon_x7(i)= lon_x7(i) + 2*pi;
    end
end

lat_x8= asin(x8(3,:) ./ a_trans5 ); % Transfer Orbit 5
length_x8= length(lat_x8);
lon_x8=zeros([1,length_x8]);
for i= 1:length_x8
    lon_x8(i)= atan2(x8(2,i),x8(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1 + P_trans2 + P_trans3 + P_trans4);
    if lon_x8(i) < -pi
        lon_x8(i)= lon_x8(i) + 2*pi;
    end
    if lon_x8(i) < -pi
        lon_x8(i)= lon_x8(i) + 2*pi;
    end
end

lat_x9= asin(x9(3,:) ./ a_trans6 ); % Transfer Orbit 6
length_x9= length(lat_x9);
lon_x9=zeros([1,length_x9]);
for i= 1:length_x9
    lon_x9(i)= atan2(x9(2,i),x9(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1 + P_trans2 + P_trans3 + P_trans4 + .5*P_trans5);
    if lon_x9(i) < -pi
        lon_x9(i)= lon_x9(i) + 2*pi;
    end
    if lon_x9(i) < -pi
        lon_x9(i)= lon_x9(i) + 2*pi;
    end
    if lon_x9(i) < -pi
        lon_x9(i)= lon_x9(i) + 2*pi;
    end
end

lat_x11= asin(x11(3,:) ./ a_trans7 ); % Transfer Orbit 7
length_x11= length(lat_x11);
lon_x11=zeros([1,length_x11]);
for i= 1:length_x11
    lon_x11(i)= atan2(x11(2,i),x11(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1 + P_trans2 + P_trans3 + P_trans4 + .5*P_trans5 + P_trans6);
    if lon_x11(i) < -pi
        lon_x11(i)= lon_x11(i) + 2*pi;
    end
    if lon_x11(i) < -pi
        lon_x11(i)= lon_x11(i) + 2*pi;
    end
    if lon_x11(i) < -pi
        lon_x11(i)= lon_x11(i) + 2*pi;
    end
    if lon_x11(i) < -pi
        lon_x11(i)= lon_x11(i) + 2*pi;
    end
end

lat_x12= asin(x12(3,:) ./ a_trans8 ); % Transfer Orbit 8
length_x12= length(lat_x12);
lon_x12=zeros([1,length_x12]);
for i= 1:length_x12
    lon_x12(i)= atan2(x12(2,i),x12(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1 + P_trans2 + P_trans3 + P_trans4 + .5*P_trans5 + P_trans6 + P_trans7);
    if lon_x12(i) < -pi
        lon_x12(i)= lon_x12(i) + 2*pi;
    end
    if lon_x12(i) < -pi
        lon_x12(i)= lon_x12(i) + 2*pi;
    end
    if lon_x12(i) < -pi
        lon_x12(i)= lon_x12(i) + 2*pi;
    end
    if lon_x12(i) < -pi
        lon_x12(i)= lon_x12(i) + 2*pi;
    end
end

lat_x13= asin(x13(3,:) ./ a_trans9 ); % Transfer Orbit 9
length_x13= length(lat_x13);
lon_x13=zeros([1,length_x13]);
for i= 1:length_x13
    lon_x13(i)= atan2(x13(2,i),x13(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1 + P_trans2 + P_trans3 + P_trans4 + .5*P_trans5 + P_trans6 + P_trans7 + P_trans8);
    if lon_x13(i) < -pi
        lon_x13(i)= lon_x13(i) + 2*pi;
    end
    if lon_x13(i) < -pi
        lon_x13(i)= lon_x13(i) + 2*pi;
    end
    if lon_x13(i) < -pi
        lon_x13(i)= lon_x13(i) + 2*pi;
    end
    if lon_x13(i) < -pi
        lon_x13(i)= lon_x13(i) + 2*pi;
    end
    if lon_x13(i) < -pi
        lon_x13(i)= lon_x13(i) + 2*pi;
    end
end

lat_x14= asin(x14(3,:) ./ a_trans10 ); % Transfer Orbit 10
length_x14= length(lat_x14);
lon_x14=zeros([1,length_x14]);
for i= 1:length_x14
    lon_x14(i)= atan2(x14(2,i),x14(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1 + P_trans2 + P_trans3 + P_trans4 + .5*P_trans5 + P_trans6 + P_trans7 + P_trans8 + P_trans9);
    if lon_x14(i) < -pi
        lon_x14(i)= lon_x14(i) + 2*pi;
    end
    if lon_x14(i) < -pi
        lon_x14(i)= lon_x14(i) + 2*pi;
    end
    if lon_x14(i) < -pi
        lon_x14(i)= lon_x14(i) + 2*pi;
    end
    if lon_x14(i) < -pi
        lon_x14(i)= lon_x14(i) + 2*pi;
    end
    if lon_x14(i) < -pi
        lon_x14(i)= lon_x14(i) + 2*pi;
    end
    if lon_x14(i) < -pi
        lon_x14(i)= lon_x14(i) + 2*pi;
    end
end

lat_x15= asin(x15(3,:) ./ a_trans11 ); % Transfer Orbit 11
length_x15= length(lat_x15);
lon_x15=zeros([1,length_x15]);
for i= 1:length_x15
    lon_x15(i)= atan2(x15(2,i),x15(1,i)) - EarthRoationalSpeed*(i + 7*P_park + P_trans1 + P_trans2 + P_trans3 + P_trans4 + .5*P_trans5 + P_trans6 + P_trans7 + P_trans8 + P_trans9 + P_trans10);
    if lon_x15(i) < -pi
        lon_x15(i)= lon_x15(i) + 2*pi;
    end
    if lon_x15(i) < -pi
        lon_x15(i)= lon_x15(i) + 2*pi;
    end
    if lon_x15(i) < -pi
        lon_x15(i)= lon_x15(i) + 2*pi;
    end
    if lon_x15(i) < -pi
        lon_x15(i)= lon_x15(i) + 2*pi;
    end
    if lon_x15(i) < -pi
        lon_x15(i)= lon_x15(i) + 2*pi;
    end
    if lon_x15(i) < -pi
        lon_x15(i)= lon_x15(i) + 2*pi;
    end
    if lon_x15(i) < -pi
        lon_x15(i)= lon_x15(i) + 2*pi;
    end
end

lat_x2= asin(x2(3,:) ./ a_tgt ); % yellowhammer Orbit
length_x2= length(lat_x2);
lon_x2=zeros([1,length_x2]);
for i= 1:length_x2
    lon_x2(i)= atan2(x2(2,i),x2(1,i)) - EarthRoationalSpeed*(i + 0);
    if lon_x2(i) < -pi
        lon_x2(i)= lon_x2(i) + 2*pi;
    end
end


figure(2)
geoplot(0,180,'.y','MarkerSize',20) % Yellowhammer Orbit
hold on
geoplot(rad2deg(lat_x1),rad2deg(lon_x1),'.','Color','[0 0 1]','MarkerSize',1) % Parking Orbit
geoplot(rad2deg(lat_x4),rad2deg(lon_x4),'.','Color','[0 1 1]','MarkerSize',1) % Transfer Orbit 1
geoplot(rad2deg(lat_x5),rad2deg(lon_x5),'.','Color','[0 1 0]','MarkerSize',1) % Transfer Orbit 2
geoplot(rad2deg(lat_x6),rad2deg(lon_x6),'.','Color','[1 1 0]','MarkerSize',1) % Transfer Orbit 3
geoplot(rad2deg(lat_x7),rad2deg(lon_x7),'.','Color','[1 0 0]','MarkerSize',1) % Transfer Orbit 4
geoplot(rad2deg(lat_x8(1:19195)),rad2deg(lon_x8(1:19195)),'.','Color','[1 0 1]','MarkerSize',1) % Transfer Orbit 5
geoplot(rad2deg(lat_x9),rad2deg(lon_x9),'.','Color','[0 0 1]','MarkerSize',1) % Transfer Orbit 6
geoplot(rad2deg(lat_x11),rad2deg(lon_x11),'.','Color','[0 1 1]','MarkerSize',1) % Transfer Orbit 7
geoplot(rad2deg(lat_x12),rad2deg(lon_x12),'.','Color','[0 1 0]','MarkerSize',1) % Transfer Orbit 8
geoplot(rad2deg(lat_x13),rad2deg(lon_x13),'.','Color','[1 1 0]','MarkerSize',1) % Transfer Orbit 9
geoplot(rad2deg(lat_x14),rad2deg(lon_x14),'.','Color','[1 0 0]','MarkerSize',1) % Transfer Orbit 10
geoplot(rad2deg(lat_x15),rad2deg(lon_x15),'.','Color','[1 0 1]','MarkerSize',1) % Transfer Orbit 11

geolimits([-45 45],[-180 180]);
legend('Yellowhammer Orbit','Parking Orbit','Transfer Orbit 1','Transfer Orbit 2','Transfer Orbit 3','Transfer Orbit 4','Transfer Orbit 5','Transfer Orbit 6','Transfer Orbit 7','Transfer Orbit 8','Transfer Orbit 9','Transfer Orbit 10','Transfer Orbit 11',Location='bestoutside')
hold off
%% Relative Motion Plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test plot - Yellowhammer Position progogated started from t1 point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - x16
% r_N_apo_trans9 is the starting point of x13, which occurs at same instant
% as starting point for this yellowhammer propogation
r= r_N_tgt_t1;                  
rdot= v_N_tgt_t1;                    
tf = (2*pi*a_tgt^(3/2))/sqrt(mu);       
dt = 1;                              
nt = ceil(tf/dt);                     
x16 = zeros(6,nt);                
t = 0;                             
x16(1:3,1) = r;                       
x16(4:6,1) = rdot;                     
for k = 1:nt-1
    x_k = x16(:,k); 
    k1 = twoB_EOM(t, x_k, mu);
    k2 = twoB_EOM(t + dt/2, x_k + k1*(dt/2), mu);
    k3 = twoB_EOM(t + dt/2, x_k + k2*(dt/2), mu);
    k4 = twoB_EOM(t + dt  , x_k + k3*dt,     mu);
    x_kp1 = x_k + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x16(:,k+1) = x_kp1; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%% The rest of the RIC math
X_chase= [x13,x14,x15];
X_target=[x16,x16,x16(1:6,1:86128)];

X_relative= X_chase - X_target;

dumbVector=linspace(1,258456,258456);

nu_closeenough = sqrt(mu / a_tgt^3).*dumbVector;

X_relative = [ X_relative(1,:).*cos(nu_closeenough) + X_relative(2,:).*sin(nu_closeenough) ;...
         -X_relative(1,:).*sin(nu_closeenough) + X_relative(2,:).*cos(nu_closeenough)];

%%%%%%%%%%%%%%%%%%%%%%%%%% RIC plotting
figure(3);
plot(X_relative(2,:),X_relative(1,:),'k','LineWidth',2);
hold on
plot(X_relative(2,1),X_relative(1,1),'ro','LineWidth',2);
plot(X_relative(2,end),X_relative(1,end),'rx','LineWidth',2);
plot(0,0,'b+','LineWidth',2);
hold off
grid('minor');
grid on;
axis equal
set(gca,'XDir','reverse');
xlabel('In-Track [km]');
ylabel('Radial [km]');
legend('InSat-1 Path','InSat-1 Start','InSat-1 End','Yellowhammer');







