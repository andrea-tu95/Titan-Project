%project code
%team members: MIRABELLI MARCO, ROMANO AGOSTINO, TURCHI ANDREA
%author: ANDREA TURCHI
%latest update: 09-06-21

close all
clear
clc
format long
set(0,'DefaultTextInterpreter','Latex')

%% Ephemerides & Problem Constants

%Constants

mu_Sun = 132712440018;                    %Sun gravitational parameter
%[km^3/s^2] 
mu_Earth = 398600.4418;                   %Earth gravitational parameter 
%[km^3/s^2]
mu_Jupiter = 126686534;                   %Jupiter gravitational parameter
%[km^3/s^2]
mu_Saturn = 37931187;                     %Saturn gravitational parameter 
%[km^3/s^2]
mu_Titan = 8978.1383;                     %Titan gravitational parameter
%[km^3/s^2]

R_Earth = 6371;                           %mean Earth radius [km]
R_Jupiter = 69911;                        %mean Jupiter radius [km]
R_Saturn = 58232;                         %mean Saturn radius [km]
R_Titan = 2576;                           %mean Titan radius [km]

JD = 86400;                               %Julian day [s]
JY = 365.25;                              %Julian year [day]
au = 149597870.700;                       %astronomical unit [km]
g0 = 9.80665;                             %reference gravity acceleration
%[m/s^2]
T_Titan = 15 + 22/24;                     %Titan rotation period [days]
                  
min_perigee = (R_Earth + 500)/au;         %minimum allowed perigee [au]
min_perigiovio = (R_Jupiter + 500)/au;    %minimum allowed perigiovio [au]
       
max_perigee = (R_Earth + 1e4)/au;         %maximum allowed perigee [au]
max_perigiovio = (R_Jupiter + 1e6)/au;    %maximum allowed perigiovio [au]

conv_v = JY*JD/(au*2*pi);                 %conversion km/s-EOS (Earth Orbit
%Speed)
conv_time = JY*JD/(2*pi);                 %conversion year*(2*pi)-s

%Optimization Setup

%variables lower bounds 
lb1 = [ 15*2*pi; 1.2*2*pi;  -0.06; -0.05; 1*2*pi; 1*2*pi; min_perigee;...
      min_perigiovio ]; 
%variables upper bounds
ub1 = [ 50*2*pi; 1.8*2*pi; -0.005;  0.05; 3*2*pi; 9*2*pi; max_perigee;...
      max_perigiovio ];

%best solution as a row in the initial matrix to faster solve the 1st
%optimization problem
X_best1 = [ 1.562348817019730e+02, 11.309733551765660,...
          -0.018355872559017, -2.759065985047989e-05,...
           13.486816633031868, 42.171658343040260,...
           8.289535726609700e-05, 0.007151913292573 ];

%randomically generate the remainig rows of the matrix 
Random_Matrix = ub1'.*(randn(79)*ones(79,1));
Swarm_init1 = [ X_best1; Random_Matrix ];             

%variables lower bounds 
lb2 = [ -pi/2; -pi; 1e-6 ]; 
%variables upper bounds
ub2 = [  pi/2;   0; 2*pi ];

%best solution as a row in the initial matrix to faster solve the 2nd
%optimization problem
X_best2 = [ 0.166206760833434, -1.024482648628141, 3.216563133724658 ];

%randomically generate the remainig rows of the matrix 
Random_Matrix = ub2'.*(randn(11)*ones(11,1));
Swarm_init2 = [ X_best2; Random_Matrix ];  

%optimization & numerical integration options
options1 = optimoptions('particleswarm', 'Display', 'iter',...
           'UseParallel', true, 'SwarmSize', 80,...
           'TolFun', 1e-18, 'MaxIter', 10000, 'InitialSwarmMatrix',...
           Swarm_init1);
options2 = odeset('RelTol', 1.0e-13, 'AbsTol', 1.0e-13);                                
options3 = optimoptions('particleswarm', 'Display', 'iter',...
           'UseParallel', true, 'SwarmSize', 12,...
           'TolFun', 1e-18, 'MaxIter', 10000, 'InitialSwarmMatrix',...
           Swarm_init2);

%Normalization

mu_Earth = mu_Earth/mu_Sun;               %adimensional Earth gravitational
%parameter
mu_Jupiter = mu_Jupiter/mu_Sun;           %adimensional Jupiter 
%gravitational parameter
mu_Saturn = mu_Saturn/mu_Sun;             %adimensional Saturn 
%gravitational parameter
mu_Titan = mu_Titan/mu_Sun;               %adimensional Titan gravitational
%parameter
mu_Sun = mu_Sun/mu_Sun;                   %adimensional Sun gravitational 
%parameter

%Ephemerides (01-05-21 00:00:00 from JPL HORIZONS)

%Earth

a_Earth = 9.994105054783823e-1;           %semi-major axis [au]
e_Earth = 1.726861984783315e-2;           %eccentricity
i_Earth = deg2rad(0);                     %inclination w.r.t. ecliptic 
%plane [rad]
RAAN_Earth = deg2rad(1.951051169092126e2);%right ascension of ascending 
%node [rad]
PA_Earth = deg2rad(2.666147990815454e2);  %perihelion argument [rad]
TA_Earth = deg2rad(1.188386647448245e2);  %true anomaly [rad]

p_Earth = a_Earth*(1 - e_Earth^2);        %semilatus rectum [au]

%Jupiter 

a_Jupiter = 5.204288073499220;            %semi-major axis [au]
e_Jupiter = 4.868543514095491e-2;         %eccentricity
i_Jupiter = deg2rad(1.303372037930412);   %inclination w.r.t. ecliptic 
%plane [rad]
RAAN_Jupiter = ...
            deg2rad(1.005101036712953e2); %right ascension of ascending 
%node [rad]
PA_Jupiter = deg2rad(2.733337349490570e2);%perihelion argument [rad]
TA_Jupiter = deg2rad(3.032743492365449e2);%true anomaly [rad]

p_Jupiter = a_Jupiter*(1 - e_Jupiter^2);  %semilatus rectum [au]

%Saturn

a_Saturn = 9.579345347178357;             %semi-major axis [au]
e_Saturn = 5.206451200189190e-2;          %eccentricity
i_Saturn = deg2rad(2.482597531912353);    %inclination w.r.t. ecliptic 
%plane [rad]
RAAN_Saturn = ...
           deg2rad(1.135750200001832e2);  %right ascension of ascending 
%node [rad]
PA_Saturn = deg2rad(3.363133626036928e2); %perihelion argument [rad]
TA_Saturn = deg2rad(2.171447899637467e2); %true anomaly [rad]

p_Saturn = a_Saturn*(1 - e_Saturn^2);     %semilatus rectum [au]

%Titan

a_Titan = 8.168306717577080e-3;           %semi-major axis [au]
e_Titan = 2.878216054864953e-2;           %eccentricity
i_Titan = deg2rad(2.770296377344247e1);   %inclination w.r.t. ecliptic
%plane [rad]
RAAN_Titan = deg2rad(1.690632268824040e2);%right ascension of ascending 
%node [rad]
PA_Titan = deg2rad(1.753306235678690e2);  %pericenter argument [rad]
TA_Titan = deg2rad(3.585451633038113e2);  %true anomaly [rad]

p_Titan = a_Titan*(1 - e_Titan^2);        %semilatus rectum [au]

%% Perifocal Reference Frame

%Earth

xp_Earth = p_Earth*cos(TA_Earth)/(1 + e_Earth*cos(TA_Earth));
yp_Earth = p_Earth*sin(TA_Earth)/(1 + e_Earth*cos(TA_Earth));
zp_Earth = 0;

Vxp_Earth = - (mu_Sun/p_Earth)^(1/2)*sin(TA_Earth);
Vyp_Earth = (mu_Sun/p_Earth)^(1/2)*(e_Earth + cos(TA_Earth));
Vzp_Earth = 0;

%Jupiter

xp_Jupiter = p_Jupiter*cos(TA_Jupiter)/(1 + e_Jupiter*cos(TA_Jupiter));
yp_Jupiter = p_Jupiter*sin(TA_Jupiter)/(1 + e_Jupiter*cos(TA_Jupiter));
zp_Jupiter = 0;

Vxp_Jupiter = - (mu_Sun/p_Jupiter)^(1/2)*sin(TA_Jupiter);
Vyp_Jupiter = (mu_Sun/p_Jupiter)^(1/2)*(e_Jupiter + cos(TA_Jupiter));
Vzp_Jupiter = 0;

%Saturn

xp_Saturn = p_Saturn*cos(TA_Saturn)/(1 + e_Saturn*cos(TA_Saturn));
yp_Saturn = p_Saturn*sin(TA_Saturn)/(1 + e_Saturn*cos(TA_Saturn));
zp_Saturn = 0;

Vxp_Saturn = - (mu_Sun/p_Saturn)^(1/2)*sin(TA_Saturn);
Vyp_Saturn = (mu_Sun/p_Saturn)^(1/2)*(e_Saturn + cos(TA_Saturn));
Vzp_Saturn = 0;

%Titan 

xp_Titan = p_Titan*cos(TA_Titan)/(1 + e_Titan*cos(TA_Titan));
yp_Titan = p_Titan*sin(TA_Titan)/(1 + e_Titan*cos(TA_Titan));
zp_Titan = 0;

Vxp_Titan = - (mu_Saturn/p_Titan)^(1/2)*sin(TA_Titan);
Vyp_Titan = (mu_Saturn/p_Titan)^(1/2)*(e_Titan + cos(TA_Titan));
Vzp_Titan = 0;

%% Heliocentric Inertial Frame (HIF)

%Earth
          
[ ROT_Earth ] = ROT_p2HIF (RAAN_Earth, PA_Earth, i_Earth);

R_HIF_Earth = ROT_Earth*[ xp_Earth; yp_Earth; zp_Earth ];
V_HIF_Earth = ROT_Earth*[ Vxp_Earth; Vyp_Earth; Vzp_Earth ];

%Jupiter
            
[ ROT_Jupiter ] = ROT_p2HIF (RAAN_Jupiter, PA_Jupiter, i_Jupiter);
           
R_HIF_Jupiter = ROT_Jupiter*[ xp_Jupiter; yp_Jupiter; zp_Jupiter ];
V_HIF_Jupiter = ROT_Jupiter*[ Vxp_Jupiter; Vyp_Jupiter; Vzp_Jupiter ];

%Saturn
           
[ ROT_Saturn ] = ROT_p2HIF (RAAN_Saturn, PA_Saturn, i_Saturn);
           
R_HIF_Saturn = ROT_Saturn*[ xp_Saturn; yp_Saturn; zp_Saturn ];
V_HIF_Saturn = ROT_Saturn*[ Vxp_Saturn; Vyp_Saturn; Vzp_Saturn ];

%Titan 

[ ROT_Titan ] = ROT_p2HIF (RAAN_Titan, PA_Titan, i_Titan);

R_HIF_Titan = ROT_Titan*[ xp_Titan; yp_Titan; zp_Titan ] + R_HIF_Saturn;
V_HIF_Titan = ROT_Titan*[ Vxp_Titan; Vyp_Titan; Vzp_Titan ] + V_HIF_Saturn;

%% Initial Conditions

X0 = [ R_HIF_Earth; V_HIF_Earth; R_HIF_Jupiter; V_HIF_Jupiter;...
     R_HIF_Saturn; V_HIF_Saturn ];

%3:1 DeltaV-EGA Maneuver

T_ell = 3*2*pi;                           %elliptic orbit period 
%[year*(2*pi)]
a_ell = (mu_Sun*(T_ell/(2*pi))^2)^(1/3);  %elliptic orbit semi-major axis 
%[au]

%Spacecraft Features

%hypergolic bipropellant engine (MMH-N2O4)
Isp = 336;                                %specific impulse [s]

Cd = 2;                                   %drag coefficient
d = 4;                                    %diameter [m] (Cassini's width)

%consider the spacecraft as a cylinder 4-meters wide
S = pi*d^2/4;                             %front surface [m^2]

%Weight Constants

c1 = 1e3;                                 %deltaV weight
c2 = 1e8;                                 %position and perigiovio error
%weight
c3 = 1e6;                                 %velocity error weight
c4 = 1e10;                                %perigee error weight

c5 = 1e3;                                 %inclination error weight
c6 = 1e10;                                %aerocapture pericenter error 
%weight
c7 = 1e4;                                 %eccentricity error weight

%Sphere Of Influence (SOI)

SOI_E = a_Earth*(mu_Earth/mu_Sun)^(2/5);  %Earth SOI radius [au]
SOI_J = a_Jupiter*(mu_Jupiter/mu_Sun)^...
        (2/5);                            %Jupiter SOI radius [au]
SOI_S = a_Saturn*(mu_Saturn/mu_Sun)^(2/5);%Saturn SOI radius [au]
SOI_T = a_Titan*(mu_Titan/mu_Saturn)^...
        (2/5);                            %Titan SOI radius [au]

%Saturn System TSI

TSI_E = 1361;                             %total solar irradiance at 1 au 
%[W/m^2]
TSI_S = TSI_E*(a_Earth/a_Saturn)^2;       %TSI at Saturn System distance
%[W/m^2]

%Titan features

T_Titan = T_Titan*JD;                     %Titan rotation period [s]
omega_T = 2*pi/T_Titan;                   %Titan angular velocity [rad/s]

J2_Titan = 33.089e-6;                     %J2 Titan
J3_Titan = -0.179e-6;                     %J3 Titan
J4_Titan = -1.077e-6;                     %J4 Titan

%% Interplanetary Transfer

%1st Optimization

rng default
if max(size(gcp)) == 0
    parpool('local')
end

nvar1 = 8;                                %problem number of variables                    

% X_opt1 = X_best1;                       %bypass the optimization process
[ X_opt1, J1 ] = particleswarm (@(X_try)cost_transfer(X_try, X0,...
                 mu_Sun, mu_Earth, mu_Jupiter, a_ell, c1, c2, c3, c4,...
                 options2), nvar1, lb1, ub1, options1);
           
dt_wait = X_opt1(1);                      %ephemerides-mission start 
%interval [year*(2*pi)]
dt1 = X_opt1(2);                          %1st trajectory arc flight time 
%[year*(2*pi)]
DV_v = X_opt1(3);                         %impulse in velocity direction 
%[EOS]
DV_n = X_opt1(4);                         %impulse in normal direction 
%[EOS]
dt2 = X_opt1(5);                          %2nd trajectory arc flight time
%[year*(2*pi)]
dt3 = X_opt1(6);                          %3rd trajectory arc flight time 
%[year*(2*pi)]
perigee = X_opt1(7);                      %target perigee [au]
perigiovio = X_opt1(8);                   %target perigiovio [au]

%propagate the planets' orbits till the mission starting epoch
[ ~, dX_planets ] = ode113 (@Prop_Orbits, [ 0 dt_wait ], X0, options2,...
                    mu_Sun);

r0_Earth = dX_planets(end, 1:3);          %Earth position at the mission 
%starting epoch [au]
v0_Earth = dX_planets(end, 4:6);          %Earth velocity at the mission 
%starting epoch [EOS]

r0_Jupiter = dX_planets(end, 7:9);        %Jupiter position at the mission 
%starting epoch [au]
v0_Jupiter = dX_planets(end, 10:12);      %Jupiter velocity at the mission
%starting epoch [EOS]

r0_Saturn = dX_planets(end, 13:15);       %Saturn position at the mission 
%starting epoch [au]
v0_Saturn = dX_planets(end, 16:18);       %Saturn velocity at the mission 
%starting epoch [EOS]

X_start_planets = [ r0_Earth'; v0_Earth'; r0_Jupiter'; v0_Jupiter';...
                  r0_Saturn'; v0_Saturn' ];

%determine spacecraft starting conditions
r0_HIF_sc = [ r0_Earth(1);
              r0_Earth(2);
              r0_Earth(3) ];              %spacecraft starting position in 
%HIF [au]
v0_HIF_sc = [ v0_Earth(1);
              v0_Earth(2);
              v0_Earth(3) ];              %spacecraft starting velocity in
%HIF [EOS]

v_ell = sqrt(2*(mu_Sun/norm(r0_HIF_sc) ...
        - mu_Sun/(2*a_ell)));             %required velocity to enter the 
%elliptic orbit [EOS]
v_inf = v_ell - norm(v0_HIF_sc);          %required hyperbolic excess 
%velocity [EOS]

v_inf_VNB = [ v_inf; 0; 0 ];              %hyperbolic excess in the VNB 
%reference frame [EOS]
%VNB reference frame is centered in the spacecraft center of mass, with:
%x-axis along the spacecraft velocity (V)
%y-axis along the spacecraft orbit normal, i.e. the angular momentum (N)
%z-axis along the spacecraft orbit binormal (B)

%define the hyperbolic excess velocity in HIF
r_vers = r0_HIF_sc/norm(r0_HIF_sc);       %spacecraft position versor
v_vers = v0_HIF_sc/norm(v0_HIF_sc);       %spacecraft velocity versor
n_vers = cross(r_vers, v_vers);           %normal versor
b_vers = cross(v_vers, n_vers);           %binormal versor

ROT_VNB2HIF = [ v_vers, n_vers, b_vers ]; %rotation matrix VNB-HIF

v_inf_HIF = ROT_VNB2HIF*v_inf_VNB;        %hyperbolic excess in the HIF 
%[EOS]

%update spacecraft starting velocity
v0_HIF_sc = v0_HIF_sc + v_inf_HIF;        %new spacecraft starting velocity
%[EOS]

%full mission starting conditions (Earth + spacecraft)
X_start1 = [ X_start_planets(1:6); r0_HIF_sc; v0_HIF_sc ];

%1st trajectory arc (till the DSM, "deep space maneuver")
[ t1, dX1 ] = ode113 (@Prop_Arc1, [ 0 dt1 ], X_start1, options2, mu_Sun);

r1_Earth = dX1(:, 1:3);                   %Earth positions [au]
v1_Earth = dX1(:, 4:6);                   %Earth velocities [EOS]

r1_sc = dX1(:, 7:9);                      %spacecraft positions [au]
v1_sc = dX1(:, 10:12);                    %spacecraft velocities [EOS]

%Deep Space Maneuver (DSM)
DV2_VNB = [ DV_v; DV_n; 0 ];              %DSM's deltaV in the VNB 
%reference frame [EOS]

%define the DSM deltaV in HIF
r_vers = r1_sc(end, :)'/...
         norm(r1_sc(end, :));             %spacecraft position versor
v_vers = v1_sc(end, :)'/...
         norm(v1_sc(end, :));             %spacecraft velocity versor
n_vers = cross(r_vers, v_vers);           %normal versor
b_vers = cross(v_vers, n_vers);           %binormal versor

ROT_VNB2HIF = [ v_vers, n_vers, b_vers ]; %rotation matrix VNB-HIF

DV2_HIF = ROT_VNB2HIF*DV2_VNB;            %DSM's deltaV in the HIF [EOS]

%update spacecraft velocity after DSM
v1_sc = v1_sc(end, :)' + DV2_HIF;         %spacecraft velocity after the 
%maneuver [EOS]

%2nd trajectory arc starting conditions (Earth + spacecraft)
X_start2 = [ r1_Earth(end, :)'; v1_Earth(end, :)'; r1_sc(end, :)'; v1_sc ];

%2nd trajectory arc (till after the Earth Gravity Assist)
[ t2, dX2 ] = ode113 (@Prop_Perturbed, [ 0 dt2 ], X_start2, options2,...
              mu_Sun, mu_Earth);

r2_Earth = dX2(:, 1:3);                   %Earth positions [au]
v2_Earth = dX2(:, 4:6);                   %Earth velocities [EOS]

r2_sc = dX2(:, 7:9);                      %spacecraft positions [au]
v2_sc = dX2(:, 10:12);                    %spacecraft velocities [EOS]

%preallocation
r_sc2E = zeros(length(t2), 3);
v_sc2E = zeros(length(t2), 3);
dist_sc2E = zeros(length(t2), 1);

for j = 1:length(t2)
    r_sc2E(j, :) = r2_sc(j, :) -...
                   r2_Earth(j, :);        %spacecraft positions w.r.t. 
    %Earth [au]
    v_sc2E(j, :) = v2_sc(j, :) -...
                   v2_Earth(j, :);        %spacecraft velocities w.r.t. 
    %Earth [EOS]
    dist_sc2E(j) = norm(r_sc2E(j, :));    %spacecraft distance from Earth 
    %C.M. [au]
end

rp_E = min(dist_sc2E);                    %perigee radius during Earth 
%gravity assist [au]
ind_E = find(dist_sc2E == rp_E);          %perigee index

%propagate Jupiter orbit before the 3rd trajectory arc
[ ~, dX_Jupiter ] = ode113 (@Prop_Keplerian, [ 0 dt1+dt2 ],...
                    X_start_planets(7:12), options2, mu_Sun);

r_Jupiter = dX_Jupiter(:, 1:3);           %Jupiter positions [au]
v_Jupiter = dX_Jupiter(:, 4:6);           %Jupiter velocities [EOS]

%3rd trajectory arc starting conditions (Jupiter + spacecraft)
X_start3 = [ r_Jupiter(end, :)'; v_Jupiter(end, :)'; r2_sc(end, :)'; ...
           v2_sc(end, :)' ];

%3rd trajectory arc (till after the Jupiter Gravity Assist)
[ t3, dX3 ] = ode113 (@Prop_Perturbed, [ 0 dt3 ], X_start3, options2,...
              mu_Sun, mu_Jupiter);

r3_Jupiter = dX3(:, 1:3);                 %Jupiter positions [au]
v3_Jupiter = dX3(:, 4:6);                 %Jupiter velocities [EOS]

r3_sc = dX3(:, 7:9);                      %spacecraft positions [au]
v3_sc = dX3(:, 10:12);                    %spacecraft velocities [EOS]

%preallocation
r_sc2J = zeros(length(t3), 3);
v_sc2J = zeros(length(t3), 3);
dist_sc2J = zeros(length(t3), 1);

for j = 1:length(t3)
    r_sc2J(j, :) = r3_sc(j, :) -...
                   r3_Jupiter(j, :);      %spacecraft positions w.r.t. 
    %Jupiter [au]
    v_sc2J(j, :) = v3_sc(j, :) -...
                   v3_Jupiter(j, :);      %spacecraft velocities w.r.t. 
    %Jupiter [EOS]
    dist_sc2J(j) = norm(r_sc2J(j, :));    %spacecraft distance from Jupiter
    %C.M. [au]
end

rp_J = min(dist_sc2J);                    %perigiovio radius during Jupiter
%gravity assist [au]
ind_J = find(dist_sc2J == rp_J);          %perigiovio index

%propagate Saturn orbit to set target position and velocity
[ ~, dX_Saturn ] = ode113 (@Prop_Keplerian, [ 0 dt1+dt2+dt3 ],...
                   X_start_planets(13:18), options2, mu_Sun);

r_Saturn = dX_Saturn(:, 1:3);             %Saturn positions [au]
v_Saturn = dX_Saturn(:, 4:6);             %Saturn velocities [EOS]

r_HIF_sc = [ r1_sc; r2_sc; r3_sc ];       %all spacecraft transfer 
%positions in HIF [au]
r_Earth = [ r1_Earth; r2_Earth ];         %all Earth positions [au]
r_Jupiter = [ r_Jupiter; r3_Jupiter ];    %all Jupiter positions [au]
t = [ t1; t2; t3 ];                       %transfer trajectory time steps 
%[year*(2*pi)]

vel_err = v3_sc(end, :) - ...
          v_Saturn(end, :);               %velocity errors [EOS]

%% Departure Phase

h_cpo = 185;                              %circular parking orbit height 
%[km]
r_cpo = (R_Earth + h_cpo)/au;             %circular parking orbit radius 
%[au]
v_cpo = sqrt(mu_Earth/r_cpo);             %circular parking orbit velocity
%[EOS]

a_hyp = - mu_Earth/(v_inf^2);             %hyperbola semi-major axis [au]
e_hyp = 1 - r_cpo/a_hyp;                  %hyperbola eccentricity

F2 = acosh((a_hyp - SOI_E)/(a_hyp*e_hyp));%hyperbolic anomaly at Earth SOI
%exit [rad]
F1 = 0;                                   %starting hyperbolic anomaly (at
%the perigee) [rad]
dt_dep = sqrt(- a_hyp^3/mu_Earth)*((e_hyp*sinh(F2) - ...
         F2) - (e_hyp*sinh(F1) - F1));    %time to reach the SOI border 
%[year*(2*pi)]

vp_hyp = sqrt(v_inf^2 + 2*mu_Earth/r_cpo);%hyperbola perigee velocity [EOS]
DV1_VNB = vp_hyp - v_cpo;                 %departure deltaV in the VNB 
%frame [EOS]

theta = atan2(r0_Earth(2), r0_Earth(1)) - ...
        pi/2 + acos(1/e_hyp);             %departure maneuver angle [rad]

r0_cpo = [ r_cpo*cos(theta);
           r_cpo*sin(theta);
           0 ];                           %spacecraft geocentric starting 
%position [au]
v0_cpo = [ - v_cpo*sin(theta);
           v_cpo*cos(theta);
           0 ];                           %spacecraft geocentric starting
%velocity [EOS]
       
T_cpo = 2*pi*sqrt(r_cpo^3/mu_Earth);      %circular parking orbit period 
%[year*(2*pi)]

[ ~, dX_cpo ] = ode113 (@Prop_Keplerian, [ 0 T_cpo ], [ r0_cpo;...
                v0_cpo ], options2, mu_Earth);

r_cpo = dX_cpo(:, 1:3);                   %circular parking orbit positions
%[au]
v_cpo = dX_cpo(:, 4:6);                   %circular parking orbit 
%velocities [EOS]

%define deltaV in geocentric frame
r_vers = r0_cpo/norm(r0_cpo);             %spacecraft position versor
v_vers = v0_cpo/norm(v0_cpo);             %spacecraft velocity versor
n_vers = cross(r_vers, v_vers);           %normal versor
b_vers = cross(v_vers, n_vers);           %binormal versor

ROT_VNB2geo = [ v_vers, n_vers, b_vers ]; %rotation matrix VNB-geocentric 
%frame

DV1_geo = ROT_VNB2geo*[ DV1_VNB; 0; 0 ];  %departure deltaV in geocentric 
%frame [EOS]

v0_hyp = v0_cpo + DV1_geo;                %hyperbola starting velocity 
%[EOS]

[ ~, dX_hyp ] = ode113 (@Prop_Keplerian, [ 0 1.5*T_cpo ], [ r0_cpo;...
                v0_hyp ], options2, mu_Earth);

r_hyp = dX_hyp(:, 1:3);                   %departure hyperbola positions 
%[au]
v_hyp = dX_hyp(:, 4:6);                   %departure hyperbola velocities
%[EOS]

%estimating the total mass of the spacecraft launched by Delta IV Heavy
[ M0 ] = Mass_sc (mu_Earth, JD/conv_time, norm(r0_cpo), ...
         norm(DV1_geo));                  %launcher total capability in 
%terms of sc mass [kg] 
     
%class 2 mass contingency in proposal stage
margin = M0*18/100;                       %mass margin [kg]
M0_est = M0 - margin;                     %best estimate of spacecraft 
%initial mass [kg]

%% Capture Phase

%determine spacecraft mass at the beginning of the capture phase
DV2 = norm(DV2_HIF)/conv_v;               %DSM deltaV [km/s]
Mp_DSM = M0_est*(1 - exp(- DV2*1e3/...
         (Isp*g0)));                      %propellant mass needed for DSM
%[kg]

M0_DSM = M0_est - Mp_DSM;                 %spacecraft total mass after DSM
%[kg]

%propagate Titan orbit till the arrival epoch
[ ~, dX_Titan ] = ode113 (@Prop_Perturbed, [ 0 dt_wait+dt1+dt2+dt3 ],...
                  [ X0(13:18); R_HIF_Titan; V_HIF_Titan ], options2, ...
                  mu_Sun, mu_Saturn);
              
r0_Titan = dX_Titan(end, 7:9);            %Titan position in HIF at the 
%arrival epcoh [au]
v0_Titan = dX_Titan(end, 10:12);          %Titan velocity in HIF at the 
%arrival epoch [EOS]

r0_Titan = r0_Titan - r_Saturn(end, :);   %Titan position in planetocentric
%frame [au]
v0_Titan = v0_Titan - v_Saturn(end, :);   %Titan velocity in planetocentric
%frame [EOS]

X_start_Titan = [ r0_Titan'; v0_Titan' ];

%2nd Optimization

nvar2 = 3;                                %problem number of variables                    

% X_opt2 = X_best2;                       %bypass the optimization process
[ X_opt2, J2 ] = particleswarm (@(X_try)cost_aerocapture(X_try, vel_err,...
                 X_start_Titan, SOI_S, mu_Saturn, mu_Titan, R_Titan, au,...
                 Cd, S, M0_DSM, i_Titan, c5, c6, c7, options2), nvar2,...
                 lb2, ub2, options3);

phi = X_opt2(1);                          %latitude [rad]
lambda = X_opt2(2);                       %longitude [rad]
dt4 = X_opt2(3);                          %time of flight [year*(2*pi)]

%starting spacecraft position on the Saturn SOI
x0_sc = SOI_S*cos(phi)*cos(lambda);
y0_sc = SOI_S*cos(phi)*sin(lambda);
z0_sc = SOI_S*sin(phi);

X_start4 = [ X_start_Titan; x0_sc; y0_sc; z0_sc; vel_err' ];

[ t4, dX4 ] = ode113 (@Prop_Capture, [ 0 dt4 ], X_start4, options2,...
              mu_Saturn, mu_Titan, R_Titan, au, Cd, S, M0_DSM);
                
r4_Titan = dX4(:, 1:3);                   %Titan positions [au]
v4_Titan = dX4(:, 4:6);                   %Titan velocities [EOS]

r4_sc = dX4(:, 7:9);                      %spacecraft positions [au]
v4_sc = dX4(:, 10:12);                    %spacecraft velocities [EOS]

%preallocation
r4_sc2T = zeros(length(t4), 3);
v4_sc2T = zeros(length(t4), 3);
dist_sc2T = zeros(length(t4), 1);
vel_sc2T = zeros(length(t4), 1);

for j = 1:length(t4)
    r4_sc2T(j, :) = r4_sc(j, :) -...
                    r4_Titan(j, :);       %spacecraft positions w.r.t. 
    %Titan [au]
    v4_sc2T(j, :) = v4_sc(j, :) -...
                    v4_Titan(j, :);       %spacecraft velocities w.r.t. 
    %Titan [EOS]
    dist_sc2T(j) = norm(r4_sc2T(j, :));   %spacecraft distance from Titan 
    %C.M. [au]
    vel_sc2T(j) = norm(v4_sc2T(j, :));    %spacecraft speed w.r.t. Titan 
    %[EOS]
end

rp_T = min(dist_sc2T);                    %pericenter radius w.r.t. Titan
%[au]

ind_T = find(dist_sc2T == rp_T);          %pericenter index
ind_enter = find(dist_sc2T <= (1300 + R_Titan)/au); 
ind_enter = ind_enter(1);                 %atmospheric entry index

[ a_c, e_c, i_c, ~, ~, TA1 ] = rv2coe (r4_sc2T(end, :)', ...
                               v4_sc2T(end, :)', mu_Titan);
i_c = rad2deg(i_c + i_Titan);             %capture orbit inclination [deg]
EA1 = 2*atan(sqrt((1 - e_c)/(1 +...
      e_c))*tan(TA1/2));                  %initial eccentric anomaly [rad]

%propagate the elliptic capture orbit till the apocenter, where the
%maneuver will be performed
TA2 = pi;                                 %final true anomaly [rad] 
EA2 = 2*atan(sqrt((1 - e_c)/(1 +...
      e_c))*tan(TA2/2));                  %final eccentric anomaly [rad]

dt5 = sqrt(a_c^3/mu_Titan)*((EA2 - e_c*sin(EA2)) - (EA1 -...
      e_c*sin(EA1)));                     %flight time till the ellipse 
%apocenter [year*(2*pi)]

X_start5 = [ r4_Titan(end, :)'; v4_Titan(end, :)'; r4_sc(end, :)'; ...
           v4_sc(end, :)' ];

[ t5, dX5 ] = ode113 (@Prop_Capture, [ t4(end) t4(end)+dt5 ], X_start5,...
              options2, mu_Saturn, mu_Titan, R_Titan, au, Cd, S, M0_DSM);

r5_Titan = dX5(:, 1:3);                   %Titan positions [au]
v5_Titan = dX5(:, 4:6);                   %Titan velocities [EOS]                  
                  
r5_sc = dX5(:, 7:9);                      %spacecraft positions [au]
v5_sc = dX5(:, 10:12);                    %spacecraft velocities [EOS]

%preallocation
r5_sc2T = zeros(length(t5), 3);
v5_sc2T = zeros(length(t5), 3);

for j = 1:length(t5)
    r5_sc2T(j, :) = r5_sc(j, :) -...
                    r5_Titan(j, :);       %spacecraft positions w.r.t. 
    %Titan [au]
    v5_sc2T(j, :) = v5_sc(j, :) -...
                    v5_Titan(j, :);       %spacecraft velocities w.r.t.
    %Titan [EOS]
end

%Hohmann transfer to reach the final operation orbit
rf_sc2T = (R_Titan + 1400)/au;            %final operation orbit radius
%[au]
vf_sc2T = sqrt(mu_Titan/rf_sc2T);         %final operation orbit velocity 
%[EOS]

a_H = (rf_sc2T + norm(r5_sc2T(end, :)))/2;%Hohmann ellipse semi-major axis
%[au]

va_H = sqrt(2*mu_Titan*(1/norm(r5_sc2T(end, :)) -...
       1/(2*a_H)));                       %Hohmann ellipse apocenter 
%velocity [EOS]
vp_H = sqrt(2*mu_Titan*(1/rf_sc2T -...
       1/(2*a_H)));                       %Hohmann ellipse pericenter
%velocity [EOS]

DV1_H = va_H - norm(v5_sc2T(end, :));     %first Hohmann transfer impulse 
%[EOS]
DV3 = norm(DV1_H)/conv_v;                 %first Hohmann transfer impulse 
%[km/s]

Mp_H1 = M0_DSM*(1 - exp(- DV3*1e3/...
        (Isp*g0)));                       %propellant mass needed for DV1_H
%[kg]
M0_H1 = M0_DSM - Mp_H1;                   %spacecraft total mass after
%DV1_H [kg]

DV2_H = vf_sc2T - vp_H;                   %second Hohmann transfer impulse
%[EOS]
DV4 = norm(DV2_H)/conv_v;                 %second Hohmann transfer impulse 
%[km/s]

Mp_H2 = M0_H1*(1 - exp(- DV4*1e3/...
        (Isp*g0)));                       %propellant mass needed for DV2_H
%[kg]
M0_ins = M0_H1 - Mp_H2;                   %spacecraft total mass after
%DV2_H [kg]
%---> M0_ins is the total mass inserted in the operation orbit

%1st impulse and orbit propagation
DV1_H_VNB = [ DV1_H; 0; 0 ];              %1st Hohmann transfer impulse in
%VNB [EOS]

r_vers = r5_sc2T(end, :)'/...
         norm(r5_sc2T(end, :));           %position versor
v_vers = v5_sc2T(end, :)'/...
         norm(v5_sc2T(end, :));           %velocity versor
n_vers = cross(r_vers, v_vers);           %normal versor
b_vers = cross(v_vers, n_vers);           %binormal versor

ROT_VNB2pf = [ v_vers, n_vers, b_vers ];  %rotation matrix 
%VNB-planetocentric frame
DV1_H_pf = ROT_VNB2pf*DV1_H_VNB;          %1st Hohmann transfer impulse in 
%planetocentric frame

dt6 = pi*sqrt(a_H^3/mu_Titan);            %Hohmann transfer half-period 
%[year*(2*pi)]
X_start6 = [ r5_Titan(end, :)'; v5_Titan(end, :)'; r5_sc(end, :)';...
             v5_sc(end, :)' + DV1_H_pf ];

[ t6, dX6 ] = ode113 (@Prop_Capture, [ t5(end) t5(end)+dt6 ], X_start6,...
              options2, mu_Saturn, mu_Titan, R_Titan, au, Cd, S, M0_H1);
                  
r6_Titan = dX6(:, 1:3);                   %Titan positions [au]
v6_Titan = dX6(:, 4:6);                   %Titan velocities [EOS]
                   
r6_sc = dX6(:, 7:9);                      %spacecraft positions [au]
v6_sc = dX6(:, 10:12);                    %spacecraft velocities [EOS]

%preallocation
r6_sc2T = zeros(length(t6), 3);
v6_sc2T = zeros(length(t6), 3);

for j = 1:length(t6)
    r6_sc2T(j, :) = r6_sc(j, :) -...
                    r6_Titan(j, :);       %spacecraft positions w.r.t. 
    %Titan [au]
    v6_sc2T(j, :) = v6_sc(j, :) -...
                    v6_Titan(j, :);       %spacecraft velocities w.r.t.
    %Titan [EOS]
end

%2nd impulse and orbit propagation
DV2_H_VNB = [ DV2_H; 0; 0 ];              %2nd Hohmann transfer impulse in 
%VNB [EOS]

r_vers = r6_sc2T(end, :)'/...
         norm(r6_sc2T(end, :));           %position versor
v_vers = v6_sc2T(end, :)'/...
         norm(v6_sc2T(end, :));           %velocity versor
n_vers = cross(r_vers, v_vers);           %normal versor
b_vers = cross(v_vers, n_vers);           %binormal versor

ROT_VNB2pf = [ v_vers, n_vers, b_vers ];  %rotation matrix 
%VNB-planetocentric frame
DV2_H_pf = ROT_VNB2pf*DV2_H_VNB;          %2nd Hohmann transfer impulse in
%planetocentric frame

dt7 = 2*pi*sqrt(rf_sc2T^3/mu_Titan);      %final operation orbit period 
%[year*(2*pi)]
X_start7 = [ r6_Titan(end, :)'; v6_Titan(end, :)'; r6_sc(end, :)';...
           v6_sc(end, :)' + DV2_H_pf ];

[ t7, dX7 ] = ode113 (@Prop_Capture, [ t6(end) t6(end)+dt7 ], X_start7,...
              options2, mu_Saturn, mu_Titan, R_Titan, au, Cd, S, M0_ins);
                  
r7_Titan = dX7(:, 1:3);                   %Titan positions [au]
v7_Titan = dX7(:, 4:6);                   %Titan velocities [EOS]
                  
r7_sc = dX7(:, 7:9);                      %spacecraft positions [au]
v7_sc = dX7(:, 10:12);                    %spacecraft velocities [EOS]

%preallocation
r7_sc2T = zeros(length(t7), 3);
v7_sc2T = zeros(length(t7), 3);

for j = 1:length(t7)
    r7_sc2T(j, :) = r7_sc(j, :) -...
                    r7_Titan(j, :);       %spacecraft positions w.r.t. 
    %Titan [au]
    v7_sc2T(j, :) = v7_sc(j, :) -...
                    v7_Titan(j, :);       %spacecraft velocities w.r.t. 
    %Titan [EOS]
end

%spacecraft ground track during capture phase
t_capture = [ t5; t6; t7 ]*conv_time;     %times [s]
r_capture = [ r5_sc2T(:, :); r6_sc2T(:, :);...
           r7_sc2T(:, :) ];               %positions [au]
v_capture = [ v5_sc2T(:, :); v6_sc2T(:, :);...
           v7_sc2T(:, :) ];               %velocities [EOS]

[ lat_capture, long_capture ] = GroundTrack (t_capture, r_capture,...
                                v_capture, mu_Titan, i_Titan, omega_T);
lat_capture = rad2deg(lat_capture);       %latitude [deg]
long_capture = rad2deg(long_capture);     %longitude [deg]

%% Operational Phase

[ t_op, dX_op ] = ode113 (@Prop_OP, [ 0 .5*2*pi ], X_start7, options2,...
                  mu_Saturn, mu_Titan, R_Titan, au, Cd, S, M0_ins, ...
                  J2_Titan, J3_Titan, J4_Titan);
                  
r_op_Titan = dX_op(:, 1:3);               %Titan positions [au]
v_op_Titan = dX_op(:, 4:6);               %Titan velocities [EOS]

r_op_sc = dX_op(:, 7:9);                  %spacecraft positions [au]
v_op_sc = dX_op(:, 10:12);                %spacecraft velocities [EOS]

%preallocation
r_op_sc2T = zeros(length(t_op), 3);
v_op_sc2T = zeros(length(t_op), 3);
a_op = zeros(length(t_op), 1);
e_op = zeros(length(t_op), 1);
i_op = zeros(length(t_op), 1);
RAAN_op = zeros(length(t_op), 1);
PA_op = zeros(length(t_op), 1);
TA_op = zeros(length(t_op), 1);

for j = 1:length(t_op)
    r_op_sc2T(j, :) = r_op_sc(j, :) -...
                      r_op_Titan(j, :);   %spacecraft positions w.r.t.
    %Titan [au]
    v_op_sc2T(j, :) = v_op_sc(j, :) -...
                      v_op_Titan(j, :);   %spacecraft velocities w.r.t. 
    %Titan [EOS]
    [ a_op(j), e_op(j), i_op(j), RAAN_op(j), PA_op(j), TA_op(j) ] = ...
    rv2coe (r_op_sc2T(j, :)', v_op_sc2T(j, :)', mu_Titan);
end

i_op = i_op + i_Titan;                    %inclination w.r.t. Titan 
%rotation axis [rad]

%spacecraft ground track during the operational phase 
[ lat_op, long_op ] = GroundTrack (t_op*conv_time, r_op_sc2T(:, :),...
                      v_op_sc2T(:, :), mu_Titan, i_Titan, omega_T);
lat_op = rad2deg(lat_op);                 %latitude [deg]
long_op = rad2deg(long_op);               %longitude [deg]

ind_week = find(t_op <= 2*pi/12/30*7);
ind_week = ind_week(end);                 %first week index

%Station Keeping Maneuvers

p_sk = a_op(end)*(1 - e_op(end)^2);       %semilatus rectum after 6 months 
%[au]
TA_sk = acos((p_sk/rf_sc2T - 1)*...
        1/e_op(end));                     %true anomaly at which perform 
%the maneuver [rad]
V_sk = [ - (mu_Titan/p_sk)^(1/2)*sin(TA_sk);
         (mu_Titan/p_sk)^(1/2)*(e_op(end) + cos(TA_sk));
         0 ];                             %velocity vector in perifocal 
%frame at TA_sk [EOS]

%change in semi-major axis and eccentricity
par = rf_sc2T*norm(V_sk)/(rf_sc2T*...
      vf_sc2T);                           %cosine of the angle between
%velocity vectors of final 
                                          %and initial orbit
DV_ae = sqrt(vf_sc2T^2 + norm(V_sk)^2 - 2*vf_sc2T*...
        norm(V_sk)*par);                  %first required deltaV [EOS]

%change in inclination
DV_i = 2*norm(V_sk)*sin((i_op(1) -...
       i_op(end))/2);                     %second required deltaV [EOS]
     
DV_SK = norm(DV_ae) + norm(DV_i);         %total deltaV for station keeping
%[EOS]
DV_SK = DV_SK/conv_v;                     %total deltaV for station keeping
%[km/s]

%for 3-years mission it is necessary to perform such maneuvers 6 times
%(every 6 months) and this allows to estimate the spacecraft dry mass
Mp_SK1 = M0_ins*(1 - exp(- DV_SK*1e3/...
         (Isp*g0)));                      %propellant mass needed for 1st
%station keeping [kg]
M0_SK1 = M0_ins - Mp_SK1;                 %s/c total mass after 1st station
%keeping [kg]

Mp_SK2 = M0_SK1*(1 - exp(- DV_SK*1e3/...
         (Isp*g0)));                      %propellant mass needed for 2nd 
%station keeping [kg]
M0_SK2 = M0_SK1 - Mp_SK2;                 %s/c total mass after 2nd station
%keeping [kg]

Mp_SK3 = M0_SK2*(1 - exp(- DV_SK*1e3/...
         (Isp*g0)));                      %propellant mass needed for 3rd 
%station keeping [kg]
M0_SK3 = M0_SK2 - Mp_SK3;                 %s/c total mass after 3rd station
%keeping [kg]

Mp_SK4 = M0_SK3*(1 - exp(- DV_SK*1e3/...
         (Isp*g0)));                      %propellant mass needed for 4th 
%station keeping [kg]
M0_SK4 = M0_SK3 - Mp_SK4;                 %s/c total mass after 4th station
%keeping [kg]

Mp_SK5 = M0_SK4*(1 - exp(- DV_SK*1e3/...
         (Isp*g0)));                      %propellant mass needed for 5th
%station keeping [kg]
M0_SK5 = M0_SK4 - Mp_SK5;                 %s/c total mass after 5th station
%keeping [kg]

Mp_SK6 = M0_SK5*(1 - exp(- DV_SK*1e3/...
         (Isp*g0)));                      %propellant mass needed for 6th 
%station keeping [kg]
Mf = M0_SK5 - Mp_SK6;                     %spacecraft final mass [kg]

%total propellant mass
Mp_tot = Mp_DSM + Mp_H1 + Mp_H2 + Mp_SK1 + Mp_SK2 + ...
         Mp_SK3 + Mp_SK4 + Mp_SK5 + Mp_SK6;        

%% Results

rpE_perc = rp_E/SOI_E*1e2;                %flyby perigee in percentage 
%w.r.t. SOI
rpJ_perc = rp_J/SOI_J*1e2;                %flyby perigiovio in percentage 
%w.r.t. SOI

rp_E = rp_E*au;                           %Earth flyby pericenter radius 
%[km]
rp_J = rp_J*au;                           %Jupiter flyby pericenter radius 
%[km]
rp_T = rp_T*au;                           %Titan aerocapture pericenter 
%radius [km]

pos_err = (r3_sc(end, :) - r_Saturn(end, :))*...
          au;                             %position errors [km]

DV1 = norm(DV1_geo)/conv_v;               %departure deltaV [km/s]

T_Eph = datenum('01-May-2021 00:00:00');
Date_Start = addtodate(T_Eph, dt_wait*conv_time, 'second');
Date_Launch = addtodate(Date_Start, - dt_dep*conv_time, 'second');
Date_DSM = addtodate(Date_Start, t1(end)*conv_time, 'second');
Date_EGA = addtodate(Date_DSM, t2(ind_E)*conv_time, 'second');
Date_JGA = addtodate(Date_EGA, (t2(end) - t2(ind_E) + t3(ind_J))*...
           conv_time, 'second');
Date_Arrival = addtodate(Date_Start, t(end)*conv_time, 'second');
Date_Capture = addtodate(Date_Arrival, dt4*conv_time, 'second');
Date_Op = addtodate(Date_Capture, (dt5 + dt6)*conv_time, 'second');

%print results
disp('------------------------------------------------------------------');
disp('INTERPLANETARY TRANSFER ERRORS');
fprintf('Position Error: %.4f km \n', norm(pos_err));
fprintf('X-axis Velocity Error: %.4f km/s \n', vel_err(1)/conv_v);
fprintf('Y-axis Velocity Error: %.4f km/s \n', vel_err(2)/conv_v);
fprintf('Z-axis Velocity Error: %.4f km/s \n', vel_err(3)/conv_v);
disp('------------------------------------------------------------------');
disp('LAUNCH DATE')
disp(datestr(Date_Launch));
disp('DEPARTURE DATE')
disp(datestr(Date_Start));
disp('DEEP SPACE MANEUVER DATE')
disp(datestr(Date_DSM));
disp('EARTH GRAVITY ASSIST DATE')
disp(datestr(Date_EGA));
disp('JUPITER GRAVITY ASSIST DATE')
disp(datestr(Date_JGA));
disp('SATURN SYSTEM ARRIVAL DATE')
disp(datestr(Date_Arrival));
disp('ELLIPTIC CAPTURE ORBIT INJECTION DATE')
disp(datestr(Date_Capture));
disp('OPERATION ORBIT INJECTION DATE')
disp(datestr(Date_Op));
disp('------------------------------------------------------------------');
disp('CLOSEST APPROACHES');
fprintf('Perigee Height: %.4f km \n', rp_E - R_Earth);
fprintf('---> %.4f percent SOI \n', rpE_perc);
fprintf('Perigiovio Height: %.4f km \n', rp_J - R_Jupiter);
fprintf('---> %.4f percent SOI \n', rpJ_perc);
fprintf('Aerocapture Height: %.4f km \n', rp_T - R_Titan);
disp('------------------------------------------------------------------');
disp('DV BUDGET');
fprintf('Earth departure: %.4f km/s \n', DV1);
fprintf('Deep Space Maneuver: %.4f km/s \n', DV2);
fprintf('1st Hohmann transfer impulse: %.4f km/s \n', DV3);
fprintf('2nd Hohmann transfer impulse: %.4f km/s \n', DV4);
fprintf('Station Keeping Maneuver: %.4f km/s \n', DV_SK);
disp('------------------------------------------------------------------');
disp('MASS BUDGET');
fprintf('Mass Margin: %.4f kg \n', margin);
fprintf('Best Estimate S/C Mass: %.4f kg \n', M0_est);
fprintf('Total Propellant Mass: %.4f kg \n', Mp_tot);
fprintf('Total Dry Mass + Attitude Control Propellant: %.4f kg \n', Mf);
disp('------------------------------------------------------------------');
fprintf('TSI at Saturn System: %.4f W/m^2 \n', TSI_S);
disp('It is necessary a RTG for the onboard power generation');
disp('------------------------------------------------------------------');

%% Plot

%heliocentric trajectory
figure
plot3(r_HIF_sc(: ,1), r_HIF_sc(:, 2), r_HIF_sc(:, 3), 'r-', 'LineWidth',...
     1.2)
hold on
plot3(r_Earth(:, 1), r_Earth(:, 2), r_Earth(:, 3), 'b-', 'LineWidth', 1.2)
hold on
plot3(r_Jupiter(:, 1), r_Jupiter(:, 2), r_Jupiter(:, 3), 'y-',...
     'LineWidth', 1.2)
hold on
plot3(r_Saturn(:, 1), r_Saturn(:, 2), r_Saturn(:, 3), 'g-', 'LineWidth',...
     1.2)
hold on
quiver3(r1_sc(end, 1), r1_sc(end, 2), r1_sc(end, 3), DV2_HIF(1)*80,...
       DV2_HIF(2)*80, DV2_HIF(3)*80, 'k', 'LineWidth', 1.2)
axis equal
grid on
set(gca,'TickLabelInterpreter','Latex')
xlabel('x [au]')
ylabel('y [au]')
zlabel('z [au]')
title('Heliocentric Transfer Trajectory')
legend('S/C', 'Earth', 'Jupiter', 'Saturn', '\DeltaV DSM', 'Location',...
      'Best')
hold off

%Earth gravity assist trajectory
figure 
plot3(r_sc2E(ind_E-70:ind_E+100, 1)*au, r_sc2E(ind_E-70:ind_E+100, 2)...
     *au, r_sc2E(ind_E-70:ind_E+100, 3)*au, 'r-', 'LineWidth', 1.2)
hold on
quiver3(0, 0, 0, v2_Earth(ind_E, 1)/conv_v*1e3, v2_Earth(ind_E, 2)/...
       conv_v*1e3, v2_Earth(ind_E, 3)/conv_v*1e3, 'b', 'LineWidth', 1.2)
hold on
quiver3(r_sc2E(ind_E, 1)*au, r_sc2E(ind_E, 2)*au, r_sc2E(ind_E, 3)*au,...
       v_sc2E(ind_E, 1)/conv_v*1e3, v_sc2E(ind_E, 2)/conv_v*1e3, ...
       v_sc2E(ind_E, 3)/conv_v*1e3, 'r', 'LineWidth', 1.2)
hold on
Earth = imread('Earth.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Earth;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Earth, R_Earth,...
          R_Earth, 50);
h_E = surface(-X,-Y,-Z,props);
set(gca,'TickLabelInterpreter','Latex')
view(3)
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Earth Gravity Assist')
legend('S/C', 'V Earth', 'S/C Direction', 'Location', 'Best')
hold off

%Jupiter gravity assist trajectory
figure 
plot3(r_sc2J(ind_J-70:ind_J+70, 1)*au, r_sc2J(ind_J-70:ind_J+70, 2)*au,...
     r_sc2J(ind_J-70:ind_J+70, 3)*au, 'r-', 'LineWidth', 1.2)
hold on
quiver3(0, 0, 0, v3_Jupiter(ind_J, 1)/conv_v*3e5, v3_Jupiter(ind_J, 2)/...
       conv_v*3e5, v3_Jupiter(ind_J, 3)/conv_v*3e5, 'y', 'LineWidth', 1.2)
hold on
quiver3(r_sc2J(ind_J, 1)*au, r_sc2J(ind_J, 2)*au, r_sc2J(ind_J, 3)*au,...
       v_sc2J(ind_J, 1)/conv_v*3e5, v_sc2J(ind_J, 2)/conv_v*3e5, ...
       v_sc2J(ind_J, 3)/conv_v*3e5, 'r', 'LineWidth', 1.2)
hold on
Jupiter = imread('Jupiter.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Jupiter;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Jupiter,...
          R_Jupiter, R_Jupiter, 50);
h_J = surface(-X,-Y,-Z,props);
set(gca,'TickLabelInterpreter','Latex')
view(3)
axis equal
axis([ -inf inf -1e7 1e7 -inf inf ])
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Jupiter Gravity Assist')
legend('S/C', 'V Jupiter', 'S/C Direction', 'Location', 'Best')
hold off

%geocentric departure trajectory
figure
plot3(r_hyp(:, 1)*au, r_hyp(:, 2)*au, r_hyp(:, 3)*au, 'r-', 'LineWidth',...
     1.2)
hold on
quiver3(0, 0, 0, v0_Earth(1)/conv_v*1e3, v0_Earth(2)/conv_v*1e3,...
       v0_Earth(3)/conv_v*1e3, 'b', 'LineWidth', 1.2)
hold on
quiver3(r0_cpo(1)*au, r0_cpo(2)*au, r0_cpo(3)*au, v0_hyp(1)/conv_v*1e3,...
       v0_hyp(2)/conv_v*1e3, v0_hyp(3)/conv_v*1e3, 'k', 'LineWidth', 1.2)
hold on
plot3(r_cpo(:, 1)*au, r_cpo(:, 2)*au, r_cpo(:, 3)*au, 'r-', 'LineWidth',...
     1.2)
hold on
Earth = imread('Earth.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Earth;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Earth, R_Earth,...
          R_Earth, 50);
h_geo = surface(-X,-Y,-Z,props);
set(gca,'TickLabelInterpreter','Latex')
view(3)
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Departure Trajectory')
legend('S/C', 'V Earth', '\DeltaV Escape', 'Location', 'Best')
hold off

%rendezvous trajectory
figure
plot3(r4_sc(:, 1)*au, r4_sc(:, 2)*au, r4_sc(:, 3)*au, 'r-', 'LineWidth',...
     1.2)
hold on
plot3(r4_Titan(:, 1)*au, r4_Titan(:, 2)*au, r4_Titan(:, 3)*au, 'color',...
     [ 1 0.64453125 0 ], 'LineWidth',1.2) 
hold on
quiver3(0, 0, 0, v_Saturn(end, 1)/conv_v*1e5, v_Saturn(end, 2)/conv_v*...
       1e5, v_Saturn(end, 3)/conv_v*1e5, 'g', 'LineWidth', 1.2)
hold on
plotSaturn (R_Saturn)
set(gca,'TickLabelInterpreter','Latex')
grid on
axis equal
axis([ -1.3e6 1.5e6 -1.12e6 1.1e6 -0.6e6 2.2e6 ])
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Rendezvous with Titan')
legend('S/C', 'Titan', 'V Saturn', 'Location', 'NorthEast')
hold off

%Titan atmosphere density 
h = linspace(0, 1300, 2e2);
rho = zeros(length(h), 1);
for j = 1:length(h)
    rho(j) = YelleModel (h(j));
end

figure
semilogx(rho, h, 'b-', 'LineWidth', 1.2)
set(gca,'TickLabelInterpreter','Latex')
grid on
xlabel('$\rho$ [kg/m$^3$]')
ylabel('Height [km]')
title('Yelle Model: Titan Average Atmospheric Density Profile')

%aerocapture sequence
%phase 1 trajectory
figure 
plot3(r4_sc2T(:, 1)*au, r4_sc2T(:, 2)*au, r4_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
quiver3(0, 0, 0, v4_Titan(ind_T, 1)/conv_v*3e3, v4_Titan(ind_T, 2)/...
       conv_v*3e3, v4_Titan(ind_T, 3)/conv_v*3e3, 'color', [ 1 ...
       0.64453125 0 ], 'LineWidth', 1.2)
hold on
Titan = imread('Titan.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Titan;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Titan, R_Titan, ...
          R_Titan, 50);
h_T1 = surface(-X,-Y,-Z,props);
set(gca,'TickLabelInterpreter','Latex')
view(3)
grid on
axis equal
axis([ -1.5e4 SOI_T*au/5 -SOI_T*au/4 SOI_T*au/4 -SOI_T*au/5 SOI_T*au/4 ])
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Phase 1')
legend('S/C', 'V Titan', 'Location', 'Best')
hold off
%phase 2 trajectory
figure 
plot3(r4_sc2T(:, 1)*au, r4_sc2T(:, 2)*au, r4_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
quiver3(0, 0, 0, v5_Titan(end, 1)/conv_v*3e3, v5_Titan(end, 2)/...
       conv_v*3e3, v5_Titan(end, 3)/conv_v*3e3, 'color', [ 1 ...
       0.64453125 0 ], 'LineWidth', 1.2)
hold on
quiver3(r5_sc2T(end, 1)*au, r5_sc2T(end, 2)*au, r5_sc2T(end, 3)*au,...
       DV1_H_pf(1)/conv_v*1e5, DV1_H_pf(2)/conv_v*1e5, DV1_H_pf(3)/...
       conv_v*1e5, 'k', 'LineWidth', 1.2)
hold on
plot3(r5_sc2T(:, 1)*au, r5_sc2T(:, 2)*au, r5_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
Titan = imread('Titan.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Titan;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Titan, R_Titan,...
          R_Titan, 50);
h_T2 = surface(-X,-Y,-Z,props);
set(gca,'TickLabelInterpreter','Latex')
view(3)
grid on
axis equal
axis([ -1.5e4 SOI_T*au/5 -SOI_T*au/4 SOI_T*au/4 -SOI_T*au/5 SOI_T*au/4 ])
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Phase 2')
legend('S/C', 'V Titan', '\DeltaV_1 Hohmann', 'Location', 'Best')
hold off
%phase 3 trajectory
figure 
plot3(r4_sc2T(:, 1)*au, r4_sc2T(:, 2)*au, r4_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
quiver3(0, 0, 0, v6_Titan(end, 1)/conv_v*3e3, v6_Titan(end, 2)/...
       conv_v*3e3, v6_Titan(end, 3)/conv_v*3e3, 'color', [ 1 ...
       0.64453125 0 ], 'LineWidth', 1.2)
hold on
quiver3(r6_sc2T(end, 1)*au, r6_sc2T(end, 2)*au, r6_sc2T(end, 3)*au, ...
       DV2_H_pf(1)/conv_v*4e4, DV2_H_pf(2)/conv_v*4e4, DV2_H_pf(3)/...
       conv_v*4e4, 'color', 'k', 'LineWidth', 1.2)
hold on
plot3(r5_sc2T(:, 1)*au, r5_sc2T(:, 2)*au, r5_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
plot3(r6_sc2T(:, 1)*au, r6_sc2T(:, 2)*au, r6_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
Titan = imread('Titan.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Titan;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Titan, R_Titan, ...
          R_Titan, 50);
h_T3 = surface(-X,-Y,-Z,props);
set(gca,'TickLabelInterpreter','Latex')
view(3)
grid on
axis equal
axis([ -1.5e4 SOI_T*au/5 -SOI_T*au/4 SOI_T*au/4 -SOI_T*au/5 SOI_T*au/4 ])
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Phase 3')
legend('S/C', 'V Titan', '\DeltaV_2 Hohmann', 'Location', 'Best')
hold off
%phase 4 trajectory
figure 
plot3(r4_sc2T(:, 1)*au, r4_sc2T(:, 2)*au, r4_sc2T(:, 3)*au, 'r-', ....
     'LineWidth', 1.2)
hold on
quiver3(0, 0, 0, v7_Titan(end, 1)/conv_v*3e3, v7_Titan(end, 2)/...
       conv_v*3e3, v7_Titan(end, 3)/conv_v*3e3, 'color', [ 1 ...
       0.64453125 0 ], 'LineWidth', 1.2)
hold on
plot3(r5_sc2T(:, 1)*au, r5_sc2T(:, 2)*au, r5_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
plot3(r6_sc2T(:, 1)*au, r6_sc2T(:, 2)*au, r6_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
plot3(r7_sc2T(:, 1)*au, r7_sc2T(:, 2)*au, r7_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
Titan = imread('Titan.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Titan;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Titan, R_Titan,...
          R_Titan, 50);
h_T4 = surface(-X,-Y,-Z,props);
set(gca,'TickLabelInterpreter','Latex')
view(3)
grid on
axis equal
axis([ -1.5e4 SOI_T*au/5 -SOI_T*au/4 SOI_T*au/4 -SOI_T*au/5 SOI_T*au/4 ])
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Phase 4')
legend('S/C', 'V Titan', 'Location', 'Best')
hold off
%total trajectory
figure 
plot3(r4_sc2T(:, 1)*au, r4_sc2T(:, 2)*au, r4_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
quiver3(0, 0, 0, v7_Titan(end, 1)/conv_v*3e3, v7_Titan(end, 2)/...
       conv_v*3e3, v7_Titan(end, 3)/conv_v*3e3, 'color', [ 1 ...
       0.64453125 0 ], 'LineWidth', 1.2)
hold on
quiver3(r5_sc2T(end, 1)*au, r5_sc2T(end, 2)*au, r5_sc2T(end, 3)*au, ...
       DV1_H_pf(1)/conv_v*1e5, DV1_H_pf(2)/conv_v*1e5, DV1_H_pf(3)/...
       conv_v*1e5, 'color', 'k', 'LineWidth', 1.2)
hold on
quiver3(r6_sc2T(end, 1)*au, r6_sc2T(end, 2)*au, r6_sc2T(end, 3)*au,...
       DV2_H_pf(1)/conv_v*4e4, DV2_H_pf(2)/conv_v*4e4, DV2_H_pf(3)/...
       conv_v*4e4, 'color', 'k', 'LineWidth', 1.2)
hold on
plot3(r5_sc2T(:, 1)*au, r5_sc2T(:, 2)*au, r5_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
plot3(r6_sc2T(:, 1)*au, r6_sc2T(:, 2)*au, r6_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
plot3(r7_sc2T(:, 1)*au, r7_sc2T(:, 2)*au, r7_sc2T(:, 3)*au, 'r-', ...
     'LineWidth', 1.2)
hold on
Titan = imread('Titan.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Titan;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Titan, R_Titan, ...
          R_Titan, 50);
h_T5 = surface(-X,-Y,-Z,props);
set(gca,'TickLabelInterpreter','Latex')
view(3)
grid on
axis equal
axis([ -1.5e4 SOI_T*au/5 -SOI_T*au/4 SOI_T*au/4 -SOI_T*au/5 SOI_T*au/4 ])
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Aerocapture')
legend('S/C', 'V Titan', '\DeltaV Hohmann', 'Location', 'Best')
hold off

%aerobraking effect on spacecraft velocity
figure
plot(vel_sc2T/conv_v, dist_sc2T*au - R_Titan, 'b-', 'LineWidth', 1.2)
hold on
plot(linspace(0, 12, 10), 1300*ones(10, 1), 'k--', 'LineWidth', 1.2)
set(gca,'TickLabelInterpreter','Latex')
grid on
axis([ 0 8 0 1450 ])
xlabel('Velocity [km/s]')
ylabel('Height [km]')
title('Aerobraking')
legend('S/C Velocity', 'Entry Height', 'Location', 'Best')
hold off

%operation orbit evolution
figure 
plot3(r_op_sc2T(:, 1)*au, r_op_sc2T(:, 2)*au, r_op_sc2T(:, 3)*au, 'r-',...
     'LineWidth', 1.2)
hold on
Titan = imread('Titan.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Titan;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Titan, R_Titan,...
          R_Titan, 50);
h_T6 = surface(-X,-Y,-Z,props);
set(gca,'TickLabelInterpreter','Latex')
view(3)
grid on
axis equal
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Operation Orbit (6 months)')
legend('S/C', 'Location', 'Best')
hold off

%capture phase spacecraft ground track
figure
Titan = imread('Titan.jpg','jpg');
image([-180 180], [90 -90], Titan)
hold on
geoshow(lat_capture, long_capture, 'Color', 'r', 'LineWidth', 1.5)
set(gca, 'ydir', 'normal', 'TickLabelInterpreter', 'Latex' )
title('Capture Phase Ground Track');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
hold off

%operational phase ground track (one month)
figure
Titan = imread('Titan.jpg','jpg');
image([-180 180], [90 -90], Titan)
hold on
geoshow(lat_op(1:ind_week), long_op(1:ind_week), 'Color', 'r', ...
        'LineWidth', 1.5)
set(gca, 'ydir', 'normal', 'TickLabelInterpreter', 'Latex' )
title('Operational Phase Ground Track (1 week)');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
hold off

%operational phase semi-major axis evolution
figure
plot(t_op/(2*pi)*12, a_op*au, 'b-', 'LineWidth', 1.2)
set(gca,'TickLabelInterpreter','Latex')
grid on
xlabel('t [months]')
ylabel('Semi-major Axis [km]')
title('Semi-major Axis Evolution')

%operational phase eccentricity evolution
figure
plot(t_op/(2*pi)*12, e_op, 'b-', 'LineWidth', 1.2)
set(gca,'TickLabelInterpreter','Latex')
grid on
xlabel('t [months]')
ylabel('Eccentricity')
title('Eccentricity Evolution')

%operational phase inclination evolution
figure
plot(t_op/(2*pi)*12, rad2deg(i_op), 'b-', 'LineWidth', 1.2)
set(gca,'TickLabelInterpreter','Latex')
grid on
xlabel('t [months]')
ylabel('Inclination [deg]')
title('Inclination Evolution')

%operational phase right ascension of ascending node evolution
figure
plot(t_op/(2*pi)*12, rad2deg(RAAN_op), 'b-', 'LineWidth', 1.2)
set(gca,'TickLabelInterpreter','Latex')
grid on
xlabel('t [months]')
ylabel('$\Omega$ [deg]')
title('Right Ascension of Ascending Node Evolution')

%operational phase pericenter argument evolution
figure
plot(t_op/(2*pi)*12, rad2deg(PA_op), 'b-', 'LineWidth', 1.2)
set(gca,'TickLabelInterpreter','Latex')
grid on
xlabel('t [months]')
ylabel('$\omega$ [deg]')
title('Pericenter Argument Evolution')

%% Auxiliary Functions

function [ ROT ] = ROT_p2HIF (RAAN, PA, inc)
%define the rotation matrix from perifocal to heliocentric inertial frame

ROT = [ cos(RAAN)*cos(PA)-sin(RAAN)*sin(PA)*cos(inc) ...
        sin(RAAN)*cos(PA)+cos(RAAN)*sin(PA)*cos(inc) ...
        sin(PA)*sin(inc); ...
        -cos(RAAN)*sin(PA)-sin(RAAN)*cos(PA)*cos(inc) ...
        -sin(RAAN)*sin(PA)+cos(RAAN)*cos(PA)*cos(inc) ...
        cos(PA)*sin(inc); ...
        sin(RAAN)*sin(inc) ...
        -cos(RAAN)*sin(inc) ...
        cos(inc) ];                       %rotation matrix perifocal-HIF
            
end

function [ J ] = cost_transfer (X_try, X0, mu_Sun, mu_Earth, ...
                 mu_Jupiter, a_ell, c1, c2, c3, c4, options2)
%define the performance index to minimize through the optimization process

dt_wait = X_try(1);
dt1 = X_try(2);
DV_v = X_try(3);
DV_n = X_try(4);
dt2 = X_try(5);
dt3 = X_try(6);
perigee = X_try(7);
perigiovio = X_try(8);

%propagate the planets' orbits till the mission starting epoch
[ ~, dX_planets ] = ode113 (@Prop_Orbits, [ 0 dt_wait ], X0, options2,...
                    mu_Sun);

r0_Earth = dX_planets(end, 1:3);          %Earth position at the mission 
%starting epoch [au]
v0_Earth = dX_planets(end, 4:6);          %Earth velocity at the mission
%starting epoch [EOS]

r0_Jupiter = dX_planets(end, 7:9);        %Jupiter position at the mission
%starting epoch [au]
v0_Jupiter = dX_planets(end, 10:12);      %Jupiter velocity at the mission
%starting epoch [EOS]

r0_Saturn = dX_planets(end, 13:15);       %Saturn position at the mission
%starting epoch [au]
v0_Saturn = dX_planets(end, 16:18);       %Saturn velocity at the mission
%starting epoch [EOS]

X_start_planets = [ r0_Earth'; v0_Earth'; r0_Jupiter'; v0_Jupiter'; ...
                  r0_Saturn'; v0_Saturn' ];

%spacecraft coincident with the Earth C.M. (patched-conics method)
r_HIF_sc = [ r0_Earth(1);
             r0_Earth(2);
             r0_Earth(3) ];               %spacecraft position in HIF [au]
v_HIF_sc = [ v0_Earth(1);
             v0_Earth(2);
             v0_Earth(3) ];               %spacecraft velocity in HIF [EOS]

v_ell = sqrt(2*(mu_Sun/norm(r_HIF_sc) ...
        - mu_Sun/(2*a_ell)));             %required velocity to enter the
%elliptic orbit [EOS]
v_inf = v_ell - norm(v_HIF_sc);           %required hyperbolic excess
%velocity [EOS]

v_inf_VNB = [ v_inf; 0; 0 ];              %hyperbolic excess in the VNB
%reference frame [EOS]
%VNB reference frame is centered in the spacecraft center of mass, with:
%x-axis along the spacecraft velocity (V)
%y-axis along the spacecraft orbit normal, i.e. the angular momentum (N)
%z-axis along the spacecraft orbit binormal (B)

%define the hyperbolic excess velocity in HIF
r_vers = r_HIF_sc/norm(r_HIF_sc);         %spacecraft position versor
v_vers = v_HIF_sc/norm(v_HIF_sc);         %spacecraft velocity versor
n_vers = cross(r_vers, v_vers);           %normal versor
b_vers = cross(v_vers, n_vers);           %binormal versor

ROT_VNB2HIF = [ v_vers, n_vers, b_vers ]; %rotation matrix VNB-HIF

v_inf_HIF = ROT_VNB2HIF*v_inf_VNB;        %hyperbolic excess in the HIF
%[EOS]

%update spacecraft starting velocity
v_HIF_sc = v_HIF_sc + v_inf_HIF;          %new spacecraft starting velocity
%[EOS]

%full mission starting conditions (Earth + spacecraft)
X_start1 = [ X_start_planets(1:6); r_HIF_sc; v_HIF_sc ];

%1st trajectory arc (till the DSM, "deep space maneuver")
[ ~, dX1 ] = ode113 (@Prop_Arc1, [ 0 dt1 ], X_start1, options2, mu_Sun);

r1_Earth = dX1(end, 1:3);                 %Earth position at the DSM [au]
v1_Earth = dX1(end, 4:6);                 %Earth velocity at the DSM [EOS]

r1_sc = dX1(end, 7:9);                    %spacecraft position at the DSM 
%[au]
v1_sc = dX1(end, 10:12);                  %spacecraft velocity at the DSM
%[EOS]

%Deep Space Maneuver (DSM)
DV_VNB = [ DV_v; DV_n; 0 ];               %DSM's deltaV in the VNB 
%reference frame [EOS]

%define the DSM deltaV in HIF
r_vers = r1_sc'/norm(r1_sc);              %spacecraft position versor
v_vers = v1_sc'/norm(v1_sc);              %spacecraft velocity versor
n_vers = cross(r_vers, v_vers);           %normal versor
b_vers = cross(v_vers, n_vers);           %binormal versor

ROT_VNB2HIF = [ v_vers, n_vers, b_vers ]; %rotation matrix VNB-HIF

DV_HIF = ROT_VNB2HIF*DV_VNB;              %DSM's deltaV in the HIF [EOS]

%update spacecraft velocity after DSM
v1_sc = v1_sc' + DV_HIF;                  %spacecraft velocity after the 
%maneuver [EOS]

%2nd trajectory arc starting conditions (Earth + spacecraft)
X_start2 = [ r1_Earth'; v1_Earth'; r1_sc'; v1_sc ];

%2nd trajectory arc (till after the Earth Gravity Assist)
[ t2, dX2 ] = ode113 (@Prop_Perturbed, [ 0 dt2 ], X_start2, options2,...
              mu_Sun, mu_Earth);

r2_Earth = dX2(:, 1:3);                   %Earth positions [au]

r2_sc = dX2(:, 7:9);                      %spacecraft positions [au]
v2_sc = dX2(:, 10:12);                    %spacecraft velocities [EOS]

%preallocation
r_sc2E = zeros(length(t2), 3);
dist_sc2E = zeros(length(t2), 1);

for j = 1:length(t2)
    r_sc2E(j, :) = r2_sc(j, :) -...
                   r2_Earth(j, :);        %spacecraft positions w.r.t.
    %Earth [au]
    dist_sc2E(j) = norm(r_sc2E(j, :));    %spacecraft distance from Earth
    %C.M. [au]
end

rp_E = min(dist_sc2E);                    %perigee radius during Earth 
%gravity assist [au]

%condition to include the Earth flyby inside the 2nd trajectory arc, and to
%reach a certain distance such that the Earth's gravitational perturbation
%can be neglected within the next arc
if norm(r2_sc(end, :)) < 2 || norm(r2_sc(end, :)) > 4
    J = inf;
    return
end

%propagate Jupiter orbit before the 3rd trajectory arc
[ ~, dX_Jupiter ] = ode113 (@Prop_Keplerian, [ 0 dt1+dt2 ],...
                    X_start_planets(7:12), options2, mu_Sun);

r_Jupiter = dX_Jupiter(end, 1:3);
v_Jupiter = dX_Jupiter(end, 4:6);

%3rd trajectory arc starting conditions (Jupiter + spacecraft)
X_start3 = [ r_Jupiter'; v_Jupiter'; r2_sc(end, :)'; v2_sc(end, :)' ];

%3rd trajectory arc (till after the Jupiter Gravity Assist)
[ t3, dX3 ] = ode113 (@Prop_Perturbed, [ 0 dt3 ], X_start3, options2,...
              mu_Sun, mu_Jupiter);

r3_Jupiter = dX3(:, 1:3);                 %Jupiter positions [au]

r3_sc = dX3(:, 7:9);                      %spacecraft positions [au]
v3_sc = dX3(:, 10:12);                    %spacecraft velocities [EOS]

%preallocation
r_sc2J = zeros(length(t3), 3);
dist_sc2J = zeros(length(t3), 1);

for j = 1:length(t3)
    r_sc2J(j, :) = r3_sc(j, :) -...
                   r3_Jupiter(j, :);      %spacecraft positions w.r.t. 
    %Jupiter [au]
    dist_sc2J(j) = norm(r_sc2J(j, :));    %spacecraft distance from Jupiter
    %C.M. [au]
end

rp_J = min(dist_sc2J);                    %perigiovio radius during Jupiter
%gravity assist [au]

%propagate Saturn orbit to set target position and velocity
[ ~, dX_Saturn ] = ode113 (@Prop_Keplerian, [ 0 dt1+dt2+dt3 ], ...
                   X_start_planets(13:18), options2, mu_Sun);

r_target = dX_Saturn(end, 1:3);           %Saturn target position [au]
v_target = dX_Saturn(end, 4:6);           %Saturn target velocity [EOS]

%determine cost function value
J = (dt_wait + dt1 + dt2 + dt3) + c1*(norm(v_inf_HIF) + norm(DV_HIF)) +...
    c2*abs(norm(r3_sc(end, :) - r_target)) + c3*abs(norm(v3_sc(end, :) ...
    - v_target)) + c4*abs(rp_E - perigee) + c2* abs(rp_J - perigiovio);

end

function [ dX ] = Prop_Orbits (~, X, mu_Sun)
%propagate the planets' orbits considering pure Keplerian motion, i.e. 
%taking into account only the Sun gravitational attraction

%Input

r_Earth = [ X(1); X(2); X(3) ];
v_Earth = [ X(4); X(5); X(6) ];

r_Jupiter = [ X(7); X(8); X(9) ];
v_Jupiter = [ X(10); X(11); X(12) ];

r_Saturn = [ X(13); X(14); X(15) ];
v_Saturn = [ X(16); X(17); X(18) ];

%Motion Equations

dr_Earth = v_Earth;
dv_Earth = - mu_Sun*r_Earth/(norm(r_Earth))^3;

dr_Jupiter = v_Jupiter;
dv_Jupiter = - mu_Sun*r_Jupiter/(norm(r_Jupiter))^3;

dr_Saturn = v_Saturn;
dv_Saturn = - mu_Sun*r_Saturn/(norm(r_Saturn))^3;

%Output

dX = [ dr_Earth; dv_Earth; dr_Jupiter; dv_Jupiter; dr_Saturn; dv_Saturn ];

end

function [ dX ] = Prop_Arc1 (~, X, mu_Sun)
%propagate considering pure Keplerian motion both the Earth and the
%spacecraft, which is not perturbed by the planet's gravitational
%attraction

%Input

r_Earth = [ X(1); X(2); X(3) ];
v_Earth = [ X(4); X(5); X(6) ];

r_sc = [ X(7); X(8); X(9) ];
v_sc = [ X(10); X(11); X(12) ];

%Motion Equations

dr_Earth = v_Earth;
dv_Earth = - mu_Sun*r_Earth/(norm(r_Earth))^3;

dr_sc = v_sc;
dv_sc = - mu_Sun*r_sc/(norm(r_sc))^3;

%Output

dX = [ dr_Earth; dv_Earth; dr_sc; dv_sc ];

end

function [ dX ] = Prop_Perturbed (~, X, mu1, mu2)
%consider one planet gravitational perturbation acting on the body motion
%body mass is neglected

%Input

r_planet = [ X(1); X(2); X(3) ];
v_planet = [ X(4); X(5); X(6) ];

r_b = [ X(7); X(8); X(9) ];
v_b = [ X(10); X(11); X(12) ];

%Motion Equations

dr_planet = v_planet;
dv_planet = - mu1*r_planet/(norm(r_planet))^3;

dr_b = v_b;
dv_b = - mu1*r_b/(norm(r_b))^3 - mu2*(r_b - r_planet)/(norm(r_b -...
       r_planet))^3;

%Output

dX = [ dr_planet; dv_planet; dr_b; dv_b ];

end

function [ dX ] = Prop_Keplerian (~, X, mu)
%propagate the orbit of a body considering pure Keplerian motion

%Input

r = [ X(1); X(2); X(3) ];
v = [ X(4); X(5); X(6) ];

%Motion Equations

dr = v;
dv = - mu*r/(norm(r))^3;

%Output

dX = [ dr; dv ];

end

function [ Mi_sc ] = Mass_sc (mu_Earth, T_GEO, r_cpo, DV)
%estimate the spacecraft initial mass launched by Delta IV Heavy

Mi_GTO = 14.22e3;                         %payload mass in GTO [kg]

a_GEO = (mu_Earth*(T_GEO/(2*pi))^2)^(1/3);%geostationary orbit semi-major 
%axis [au]
r_GEO = a_GEO;                            %geostationary orbit radius [au]

a_GTO = (r_GEO + r_cpo)/2;                %GTO semi-major axis [au]
vp_GTO = sqrt(2*mu_Earth*(1/r_cpo -...
         1/(2*a_GTO)));                   %GTO perigee velocity [EOS]

v_cpo = sqrt(mu_Earth/r_cpo);             %circular parking orbit velocity
%[EOS]
DV_GTO = vp_GTO - v_cpo;                  %deltaV for GTO insertion [EOS]

Mi_sc = Mi_GTO*DV_GTO/DV;                 %spacecraft initial mass estimate
%[kg]

end

function plotSaturn (R_Saturn)
%plot Saturn 3D Image 

Saturn = imread('Saturn.jpg','jpg');
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Saturn;
center = [0; 0; 0];
[X,Y,Z] = ellipsoid(center(1), center(2), center(3), R_Saturn, R_Saturn,...
          R_Saturn, 50);
surface(-X,-Y,-Z,props);
view(3)
axis equal

%Ring D

M = 1e2;
N = 1e2;
R1 = 66900;                               %inner radius 
R2 = 74510;                               %outer radius
nR = linspace(R1,R2,M);
nT = linspace(0,2*pi,N);
[R,T] = meshgrid(nR,nT) ;
%convert grid to cartesian coordinates
X = R.*cos(T); 
Y = R.*sin(T);
[m,n] = size(X);
%plot grid
hold on
axis equal
%plot internal grid lines
for i=1:m
    plot(X(i,:),Y(i,:),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end
for j=1:n
    plot(X(:,j),Y(:,j),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end

%Ring C

M = 1e2;
N = 1e2;
R1 = 74658;                               %inner radius 
R2 = 91975;                               %outer radius
nR = linspace(R1,R2,M);
nT = linspace(0,2*pi,N);
[R,T] = meshgrid(nR,nT) ;
%convert grid to cartesian coordinates
X = R.*cos(T); 
Y = R.*sin(T);
[m,n] = size(X);
%plot grid
hold on
axis equal
%plot internal grid lines
for i=1:m
    plot(X(i,:),Y(i,:),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end
for j=1:n
    plot(X(:,j),Y(:,j),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end

%Ring B

M = 1e2;
N = 1e2;
R1 = 91975;                               %inner radius 
R2 = 117507;                              %outer radius
nR = linspace(R1,R2,M);
nT = linspace(0,2*pi,N);
[R,T] = meshgrid(nR,nT) ;
%convert grid to cartesian coordinates
X = R.*cos(T); 
Y = R.*sin(T);
[m,n] = size(X);
%plot grid
hold on
axis equal
%plot internal grid lines
for i=1:m
    plot(X(i,:),Y(i,:),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end
for j=1:n
    plot(X(:,j),Y(:,j),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end

%Ring A

M = 1e2;
N = 1e2;
R1 = 122340;                              %inner radius 
R2 = 136780;                              %outer radius
nR = linspace(R1,R2,M);
nT = linspace(0,2*pi,N);
[R,T] = meshgrid(nR,nT) ;
%convert grid to cartesian coordinates
X = R.*cos(T); 
Y = R.*sin(T);
[m,n] = size(X);
% Plot grid
hold on
axis equal
%plot internal grid lines
for i=1:m
    plot(X(i,:),Y(i,:),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end
for j=1:n
    plot(X(:,j),Y(:,j),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end

%Ring F

M = 1e2;
N = 1e2;
R1 = 139826;                              %inner radius 
R2 = R1;                                  %outer radius
nR = linspace(R1,R2,M);
nT = linspace(0,2*pi,N);
[R,T] = meshgrid(nR,nT) ;
%convert grid to cartesian coordinates
X = R.*cos(T); 
Y = R.*sin(T);
[m,n] = size(X);
%plot grid
hold on
axis equal
%plot internal grid lines
for i=1:m
    plot(X(i,:),Y(i,:),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end
for j=1:n
    plot(X(:,j),Y(:,j),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end

%Ring G

M = 1e2;
N = 1e2;
R1 = 166000;                              %inner radius 
R2 = 173000;                              %outer radius
nR = linspace(R1,R2,M);
nT = linspace(0,2*pi,N);
[R,T] = meshgrid(nR,nT) ;
%convert grid to cartesian coordinates
X = R.*cos(T); 
Y = R.*sin(T);
[m,n] = size(X);
%plot grid
hold on
axis equal
%plot internal grid lines
for i=1:m
    plot(X(i,:),Y(i,:),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end
for j=1:n
    plot(X(:,j),Y(:,j),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
end

%Ring E (don't plot because it's very thin and faint, but we consider it
%for the trajectory in order to not cross it)

% M = 1e3;
% N = 1e2;
% R1 = 180000;                            %inner radius 
% R2 = 480000;                            %outer radius
% nR = linspace(R1,R2,M);
% nT = linspace(0,2*pi,N);
% [R,T] = meshgrid(nR,nT) ;
% %convert grid to cartesian coordinates
% X = R.*cos(T); 
% Y = R.*sin(T);
% [m,n] = size(X);
% %plot grid
% hold on
% axis equal
% %plot internal grid lines
% for i=1:m
%     plot(X(i,:),Y(i,:),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
% end
% for j=1:n
%     plot(X(:,j),Y(:,j),'color',[0.9375 0.8984 0.5469],'linewidth',1.2); 
% end

end

function [ rho ] = YelleModel (h)
%interpolate data of the Yelle model for Titan atmosphere (average)

%height over Titan surface in km
height = [ 0, 10, 30, 50, 70, 90, 110, 130, 154, 174, 194, 235, 285,...
           335, 385, 435, 490, 540, 590, 640, 690, 740, 790, 840, 900, ...
           950, 1000, 1050, 1100, 1150, 1200, 1250, 1300 ];
%atmospheric density in kg/m^3
density = [ 5.44, 3.56, 1.23, 0.358, 0.827e-1, 0.317e-1, 0.160e-1, ...
            0.875e-2, 0.456e-2, 0.274e-2, 0.168e-2, 0.650e-3, 0.222e-3,...
            0.812e-4, 0.310e-4, 0.119e-4, 0.399e-5, 0.138e-5, 0.450e-6,...
            0.149e-6, 0.527e-7, 0.203e-7, 0.845e-8, 0.377e-8, 0.154e-8,...
            0.767e-9, 0.390e-9, 0.202e-9, 0.107e-9, 0.576e-10, ...
            0.316e-10, 0.177e-10, 0.100e-10 ];
        
rho = spline(height, density, h);         %atmospheric density at the
%desired height [kg/m^3]  

end

function [ J ] = cost_aerocapture (X_try, v0_sc, X0_Titan, SOI_S,...
                 mu_Saturn, mu_Titan, R_Titan, au, Cd, S, M, i_Titan, ...
                 c1, c2, c3, options)

phi = X_try(1);                           %latitude [rad]
lambda = X_try(2);                        %longitude [rad]
dt = X_try(3);                            %time of flight [year*(2*pi)]

%starting spacecraft position on the Saturn SOI
x0_sc = SOI_S*cos(phi)*cos(lambda);
y0_sc = SOI_S*cos(phi)*sin(lambda);
z0_sc = SOI_S*sin(phi);

X_start = [ X0_Titan; x0_sc; y0_sc; z0_sc; v0_sc' ];

[ t, dX_capture ] = ode113 (@Prop_Capture, [ 0 dt ], X_start, options,...
                    mu_Saturn, mu_Titan, R_Titan, au, Cd, S, M);
                
r_Titan = dX_capture(:, 1:3);             %Titan positions [au]
v_Titan = dX_capture(:, 4:6);             %Titan velocities [EOS]

r_sc = dX_capture(:, 7:9);                %spacecraft positions [au]
v_sc = dX_capture(:, 10:12);              %spacecraft velocities [EOS]

%preallocation
r_sc2T = zeros(length(t), 3);
v_sc2T = zeros(length(t), 3);
dist_sc2T = zeros(length(t), 1);
dist_sc2S = zeros(length(t), 1);

for j = 1:length(t)
    r_sc2T(j, :) = r_sc(j, :) -...
                   r_Titan(j, :);         %spacecraft positions w.r.t. 
    %Titan [au]
    v_sc2T(j, :) = v_sc(j, :) -...
        v_Titan(j, :);                    %spacecraft velocities w.r.t.
    %Titan [EOS]
    dist_sc2T(j) = norm(r_sc2T(j, :));    %spacecraft distance from Titan
    %C.M. [au]
    dist_sc2S(j) = norm(r_sc(j, :));      %spacecraft distance from Saturn
    %C.M. [au]
end

rp_T = min(dist_sc2T);                    %pericenter radius w.r.t. Titan 
%[au]
rp_S = min(dist_sc2S);                    %pericenter radius w.r.t. Saturn
%[au]

rp_target = (250 + R_Titan)/au;           %target pericenter radius [au]
r_out = (1310 + R_Titan)/au;              %target exit from Titan 
%atmosphere [au]

[ ~, e, i, ~, ~, ~ ] = rv2coe (r_sc2T(end, :)', v_sc2T(end, :)',...
                       mu_Titan);         %classical orbital parameters 
%after the aerocapture

%spacecraft must not go through Saturn's rings (not even the E one)
if rp_S <= 480000/au
    J = inf;
    return
end

%pericenter must be included inside the aerocapture trajectory arc
if rp_T == norm(r_sc2T(end, :))
    J = inf;
    return
end

%control to accelerate the optimization avoiding the calculation of
%unacceptable attempts
if any(dist_sc2T < (R_Titan + 100)/au)
    J = inf;
    return
end

J = c1*abs(abs(i) - pi/2 + i_Titan) + c2*abs(rp_T - rp_target) + ...
    c2*abs(norm(r_sc2T(end, :)) - r_out) + c3*...
    abs(e - 0.5);                         %cost function value

end

function [ dX ] = Prop_Capture (~, X, mu_Saturn, mu_Titan, R_Titan, au,...
                  Cd, S, M)
%consider one planet gravitational perturbation acting on the body motion
%body mass is used in the drag equation
%the planet is considered non-rotating, therefore the velocity relative to 
%the atmosphere is equal to the spacecraft one w.r.t. the planet itself

%Input

r_Titan = [ X(1); X(2); X(3) ];
v_Titan = [ X(4); X(5); X(6) ];

r_sc = [ X(7); X(8); X(9) ];
v_sc = [ X(10); X(11); X(12) ];

r_sc2T = r_sc - r_Titan;

%Motion Equations

dr_Titan = v_Titan;
dv_Titan = - mu_Saturn*r_Titan/(norm(r_Titan))^3;

dr_sc = v_sc;
if norm(r_sc2T) <= (R_Titan + 1300)/au
    rho = YelleModel (norm(r_sc2T)*au - R_Titan);
    dv_sc = - mu_Saturn*r_sc/(norm(r_sc))^3 - mu_Titan*(r_sc - ...
            r_Titan)/(norm(r_sc - r_Titan))^3 - 0.5*rho*Cd*S/M*1e3*au*...
            norm(v_sc - v_Titan)*(v_sc - v_Titan);
else
    dv_sc = - mu_Saturn*r_sc/(norm(r_sc))^3 - mu_Titan*(r_sc -...
            r_Titan)/(norm(r_sc - r_Titan))^3;
end

%Output

dX = [ dr_Titan; dv_Titan; dr_sc; dv_sc ];

end

function [ a, e, i, RAAN, PA, TA ] = rv2coe (rVect, vVect ,mu)
%classical orbital elements from position and velocity vectors

iVers = [ 1 0 0 ]; 
jVers = [ 0 1 0 ];
kVers = [ 0 0 1 ];

r = norm(rVect);                          %position [au]
v = norm(vVect);                          %speed [EOS]
E = (v^2)/2 - mu/r;                       %energy 
a = - mu/(2*E);                           %semi-major axis [au]
hVect = cross(rVect, vVect);              %angular momentum 
h = norm(hVect);                          %angular momentum module
eccVect = (cross(vVect, hVect))/mu -...
          rVect/r;                        %eccentricity vector
e = norm(eccVect);                        %eccentricity
i = acos(hVect(3)/h);                     %inclination [rad]
nVers = (cross(kVers, hVect))/(norm(cross...
        (kVers, hVect)));                 %nodes line versor
%right ascension of ascending node [rad]
if dot(nVers, jVers) >= 0
    RAAN = acos(dot(iVers, nVers));                  
else
    RAAN = 2*pi - acos(dot(iVers, nVers));
end
eccVers = eccVect/e;                      %eccentricity versor
%pericenter argument [rad]
if dot(eccVect, kVers) >= 0
    PA = acos(dot(nVers, eccVers));
else
    PA = 2*pi - acos(dot(nVers, eccVers));
end
rVers = rVect/r;                          %position versor
%true anomaly [rad]
if dot(rVect, vVect) >= 0
    TA = acos(dot(eccVers, rVers));
else
    TA = 2*pi - acos(dot(eccVers, rVers));
end
    
end

function [ dX ] = Prop_OP (~, X, mu_Saturn, mu_Titan, R_Titan, au, Cd,...
                  S, M, J2, J3, J4)
%consider Saturn gravitational perturbation acting on the spacecraft motion
%spacecraft mass is used in the drag equation
%Titan is considered non-rotating, therefore the velocity relative to 
%the atmosphere is equal to the spacecraft one w.r.t. Titan itself
%for the final operation orbit evolution also Titan J2, J3 and J4 are
%taking into account (drag + geopotential + third-body)

%Input

x_Titan = X(1);
y_Titan = X(2);
z_Titan = X(3);

Vx_Titan = X(4);
Vy_Titan = X(5);
Vz_Titan = X(6);

r_Titan = [ x_Titan; y_Titan; z_Titan ];
v_Titan = [ Vx_Titan; Vy_Titan; Vz_Titan ];

x_sc = X(7);
y_sc = X(8);
z_sc = X(9);

Vx_sc = X(10);
Vy_sc = X(11);
Vz_sc = X(12);

r_sc = [ x_sc; y_sc; z_sc ];
v_sc = [ Vx_sc; Vy_sc; Vz_sc ];

r_sc2T = r_sc - r_Titan;

%Motion Equations

dx_Titan = Vx_Titan;
dy_Titan = Vy_Titan;
dz_Titan = Vz_Titan;

dVx_Titan = - mu_Saturn*x_Titan/(norm(r_Titan))^3;
dVy_Titan = - mu_Saturn*y_Titan/(norm(r_Titan))^3;
dVz_Titan = - mu_Saturn*z_Titan/(norm(r_Titan))^3;

dx_sc = Vx_sc;
dy_sc = Vy_sc;
dz_sc = Vz_sc;
if norm(r_sc2T) <= (R_Titan + 1300)/au
    rho = YelleModel (norm(r_sc2T)*au - R_Titan);
    dVx_sc = - mu_Saturn*x_sc/(norm(r_sc))^3 ... 
             - 0.5*rho*Cd*S/M*1e3*au*norm(v_sc - v_Titan)*(Vx_sc -...
               Vx_Titan) ...
             - mu_Titan*(x_sc - x_Titan)/(norm(r_sc - r_Titan))^3 ...
             + 3*mu_Titan*J2*(R_Titan/au)^2/(2*(norm(r_sc - r_Titan))^...
               5)*(5*(z_sc - z_Titan)^2/(norm(r_sc - r_Titan))^2 -...
               1)*(x_sc - x_Titan) ...
             + 5*mu_Titan*J3*(R_Titan/au)^3*(x_sc - x_Titan)*(z_sc -...
               z_Titan)/(2*(norm(r_sc - r_Titan))^7)*(7*(z_sc - ...
               z_Titan)^2/(norm(r_sc - r_Titan))^2 - 3) ...
             + 15*mu_Titan*J4*(R_Titan/au)^4*(x_sc - x_Titan)/(8*...
               (norm(r_sc - r_Titan))^7)*(1 - 14*(z_sc - z_Titan)^2/...
               (norm(r_sc - r_Titan))^2 + 21*(z_sc - z_Titan)^4/...
               (norm(r_sc - r_Titan))^4);
    dVy_sc = - mu_Saturn*y_sc/(norm(r_sc))^3 ... 
             - 0.5*rho*Cd*S/M*1e3*au*norm(v_sc - v_Titan)*(Vy_sc -...
               Vy_Titan) ...
             - mu_Titan*(y_sc - y_Titan)/(norm(r_sc - r_Titan))^3 ...
             + 3*mu_Titan*J2*(R_Titan/au)^2/(2*(norm(r_sc - r_Titan))^...
               5)*(5*(z_sc - z_Titan)^2/(norm(r_sc - r_Titan))^2 -...
               1)*(y_sc - y_Titan) ...
             + 5*mu_Titan*J3*(R_Titan/au)^3*(y_sc - y_Titan)*(z_sc -...
               z_Titan)/(2*(norm(r_sc - r_Titan))^7)*(7*(z_sc -...
               z_Titan)^2/(norm(r_sc - r_Titan))^2 - 3) ...
             + 15*mu_Titan*J4*(R_Titan/au)^4*(y_sc - y_Titan)/(8*...
               (norm(r_sc - r_Titan))^7)*(1 - 14*(z_sc - z_Titan)^2/...
               (norm(r_sc - r_Titan))^2 + 21*(z_sc - z_Titan)^4/...
               (norm(r_sc - r_Titan))^4);
    dVz_sc = - mu_Saturn*z_sc/(norm(r_sc))^3 ...
             - 0.5*rho*Cd*S/M*1e3*au*norm(v_sc - v_Titan)*(Vz_sc -...
               Vz_Titan) ...
             - mu_Titan*(z_sc - z_Titan)/(norm(r_sc - r_Titan))^3 ...
             + 3*mu_Titan*J2*(R_Titan/au)^2/(2*(norm(r_sc - r_Titan))^...
               5)*(5*(z_sc - z_Titan)^2/(norm(r_sc - r_Titan))^2 -...
               3)*(z_sc - z_Titan) ...
             + 5*mu_Titan*J3*(R_Titan/au)^3/(2*(norm(r_sc - r_Titan))^...
               5)*(3/5 - 6*(z_sc - z_Titan)^2/(norm(r_sc - r_Titan))^...
               2 + 7*(z_sc - z_Titan)^4/(norm(r_sc - r_Titan))^4) ...
             + 15*mu_Titan*J4*(R_Titan/au)^4*(z_sc - z_Titan)/(8*...
               (norm(r_sc - r_Titan))^7)*(5 - 70*(z_sc - z_Titan)^2/(3*...
               (norm(r_sc - r_Titan))^2) + 21*(z_sc - z_Titan)^4/...
               (norm(r_sc - r_Titan))^4);
else
    dVx_sc = - mu_Saturn*x_sc/(norm(r_sc))^3 ... 
             - mu_Titan*(x_sc - x_Titan)/(norm(r_sc - r_Titan))^3 ...
             + 3*mu_Titan*J2*(R_Titan/au)^2/(2*(norm(r_sc - r_Titan))^...
               5)*(5*(z_sc - z_Titan)^2/(norm(r_sc - r_Titan))^2 -...
               1)*(x_sc - x_Titan) ...
             + 5*mu_Titan*J3*(R_Titan/au)^3*(x_sc - x_Titan)*(z_sc -...
               z_Titan)/(2*(norm(r_sc - r_Titan))^7)*(7*(z_sc -...
               z_Titan)^2/(norm(r_sc - r_Titan))^2 - 3) ...
             + 15*mu_Titan*J4*(R_Titan/au)^4*(x_sc - x_Titan)/(8*...
               (norm(r_sc - r_Titan))^7)*(1 - 14*(z_sc - z_Titan)^2/...
               (norm(r_sc - r_Titan))^2 + 21*(z_sc - z_Titan)^4/...
               (norm(r_sc - r_Titan))^4);
    dVy_sc = - mu_Saturn*y_sc/(norm(r_sc))^3 ... 
             - mu_Titan*(y_sc - y_Titan)/(norm(r_sc - r_Titan))^3 ...
             + 3*mu_Titan*J2*(R_Titan/au)^2/(2*(norm(r_sc - r_Titan))^...
               5)*(5*(z_sc - z_Titan)^2/(norm(r_sc - r_Titan))^2 -...
               1)*(y_sc - y_Titan) ...
             + 5*mu_Titan*J3*(R_Titan/au)^3*(y_sc - y_Titan)*(z_sc -...
               z_Titan)/(2*(norm(r_sc - r_Titan))^7)*(7*(z_sc -...
               z_Titan)^2/(norm(r_sc - r_Titan))^2 - 3) ...
             + 15*mu_Titan*J4*(R_Titan/au)^4*(y_sc - y_Titan)/(8*...
               (norm(r_sc - r_Titan))^7)*(1 - 14*(z_sc - z_Titan)^2/...
               (norm(r_sc - r_Titan))^2 + 21*(z_sc - z_Titan)^4/...
               (norm(r_sc - r_Titan))^4);
    dVz_sc = - mu_Saturn*z_sc/(norm(r_sc))^3 ...
             - mu_Titan*(z_sc - z_Titan)/(norm(r_sc - r_Titan))^3 ...
             + 3*mu_Titan*J2*(R_Titan/au)^2/(2*(norm(r_sc - r_Titan))^...
               5)*(5*(z_sc - z_Titan)^2/(norm(r_sc - r_Titan))^2 - 3)*...
               (z_sc - z_Titan) ...
             + 5*mu_Titan*J3*(R_Titan/au)^3/(2*(norm(r_sc - r_Titan))^...
               5)*(3/5 - 6*(z_sc - z_Titan)^2/(norm(r_sc - r_Titan))^...
               2 + 7*(z_sc - z_Titan)^4/(norm(r_sc - r_Titan))^4) ...
             + 15*mu_Titan*J4*(R_Titan/au)^4*(z_sc - z_Titan)/(8*...
               (norm(r_sc - r_Titan))^7)*(5 - 70*(z_sc - z_Titan)^2/...
               (3*(norm(r_sc - r_Titan))^2) + 21*(z_sc - z_Titan)^4/...
               (norm(r_sc - r_Titan))^4);
end

%Output

dX = [ dx_Titan; dy_Titan; dz_Titan; dVx_Titan; dVy_Titan; dVz_Titan;...
       dx_sc; dy_sc; dz_sc; dVx_sc; dVy_sc; dVz_sc ];

end

function [ lat, long ] = GroundTrack (t, r, v, mu_Titan, i_Titan, omega_T)
%determine the satellite latitude and longitude (t must be in seconds!)

%preallocation
i = zeros(length(t), 1);
RAAN = zeros(length(t), 1);
PA = zeros(length(t), 1);
TA = zeros(length(t), 1);
Xi = zeros(length(t), 1);
omega = zeros(length(t), 1);
lat = zeros(length(t), 1);
long = zeros(length(t), 1);

for j = 1:length(t)
    [ ~, ~, i(j), RAAN(j), PA(j), TA(j) ] = rv2coe (r(j, :)', ...
                                            v(j, :)', mu_Titan);
    i(j) = i(j) + i_Titan;                %inclination w.r.t. Titan 
    %rotation axis [rad]
    Xi(j) = PA(j) + TA(j);                %latitude argument [rad]
    omega(j) = RAAN(j) - omega_T*t(j);    %geographical longitude of
    %ascending node [rad]
    lat(j) = asin(sin(i(j))*sin(Xi(j)));  %latitude [rad]
    long(j) = atan2((sin(omega(j))*cos(Xi(j)) + cos(omega(j))*...
              sin(Xi(j))*cos(i(j))), (cos(omega(j))*cos(Xi(j)) -...
              sin(omega(j))*sin(Xi(j))*...
              cos(i(j))));                %longitude [rad]
end

end
