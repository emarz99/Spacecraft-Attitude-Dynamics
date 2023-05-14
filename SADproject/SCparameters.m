 %% Space Attitude Dynamics
clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% Relevant Data
muE  = astroConstants(13);     % Earth planetary constant
R_E  = astroConstants(23);     % Earth mean radius
i_E  = astroConstants(8);      % Earth Obliquity
c    = astroConstants(5);      % Speed of light
T_E  = 365.256*24*3600;        % Earth orbit period 
AU   = astroConstants(2);      % Astronomical unit 
i_m  = deg2rad(11.5);          % Magnetic axis inclination
om_E = 7.292*1e-5;             % Earth angular speed [rad/s]

%% SRP distubance
Fe_Sun       = astroConstants(31);     % Direct solar radiation 
Fe_reflected = 500;                    % Radiation reflected by the Earth
Fe_E         = 117;                    % Earth radiation

n_Sun = 2*pi/T_E;

%% Magnetic distubance
m    = [0.01 0.05 0.01]';      % Residual magnetic induction

g_01 = -29615*1e-9;      g_11 = -1728*1e-9;       h_11 = 5186*1e-9;

H0 = (g_01^2 + g_11^2 + h_11^2)^(1/2);

%% Orbit Characterisation

% Sun-Synchronous Polar Orbit (Circular R=a)
altittude = 1200;       %[km]

a   = R_E+altittude;               %[km]  
e   = 0;
i   = acos(-(a/12352)^(7/2));      %[rad]
OM  = deg2rad(0);                  %[rad]
om  = deg2rad(0);                  %[rad]
th0 = deg2rad(0);                  %[rad]

i_deg = rad2deg(i);

Ai_LN = [1    0      0   ;
         0  cos(i) sin(i);
         0 -sin(i) cos(i)];

if a/R_E > csc(i-i_E)
    disp('Orbit without eclipse phase')
end

c_m = (R_E^3*H0)/a^3;
n   = sqrt(muE/a^3);        % Angular speed
Tsc = 2*pi*sqrt(a^3/muE);   % Orbit period

%% Initial condition (Body frame)

om0 = [1 1 1]';

%% S/C Design (6U Cubesat)

% Main body
n1 = [1 0 0];     n4 = -n1;
n2 = [0 1 0];     n5 = -n2;
n3 = [0 0 1];     n6 = -n3;

X = 0.2;    %[m]
Y = 0.3;    %[m]
Z = 0.1;    %[m]

% Solar panels
i_p = deg2rad(30);
n7 = [0  cos(i_p) sin(i_p)];     n8  = -n7;
n9 = [0 -cos(i_p) sin(i_p)];     n10 = -n9;

x = 0.2;       %[m]
y = 0.3;    %[m]  thickness
z = 0.0025;       %[m] 

N_Bi = [n1; n2; n3; n4; n5; n6; n7; n8; n9; n10];

A1 = Y*Z;   A4 = A1;    
A2 = X*Z;   A5 = A2;    
A3 = X*Y;   A6 = A3;   A7 = A3;   A8 = A3;   A9 = A3;   A10 = A3;

Ai = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10];

b0 = zeros(1,3);
b1 = [-X/2*(1+cos(i_p))   0   Z/2*(1+sin(i_p))];
b2 = [-X/2*(1+cos(i_p))   0  -Z/2*(1+sin(i_p))];

Ri = [b0; b0; b0; b0; b0; b0; b1; b1; b2; b2];

% Inertia Matrix
M_bus    = 6*1.33;     %[kg]   % 1.33kg per U
M_panels = 12*50e-3;   %[kg]   % 50g per U

I_body  = (M_bus/3)*diag([Y^2+Z^2, X^2+Z^2, X^2+Y^2]); 

I_panel1 = (M_panels/3)*diag([y^2+z^2, x^2+z^2, x^2+y^2]) + M_panels*(dot(b1,b1)*eye(3)-b1'*b1); 
I_panel2 = (M_panels/3)*diag([y^2+z^2, x^2+z^2, x^2+y^2]) + M_panels*(dot(b2,b2)*eye(3)-b2'*b2); 

J = I_body + I_panel1 + I_panel2;

%% Thrusters (8 thrusters configuration)
T = 0.005;       %[N]

i_t1 = deg2rad(0);
i_t2 = deg2rad(90);

b_t1 = diag([Y/2*sin(i_t1)+Z/2*cos(i_t1) X/2*sin(i_t1) X/2*cos(i_t1)]);
b_t2 = diag([Y/2*sin(i_t2)+Z/2*cos(i_t2) X/2*sin(i_t2) X/2*cos(i_t2)]);

R_t1 = [-1   0   1;
         1   0  -1;
        -1   0  -1;
         1   0   1]*b_t1;

R_t2 = [ 1  -1   0;
        -1  -1   0;
         1   1   0;
        -1   1   0]*b_t2;


R_t = [R_t1; R_t2];

%% Inertia Wheels (pyramid configuration) 

maxT_iw = 0.1*1e-3;       %[Nm]
maxL_iw = 6*1e-3;         %[Nms]

c1_iw = 1/sqrt(3);
c2_iw = sqrt(3)/4;

A_iw     = [-1  1  1 -1;
            -1 -1  1  1;
             1  1  1  1]*c1_iw;

pinvA_iw = [-1 -1  1;
             1 -1  1;
             1  1  1;
            -1  1  1]*c2_iw;

%% Observer and Linear control 

% State obeserver matrices
Kx = (J(3,3) - J(2,2))/J(1,1);
Ky = (J(3,3) - J(1,1))/J(2,2);
Kz = (J(2,2) - J(1,1))/J(3,3);

A_lin = [     zeros(3)       ,        eye(3)     ;
         -Kx*n^2     0    0     0    (1-Kx)*n   0;
             0   -Ky*n^2  0  (Ky-1)*n    0      0;
             0       0    0     0        0      0];

B_lin = [zeros(3); diag(1./diag(J))];
C_lin = [eye(3), zeros(3)];

%tuning parameters (Q weigthing matrix)
W = diag(ones(3,1)); 
Q = C_lin'*W*C_lin ; 
Q(4:6, 4:6) = J;

%tuning parameters (R weigthing matrix)
M_max = A_iw(1,:)*(maxT_iw*[0 1 1 0]')*ones(3,1);

R = diag(1./(M_max));

S = zeros(6,3); 

[K,P,CLP] = lqr(A_lin,B_lin,Q,R,S);     % CLP: eig(A-B*K) closed loop poles

L = place(A_lin',C_lin',CLP)';

%% Sensors Accuracy

bias_gyro = 0.3;   %[deg/h]
ARW_gyro  = 0.15;  %[deg/sqrt(h)]

err_SS = 0.1;  %[deg]

err_HS = 0.1;  %[deg]
