% This file is used for stability analyisis for both the full and the reduced set of 
% equations of motion and it is based on the cit2s.m file
%
% Calculation of state matrix and input matrix for calculation
% of symmetric aircraft response to atmospheric turbulence.
% The system model is in the form
%
%       .
%       x = Ax + Bu
%       -    -    -
%
% with 
%     x = [u/V alpha theta qc/V u_g/V alpha_g alpha_g*]' 
% and 
%     u = [delta_e w_1 w_3]'.
%
% The turbulence filters are according to Dryden.

%
%  Aircraft: Cessna Citation I Ce-500
%

% INPUT TURBULENCE- AND AIRCRAFT PARAMETERS

% AIRCRAFT FLIGHT CONDITION 'CRUISE'.
V     = 121.3;
m     = 5445;
twmuc = 2*209;
KY2   = 0.950;
c     = 2.022;
S     = 24.2;
lh    = 5.5;
g     = 9.80665;

% TURBULENCE PARAMETERS
sigma = 1;   %m/sec
Lg    = 150; %m

sigmaug_V = sigma/V;
sigmaag   = sigma/V;

% AIRCRAFT SYMMETRIC AERODYNAMIC DERIVATIVES : 
CX0 = 0.0000;     CZ0  =-0.5640;     Cm0  =  0.0000;
CXu =-0.0850;     CZu  =-1.1280;     Cmu  =  0.0000;
CXa = 0.2140;     CZa  =-5.4300;     Cma  = -0.5180;
CXq = 0.0000;     CZq  =-4.0700;     Cmq  = -7.3500;
CXd = 0.0000;     CZd  =-0.5798;     Cmd  = -1.4440;
CXfa= 0.0000;     CZfa =-1.6200;     Cmfa = -4.2000;
                  CZfug= 0.0000;     Cmfug= -Cm0*lh/c;
                  CZfag= CZfa-CZq;   Cmfag=  Cmfa-Cmq;

% CALCULATION OF AIRCRAFT SYMMETRIC STABILITY DERIVATIVES
xu   = (V/c)*(CXu/twmuc);
xa   = (V/c)*(CXa/twmuc);
xt   = (V/c)*(CZ0/twmuc);
xq   = 0;
xd   = (V/c)*(CXd/twmuc);
xug  = xu;
xfug = 0;
xag  = xa;
xfag = 0;

zu   = (V/c)*( CZu/(twmuc-CZfa));
za   = (V/c)*( CZa/(twmuc-CZfa));
zt   = (V/c)*(-CX0/(twmuc-CZfa));
zq   = (V/c)*((CZq+twmuc)/(twmuc-CZfa));
zd   = (V/c)*( CZd/(twmuc-CZfa));
zug  = zu;
zfug = (V/c)*( CZfug/(twmuc-CZfa));
zag  = za;
zfag = (V/c)*( CZfag/(twmuc-CZfa));

mu   = (V/c)*(( Cmu+CZu*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
ma   = (V/c)*(( Cma+CZa*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mt   = (V/c)*((-CX0*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mq   = (V/c)*(Cmq+Cmfa*(twmuc+CZq)/(twmuc-CZfa))/(twmuc*KY2);
md   = (V/c)*((Cmd+CZd*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mug  = mu;
mfug = (V/c)*(Cmfug+CZfug*Cmfa/(twmuc-CZfa))/(twmuc*KY2);
mag  = ma;
mfag = (V/c)*(Cmfag+CZfag*Cmfa/(twmuc-CZfa))/(twmuc*KY2);

%------------------------------------------------------------------------------
%----------------------- STABILITY ANALYSIS -----------------------------------
%------------------------------------------------------------------------------


%% State space matrices calculation for the complete set of equations of motion
A_complete_system=[xu xa xt 0    xug                  xag       0;
   zu za zt zq   zug-zfug*V/Lg*(c/V)  zag       zfag*(c/V);
   0  0  0  V/c  0                    0         0;
   mu ma mt mq   mug-mfug*V/Lg*(c/V)  mag       mfag*(c/V);
   0  0  0  0   -V/Lg                 0         0;
   0  0  0  0    0                    0         1;
   0  0  0  0    0                   -(V/Lg)^2 -2*V/Lg];

B_complete_system=...
 [xd 0                                 0;
  zd zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) zfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  0                                 0;
  md mfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) mfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  sigmaug_V*sqrt(2*V/Lg)            0;
  0  0                                 sigmaag*sqrt(3*V/Lg);
  0  0                                 (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3)];


%C_complete_system = eye(7); D_complete_system = zeros(7,3);


nz_row = (V/g)*(A_complete_system(3,:)-A_complete_system(2,:));
C_complete_system = [eye(4) zeros(4,3);
                     nz_row]; 
D_complete_system = [zeros(4,3);
                     (V/g)*(B_complete_system(3,:)-B_complete_system(2,:))];
                 
                     %-V/g*[zd, zfug*c/V*sigmaug_V*sqrt(2*V/Lg), zfag*(c/V)*sigmaag*sqrt(3*V/Lg)]]
                

complete_system= ss(A_complete_system, B_complete_system, C_complete_system, D_complete_system);

complete_system.StateName    = {'u','\alpha','\theta','qc/V','u_g','alpha_g','alpha_g_dot'};
complete_system.InputName    = {'\delta_{el}','w1','w3'};
complete_system.OutputName   = {'u','\alpha','\theta','qc/V','n_z'};

poles_complete_system = pole(complete_system); %this extracts the eigenvaues that are used in stability analysis

%[wn,zeta, poles] = damp(complete_system); %this line also extracts the natural frequency and damping of the poles if needed

pzmap(complete_system);
grid on;
hm = findobj(gca, 'Type', 'Line');          % Handle To 'Line' Objects
hm(2).MarkerSize = 15;                      % ‘Zero’ Marker
hm(3).MarkerSize = 15; 

%% State space matrices calculation for the reduced set of equations of motion,
%% for the short-period assumptions 


A_sp=[ za 0 zq   zug-zfug*V/Lg*(c/V)  zag       zfag*(c/V);
       0  0  V/c  0                    0         0;
       ma 0 mq   mug-mfug*V/Lg*(c/V)  mag       mfag*(c/V);
       0  0  0   -V/Lg                 0         0;
       0  0  0    0                    0         1;
       0  0  0    0                   -(V/Lg)^2 -2*V/Lg];

B_sp= [zd zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) zfag*(c/V)*sigmaag*sqrt(3*V/Lg);
       0  0                                 0;
       md mfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) mfag*(c/V)*sigmaag*sqrt(3*V/Lg);
       0  sigmaug_V*sqrt(2*V/Lg)            0;
       0  0                                 sigmaag*sqrt(3*V/Lg);
       0  0                                 (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3)];

 
nz_row_sp = (V/g)*(A_sp(2,:)-A_sp(1,:));

C_sp = [eye(3) zeros(3,3);
                 nz_row_sp];
                 
D_sp = [zeros(3,3);
        (V/g)*(B_sp(2,:)-B_sp(1,:))];

sp_system= ss(A_sp, B_sp, C_sp, D_sp);
poles_sp = pole(sp_system);

% pzmap(sp_system);
% grid on;
% hm = findobj(gca, 'Type', 'Line');          % Handle To 'Line' Objects
% hm(2).MarkerSize = 15;                      % ‘Zero’ Marker
% hm(3).MarkerSize = 15;   

sp_system.StateName    = {'\alpha','\theta','qc/V','u_g','alpha_g','alpha_g_dot'};
sp_system.InputName    = {'\delta_{el}','w1','w3'};
sp_system.OutputName   = {'\alpha','\theta','qc/V','n_z'};

save complete_system 
save sp_system

