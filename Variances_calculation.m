%In this file, the variances for both the complete and the
%reduced aircraft model are calculated

Spectral_analysis

%--------------------------------------------------------------------------
%% ---------------- USING THE ANALYTICAL POWER SPECTRA --------------------
%--------------------------------------------------------------------------

%% ----------- For the complete model -----------
%NUMERICAL INTEGRATION OF PSD ’s
do = diff(omega)' ; % compute " difference vector" in omega
% i.e., omega(k+1) - omega(k);
% then perform ( very crude) integration
var_complete(1) = sum(do.*Suu_1(1:Nf-1)) ;
var_complete(2) = sum(do.*Saa_1(1:Nf-1)) ;
var_complete(3) = sum(do.*Stt_1(1:Nf-1)) ;
var_complete(4) = sum(do.*Sqq_1(1:Nf-1)) ;
var_complete(5) = sum(do.*Snn_1(1:Nf-1)) ;
disp (  'Complete model: Numerical integration of PSD yields the variances :' )
var_analytical = var_complete/pi;
disp(var_analytical);

%% ----------- For the reduced model -----------
%NUMERICAL INTEGRATION OF PSD ’s
do = diff(omega)' ; % compute " difference vector" in omega
% i.e., omega(k+1) - omega(k);
% then perform ( very crude) integration
var_sp(1) = sum(do.*Saa_1_sp(1:Nf-1)) ;
var_sp(2) = sum(do.*Stt_1_sp(1:Nf-1)) ;
var_sp(3) = sum(do.*Sqq_1_sp(1:Nf-1)) ;
var_sp(4) = sum(do.*Snn_1_sp(1:Nf-1)) ;
disp (  'Reduced model: Numerical integration of PSD yields the variances :' )
var_analytical_sp = var_sp/pi;
disp(var_analytical_sp);

%--------------------------------------------------------------------------
%% ---------------- USING THE IMPULSE RESPONSE METHOD --------------------
%--------------------------------------------------------------------------
%based on exampl72.m

dt = 0.01; 
T  = 200; 
t = [0:dt:T];
N = length(t);

%% ----------- For the complete model -----------


% ZERO INPUT, INITIAL CONDITION EQUALS B (for input 2 and 3 in this case)
u = zeros(3,N); x0=B_complete_system(:,2)+B_complete_system(:,3);

% CALCULATION OF IMPULSE RESPONSES
h = lsim(A_complete_system,B_complete_system,C_complete_system,D_complete_system,u,t,x0);	

% PLOT IMPULSE RESPONSE
subplot(2,1,1)
plot(t,h(:,1)); xlabel('Time [sec]'); ylabel('h u w3(t)');
subplot(2,1,2)
plot(t,h(:,2)); xlabel('Time [sec]'); ylabel('h alpha w3(t)');


% CALCULATION OF PRODUCT MATRIX OF IMPULSE RESPONSES
h11=h(:,1).*h(:,1);h12=h(:,1).*h(:,2);h13=h(:,1).*h(:,3);h14=h(:,1).*h(:,4);h15=h(:,1).*h(:,5);
                   h22=h(:,2).*h(:,2);h23=h(:,2).*h(:,3);h24=h(:,2).*h(:,4);h25=h(:,2).*h(:,5);
                                      h33=h(:,3).*h(:,3);h34=h(:,3).*h(:,4);h35=h(:,3).*h(:,5);
                                                         h44=h(:,4).*h(:,4);h45=h(:,4).*h(:,5);
                                                                            h55=h(:,5).*h(:,5);
 % PLOT (CROSS) PRODUCTS OF IMPULSE RESPONSES
plot(t,h11); xlabel('Time [sec]'); ylabel('h1*h1(t)');


% INTEGRATION OF PRODUCT MATRIX OF IMPULSE RESPONSES
var11(1)=0; var12(1)=0; var13(1)=0; var14(1)=0;
            var22(1)=0; var23(1)=0; var24(1)=0;
                        var33(1)=0; var34(1)=0;
                                    var44(1)=0; var55(1)=0;
                                    
dth11 = dt*h11; dth12 = dt*h12; dth13 = dt*h13; dth14 = dt*h14;
                dth22 = dt*h22; dth23 = dt*h23; dth24 = dt*h24;
                                dth33 = dt*h33; dth34 = dt*h34;
                                                dth44 = dt*h44;
                                                dth55 = dt*h55;
                                                
 % Integral loop:                                                               
for i=1:N-1
    
    var11(i+1) = var11(i) + dth11(i);
    var12(i+1) = var12(i) + dth12(i);
    var13(i+1) = var13(i) + dth13(i);
    var14(i+1) = var14(i) + dth14(i);
    var22(i+1) = var22(i) + dth22(i);
    var23(i+1) = var23(i) + dth23(i);
    var24(i+1) = var24(i) + dth24(i);
    var33(i+1) = var33(i) + dth33(i);
    var34(i+1) = var34(i) + dth34(i);
    var44(i+1) = var44(i) + dth44(i);
    var55(i+1) = var55(i) + dth55(i);
end

disp (  'Complete model: Using impulse response method yields the variances :' )
var_impulse = [var11(end) var22(end) var33(end) var44(end) var55(end)];
disp(var_impulse)

%% ----------- For the reduced model -----------

% ZERO INPUT, INITIAL CONDITION EQUALS B (for input 2 and 3 in this case)
u = zeros(3,N); x0_sp=B_sp(:,2)+B_sp(:,3);

% CALCULATION OF IMPULSE RESPONSES
h_sp = lsim(A_sp,B_sp,C_sp,D_sp,u,t,x0_sp);	

% PLOT IMPULSE RESPONSE
subplot(2,1,1)
plot(t,h_sp(:,1)); xlabel('Time [sec]'); ylabel('h u w3(t)');
subplot(2,1,2)
plot(t,h_sp(:,2)); xlabel('Time [sec]'); ylabel('h alpha w3(t)');


% CALCULATION OF PRODUCT MATRIX OF IMPULSE RESPONSES
h11_sp=h_sp(:,1).*h_sp(:,1);h12_sp=h_sp(:,1).*h_sp(:,2);h13_sp=h_sp(:,1).*h_sp(:,3);h14_sp=h_sp(:,1).*h_sp(:,4);
                            h22_sp=h_sp(:,2).*h_sp(:,2);h23_sp=h_sp(:,2).*h_sp(:,3);h24_sp=h_sp(:,2).*h_sp(:,4);
                                                        h33_sp=h_sp(:,3).*h_sp(:,3);h34_sp=h_sp(:,3).*h_sp(:,4);
                                                                                    h44_sp=h_sp(:,4).*h_sp(:,4);
                                                                                                               
 % PLOT (CROSS) PRODUCTS OF IMPULSE RESPONSES
plot(t,h11_sp); xlabel('Time [sec]'); ylabel('h1*h1(t)');


% INTEGRATION OF PRODUCT MATRIX OF IMPULSE RESPONSES
var11_sp(1)=0; var12_sp(1)=0; var13_sp(1)=0; var14_sp(1)=0;
               var22_sp(1)=0; var23_sp(1)=0; var24_sp(1)=0;
                              var33_sp(1)=0; var34_sp(1)=0;
                                             var44_sp(1)=0; 
                                    
dth11_sp = dt*h11_sp; dth12_sp = dt*h12_sp; dth13_sp = dt*h13_sp; dth14_sp = dt*h14_sp;
                      dth22_sp = dt*h22_sp; dth23_sp = dt*h23_sp; dth24_sp = dt*h24_sp;
                                            dth33_sp = dt*h33_sp; dth34_sp = dt*h34_sp;
                                                                  dth44_sp = dt*h44_sp;
                                                                  
                                                
 % Integral loop:                                                               
for i=1:N-1
    
    var11_sp(i+1) = var11_sp(i) + dth11_sp(i);
    var12_sp(i+1) = var12_sp(i) + dth12_sp(i);
    var13_sp(i+1) = var13_sp(i) + dth13_sp(i);
    var14_sp(i+1) = var14_sp(i) + dth14_sp(i);
    var22_sp(i+1) = var22_sp(i) + dth22_sp(i);
    var23_sp(i+1) = var23_sp(i) + dth23_sp(i);
    var24_sp(i+1) = var24_sp(i) + dth24_sp(i);
    var33_sp(i+1) = var33_sp(i) + dth33_sp(i);
    var34_sp(i+1) = var34_sp(i) + dth34_sp(i);
    var44_sp(i+1) = var44_sp(i) + dth44_sp(i);
end

disp (  'Reduced model: Using impulse response method yields the variances :' )
var_impulse_sp = [var11_sp(end) var22_sp(end) var33_sp(end) var44_sp(end)];
disp(var_impulse_sp)

%--------------------------------------------------------------------------
%% ---------------- USING THE VAR.M ROUTINE -------------------------------
%--------------------------------------------------------------------------
%% ----------- For the complete model -----------

for i=1:5
  var_2(i)= var(y(:,i));
end

disp (  'Complete model: Using the var.m routine yields the variances :' )
disp(var_2)

%% ----------- For the reduced model -----------

for i=1:4
  var_2_sp(i)= var(y_sp(:,i));
end

disp (  'Reduced model: Using the var.m routine yields the variances :' )
disp(var_2_sp)






