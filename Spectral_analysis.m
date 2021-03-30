%In this file, the spectral analysis for both the complete and the
%reduced aircraft model is performed
% Based on the examp73.m file

load complete_system
load sp_system

plotting = "True";
%--------------------------------------------------------------------------
%% ----------------------- ANALYTICAL METHOD ------------------------------
%--------------------------------------------------------------------------

% TIME AXIS INPUT VECTOR DEFINITION
        dt = 0.01; %if you use 0.1 plots look better
        T  = 200; 
        t = [0:dt:T];
        N = length(t);
        Fs = 1/dt;


        
        % GET INPUT PARAMETERS
        Wc = sigma;
        Nf = N;

        % DEFINE NOISE INTENSITY
        W  = Wc/dt;    % discrete time covariance, remember?
      
        % DEFINE FREQUENCY AXIS
        omega = logspace(-2,2,Nf);
        
        %------------------------------------------------------------------
        
        %--------------------- For the complete model ---------------------
        
        %------------------------------------------------------------------

        % COMPUTE FREQUENCY RESPONSE
        mag_1 = bode(A_complete_system,B_complete_system,C_complete_system,D_complete_system,2,omega);
        mag_2 = bode(A_complete_system,B_complete_system,C_complete_system,D_complete_system,3,omega);
        mag = mag_1 + mag_2;

        % COMPUTE POWER SPECTRA
        Suu_1 = mag(:,1).^2;
        Saa_1 = mag(:,2).^2;
        Stt_1 = mag(:,3).^2;
        Sqq_1 = mag(:,4).^2;
        Snn_1 = mag(:,5).^2;

        if plotting=="True"
        % PLOT POWER SPECTRA
            figure(1)
            subplot(2,3,1); 
            loglog(omega,Suu_1); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\hat{u}\hat{u}}$ [$\frac{1}{rad/s}$]','interpreter','latex');
            grid on
            subplot(2,3,2); 
            loglog(omega,Saa_1); xlabel('$\omega$ [rad/sec]','interpreter','latex');  ylabel('$S_{\alpha\alpha}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on
            subplot(2,3,3); 
            loglog(omega,Stt_1); xlabel('$\omega$ [rad/sec]','interpreter','latex');  ylabel('$S_{\theta\theta}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on
            subplot(2,3,4); 
            loglog(omega,Sqq_1);xlabel('$\omega$ [rad/sec]','interpreter','latex');  ylabel('$S_{\frac{qc}{V}\frac{qc}{V}}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on
            subplot(2,3,5); 
            loglog(omega,Snn_1); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{n_z n_z}$ [$\frac{1}{rad/s}$]','interpreter','latex');
            grid on
        end
        
        %------------------------------------------------------------------
        
        %--------------------- For the reduced model ---------------------
        
        %------------------------------------------------------------------

        % COMPUTE FREQUENCY RESPONSE
        mag_sp_1 = bode(A_sp,B_sp,C_sp,D_sp,2,omega);
        mag_sp_2 = bode(A_sp,B_sp,C_sp,D_sp,3,omega);
        mag_sp = mag_sp_1 + mag_sp_2;

        % COMPUTE POWER SPECTRA
        Saa_1_sp = mag_sp(:,1).^2;
        Stt_1_sp = mag_sp(:,2).^2;
        Sqq_1_sp = mag_sp(:,3).^2;
        Snn_1_sp = mag_sp(:,4).^2;
        
        if plotting=="True" 
   
            % PLOT POWER SPECTRA
            figure(2)
            subplot(2,2,1); 
            loglog(omega,Saa_1_sp); xlabel('$\omega$ [rad/sec]','interpreter','latex');  ylabel('$S_{\alpha\alpha}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on
            subplot(2,2,2); 
            loglog(omega,Stt_1_sp); xlabel('$\omega$ [rad/sec]','interpreter','latex');  ylabel('$S_{\theta\theta}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on
            subplot(2,2,3); 
            loglog(omega,Sqq_1_sp);xlabel('$\omega$ [rad/sec]','interpreter','latex');  ylabel('$S_{\frac{qc}{V}\frac{qc}{V}}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on
            subplot(2,2,4); 
            loglog(omega,Snn_1_sp); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{n_z n_z}$ [$\frac{1}{rad/s}$]','interpreter','latex');
            grid on
        end
        
%--------------------------------------------------------------------------
%% ----------------------- EXPERIMENTAL METHOD FFT.m ----------------------
%-------------------------------------------------------------------------- 


%------------------------------------------------------------------
        
%--------------------- For the complete model ---------------------
        
%------------------------------------------------------------------
        
nn = zeros(1,N);             % zero input elevator
u_g = randn(1,N)/sqrt(dt);    % scaled input long. turbulence,
                             % note the sqrt(dt) because of lsim
w_g = randn(1,N)/sqrt(dt);    % scaled input vert. turbulence,
                             % note the sqrt(dt) because of lsim
u1  = [nn' u_g' nn'];          % input vector definition (combined
                             % longitudinal and vertical turbulence).
u2  = [nn' nn' w_g'];
        
y1 = lsim(A_complete_system,B_complete_system,C_complete_system,D_complete_system,u1,t);
y2 = lsim(A_complete_system,B_complete_system,C_complete_system,D_complete_system,u2,t);
        
y= y1 + y2;

% FFT ALL SIGNALS
U      = dt*fft(y(:,1));
ALPHA  = dt*fft(y(:,2));
THETA  = dt*fft(y(:,3));
QCV    = dt*fft(y(:,4));
NZ     = dt*fft(y(:,5));

% COMPUTE PSDs
Suu_2  = (1/T)*     U.*conj(U);
Saa_2  = (1/T)* ALPHA.*conj(ALPHA);
Stt_2  = (1/T)* THETA.*conj(THETA);
Sqq_2  = (1/T)*   QCV.*conj(QCV);
Snn_2  = (1/T)*    NZ.*conj(NZ);


% DEFINE FREQUENCY VECTOR FOR PLOTTING
omega_2 = 2*pi*Fs*(0:(N/2)-1)/N;
        

if plotting=="True"
    % PLOT POWER SPECTRA
    figure(3)

    subplot(2,3,1); 
    loglog(omega,Suu_1); 
    hold on
    loglog(omega_2,Suu_2(1:N/2)); 
    xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\hat{u}\hat{u}}$ [$\frac{1}{rad/s}$]','interpreter','latex');
    grid on

    subplot(2,3,2); 
    loglog(omega,Saa_1);
    hold on
    loglog(omega_2,Saa_2(1:N/2));xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\alpha\alpha}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
    grid on

    subplot(2,3,3); 
    loglog(omega,Stt_1);
    hold on
    loglog(omega_2,Stt_2(1:N/2));xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\theta\theta}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
    grid on

    subplot(2,3,4);
    loglog(omega,Sqq_1);
    hold on
    loglog(omega_2,Sqq_2(1:N/2)); xlabel('$\omega$ [rad/sec]','interpreter','latex');ylabel('$S_{\frac{qc}{V}\frac{qc}{V}}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');

    grid on

    subplot(2,3,5); 
    loglog(omega,Snn_1);
    hold on
    loglog(omega_2,Snn_2(1:N/2)); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{n_z n_z}$ [$\frac{1}{rad/s}$]','interpreter','latex');
    grid on
end
%------------------------------------------------------------------
        
%--------------------- For the reduced model ---------------------
        
%------------------------------------------------------------------


dt_sp = dt;
T_sp = T;
t_sp = 0:dt_sp:T_sp;
N_sp = length(t_sp);

nn_sp = zeros(1,N_sp);                % zero input elevator
u_g_sp = randn(1,N_sp)/sqrt(dt_sp);    % scaled input long. turbulence,
                                      % note the sqrt(dt) because of lsim
w_g_sp = randn(1,N_sp)/sqrt(dt_sp);    % scaled input vert. turbulence,
                                      % note the sqrt(dt) because of lsim
u_sp_1  = [nn_sp' u_g_sp' nn_sp'];     % input vector definition (combined
u_sp_2  = [nn_sp' nn_sp' w_g_sp'];     % longitudinal and vertical turbulence).

% SIMULATION OF MOTION VARIABLES
y1_sp = lsim(A_sp,B_sp,C_sp,D_sp,u_sp_1,t_sp);
y2_sp = lsim(A_sp,B_sp,C_sp,D_sp,u_sp_2,t_sp);

y_sp = y1_sp + y2_sp;

% FFT ALL SIGNALS
ALPHA_sp  = dt*fft(y_sp(:,1));
THETA_sp  = dt*fft(y_sp(:,2));
QCV_sp    = dt*fft(y_sp(:,3));
NZ_sp     = dt*fft(y_sp(:,4));

% COMPUTE PSDs
Saa_2_sp  = (1/T_sp)* ALPHA_sp.*conj(ALPHA_sp);
Stt_2_sp  = (1/T_sp)* THETA_sp.*conj(THETA_sp);
Sqq_2_sp  = (1/T_sp)*   QCV_sp.*conj(QCV_sp);
Snn_2_sp  = (1/T_sp)*    NZ_sp.*conj(NZ_sp);

% DEFINE FREQUENCY VECTOR FOR PLOTTING
omega_2_sp = 2*pi*Fs*(0:(N_sp/2)-1)/N_sp;
        
        if plotting=="True"
            % PLOT POWER SPECTRA
            figure(4)

            subplot(2,2,1); 
            loglog(omega_2_sp,Saa_2_sp(1:N/2));xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\alpha\alpha}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            hold on
            loglog(omega,Saa_1_sp);
            grid on

            subplot(2,2,2); 
            loglog(omega,Stt_1_sp);
            hold on
            loglog(omega_2_sp,Stt_2_sp(1:N/2));xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\theta\theta}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on

            subplot(2,2,3);
            loglog(omega,Sqq_1_sp);
            hold on
            loglog(omega_2_sp,Sqq_2_sp(1:N/2)); xlabel('$\omega$ [rad/sec]','interpreter','latex');ylabel('$S_{\frac{qc}{V}\frac{qc}{V}}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');

            grid on

            subplot(2,2,4); 
            loglog(omega,Snn_1_sp);
            hold on
            loglog(omega_2_sp,Snn_2_sp(1:N/2)); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{n_z n_z}$ [$\frac{1}{rad/s}$]','interpreter','latex');
            grid on

        end 
%--------------------------------------------------------------------------
%% ----------------------- EXPERIMENTAL METHOD ----------------------------
%% ------------------------ SMOOTHING FILTER ------------------------------
%-------------------------------------------------------------------------- 


%------------------------------------------------------------------
        
%--------------------- For the complete model ---------------------
        
%------------------------------------------------------------------

Sxx_2 =[Suu_2(1:N/2) Saa_2(1:N/2) Stt_2(1:N/2) Sqq_2(1:N/2) Snn_2(1:N/2)];
Sxx_filt = [zeros(length(Suu_2(1:N/2)),1) zeros(length(Saa_2(1:N/2)),1) zeros(length(Stt_2(1:N/2)),1) zeros(length(Sqq_2(1:N/2)),1) zeros(length(Snn_2(1:N/2)),1)];

for i=1:5
    
    Sxx_filt(1,i) = Sxx_2(1,i);
    
    for j=2:(length(Sxx_2(:,i))-1)
        Sxx_filt(j,i)= 0.25*Sxx_2(j-1,i) + 0.5*Sxx_2(j,i) + 0.25*Sxx_2(j+1,i);
    end
        
    Sxx_filt(length(Sxx_2(:,i)),i)= Sxx_2(length(Sxx_2(:,i)),i);

end   
        

        if plotting=="True"
            figure(5)

            subplot(3,2,1); 
            loglog(omega,Suu_1);
            hold on
            loglog(omega_2,Suu_2(1:N/2)); 
            hold on
            loglog(omega_2, Sxx_filt(:,1),'color',[0 0.5 0]); 
            xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\hat{u}\hat{u}}$ [$\frac{1}{rad/s}$]','interpreter','latex');
            grid on

            subplot(3,2,2); 
            loglog(omega,Saa_1);
            hold on
            loglog(omega_2,Saa_2(1:N/2)); 
            hold on
            loglog(omega_2, Sxx_filt(:,2),'color',[0 0.5 0]); 
            xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\alpha\alpha}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on

            subplot(3,2,3); 
            loglog(omega,Stt_1);
            hold on
            loglog(omega_2,Stt_2(1:N/2)); 
            hold on
            loglog(omega_2, Sxx_filt(:,3),'color',[0 0.5 0]); 
            xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\theta\theta}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on 

            subplot(3,2,4); 
            loglog(omega,Sqq_1);
            hold on
            loglog(omega_2,Sqq_2(1:N/2)); 
            hold on
            loglog(omega_2, Sxx_filt(:,4),'color',[0 0.5 0]); 
            xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\frac{qc}{V}\frac{qc}{V}}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on

            subplot(3,2,5); 
            loglog(omega,Snn_1);
            hold on
            loglog(omega_2,Snn_2(1:N/2)); 
            hold on
            loglog(omega_2, Sxx_filt(:,5),'color',[0 0.5 0]); 
            xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{n_z n_z}$ [$\frac{1}{rad/s}$]','interpreter','latex');
            grid on
        end

%------------------------------------------------------------------
        
%--------------------- For the reduced model ---------------------
        
%------------------------------------------------------------------

Sxx_2_sp =[Saa_2_sp(1:N/2) Stt_2_sp(1:N/2) Sqq_2_sp(1:N/2) Snn_2_sp(1:N/2)];
Sxx_filt_sp = [zeros(length(Saa_2_sp(1:N/2)),1) zeros(length(Stt_2_sp(1:N/2)),1) zeros(length(Sqq_2_sp(1:N/2)),1) zeros(length(Snn_2_sp(1:N/2)),1)];

for i=1:4
    
    Sxx_filt_sp(1,i) = Sxx_2_sp(1,i);
    
    for j=2:(length(Sxx_2_sp(:,i))-1)
        Sxx_filt_sp(j,i)= 0.25*Sxx_2_sp(j-1,i) + 0.5*Sxx_2_sp(j,i) + 0.25*Sxx_2_sp(j+1,i);
    end
        
    Sxx_filt_sp(length(Sxx_2_sp(:,i)),i)= Sxx_2_sp(length(Sxx_2_sp(:,i)),i);

end   

        if plotting=="True"
            figure(6)

            subplot(2,2,1); 
            loglog(omega,Saa_1_sp);
            hold on
            loglog(omega_2_sp,Saa_2_sp(1:N/2)); 
            hold on
            loglog(omega_2_sp, Sxx_filt_sp(:,1),'color',[0 0.5 0]); 
            xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\alpha\alpha}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on

            subplot(2,2,2); 
            loglog(omega,Stt_1_sp);
            hold on
            loglog(omega_2_sp,Stt_2_sp(1:N/2)); 
            hold on
            loglog(omega_2_sp, Sxx_filt_sp(:,2),'color',[0 0.5 0]); 
            xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\theta\theta}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on 

            subplot(2,2,3); 
            loglog(omega,Sqq_1_sp);
            hold on
            loglog(omega_2_sp,Sqq_2_sp(1:N/2));
            hold on
            loglog(omega_2_sp, Sxx_filt_sp(:,3), 'color',[0 0.5 0]); 
            xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{\frac{qc}{V}\frac{qc}{V}}$ [$\frac{rad^2}{rad/s}$]','interpreter','latex');
            grid on

            subplot(2,2,4); 
            loglog(omega,Snn_1_sp);
            hold on
            loglog(omega_2_sp,Snn_2_sp(1:N/2)); 
            hold on
            loglog(omega_2_sp, Sxx_filt_sp(:,4), 'color',[0 0.5 0]); 
            xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('$S_{n_z n_z}$ [$\frac{1}{rad/s}$]','interpreter','latex');
            grid on

        end

