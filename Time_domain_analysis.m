%In this file, the time-domain analysis for both the complete and the
%reduced aircraft model is performed
% Based on the examp71.m file

load complete_system
load sp_system

% -----------------------------------------------------------------------
%% ---------  Time-domain analysis for the complete model ---------------
%------------------------------------------------------------------------

        % TIME AXIS INPUT VECTOR DEFINITION
        dt = 0.01;
        T  = 60; 
        t = 0:dt:T;
        N = length(t);

        % INPUT VECTOR DEFINITION
        nn = zeros(1,N);             % zero input elevator
        u_g = randn(1,N)/sqrt(dt);    % scaled input long. turbulence,
                                     % note the sqrt(dt) because of lsim
        w_g = randn(1,N)/sqrt(dt);    % scaled input vert. turbulence,
                                     % note the sqrt(dt) because of lsim
        u1  = [nn' u_g' nn'];          % input vector definition (combined
                                     % longitudinal and vertical turbulence).
        u2  = [nn' nn' w_g'];

        % SIMULATION OF MOTION VARIABLES
        y1 = lsim(A_complete_system,B_complete_system,C_complete_system,D_complete_system,u1,t);
        y2 = lsim(A_complete_system,B_complete_system,C_complete_system,D_complete_system,u2,t);
        
        y= y1 + y2;
        
        figure()
        % PLOTTING RESULTS
        subplot(5,1,1);
        plot(t,y(:,1))
        xlabel('time [s]'); ylabel('$\hat{u}$ [-]', 'FontSize',13,'interpreter','latex'); title('Airspeed Deviation', 'FontSize',12);

        subplot(5,1,2);
        plot(t,y(:,2)*180/pi)
        xlabel('time [s]'); ylabel('$\alpha$ [deg]', 'FontSize',13, 'interpreter','latex'); title('Angle of Attack', 'FontSize',12);

        subplot(5,1,3);
        plot(t,y(:,3)*180/pi)
        xlabel('time [s]'); ylabel('$\theta$ [deg]', 'FontSize',13, 'interpreter','latex'); title('Pitch Angle', 'FontSize',12);

        subplot(5,1,4);
        plot(t,y(:,4)*180/pi)
        xlabel('time [s]'); ylabel('$\frac{qc}{V}$ [deg]', 'FontSize',13, 'interpreter','latex'); title('Pitch Rate', 'FontSize',12);

        subplot(5,1,5);
        plot(t,y(:,5))
        xlabel('time [s]'); ylabel('$n_z$ [-]', 'FontSize',13,'interpreter','latex'); title('Load Factor', 'FontSize',12);
        
        hold on
        
        sgtitle('Complete model','FontSize',15)
        
        

% -----------------------------------------------------------------------
%% ---------  Time-domain analysis for the reduced model ---------------
%------------------------------------------------------------------------

        % TIME AXIS INPUT VECTOR DEFINITION
        dt_sp = 0.01;
        T_sp  = 60; 
        t_sp = 0:dt_sp:T_sp;
        N_sp = length(t_sp);

        % INPUT VECTOR DEFINITION
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

        figure()
        % PLOTTING RESULTS
        subplot(4,1,1);
        plot(t_sp,y_sp(:,1)*180/pi)
        xlabel('time [s]'); ylabel('$\alpha$ [deg]', 'FontSize',13,'interpreter','latex'); title('Angle of Attack', 'FontSize',12);

        subplot(4,1,2);
        plot(t_sp,y_sp(:,2)*180/pi)
        xlabel('time [s]'); ylabel('$\theta$ [deg]', 'FontSize',13, 'interpreter','latex'); title('Pitch Angle', 'FontSize',12);

        subplot(4,1,3);
        plot(t_sp,y_sp(:,3)*180/pi)
        xlabel('time [s]'); ylabel('$\frac{qc}{V}$ [deg]','FontSize',13, 'interpreter','latex'); title('Pitch Rate', 'FontSize',12);

        subplot(4,1,4);
        plot(t_sp,y_sp(:,4))
        xlabel('time [s]'); ylabel('$n_z$ [-]', 'FontSize',13, 'interpreter','latex'); title('Load Factor', 'FontSize',12);
        
        hold on
        
        sgtitle('Reduced model','FontSize',15);
