%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------------- %%
%                              EXAMPLE FILE                               %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Use:
%  Calculates the forward response and sensitivity distribution for a given 
%  layered half-space and loop-loop configuration. Typical characteristics 
%  of the half-space are stored in the Model structure (M) while the sensor 
%  characteristics are stored in the Sensor structure (S). This file is
%  used as example to demonstrate the fieldnames within the structures. 
%
%  Created by Daan Hanssens
%  Ghent University, Belgium
%  December, 2016
%
%  Cite:
%  Hanssens, D., Delefortrie, S., De Pue, J., Van Meirvenne, M. and P.
%  De Smedt, Practical aspects of frequency domain electromagnetic forward 
%  and sensitivity modelling of a magnetic dipole in a multi-layered 
%  half-space: Submitted to Geoscience and Remote Sensing Magazine
%

   clc; clear all; close all;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------- User-input -------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %
    % Sensor characteristics (S structure)
    %

        S.x=        4;                                                     % x-coordinate receiver (m)
        S.y=        0;                                                     % y-coordinate receiver (m)
        S.z=        -0.16;                                                 % z-coordinate receiver (m) - positive z-axis pointed down
        S.height=   0.16;                                                  % Height of transmitter (m)
        S.freq=     9000;                                                  % Frequency (Hz)
        S.mom=      1;                                                     % Transmitter moment (A.m^2)
        S.ori=      'ZZ';                                                  % Coil orientation (Transmitter(X,Y,Z),Receiver(X,Y,Z))


    %
    % Model characteristics (M structure)
    %

        M.sus=      linspace(100,100,100) .* 1e-5;                         % Susceptibility of layer(s) (-)
        M.con=      linspace(0.001,0.001,100);                             % Conductivity of layer(s) (S/m)
        M.perm=     linspace(1e-12,1e-12,100);                             % Permittivity of layer(s) (F/m)
        M.thick=    logspace(log10(0.1),log10(0.5),100);                   % Layer(s) thickness (m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------- Modelling --------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %
    % Calculate forward response (ppm)
    %

        [FWD_IP,FWD_QP] = FDEM1DFWD_RC(S,M);                               % Reflection coefficient approach
%        [FWD_IP,FWD_QP] = FDEM1DFWD_PM(S,M);                              % Propagation matrix approach


    %
    % Calculate sensitivity response (ppm)
    %
    
        par= 'con';  % ('con','sus','perm')
        [SENS_IP,SENS_QP,Err] = FDEM1DSENS_RC(S,M,par);                    % Reflection coefficient approach
%         [SENS_IP,SENS_QP,Err] = FDEM1DSENS_PM(S,M,par);                  % Propagation matrix approach
        
        
        %
        % Plot sensitivity
        % (Optionally: exclude basement layer due to increased sensitivity 
        % of infinite basement layer)
        %

            figure(); 
            subplot(2,1,1);
            plot(cumsum(M.thick(1:end-1)),SENS_QP(1:end-1),'k.'); 
            legend(['Error < ',num2str(abs(Err.QP))]); legend boxoff;
            ylabel(['Sensitivity: QP, ',par]); 
            subplot(2,1,2);
            plot(cumsum(M.thick(1:end-1)),SENS_IP(1:end-1),'k.'); 
            legend(['Error < ',num2str(abs(Err.IP))]); legend boxoff;
            xlabel('Depth (m)'); ylabel(['Sensitivity: IP, ',par]);
    
    
