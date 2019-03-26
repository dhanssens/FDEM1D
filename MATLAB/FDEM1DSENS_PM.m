%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------------- %%
%              SENSITIVITY DISTRIBUTION (PROPAGATION MATRIX)              %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [SENS_IP,SENS_QP] = FDEM1DSENS_PM(S,M,par)
%
%  Use:
%  Calculates the sensitivity distribution of a given layered soil medium 
%  and loop-loop configuration towards a certain physical property using 
%  the brute-force or perturbation method. Typical characteristics of the 
%  soil medium are stored in the Model structure (M) while the sensor 
%  characteristics are stored in the Sensor structure (S).
%
%  Input:
%  S (structure)           Sensor characteristics
%  M (structure)           Model characteristics
%  par                     Sensitivity parameter ('con','sus','perm')
%
%  Output:
%  SENS_IP                 IP sensitivity
%  SENS_QP                 QP sensitivity
%  Error (structure)       Estimated max. error (.IP and .QP)
%
%  Created by Daan Hanssens
%  UGent, Belgium
%  january 2017
%
%  Cite:
%  Hanssens, D., Delefortrie, S., De Pue, J., Van Meirvenne, M., 
%  and P. De Smedt. Frequency-Domain Electromagnetic Forward and 
%  Sensitivity Modeling: Practical Aspects of modeling a Magnetic Dipole 
%  in a Multilayered Half-Space. IEEE Geoscience and Remote Sensing 
%  Magazine, 7(1), 74-85
%

function [SENS_IP,SENS_QP,Error] = FDEM1DSENS_PM(S,M,par)
        
    %
    % Store original profile
    %
    
        op= M.(par);             
       
        
    %
    % Calculate partial derivatives
    %
    
        % Loop Model layers
        for i= 1:numel(M.(par))

            %
            % Get original response
            %
            
                if i==1; [FWD_IP_ori,FWD_QP_ori] = FDEM1DFWD_PM(S,M); end;

                
            %
            % Get altered response (forward)
            %
            
                M.(par)= op;                                               % Get original profile
                pert(i)= op(i) .* 0.01;                                    % Get relative perturbation (1%)
                M.(par)(i)= op(i) + pert(i);
                [FWD_IP_alt_p(i),FWD_QP_alt_p(i)] = FDEM1DFWD_PM(S,M);


            %
            % Get altered response (backward)
            %
            
                M.(par)(i)= op(i) - pert(i);
                [FWD_IP_alt_n(i),FWD_QP_alt_n(i)] = FDEM1DFWD_PM(S,M);

                
        end
        
    %
    % First derivative (Output)
    %
            
        SENS_QP= (FWD_QP_alt_p - FWD_QP_ori) ./ pert;
       	SENS_IP= (FWD_IP_alt_p - FWD_IP_ori) ./ pert;
          
                
    %
    % Second derivative
    %
            
    	SENS_QP_pert_sd= (FWD_QP_alt_p - 2.*FWD_QP_ori + ...
                                FWD_QP_alt_n) ./ pert .^2;
       	SENS_IP_pert_sd= (FWD_IP_alt_p - 2.*FWD_IP_ori + ...
                                FWD_IP_alt_n) ./ pert .^2;
    
    %
    % Estimate maximum error (Output)
    %

        Error.QP= max(SENS_QP_pert_sd.* pert/2);
        Error.IP= max(SENS_IP_pert_sd.* pert/2);
            
          
end
