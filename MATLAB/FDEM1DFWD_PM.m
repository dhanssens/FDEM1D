%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------------- %%
%                 FORWARD RESPONSE (PROPAGATION MATRIX)                   %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [FWD_IP,FWD_QP] = FDEM1DFWD_PM(S,M)
%
%  Use:
%  Calculates the forward response (ppm) of a given layered half-space and
%  loop-loop configuration. Calculation of the Hankel transform makes use 
%  of a digital filtering (Guptasarma and Singh filter). Typical 
%  characteristics of the half-space are stored in the Model structure (M) 
%  while the sensor characteristics are stored in the Sensor structure (S).
%
%  Input:
%  S (structure)        Sensor characteristics
%  S.x=                 x-coordinate receiver (m)
%  S.y=                 y-coordinate receiver (m)
%  S.z=                 z-coordinate receiver (m)
%  S.height=            Height of transmitter (m)
%  S.freq=              Frequency (Hz)
%  S.mom=               Transmitter moment (A.m^2)
%  S.r=                 Coil spacing (m)
%  S.omega=             Angular frequency (Rad/s)
%  S.ori=               Coil orientation (2 letter combination of X, Y, and Z)
%
%  M (structure)        Model characteristics
%  M.sus=               Susceptibility of layer(s) (SI unit)
%  M.con=               Conductivity of layer(s) (S/m)
%  M.perm=              Permittivity of layer(s) (F/m)
%  M.thick=             Layer(s) thickness (m)
%
%  Output:
%  FWD_IP              	IP response (ppm)
%  FWD_QP              	QP response (ppm)
%
%  Based on:      
%  Ward, S. H., and G. W. Hohmann, 1987, Electromagnetic theory for 
%  geophysical applications, in M. N. Nabighian, ed., Electromagnetic 
%  methods in applied geophysics, vol. 1, Theory: SEG Investigations in 
%  Geophysics 3, 131–311.
%
%  Guptasarma, D., and B. Singh, 1997, New digital linear filters for     
%  Hankel J0 and J1 transforms: Geophysical Prospection, 45, 745-762. 
%
%  Farquharson, C. G., and D. W., Oldenburg, 1996, Approximate
%  sensitivities for the electromagnetic inverse problem: Geophysics, 126,
%  235-252.
%
%  Minsley, B. J., 2011, A trans-dimensional Bayesian Markov chain Monte 
%  Carlo algorithm for model assessment using frequency-domain 
%  electromagnetic data: Geophysical Journal International, 187, 252-272.
%
%  Elwaseif, M., Robinson, J., Day-Lewis, F. D., Ntarlagiannis, D., Slater,
%  L. D., Lane, J. W., Minsley, B. J., and G. Schultz, 2017, A Matlab-
%  based frequency-domain electromagnetic inversion code (FEMIC) with
%  graphical user interface: Computers and Geosciences, 99, 61-71.
%
%  Created by Daan Hanssens
%  Ghent University, Belgium
%  January, 2017
%
%  Cite:
%  Hanssens, D., Delefortrie, S., De Pue, J., Van Meirvenne, M. and P.
%  De Smedt, Practical aspects of frequency domain electromagnetic forward 
%  and sensitivity modelling of a magnetic dipole in a multi-layered 
%  half-space: Submitted to Geoscience and Remote Sensing Magazine
%

function [FWD_IP,FWD_QP] = FDEM1DFWD_PM(S,M)

    %
    % Update parameters (optional; i.e. if not included in S structure)
    %
    
        S.r= sqrt(S.x^2 + S.y^2 + (S.z+S.height)^2);                       % Coil spacing (m)
        S.omega= 2*pi*S.freq;                                              % Angular frequency (Rad/s)

        
    %
    % Get magnetic fields
    %
    
        [H, Hn]= magneticFields(S,M);
    
        
    %
    % Normalization of magnetic field (ppm)
    %
    
        H_rsp= 1e6 .* H ./ Hn;
    
        
    %
    % Get forward response (ppm)
    %
    
        FWD_IP= real(H_rsp);
        FWD_QP= imag(H_rsp);

        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------------- %%
%                      MAGNETIC FIELD CALCULATION                         %  
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Use:
%  Calculate magnetic fields for an x-directed transmitter and different  
%  X,Y,Z-receiver orientations ('XX','XY','XZ'), Y-directed transmitter   
%  and different X,Y,Z-receiver orientations ('YX','YY','YZ') and Z-      
%  directed transmitter and different x,y,z-receiver orientations ('ZX',  
%  'ZY','ZZ').
%
%  Input:
%  S (structure)           Sensor characteristics
%  M (structure)           Model characteristics
%                                                                         
%  Based on:                                                              
%  Ward, S. H., and G. W. Hohmann, 1987, Electromagnetic theory for       
%  geophysical applications, in M. N. Nabighian, ed., Electromagnetic     
%  methods in applied geophysics, vol. 1, Theory: SEG Investigations in   
%  Geophysics 3, 131–311.  
%
%  Minsley, B. J., 2011, A trans-dimensional Bayesian Markov chain Monte 
%  Carlo algorithm for model assessment using frequency-domain 
%  electromagnetic data: Geophysical Journal International, 187, 252-272.
%
%  Created by Daan Hanssens
%  Ghent University, Belgium
%  January, 2017
%

function [H, Hn] = magneticFields(S,M)
     
    %
    % Variables
    %
    
        mu0=    4*pi*10^-7;                                                % Vacuum permeability (N/A²)
        eps0=   1/(35950207149.4727056*pi);                                % Vacuum permittivity (F/m)
    
        
    %
    % Get orientation (S.ori)
    %
    
        if strcmpi(S.ori,'ZZ')

            %
            % Calculate Hzz (HCP)
            %
            
                H= S.mom/(4*pi)*digitalFilter('rTEp2',0,S,M,mu0,eps0);
                Hn= S.mom/(4*pi)*digitalFilter('rTE02',0,S,M,mu0,eps0);

                
        elseif strcmpi(S.ori,'ZY')

            %
            % Calculate Hzy (NULL)
            %
            
                H= -S.mom/(4*pi)*S.y/S.r*...
                    digitalFilter('rTEn2',1,S,M,mu0,eps0);
                Hn= S.mom/(4*pi)*digitalFilter('rTE02',0,S,M,mu0,eps0);

                
        elseif strcmpi(S.ori,'ZX')

            %
            % Calculate Hzx (PRP)
            %
            
                H= -S.mom/(4*pi)*S.x/S.r*...
                    digitalFilter('rTEn2',1,S,M,mu0,eps0);
                Hn= S.mom/(4*pi)*digitalFilter('rTE02',0,S,M,mu0,eps0);

                
        elseif strcmpi(S.ori,'XX')

            %
            % Calculate Hxx (VCA)
            %
            
                H= -S.mom/(4*pi)*(1/S.r-2*S.x^2/S.r^3)*...
                    digitalFilter('rTEn1',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.x^2/S.r^2*...
                    digitalFilter('rTEn2',0,S,M,mu0,eps0);
                Hn= -S.mom/(4*pi)*(1/S.r-2*S.x^2/S.r^3)*...
                    digitalFilter('rTE01',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.x^2/S.r^2*...
                    digitalFilter('rTE02',0,S,M,mu0,eps0);

                
        elseif strcmpi(S.ori,'XY')

            %
            % Calculate Hxy (NULL)
            %
            
                H= S.mom/(2*pi)*S.x*S.y/S.r^3*...
                    digitalFilter('rTEn1',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.x*S.y/S.r^2*...
                    digitalFilter('rTEn2',0,S,M,mu0,eps0);
                Hn= -S.mom/(4*pi)*(1/S.r-2*S.x^2/S.r^3)*...
                    digitalFilter('rTE01',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.x^2/S.r^2*...
                    digitalFilter('rTE02',0,S,M,mu0,eps0);

                
        elseif strcmpi(S.ori,'XZ')

            %
            % Calculate Hxz (PRP)
            %
            
                H= S.mom/(4*pi)*S.x/S.r*...
                    digitalFilter('rTEp2',1,S,M,mu0,eps0);
                Hn= -S.mom/(4*pi)*(1/S.r-2*S.x^2/S.r^3)*...
                    digitalFilter('rTE01',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.x^2/S.r^2*...
                    digitalFilter('rTE02',0,S,M,mu0,eps0);

                
        elseif strcmpi(S.ori,'YX')

            %
            % Calculate Hyx (NULL)
            %
            
                H= S.mom/(2*pi)*S.x*S.y/S.r^3*...
                    digitalFilter('rTEn1',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.x*S.y/S.r^2*...
                    digitalFilter('rTEn2',0,S,M,mu0,eps0);
                Hn= -S.mom/(4*pi)*(1/S.r-2*S.y^2/S.r^3)*...
                    digitalFilter('rTE01',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.y^2/S.r^2*...
                    digitalFilter('rTE02',0,S,M,mu0,eps0);

                
        elseif strcmpi(S.ori,'YY')

            %
            % Calculate Hyy (VCP)
            %
            
                H= -S.mom/(4*pi)*(1/S.r-2*S.y^2/S.r^3)*...
                    digitalFilter('rTEn1',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.y^2/S.r^2*...
                    digitalFilter('rTEn2',0,S,M,mu0,eps0);
                Hn= -S.mom/(4*pi)*(1/S.r-2*S.y^2/S.r^3)*...
                    digitalFilter('rTE01',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.y^2/S.r^2*...
                    digitalFilter('rTE02',0,S,M,mu0,eps0);

                
        elseif strcmpi(S.ori,'YZ')

            %
            % Calculate Hyz (NULL)
            %
            
                H= S.mom/(4*pi)*S.y/S.r*...
                    digitalFilter('rTEp2',1,S,M,mu0,eps0);
                Hn= -S.mom/(4*pi)*(1/S.r-2*S.y^2/S.r^3)*...
                    digitalFilter('rTE01',1,S,M,mu0,eps0) - ...
                    S.mom/(4*pi)*S.y^2/S.r^2*...
                    digitalFilter('rTE02',0,S,M,mu0,eps0);

                
        else

            %
            % Error message if T-R orientation is unknown
            %
            
                error('Transmitter-receiver orientation unknown.');

                
        end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------------- %%
%                             DIGITAL FILTER                              %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Use:
%  Solves the Hankel Transform of the zeroth or first order by using a    
%  Guptasarma filtering routine. 
%
%  Input:
%  S (structure)           Sensor characteristics
%  M (structure)           Model characteristics
%  funname                 Sampled function
%  order                   Order of Bessel function
%                                                                         
%  Based on:                                                              
%  Guptasarma, D., and B. Singh, 1997, New digital linear filters for     
%  Hankel J0 and J1 transforms: Geophysical Prospection, 45, 745-762. 
%
%  Minsley, B. J., 2011, A trans-dimensional Bayesian Markov chain Monte 
%  Carlo algorithm for model assessment using frequency-domain 
%  electromagnetic data: Geophysical Journal International, 187, 252-272.
%
%  Elwaseif, M., Robinson, J., Day-Lewis, F. D., Ntarlagiannis, D., Slater,
%  L. D., Lane, J. W., Minsley, B. J., and G. Schultz, 2017, A Matlab-
%  based frequency-domain electromagnetic inversion code (FEMIC) with
%  graphical user interface: Computers and Geosciences, 99, 61-71.
%
%  Created by Daan Hanssens
%  Ghent University, Belgium
%  January, 2017
%

function y = digitalFilter(funname,order,S,M,mu0,eps0)
   
    %
    % Load Guptasarma and Singh filter
    %
    
        if order == 0;                                                     % Load 120-point filter
            filter.a = -8.3885;
            filter.s = 0.090422646867;
            filter.w = [9.62801364263000e-07 -5.02069203805000e-06 1.25268783953000e-05 -1.99324417376000e-05 2.29149033546000e-05 -2.04737583809000e-05 1.49952002937000e-05 -9.37502840980000e-06 5.20156955323000e-06 -2.62939890538000e-06 1.26550848081000e-06 -5.73156151923000e-07 2.76281274155000e-07 -1.09963734387000e-07 7.38038330280000e-08 -9.31614600001000e-09 3.87247135578000e-08 2.10303178461000e-08 4.10556513877000e-08 4.13077946246000e-08 5.68828741789000e-08 6.59543638130000e-08 8.40811858728000e-08 1.01532550003000e-07 1.26437360082000e-07 1.54733678097000e-07 1.91218582499000e-07 2.35008851918000e-07 2.89750329490000e-07 3.56550504341000e-07 4.39299297826000e-07 5.40794544880000e-07 6.66136379541000e-07 8.20175040653000e-07 1.01015545059000e-06 1.24384500153000e-06 1.53187399787000e-06 1.88633707689000e-06 2.32307100992000e-06 2.86067883258000e-06 3.52293208580000e-06 4.33827546442000e-06 5.34253613351000e-06 6.57906223200000e-06 8.10198829111000e-06 9.97723263578000e-06 1.22867312381000e-05 1.51305855976000e-05 1.86329431672000e-05 2.29456891669000e-05 2.82570465155000e-05 3.47973610445000e-05 4.28521099371000e-05 5.27705217882000e-05 6.49856943660000e-05 8.00269662180000e-05 9.85515408752000e-05 0.000121361571831000 0.000149454562334000 0.000184045784500000 0.000226649641428000 0.000279106748890000 0.000343716968725000 0.000423267056591000 0.000521251001943000 0.000641886194381000 0.000790483105615000 0.000973420647376000 0.00119877439042000 0.00147618560844000 0.00181794224454000 0.00223860214971000 0.00275687537633000 0.00339471308297000 0.00418062141752000 0.00514762977308000 0.00633918155348000 0.00780480111772000 0.00961064602702000 0.0118304971234000 0.0145647517743000 0.0179219149417000 0.0220527911163000 0.0271124775541000 0.0333214363101000 0.0408864842127000 0.0501074356716000 0.0612084049407000 0.0745146949048000 0.0900780900611000 0.107940155413000 0.127267746478000 0.146676027814000 0.162254276550000 0.168045766353000 0.152383204788000 0.101214136498000 -0.00244389126667000 -0.154078468398000 -0.303214415655000 -0.297674373379000 0.00793541259524000 0.426273267393000 0.100032384844000 -0.494117404043000 0.392604878741000 -0.190111691178000 0.0743654896362000 -0.0278508428343000 0.0109992061155000 -0.00469798719697000 0.00212587632706000 -0.000981986734159000 0.000444992546836000 -0.000189983519162000 7.31024164292000e-05 -2.40057837293000e-05 6.23096824846000e-06 -1.12363896552000e-06 1.04470606055000e-07];
        
        elseif order == 1;                                                 % Load 140-point filter
            filter.a = -7.91001919;
            filter.s = 0.087967143957;
            filter.w = [-6.76671159511000e-14 3.39808396836000e-13 -7.43411889153000e-13 8.93613024469000e-13 -5.47341591896000e-13 -5.84920181906000e-14 5.20780672883000e-13 -6.92656254606000e-13 6.88908045074000e-13 -6.39910528298000e-13 5.82098912530000e-13 -4.84912700478000e-13 3.54684337858000e-13 -2.10855291368000e-13 1.00452749275000e-13 5.58449957721000e-15 -5.67206735175000e-14 1.09107856853000e-13 -6.04067500756000e-14 8.84512134731000e-14 2.22321981827000e-14 8.38072239207000e-14 1.23647835900000e-13 1.44351787234000e-13 2.94276480713000e-13 3.39965995918000e-13 6.17024672340000e-13 8.25310217692000e-13 1.32560792613000e-12 1.90949961267000e-12 2.93458179767000e-12 4.33454210095000e-12 6.55863288798000e-12 9.78324910827000e-12 1.47126365223000e-11 2.20240108708000e-11 3.30577485691000e-11 4.95377381480000e-11 7.43047574433000e-11 1.11400535181000e-10 1.67052734516000e-10 2.50470107577000e-10 3.75597211630000e-10 5.63165204681000e-10 8.44458166896000e-10 1.26621795331000e-09 1.89866561359000e-09 2.84693620927000e-09 4.26886170263000e-09 6.40104325574000e-09 9.59798498616000e-09 1.43918931885000e-08 2.15798696769000e-08 3.23584600810000e-08 4.85195105813000e-08 7.27538583183000e-08 1.09090191748000e-07 1.63577866557000e-07 2.45275193920000e-07 3.67784458730000e-07 5.51470341585000e-07 8.26916206192000e-07 1.23991037294000e-06 1.85921554669000e-06 2.78777669034000e-06 4.18019870272000e-06 6.26794044911000e-06 9.39858833064000e-06 1.40925408889000e-05 2.11312291505000e-05 3.16846342900000e-05 4.75093313246000e-05 7.12354794719000e-05 0.000106810848460000 0.000160146590551000 0.000240110903628000 0.000359981158972000 0.000539658308918000 0.000808925141201000 0.00121234066243000 0.00181650387595000 0.00272068483151000 0.00407274689463000 0.00609135552241000 0.00909940027636000 0.0135660714813000 0.0201692550906000 0.0298534800308000 0.0439060697220000 0.0639211368217000 0.0916763946228000 0.128368795114000 0.173241920046000 0.219830379079000 0.251193131178000 0.232380049895000 0.117121080205000 -0.117252913088000 -0.352148528535000 -0.271162871370000 0.291134747110000 0.317192840623000 -0.493075681595000 0.311223091821000 -0.136044122543000 0.0512141261934000 -0.0190806300761000 0.00757044398633000 -0.00325432753751000 0.00149774676371000 -0.000724569558272000 0.000362792644965000 -0.000185907973641000 9.67201396593000e-05 -5.07744171678000e-05 2.67510121456000e-05 -1.40667136728000e-05 7.33363699547000e-06 -3.75638767050000e-06 1.86344211280000e-06 -8.71623576811000e-07 3.61028200288000e-07 -1.05847108097000e-07 -1.51569361490000e-08 6.67633241420000e-08 -8.33741579804000e-08 8.31065906136000e-08 -7.53457009758000e-08 6.48057680299000e-08 -5.37558016587000e-08 4.32436265303000e-08 -3.37262648712000e-08 2.53558687098000e-08 -1.81287021528000e-08 1.20228328586000e-08 -7.10898040664000e-09 3.53667004588000e-09 -1.36030600198000e-09 3.52544249042000e-10 -4.53719284366000e-11];        
        end
    
        
    %
    % Get lambda values
    %
    
        flen = numel(filter.w);
        ind = 1:flen;
        l1 = 1./S.r;
        l2 = 10.^(filter.a+(ind-1).*filter.s);
        lambda = l1' .* l2';
    
        
    %
    % Evaluate function at lambda
    %
    
        if strcmp(funname,'rTE01')
            YF = PM01(lambda,S,mu0,eps0);
        elseif strcmp(funname,'rTE02')
            YF = PM02(lambda,S,mu0,eps0);
        elseif strcmp(funname,'rTEn1')
            YF = PMn1(lambda,S,M,mu0,eps0);
        elseif strcmp(funname,'rTEn2')
            YF = PMn2(lambda,S,M,mu0,eps0); 
        elseif strcmp(funname,'rTEp2')
            YF = PMp2(lambda,S,M,mu0,eps0);
        end

        
    %
    % Calculate output, considering weights and r
    %
    
        y = dot(YF',filter.w') ./ S.r;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------------- %%
%                           PROPAGATION MATRIX                            %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Use:
%  Calculates the P(2,1)/P(1,1) ratio of the propagation matrix for a given 
%  layered half-space and lambda value.
%
%  Input:
%  S (structure)           Sensor characteristics
%  M (structure)           Model characteristics
%  lambda                  Hankel transform parameter
%
%  Based on:
%  Farquharson, C. G., and D. W., Oldenburg, 1996, Approximate
%  sensitivities for the electromagnetic inverse problem: Geophysics, 126,
%  235-252.
%
%  Created by Daan Hanssens
%  Ghent University, Belgium
%  January, 2017
%

function PP = propagationMatrix(lambda,S,M,mu0,eps0)

    %
    % Calculate mu
    %
    
        mu= mu0.*(1+M.sus);
    
        
    %
    % Calculate u0 and u
    %
    
        lambda= reshape(lambda,1,1,[]);                                    % Reshape lambda
        u0= sqrt(lambda.^2 - S.omega.^2.*mu0.*eps0);
        u= sqrt(repmat(lambda.^2,1,length(M.con),1) - ...
                repmat(S.omega.^2.*mu.*M.perm,1,1,numel(lambda)) + ...
                repmat(1i.*S.omega.*mu.*M.con,1,1,numel(lambda)));
    
            
    %
    % Calculate m(1)
    %
    
        PM =  [1./2.*(1+(mu0.*u(1,1,:))./(mu(1,1,:).*u0)) ...
                    1./2.*(1-(mu0.*u(1,1,:))./(mu(1,1,:).*u0));
                    1./2.*(1-(mu0.*u(1,1,:))./(mu(1,1,:).*u0)) ...
                    1./2.*(1+(mu0.*u(1,1,:))./(mu(1,1,:).*u0))];
  

    for i=2:length(M.con)

        %
        % Calculate m(i)
        %
        
            m =  [1./2.*(1+(mu(1,i-1,:).*u(1,i,:))./(mu(1,i,:).*u(1,i-1,:))) ...
                        1./2.*(1-(mu(1,i-1,:).*u(1,i,:))./(mu(1,i,:).*u(1,i-1,:)));
                        1./2.*(1-(mu(1,i-1,:).*u(1,i,:))./(mu(1,i,:).*u(1,i-1,:))).*exp(-2.*u(1,i-1,:).*M.thick(1,i-1,:)) ...
                        1./2.*(1+(mu(1,i-1,:).*u(1,i,:))./(mu(1,i,:).*u(1,i-1,:))).*exp(-2.*u(1,i-1,:).*M.thick(1,i-1,:))];

                    
        %
        % Update propagation matrix PM
        %
        
            for ii=1:numel(lambda)
                PM(:,:,ii) = PM(:,:,ii)*m(:,:,ii);
            end
             
    end
     
    
    %
    % Calculate P(2,1)/P(1,1) ratio
    %
        PP= reshape(PM(2,1,:)./PM(1,1,:),[],1,1);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------------------------------------------------------- %
%                            LAMBDA FUNCTIONS                             %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Use:
%  Additional lambda functions for Hankel calculation of Primary and 
%  Secondary magnetic fields.
%
%  Input:
%  S (structure)           Sensor characteristics
%  M (structure)           Model characteristics
%  lambda                  Hankel transform parameter
%
%  Based on:
%  Ward, S. H., and G. W. Hohmann, 1987, Electromagnetic theory for       
%  geophysical applications, in M. N. Nabighian, ed., Electromagnetic     
%  methods in applied geophysics, vol. 1, Theory: SEG Investigations in   
%  Geophysics 3, 131–311.  
%
%  Guptasarma, D., and B. Singh, 1997, New digital linear filters for     
%  Hankel J0 and J1 transforms: Geophysical Prospection, 45, 745-762.  
%
%  Created by Daan Hanssens
%  Ghent University, Belgium
%  January, 2017
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------- Function 1/5 ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = PM01(lambda,S,mu0,eps0)
       
    %
    % Define variable u0
    %
    
        u0= sqrt(lambda.^2-S.omega.^2.*mu0.*eps0);
    
        
    %
    % Output
    %
    
        y= (exp(-u0.*(S.z+S.height))).*lambda;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------- Function 2/5 ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = PM02(lambda,S,mu0,eps0)
    
    %
    % Define variable u0
    %
    
        u0= sqrt(lambda.^2-S.omega.^2.*mu0.*eps0);
    
        
    %
    % Output
    %
    
        y= (exp(-u0.*(S.z+S.height))).*lambda.^3./u0;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------- Function 3/5 ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = PMn1(lambda,S,M,mu0,eps0)
       
    %
    % Calculate PM
    %
    
        PM= propagationMatrix(lambda,S,M,mu0,eps0);
    
        
    %
    % Define variable u0
    %
    
        u0= sqrt(lambda.^2-S.omega.^2.*mu0.*eps0);
    
        
    %
    % Output
    %
    
        y= (-PM.*exp(u0.*(S.z-S.height))).*lambda;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------- Function 4/5 ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = PMn2(lambda,S,M,mu0,eps0)
    
    %
    % Calculate PM
    %
    
        PM= propagationMatrix(lambda,S,M,mu0,eps0);
    
        
    %
    % Define variable u0
    %
    
        u0= sqrt(lambda.^2-S.omega.^2.*mu0.*eps0);
    
        
    %
    % Output
    %
    
        y= (-PM.*exp(u0.*(S.z-S.height))).*lambda.^2;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------- Function 5/5 ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = PMp2(lambda,S,M,mu0,eps0)
       
    %
    % Calculate PM
    %
    
        PM= propagationMatrix(lambda,S,M,mu0,eps0);
    
        
    %
    % Define variable u0
    %
    
        u0= sqrt(lambda.^2-S.omega.^2.*mu0.*eps0);
    
    
    %
    % Output
    %
    
        y= (PM.*exp(u0.*(S.z-S.height))).*lambda.^3./u0;
    
end
