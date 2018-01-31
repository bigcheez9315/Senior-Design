%Pad heat transfer analysis
  %Givens
    
    %specific heat of tissue
        c_tis = 3.75; %kJ/kg
    %density of tissue
        rho_tis = 1000; %kg/m^3
    %volume of tissue
        V_tis = .000262; %m^3
    %change in temperature between center of body and final cooling temp
        T_delta = 23.9; %Celsius
    %hydraulic diameter of pad
        hydd = .017; % meters assuming 2/3 x 2/3 inches
    %kinematic viscosity of water 
        visc_w = 4.116.*(10.^(-7)); % m^2*s^-1
    %volumetric flow of pump
        flow_p = .00042; % m^3/s
    %Prandtl number assuming 35 deg F
        Pr = 13;
    %length of pad (basing off universal pad on amazon ~ 32 in.
        L = .82; % meters
    %viscosity ratio
        mu = 1; 
    %thermal conductivity of water
        k_w = .58; %W/m.*K
    %thermal conductivity of flesh (mix of fat and muscle)
        k_f = .375; %W/m.*K
    %desired 55 deg F distance to heat source depth
        depth = .05; % m
    
    %contact area = Length water travels x width of pass
    % width of pass in pad
    w = .025; % looking at universal pad, width seems to be 1 in. ~ .025 m
    A = w*L; %m^2
    
    
    
    
%Calculations
    %skin contact
        %initial energy (J)
            q_c = 1000*(c_tis.*rho_tis.*V_tis.*T_delta);
        %Reynold's number
            %flow velocity through pad
            v_flow = flow_p./(pi.*(hydd.^2)./4);
            Re = v_flow.*L/(visc_w);
        %Average Nu for turbulent flow over a flat plate: https://www.sfu.ca/~mbahrami/ENSC%20388/Notes/Forced%20Convection.pdf
            Nu = 0.037.*(Re.^0.8).*(Pr.^(1./3));

        %for laminar flow Re < 10000
            %Graetz number
           %     Gr = hydd.*Re.*Pr./L;
        %Nusselt number
          %  Nu = (3.66 + (.085.*Gr)./(1 + (.047.*Gr).^(2/3))).*(mu.^.14);

        
    
    %inside tissue
        % power generated inside of tissue  (W)
            q = k_f.*T_delta./depth;
    %time (let's set this to 5 minutes
        t = 300; %s
    %Let's estimate the power needed in order to cool
    % the skin to desired temp in 5 minutes. q_h
    q_h = (q_c/t) + q;
    
    % Let's calculate Temperature Difference, T_d for 
            %heat transfer coefficient
            h = Nu.*k_w./hydd;
        
            T_d = q_h/(A*h); % Temp. difference between surface of pad
                             % and water temp
            
     % Let's calculate the temperature difference between inlet and 
     %o outlet of the pad
     
     t_pad =  L/v_flow;   %time water is in pad
     rho_water = 1000; %density of water
     C_p = 4185.5;  %Specific heat of water J/(Kg*K)
     V_water_pad = L*(pi * hydd^2)/4; % Volume of water in pad
     m = rho_water*V_water_pad; %mass of water in pad
     T_dif_pad = (q_h*t_pad)/(m*C_p);  % Temp dif. between water inlet and outlet of pad

% Waterblock Heat Transfer Analysis:            
  % We want the temperature between inlet and oulet of waterblock 
  % to be atleast T_dif_pad. Now let's start calculating it
  
  % There are 4 zones of heat transfer
  % Zone 1: Heat Transfer from the Cooling Channels to the Water Bath
  
  % This heat transfer occurs my forced convection of the water traveling 
  % through the pad. q = hA(T_dif)
  
  % We need to calculate the area and convective heat transfer coefficient
  % inside the pad
  
  %Givens:
    %  diameter of waterblock inlet--units in meters. Specs from amazon waterblock
    d_wb = .0095; % m
    % Cross-sectional area of inlet m^2    
    A_c = (pi*d_wb^2)/4;  %m^2
    % Length of single pass of water block. To get the total length water travels,
    % we need to multiply this by number of passes water takes in water block
    L_pass = .042; % m
    % I got this number from looking at an image of the inside 
    % of a waterblock and saw it had 8 passes              
    N_pass = 8; % unitless
    % Total length water travels 
    L_wb = L_pass*N_pass;   %m
    % Perimeter of cooling channel. Same as total surface area of block channels 
    % touching the water
    P = (pi*d_wb)*L_wb; % m
    % hydraulic diameter of water block
    hydd_wb = 4*A_c/P; % m
    % speed of water through water block
    v_flow_wb = flow_p/A_c; % m/s
    % Reynold's number of flow through water block
    Re_wb = (v_flow_wb*hydd_wb)/visc_w; % Re > 10^5 so flow is turblent
    
    k_copper = 401; % W/(m*K)
    A_copperplate = 0.0018 ; % in m^2 based off dimensions of amazon water block(42mm x 42mm)
    t_copper = .00065; % thickness of copper plate  in meters 
    % Use Dittus Boelter equation to get Nusselt Number for turbulent flow
    % in pipe
    Nu_wb = .023*Re_wb^.8*Pr^.3; 
    
    % calculate convective heat transfer coefficient from Nu_wb
    h_wb = (Nu_wb*k_w)/hydd_wb; 
    
    
    % q_cwb = h_wb*P*(T_s - T_in) where T_s is surface temp of water block
    % and T_in is inlet temperature of water
    
    % For a given inlet temperature we can calculate T_cold of peltier
    % cell. We can also calculate the T_hot in a future step. We can use 
    % this temperature difference together with the applied amperage to
    % estimate Q_c from the TEC correlation chart and see if this is too
    % high. 
    
    % First let's calculate T_cold:
    
    T_initial = 23; % We will use this for sample calculation of cold side.
                    % This is ambient temperature so it will be initial
                    % temp. of water before cooling
     % Thermal resistance as heat travels from water to water block through
     % copper plate 
     R_cold =  (1/(h_wb*P))+(t_copper/(k_copper*A_copperplate));          
     T_coldside = T_initial - q_h*R_cold; % see eq. 20 in reference pdf for details
     

  % Zone 3: Heat Transfer and Heat Production in the Peltier Cell
  % Zone 4: Heat Transfer through the Copper Plate and Heat Sink
  % Here our goal is to find the hot side temperature of the copper plate
  
  % First we need to find the rate at which heat is transferred 
  % from hot side to air. This heat transfer invovles conduction through
  % the copper plate and forced convection of the air traveling through the
  % heat sink
  
  % Reynold's number for air flow:
  
  % We need volumetric flow rate to get reynold's number
  % I estimate volumetric flow rate from the specs for a fan I found on
  % amazon. Range is between 15.7-54.8 CFM(ft^3/min) = .216- .91 m^3/s
  % Here is link to fan: https://www.amazon.com/dp/B00K7809O2/ref=psdc_11036281_t2_B004PLXV8S
  % I will pick .5 as average flow rate
  flow_fan = .5 ; % m^3/s
  
  % Heat sink has dimensions 90 x 140 mm = .0126 m^2
  A_heatsink = .0126; % m^2
  % Length of the heat sink is obtained from specs: 
  L_heatsink = .051; %m
  % Spacing between fins in heat sink. Taken from estimates online
  w = .0075; % m
  hydd_heatsink = 2*w; %m
  % kinematic viscosity of air
  visc_air = 1.48*(10^-5);% m^2/s
  % Get prantl number and conductivity for air from engineering toolbox
  Pr_air = .7;
  k_air = 1.4;
  %Temperature of ambient air
  T_ambient = 22;
  Re_heatsink = ((flow_fan/A_heatsink)*hydd_heatsink)/visc_air;
  % Since Reynold's number is less than 10^5 the flow is laminar
  % For laminar flow across parallel plates, we can calculate the Nu
  
  Nu_heatsink = .664*(Re_heatsink^.5)*(Pr_air)^(1/3); 
  h_heatsink = (k_air*Nu_heatsink)/L_heatsink; 
  
  % Calculate thermal resistance as heat travels from hot side of TEC to
  % ambient air. This involves conduction through copper plate and
  % convection of air through heat sink
  
  R_Hot=(1/(h_heatsink*A_heatsink))+(t_copper/(k_copper*A_copperplate));
  
  % we can calculate the rate that heat is generated inside of the TEC by
  % using the upper limit which is calculated as P=I^2*R
  % We can calculate R from product Specs. For one of the ones we ordered
  % the V_max = 15.9 V and I_Max = 6A and V=IR --> 
  % R = V_max/I_max = 2.65 ohms. 

  % let's set an arbitrary current to be 3A
  I = 3; % A
  R = 2.65; % ohms
  q_hotside = I^2*R;
  
  % Since we now Q generated from hot side, we now can solve for T_hotside
  
  T_hotside = q_hotside*R_Hot + T_ambient;
  
 % From the TEC performance chart we can extrapolate the 
 % Q_c, the cooling rate of the water to be .68*92 = 62.56 W-- too small
 
  
  
  
  
  
  
  
  
  
  
            
