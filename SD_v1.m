%Givens
    
    %specific heat of tissue
        c_tis = 3.75; %kJ/kg
    %density of tissue
        rho_tis = 1000; %kg/m^3
    %volume of tissue
        V_tis = .000262; %m^3
    %change in temperature
        T_delta = 23.9; %Celsius
    %hydraulic diameter of pad
        hydd = .017; %meters assuming 2/3 x 2/3 inches
    %dynamic viscosity of water
        visc_w = 1.5.*(10.^(-6)); %Pa*s
    %volumetric flow of pump
        flow_p = .00042; %m^3/s
    %Prandtl number assuming 35 deg F
        Pr = 13;
    %length of pad
        L = .61; %meters
    %viscosity ratio
        mu = 1;
    %thermal conductivity of water
        k_w = .58; %W/m.*K
    %thermal conductivity of flesh (mix of fat and muscle)
        k_f = .375; 
    %desired 55 deg F distance to heat source depth
        depth = .05;
    %contact area
        A = .0103; %m^2
    %difference between water temp and skin temp
        T_d = 37;
    
    
    
%Calculations
    %skin contact
        %initial energy
            q_c = c_tis.*rho_tis.*V_tis.*T_delta;
        %Reynold's number
            %flow velocity through pad
            v_flow = flow_p./(pi.*(hydd.^2)./4);
            Re = v_flow.*hydd./visc_w;
        %Seider_Tate Correlation
            Nu = 0.023.*(Re.^0.8).*(Pr.^(1./3));

        %for laminar flow Re < 10000
            %Graetz number
           %     Gr = hydd.*Re.*Pr./L;
        %Nusselt number
          %  Nu = (3.66 + (.085.*Gr)./(1 + (.047.*Gr).^(2/3))).*(mu.^.14);
        %heat transfer coefficient
            h = Nu.*k_w./hydd;
        %q_c rate
            q_cdot = h.*A.*T_d;
        
    
    %inside tissue
        %energy
            q = k_f.*T_delta./depth;
