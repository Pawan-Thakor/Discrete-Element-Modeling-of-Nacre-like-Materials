%% Defining Properties for Shear Stress-Strain Curve of each spring interface

G=100*10^6;                    %Shear Modulus
H_1=3.5*10^6;                    %Hardening Value
tau_y=7*10^6;                %Yeild Shear Stress value
gamma_u=1.5;                     %Max Allowed Plastic Strain 
gamma_p_s=1;                 %Plastic Strain at start of Softening region

tau_s = tau_y + H_1*gamma_p_s;               %Max Shear Stress Point
G_1 = H_1*G/(G+H_1);                         %Slope of Hardening in Shear stress-strain curve 
G_2 = -tau_s/(gamma_u-(gamma_p_s+tau_s/G));  %Slope of Softening in Shear stress-strain curve
H_2 = G*G_2/(G-G_2);                         %Softening Value
gamma_p_s_max = gamma_p_s - tau_s/H_2;