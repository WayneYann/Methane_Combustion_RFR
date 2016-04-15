%Dummy variables
Cp_G = 10^3; %J/kg-K
hWG = 0;

R = 8.3144598;     %J/(K*mol)
y_G0 = 3500;                %ppmV (0.35%)
T_G0 = 298;                 %K
P0 = 1.2*10^5;                  %Pa
rho_G0 = 28.97*10^(-3)*P0/(R*T_G0); %kg/m^3
L_R = 0.5;                  %m
D_R = 0.05;                 %m
d_W = 0.0012;               %m
rho_W = 7700;               %kg/m^3
Cp_W = 500;                 %J/kg-K
k_W = 19.51;                %W/m-K
T_Pre = 673;                %K
delH_R = -802500;           %J/mol
A_he = 1.58*10^(11);        %m^3/s-m^3 cat
Ea = 1.125*10^5;            %J/mol
d_pores = 6.377*10^(-9);    %m
e_pores = 0.519;
tau_pores = 2;
L_C = 46*10^(-6);
hWS = 0.005*0.042;


%Particle Bed Catalyst
rho_PB_C = 1257;            %kg/m^3
Cp_PB_C = 836;              %J/kg-K
k_PB_C = 0.042;             %W/m-K
e_PB_C = 0.36;              

%Particle Bed Inert
rho_PB_I = 2000;            %kg/m^3
Cp_PB_I = 890;              %J/kg-K
k_PB_I = 1.13;              %W/m-K
e_PB_I = 0.36;              

%Structured Bed
rho_SB = 2500;              %kg/m^3
Cp_SB = 965;                %J/kg-K
k_SB = 2.15;                %W/m-K
e_SB = 0.66;

u0 = 0.1;           %m/s
d_P = 0.003;      %m
N_m = 400;          %cpsi
D_H = 0.00116;  %m
aG = 4/D_H;
f_BC = 0.5;         
v = u0*e_SB;                %m/s

