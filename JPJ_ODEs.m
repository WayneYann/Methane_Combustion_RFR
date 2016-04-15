function du = JPJ_ODEs(t,u,N,NE,stepsize)

JPJ_variables; %import Variables

du = zeros(NE*N,1);
%These only work for structured bed
e = e_SB;
rho_S = rho_SB;
Cp_S = Cp_SB;
k_Seff = k_SB*(1-e_SB);
aS = aG*(e/(1-e));


Tg = u(2:4:(NE*N));
Ts = u(3:4:(NE*N));
TgAVG = sum(Tg)/N;
TsAVG = sum(Ts)/N;
%Using average T for reaction
%{
kVC = A_he*exp(-1*Ea/(8.314*TsAVG));
rho_G = 28.97*10^(-3)*P0/(R*TgAVG);
CG = P0/(R*TgAVG);
D_AB = 9.86*10^(-10)*TgAVG^1.75;
Deff = D_AB+(v*D_H)^2/(196*D_AB);
k_G = 0.01679+6.073*10^(-5)*TgAVG;
mu = 7.701*10^(-6)+4.166*10^(-8)*TgAVG-7.531*10^(-12)*TgAVG^2;
Re = u0*D_H*rho_G/mu;
Sc = mu/rho_G/D_AB;
Pe = 1/(0.73*e/(Re*Sc)+0.5/(1+9.7*e/Re/Sc));
k_Geff = u0*D_H*rho_G*Cp_G/Pe;
Dk = 97*d_pores/2*(TsAVG/(28*10^-3))^0.5;
Deff_pore = Dk*e_pores/tau_pores;
thiele = L_C*(kVC/Deff_pore)^(0.5);
n = tanh(thiele)/thiele;
KG = (2+1.1*Sc^(1/3)*Re^(0.6))*D_AB/D_H;
Pr = Cp_G*mu/k_G;
Gz = D_H/0.5*Re*Pr;
NuT = 2.977*(1+3.6*(Gz)^.5*exp(-50/Gz));
%h = (2+1.1*Pr^(1/3)*Re^(0.6))*k_G/D_H(i);
h = NuT*k_G/D_H;
kheS = L_C*kVC;
%}

%Matrix for derivatives (used in finite difference)
yT02 = zeros(NE,2);
for i = 1:NE:(NE*(N-1)+1)



kVC = A_he*exp(-1*Ea/(8.314*u(i+2)));
rho_G = 28.97*10^(-3)*P0/(R*u(i+1));
CG = P0/(R*u(i+1));
D_AB = 9.86*10^(-10)*u(i+1)^1.75;
Deff = D_AB+(v*D_H)^2/(196*D_AB);
k_G = 0.01679+6.073*10^(-5)*u(i+1);
mu = 7.701*10^(-6)+4.166*10^(-8)*u(i+1)-7.531*10^(-12)*u(i+1)^2;
Re = u0*D_H*rho_G/mu;
Sc = mu/rho_G/D_AB;
Pe = 1/(0.73*e/(Re*Sc)+0.5/(1+9.7*e/Re/Sc));
k_Geff = u0*D_H*rho_G*Cp_G/Pe;
Dk = 97*d_pores/2*(u(i+2)/(28*10^-3))^0.5;
Deff_pore = Dk*e_pores/tau_pores;
thiele = L_C*(kVC/Deff_pore)^(0.5);
n = tanh(thiele)/thiele;   
KG = (2+1.1*Sc^(1/3)*Re^(0.6))*D_AB/D_H;
Pr = Cp_G*mu/k_G;
Gz = D_H/0.5*Re*Pr;
NuT = 2.977*(1+3.6*(Gz)^.5*exp(-50/Gz));
h = NuT*k_G/D_H;
kheS = L_C*kVC;


    if i > 0.25*(NE*N) && i < 0.75*(NE*N)
    ys = u(i)/(1+(n*kheS/KG));
    rheS = -1*n*kheS*CG*ys/10^6;
    elseif i >= 0.75*(NE*N)
    ys = u(i)/(1+(n*kheS/KG));
    rheS = -1*n*kheS*CG*ys/10^6*10^-15;
    else
    rheS = 0;
    ys = u(i);
    end;
  
    if i == 1
       yT02(1,1) = y_G0 + e*Deff/u0*(u(i)-y_G0)/(2*stepsize);
       yT02(2,1) = T_G0 + e*k_Geff/(u0*rho_G0*Cp_G)*(u(i+1)-T_G0)/(2*stepsize);
       yT02(3,1) = u(i+2+NE);
       yT02(4,1) = u(i+3+NE);
       
       yT02(1,2) = u(i+NE);
       yT02(2,2) = u(i+1 + NE);
       yT02(3,2) = u(i+2 + NE);
       yT02(4,2) = u(i+3 + NE); 
     
    elseif i == (NE*(N-1)+1)
    
       yT02(1,1) = u(i-NE);
       yT02(2,1) = u(i+1 - NE);
       yT02(3,1) = u(i+2 - NE);
       yT02(4,1) = u(i+3 - NE);
       
       yT02(:,2) = yT02(:,1);
        
    else
       yT02(1,1) = u(i-NE);
       yT02(2,1) = u(i+1 - NE);
       yT02(3,1) = u(i+2 - NE);
       yT02(4,1) = u(i+3 - NE);
       
       yT02(1,2) = u(i+NE);
       yT02(2,2) = u(i+1 + NE);
       yT02(3,2) = u(i+2 + NE);
       yT02(4,2) = u(i+3 + NE);  
                
    end;
    
    ddz = (yT02(:,2) - yT02(:,1))/(2*stepsize);
    d2dz2 = (yT02(:,2) - 2*[u(i); u(i+1); u(i+2); u(i+3)] + yT02(:,1))/(stepsize^2);
    
    du(i) = (-1*u0/e_SB*rho_G0/rho_G*ddz(1) + Deff*d2dz2(1)...
        -1*aG*KG*(u(i)-ys));
    
    du(i+1) = (-u0/e*(rho_G0/rho_G)*ddz(2) + k_Geff/(rho_G*Cp_G)*d2dz2(2)...
        -1*aG*h/(rho_G*Cp_G)*(u(i+1)-u(i+2)) -4*hWG/(e*D_R*rho_G*Cp_G)*(u(i+1)-u(i+3)));
    
    du(i+2) = (k_Seff/(rho_S*Cp_S)*d2dz2(3)...
        -aS*h/rho_S/Cp_S*(u(i+2)-u(i+1))+aS*rheS*delH_R/(rho_S*Cp_S) -1*4*hWS/((1-e)*D_R*rho_S*Cp_S)*(u(i+2)-u(i+3)));
    
    du(i+3) = ( k_W/(rho_W*Cp_W)*d2dz2(4)...
  +hWG*D_R*(u(i+1)-u(i+3))+hWS*D_R*(u(i+2)-u(i+3)))/(rho_W*Cp_W*d_W*(D_R+d_W));

end

end
