function du = JPJ_ODEs_No_Diffusion(t,u,N,NE,stepsize,y_G0,u0)

JPJ_variables_no_Diffusion; %import Variables

du = zeros(NE*N,1);
%These only work for structured bed
e = e_SB;
rho_S = rho_SB;
Cp_S = Cp_SB;
aS = aG*(e/(1-e));
%Matrix for derivatives (used in finite difference)
yT02 = zeros(NE,1);
for i = 1:NE:(NE*(N-1)+1)

kVC = A_he*exp(-1*Ea/(8.314*u(i+2)));
rho_G = 28.97*10^(-3)*P0/(R*u(i+1));
CG = P0/(R*u(i+1));
D_AB = 9.86*10^(-10)*u(i+1)^1.75;
k_G = 0.01679+6.073*10^(-5)*u(i+1);
mu = 7.701*10^(-6)+4.166*10^(-8)*u(i+1)-7.531*10^(-12)*u(i+1)^2;
Re = u0*D_H*rho_G/mu;
Sc = mu/rho_G/D_AB;
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

%This code adds in the effect of the inert
    if i > 0.25*(NE*N) && i < 0.75*(NE*N)
    ys = u(i)/(1+(n*kheS/KG));
    rheS = -1*n*kheS*CG*ys/10^6;
    else
    rheS = 0;
    ys = u(i);
    end;
  
   if i == 1
     
     %Use this code test temporary changes in feed concentration
     %if t < 100
      % yT02(1) = 2*y_G0;
       %else
       yT02(1) = y_G0;
      %end
       yT02(2) = T_G0;
       yT02(3) = u(i+2);
     
    else
       yT02(1) = u(i-NE);
       yT02(2) = u(i+1 - NE);
       yT02(3) = u(i+2 - NE);
                
    end;
    
    ddz = ([u(i); u(i+1); u(i+2)] - yT02)/stepsize;
    du(i) = (-1*u0/e_SB*rho_G0/rho_G*ddz(1) -1*aG*KG*(u(i)-ys));

    du(i+1) = (-u0/e*(rho_G0/rho_G)*ddz(2) -1*aG*h/(rho_G*Cp_G)*(u(i+1)-u(i+2)));
    
    du(i+2) = (-aS*h/rho_S/Cp_S*(u(i+2)-u(i+1))+aS*rheS*delH_R/(rho_S*Cp_S));
   
end

end
