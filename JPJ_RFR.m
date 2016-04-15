JPJ_variables
tswitch = 2500;
N = 300; %Number of discretized points
Total_Runs = 1;
NE = 4; %Number of Equations per point
tspan = [0 tswitch];
stepsize = L_R/N;
%concentration matrix
uinitial = ones(NE*N,1)*T_Pre;
for i = 1:NE:(NE*(N-1)+1)
        uinitial(i) = 0;
end

options = odeset('RelTol', 10^-3);
for RunNumber = 1:Total_Runs;
display(RunNumber)
tic
startcpu = cputime;
[T,Y, Eigenmatrix] =  ode15s_JPJ(@(t,u) JPJ_ODEs(t,u,N,NE,stepsize),tspan,uinitial,options);
cputimeelapsed = cputime-startcpu
toc

e = e_SB;
rho_G = 28.97*10^(-3)*P0./(R*Y(:,2));
D_AB = 9.86*10^(-10)*Y(:,2).^1.75;
Deff = D_AB+(v*D_H)^2./(196*D_AB);
mu = 7.701*10^(-6)+4.166*10^(-8)*Y(:,2)-7.531*10^(-12)*Y(:,2).^2;
Re = u0*D_H*rho_G./mu;
Sc = mu./rho_G./D_AB;
Pe = 1./(0.73*e./(Re.*Sc)+0.5./(1+9.7*e./Re./Sc));
k_Geff = u0*D_H*rho_G.*Cp_G./Pe;

T_Length = length(T);
Length_Vector = ones(T_Length,1);
Y0 = zeros(T_Length,4); %adding on values at X0 for all t
Y0(:,1) = Length_Vector*y_G0 + e*Deff/u0.*(Y(:,1)-Length_Vector*y_G0)./(2*stepsize);
Y0(:,2) = Length_Vector*T_G0 + e*k_Geff/(u0*rho_G0*Cp_G).*(Y(:,2)-Length_Vector*T_G0)./(2*stepsize);
Y0(:,3) = Y(:,3);
Y0(:,4) = Y(:,4);

Y = [Y0 Y];
xinitial = Y(end,1:(NE*N));

for i = 1:N
    uinitial(NE*(i-1)+1) = xinitial(NE*(N-i)+1);
    uinitial(NE*(i-1)+2) = xinitial(NE*(N-i)+2);
    uinitial(NE*(i-1)+3) = xinitial(NE*(N-i)+3);
    uinitial(NE*(i-1)+4) = xinitial(NE*(N-i)+4);
end

end
xmesh = 0:stepsize:L_R; %goes from x0 to xN

figure(1)
surf(xmesh,T,Y(:,1:NE:(NE*N+1)),'edgecolor','none');
xlabel('Distance (m)','fontsize',20)
ylabel('Time (s)','fontsize',20)
zlabel('Methane Concentration (ppm)','fontsize',20)
title('Concentration profile','fontsize',20)

figure(2)
surf(xmesh,T,Y(:,2:NE:(NE*N+2)),'edgecolor','none');
xlabel('Distance (m)','fontsize',20)
ylabel('Time (s)','fontsize',20)
zlabel('Gas Temperature (K)','fontsize',20)
title('Gas Temperature','fontsize',20)

figure(3)
surf(xmesh,T,Y(:,3:NE:(NE*N+3)),'edgecolor','none');
xlabel('Distance (m)','fontsize',20)
ylabel('Time (s)','fontsize',20)
zlabel('Solids Temperature (K)','fontsize',20)
title('Solids Temperature','fontsize',20)


figure(4)
surf(xmesh,T,Y(:,4:NE:(NE*N+4)),'edgecolor','none');
xlabel('Distance (m)','fontsize',20)
ylabel('Time (s)','fontsize',20)
zlabel('Wall Temperature (K)','fontsize',20)
title('Wall Temperature','fontsize',20)


tmiddle = round(tswitch/2);
tmiddleIndex = 0;
ttest =100;
tten = round(0.01*tswitch);
ttenIndex = 0;
ttestten = 100;

for i = 1:length(T)
    if abs((T(i) - tmiddle)) < ttest
        tmiddleIndex    = i;  
        ttest = abs((T(i) - tmiddle));
    end
    
     if abs((T(i) - tten)) < ttestten
        ttenIndex    = i;  
        ttestten = abs((T(i) - tmiddle));
    end
end;

figure(5);
whitebg('white')
plot(xmesh,Y(ttenIndex,1:NE:(NE*N+1)),'b-',xmesh...
    ,Y(tmiddleIndex,1:NE:(NE*N+1)),'b:',...
    xmesh,Y(length(T),1:NE:(NE*N+1)),'b--'...
,'LineWidth',2);
xlabel('Distance (m)','fontsize',20)
zlabel('Methane Concentration (ppm)','fontsize',20)
title('Concentration profile in Reactor','fontsize',20)
 set(gcf, 'color', [1 1 1])
 
 figure(6);
whitebg('white')
plot(xmesh,Y(ttenIndex,2:NE:(NE*N+2)),'b-',xmesh...
    ,Y(tmiddleIndex,2:NE:(NE*N+2)),'b:',...
    xmesh,Y(length(T),2:NE:(NE*N+2)),'b--'...
,'LineWidth',2);
xlabel('Distance (m)','fontsize',20)
zlabel('Gas Temperature (K))','fontsize',20)
title('Temperature Profile in Reactor','fontsize',20)
 set(gcf, 'color', [1 1 1])
 
  figure(7);
whitebg('white')
plot(Eigenmatrix(:,end));
xlabel('Jacobian Number','fontsize',15)
ylabel('Condition Number','fontsize',15)
title('Full Solution Condition Number: 300 points','fontsize',15)
 set(gcf, 'color', [1 1 1])




