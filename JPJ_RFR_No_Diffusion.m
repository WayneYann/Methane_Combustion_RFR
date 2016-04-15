y_G0 = 3500;
tswitch = 2200;
T_Pre = 673;
u0 = 0.1;           %m/s
N = 50; %Number of discretized points
Total_Runs = 5;

JPJ_variables_no_Diffusion
NE = 3; %Number of Equations per point
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
[T,Y,Eigenmatrix] =  ode15s_JPJ(@(t,u) JPJ_ODEs_No_Diffusion(t,u,N,NE,stepsize,y_G0,u0),tspan,uinitial,options); 
cputimeelapsed = cputime-startcpu
toc

T_Length = length(T);
Length_Vector = ones(T_Length,1);
Y0 = zeros(T_Length,3); %adding on values at X0 for all t
Y0(:,1) = Length_Vector*y_G0;
Y0(:,2) = Length_Vector*T_G0;
Y0(:,3) = Y(:,3);

Y = [Y0 Y];
xinitial = Y(end,1:(NE*N));

for i = 1:N
    uinitial(NE*(i-1)+1) = xinitial(NE*(N-i)+1);
    uinitial(NE*(i-1)+2) = xinitial(NE*(N-i)+2);
    uinitial(NE*(i-1)+3) = xinitial(NE*(N-i)+3);
end

end
xmesh = 0:stepsize:L_R; %goes from x0 to xN

figure(1)
surf(xmesh,T,Y(:,1:NE:(NE*N+1)),'edgecolor','none');
xlabel('Distance (m)','fontsize',20)
ylabel('Time (s)','fontsize',20)
zlabel('Concentration (ppm)','fontsize',20)
title('Concentration: 100 second drop','fontsize',20)

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
title('Solids Temperature Profile','fontsize',20)

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
        ttestten = abs((T(i) - tten));
    end
end;

figure(5);
whitebg('white')
plot(xmesh,Y(ttenIndex,1:NE:(NE*N+1)),'b-',xmesh...
    ,Y(tmiddleIndex,1:NE:(NE*N+1)),'b:',...
    xmesh,Y(length(T),1:NE:(NE*N+1)),'b--'...
,'LineWidth',2);
xlabel('Distance (m)','fontsize',20)
ylabel('Concentration (ppm)','fontsize',20)
title('Concentration Profile','fontsize',20)
legend('Early: 22 (s)','Middle: 1100 (s)','End: 2200 (s)')
 set(gcf, 'color', [1 1 1])
 
 figure(6);
whitebg('white')
plot(xmesh,Y(ttenIndex,3:NE:(NE*N+3)),'b-',xmesh...
    ,Y(tmiddleIndex,3:NE:(NE*N+3)),'b:',...
    xmesh,Y(length(T),3:NE:(NE*N+3)),'b--'...
,'LineWidth',2);
xlabel('Distance (m)','fontsize',20)
ylabel('Gas Temperature (K))','fontsize',20)
title('Temperature Profile','fontsize',20)
legend('Early: 22 (s)','Middle: 1100 (s)','End: 2200 (s)')
 set(gcf, 'color', [1 1 1])

   figure(7);
whitebg('white')
plot(Eigenmatrix(:,end));
xlabel('Jacobian Number','fontsize',15)
ylabel('Condition Number','fontsize',15)
title('Simplified Solution Condition Number: 300 points','fontsize',15)
 set(gcf, 'color', [1 1 1])
 
 

 

