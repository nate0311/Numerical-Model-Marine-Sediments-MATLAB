%This is a numerical model for 227Ac in deep-sea sediments (This code is also valid 
%for the behavior of other radioisotopes in marine sediments as well). 
%The only transport mechanism is diffusion, and advection has been ignored
%Equation that will be used numerical in this model -->
%dC/dt=(Ds+K*Db)*d2C/dx2 - lambda*(1+K)*C + P 

%%
lambdaAc=log(2)/21.77;                      %decay constant 227Ac in 1/yr
lambdaAc1=(log(2)/21.77)/(365.25*24*60);    %decay constant 227Ac in 1/min
 
lambdaPa=log(2)/32760;                      %decay constant 231Pa in 1/yr
lambdaPa1=(log(2)/32760)/(365.25*24*60);    %decay constant 231Pa in 1/min
 
F=0.49;                                     %fraction released of mobile Ac 
rho=2.5;                                    %density of sediments in g/cm3
 
Apa=3.09;                                   %activity of 231Pa in dpm/g
Apaa=Apa*(365.25*24*60);                    %activity of 231Pa in dpy/g
 
Dm=2.9e-6*(60*60*24*365.25);                 %molecular diffusion(cm2/yr) 
Db=0.007;                                   %bioturbation rate (cm2/yr)
 
por=0.87;                                   %porosity
Kd=15125;                                   %distribution coefficient (cc/g)  
K=(Kd*rho*(1-por))/(por);                 	 %dimensionless coefficient(no units) 
 
dx=0.1;                                     %length of interval
N=100;                                      %number of Boxes in model
 
%%
%Create arrays for porosity and K
 
%porosity array with depth starting at 0.1 cm 
%equation is the fit to porosity profile
x=[1:1:N+1];
xlength=length(x);
por1=zeros(xlength,1);

for j=1:1:xlength; 
        por1(j)=0.81501 - 0.086889*log10(j*dx); %Fit porosity profile at location of interests. This fit is from North Pacific
end
 
%K array: start by defining Kd in the upper few cm
%Kd starts at measured K value from plots: Insert starting K value above

K1=zeros(xlength,1);

for i=1:1:N+1; 
      K1(i)=(Kd*rho*(1-por1(i)))/(por1(i));
end
 
%create Ds from Dm
Ds=Dm*(por1.^2);
%%
%initial conditions:
 
por=por1(1);
porl=por1(N+1);                       	%create porosity for last interval
Plast=(F*rho*(1-porl)*Apaa)/(porl);   	%production term for last interval box
 
%P=[(F*rho*(1-por)*Apaa)/(por)];        	%production term units in atoms/cc-yr
 
C=((Plast)/(lambdaAc*(1+K1(N)))*ones(N,1)); %create initial array of all P
 
C(N+1)=((Plast)/(lambdaAc*(1+K1(N))));   	%boundary conditions P: (atoms/cc-yr)
 
dt=0.00004;                                    %time step in units of yr
sumtime=0;
 
dC=zeros(xlength,1);
 
%%
%for loop
for timestep=1:2500000;
  
CW00=3;%INITIAL CONCENTRATION AT SWI (UNITS OF dpm/m3)
%concentration in water above SWI: 3.0 dpm/m3 = 0.000003 dpm/cc = 49.6 atoms/cc
 
CW0=CW00/1000000;%convert m3 to cc
CW=CW0/lambdaAc1; %convert dpm to atoms   %initial concentration,
 
%first box dC
dC(1)=dt*[(Ds(1)+K1(1)*Db)*por1(1)*((CW - C(1))/(dx^2/2)) + (Ds(1)+K1(1)*Db)*por1(1)*((C(2)-C(1))/(dx^2))- lambdaAc*(1+K1(1))*C(1) + ((F*rho*(1-por1(1))*Apaa)/(por1(1)))];
%because the distance from middle of box 1 to overlying water is 1/2 of a box dimension.
    
for i=2:1:N;%for loop for dC 
dC(i)=dt*[((Ds(i)+K1(i)*Db)*por1(i)*(C(i-1) - 2*C(i) + C(i+1))/(dx^2)) - lambdaAc*(1+K1(i))*C(i) + ((F*rho*(1-por1(i))*Apaa)/(por1(i)))];   
end
 
for i=1:1:N;%increment boxes with dc   
    C(i) = C(i) + dC(i);    
end
 
sumtime=sumtime+dt;
sumdays=sumtime*365.25;
end
 
sumetime=sumtime+dt;

%% Plot 227Ac vs. depth in dpm/cc

C1=lambdaAc1*C; % multiple C by 227Ac decay constant in 1/min to convert to dpm (atoms*lambda=Activity)

xx=[0:dx:N*dx]; % x array for plotting

figure(1);%dpm/cc vs. depth
plot(C1,xx);
ylim([0.0 10]);
set(gca, 'YDir','reverse');
title('227Ac vs. depth','fontsize',18);
xlabel('227Ac (dpm/cc)','fontsize',18);
ylabel('Depth (cm)','fontsize',18);

%% 227Ac Flux
Ds(1);   %Ds in first interval
por1(1); %porosity in first interval

Slope=(C1(1)-CW0)/(dx/2); %units=dpm/(g*cm)
Flux3=(Slope)*Ds(1)*por1(1)*100^2; %dpm/m2-yr


%% Sum Square Values

B1=[2.33,3.35,3.05,3.11]; %measured 227Ac values (dpm/g)
x1=[0.5,1.5,2.5,7.5];     %depth of 227Ac measured values (cm)
Er=[0.09,0.15,0.11,0.18]; %measured 227Ac error values

%%mass summation (denominator)
mx=rho*(1-por1)*dx; %g/cc*(1-porosity)
mxb= reshape(mx(1:100), [], 10);
summass = sum(mxb); %sum every 10 intervals. Sum down to 10cm
summass1=sum(summass);%total vector sum, all 10cm 
 
%%mass summation (numerator)  
C3=rho*(1-por1).*C1*dx; %g/cc*(1-porosity)
C3b= reshape(C3(1:100), [], 10);
Actmass = sum(C3b); %sum every 10 intervals. Sum down to 10cm
sumActmass1=sum(Actmass);%total vector sum, all 10cm


%%take model values and average for every cm, including the mass, final units are in dpm/g

for i=1:1:length(summass); 
      Cmodel(i)=Kd*Actmass(i)/summass(i)+ (1-F)*Apa;
end

%compare measured and calculated 227Ac
sumsquare=0;
for i=1:1:3;
    sumsquare=sumsquare+(B1(i)-Cmodel(i))^2;
end
sumsquare;


