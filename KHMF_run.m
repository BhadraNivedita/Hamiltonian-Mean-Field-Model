%% Keep the file HMF.m in the same folder
tic
clear all; close all; clc;

N=100; K=1;%omega=5; % No of particles and coupling parameter

Uspan=0.1:0.1:1.5; tspan=[0 150];
combination=20; average=100;

Mav=zeros(length(Uspan),1);
Tav=zeros(length(Uspan),1);
Hav=zeros(length(Uspan),1);

fid=fopen('HMF_diff_amp.dat','w');

omega=1;

for amp=[  .01 .1 .5 1 5]

for l=1:length(Uspan)

U=Uspan(l);

for k=1:combination

%% Initial positions and momenta

pos = 0.1*(randn(N,1)-0.5); % Random initial positions
 
PE=0;
for i=1:N
     for j=1:N
         PE = PE + cos(pos(i)-pos(j));
     end
end
PE = (K/(2*N))*(N^2-PE);
 
a = 2*sqrt((6/N)*(U*N-PE));
 
mom = a*(rand(N,1)-0.5); % Initial momenta

initial=[pos mom];

%% ODE45

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

[t,x]=ode45(@(t,y) KHMF(t,y,N,omega,amp),tspan,initial,options);

L=length(t);
H=zeros(L,1);
T=zeros(L,1);
M=zeros(L,1);

%% Calculating potential energy, kinetic energy and total energy in each time step

for m=1:L
    
KE=0; PE=0;

for j=(N+1):2*N
    KE = KE + x(m,j).^2;
end

KE = KE/2;

for i=1:N
    for j=1:N
        PE = PE + cos(x(m,i)-x(m,j));
    end
end

PE = (amp*cos(omega*t(m))./(2*N))*(N^2-PE);

H(m) = KE + PE;
T(m) = 2*KE/N;

%% Calculating Magnetization (Order parameter)

Mx=0; My=0;
for j=1:N
Mx = Mx + cos(x(m,j));
My = My + sin(x(m,j));
end
M(m) = sqrt(Mx.^2+My.^2)/N;
end

Mav(l) = Mav(l) + mean(M(end-average:end));
Tav(l) = Tav(l) + mean(T(end-average:end));
Hav(l) = Hav(l) + mean(H(end-average:end));

end

Mav(l) = Mav(l)/combination;
Tav(l) = Tav(l)/combination;
Hav(l) = Hav(l)/combination;


end

%plot(Tav,'.-r'); hold on
plot(Mav,'-ob');hold on
%hold off
%toc


fprintf(fid, '%f\t%f\t%f\n', [Hav Tav Mav]');
fprintf(fid,'\n');
end
fclose(fid);
toc