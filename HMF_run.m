%% Keep the file HMF.m in the same folder
tic
clear all; close all; clc;

N=20; K=1; % No of particles and coupling parameter

Uspan=0.1:0.1:1.5; tspan=[0 100];
combination=10; average=1000;

Mav=zeros(length(Uspan),1);
Tav=zeros(length(Uspan),1);

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
 
mom = sqrt(2*(U*N-PE)/N)*(-1).^(floor(2*rand(N,1))); % Initial momenta

initial=[pos mom];

%% ODE45

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

[t,x]=ode45(@(t,y) HMF(t,y,N,K),tspan,initial,options);
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

PE = (K./(2*N))*(N^2-PE);

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

end

Mav(l) = Mav(l)/combination;
Tav(l) = Tav(l)/combination;

end


toc

fid=fopen('HMF_100_mat.dat','w');
fprintf(fid, '%f %f %f \n', [Uspan' Tav Mav]');
fclose(fid);