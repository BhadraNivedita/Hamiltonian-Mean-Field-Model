clear; clc

a=load('qss_psection.dat');

M=a(:,1);
% plot(M,'.-'); hold on
% plot(linspace(0,1000,10),mean(M)*ones(10,1))

x=[];
for i=2:length(M)-1
    if M(i)>(mean(M)-0.001) && M(i)<(mean(M)+0.001) && M(i)>M(i-1) && M(i)<M(i+1)
        x=[x,i];
    end
end

aa=a(:,2:end); np=size(aa,2)/2;

for i=1:np
    p=aa(:,2*i-1); q=mod(abs(aa(:,2*i)),2*pi); plot(q(x),p(x),'.'); hold on
end
hold off

% plot(x,mean(M)*ones(length(x),1),'o'); hold off