clear; clc;

data=load('qss_veldist.dat');

N=20;
fsize=20;

[n1,x1]=hist(data,N);
f1=diff(x1(1:2));

plot(x1,smooth(n1)/sum(n1),'-ob','Linewidth',2,'DisplayName','w=0')