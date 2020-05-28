%analysis of relaxation time 
clear;close;clc;
data1=load('w0N100.dat');data1=data1(:,2);
data2=load('w0N500.dat');data2=data2(:,2);
data3=load('w0N1000.dat');data3=data3(:,2);

data4=load('w10N100.dat');data4=data4(:,2);
data5=load('w10N500.dat');data5=data5(:,2);
data6=load('w10N1000.dat');data6=data6(:,2);
l=length(data1);
n=9999;N1=100;N2=500;N3=1000;
com=l/n;
sum1=zeros(n,1);sum2=zeros(n,1);sum3=zeros(n,1);
sum4=zeros(n,1);sum5=zeros(n,1);sum6=zeros(n,1);

for jj=1:com
    for ii=1:n
        sum1(ii)=sum1(ii)+data1((jj-1)*n+ii);
        sum2(ii)=sum2(ii)+data2((jj-1)*n+ii);
        sum3(ii)=sum3(ii)+data3((jj-1)*n+ii);
        sum4(ii)=sum4(ii)+data4((jj-1)*n+ii);
        sum5(ii)=sum5(ii)+data5((jj-1)*n+ii);
        sum6(ii)=sum6(ii)+data6((jj-1)*n+ii);
    end
end
figure()
subplot(1,2,1)
plot(sum1/com,'.r')
hold on 
plot(sum2/com,'.k')
hold on 
plot(sum3/com,'.b')
set(gca, 'XScale', 'log')
%axis([0 10000 0 0.05])
subplot(1,2,2)
plot(sum4/com,'.r')
hold on 
plot(sum5/com,'.k')
hold on 
plot(sum6/com,'.b')
set(gca, 'XScale', 'log')

print -depsc -painters relaxationw0.eps