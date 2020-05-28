%relaxation time for N=100,1000,10000

clear;close;clc;
data1=load('w10N100.dat');data1=data1(:,2);
data2=load('w10N500.dat');data2=data2(:,2);
data3=load('w10N1000.dat');data3=data3(:,2);
data4=load('w10N2000.dat');data4=data4(:,2);
data5=load('w10N5000.dat');data5=data5(:,2);
data6=load('w10N10000.dat');data6=data6(:,2);


l=length(data1);ll=length(data4);
n=9999;nn=1999;
com=l/n;com1=ll/nn;
sum1=zeros(n,1);sum2=zeros(n,1);sum3=zeros(n,1);
sum4=zeros(nn,1);sum5=zeros(nn,1);

for jj=1:com
    for ii=1:n
        sum1(ii)=sum1(ii)+data1((jj-1)*n+ii);
        sum2(ii)=sum2(ii)+data2((jj-1)*n+ii);
        sum3(ii)=sum3(ii)+data3((jj-1)*n+ii);
        
    end
end
for jj=1:com1
    for ii=1:nn
sum4(ii)=sum4(ii)+data4((jj-1)*nn+ii);
sum5(ii)=sum5(ii)+data5((jj-1)*nn+ii);
    end
end
figure()
plot(sum5/com,'.-c','Linewidth',1,'DisplayName','N = 5000')
hold on 
plot(sum4/com,'.-g','Linewidth',1,'DisplayName','N = 2000')
hold on 
plot(sum3/com,'.-r','Linewidth',1,'DisplayName','N = 1000')
hold on 
plot(sum2/com,'.-k','Linewidth',1,'DisplayName','N = 500')
hold on 
plot(sum1/com,'.-b','Linewidth',1,'DisplayName','N = 100')
axis([0 10000 0.12 0.45])
hold off
set(gca,'LineWidth',1.5)
set(gca, 'FontSize', 14)
set(gca, 'XScale', 'log')
%legend('show')

rect=[0.25, 0.25, .25, .25];
legend('show')
xlabel('Time','Interpreter','LaTex','Fontsize',20);
ylabel('T','Interpreter','LaTex','Fontsize',20)

print -depsc -painters relaxationw10.eps