% plot MQSSPT for both w=0 and w=10



data1=load('dataw0.dat');
data2=load('dataw10.dat');
MQ1=data1(:,2);M1=data1(:,1);
MQ2=data2(:,2);M2=data2(:,1);


figure()


plot(MQ1,M1,'-or','Linewidth',2)


hold on 


plot(smooth(MQ2),M2,'-ob','Linewidth',2)

hold off

print -depsc -painters QSSPT.eps


