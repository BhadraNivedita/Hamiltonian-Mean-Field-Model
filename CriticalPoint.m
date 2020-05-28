% finds the critical point of hmf caloric curve by linear fitting

clear all; close all

a=load('w15.dat');

xx=a(:,1); yy=smooth(a(:,4));
x1=xx(41:46); y1=yy(41:46);
x2=xx(54:75); y2=yy(54:75);

figure(1)
plot(xx,yy,'*-b',x1,y1,'r',x2,y2,'r','LineWidth',2)

f1=fit(x1,y1,'poly1'); x1=xx(32:50);
f2=fit(x2,y2,'poly1'); x2=xx(40:75);

Uc=(f2.p2-f1.p2)/(f1.p1-f2.p1)
xc=ones(10,1)*Uc; yc=linspace(0,0.6,10);

figure(2)
plot(xx,yy,'*-b',x1,f1(x1),'r',x2,f2(x2),'r',xc,yc,'--k','LineWidth',2)

