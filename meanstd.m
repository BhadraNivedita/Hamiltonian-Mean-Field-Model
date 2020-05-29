clear all;
close all;
clc;

w1=load('ts_w_1.dat');
m1=mean(w1); 
std1=std(w1);

w2=load('ts_w_2.dat');
m2=mean(w2); 
std2=std(w2);


w3=load('ts_w_3.dat');
m3=mean(w3);
std3=std(w3);


w4=load('ts_w_4.dat');
m4=mean(w4); 
std4=std(w4);


w5=load('ts_w_5.dat');
m5=mean(w5);
std5=std(w5);


w6=load('ts_w_6.dat');
m6=mean(w6); 
std6=std(w6);


w7=load('ts_w_7.dat');
m7=mean(w7); 
std7=std(w7);

w8=load('ts_w_8.dat');
m8=mean(w8); 
std8=std(w8);


w9=load('ts_w_9.dat');
m9=mean(w9); 
std9=std(w9);


w10=load('ts_w_10.dat');
m10=mean(w10); 
std10=std(w10);


w11=load('ts_w_11.dat');
m11=mean(w11); 
std11=std(w11);

w12=load('ts_w_12.dat');
m12=mean(w12); 
std12=std(w12);

w13=load('ts_w_13.dat');
m13=mean(w13); 
std13=std(w13);


w14=load('ts_w_14.dat');
m14=mean(w14);
std14=std(w14);

w15=load('ts_w_15.dat');
m15=mean(w15); 
std15=std(w15);

w16=load('ts_w_16.dat');
m16=mean(w16); 
std16=std(w16);

w17=load('ts_w_17.dat');
m17=mean(w17); 
std17=std(w17);

w18=load('ts_w_18.dat');
m18=mean(w18); 
std18=std(w18);

w19=load('ts_w_19.dat');
m19=mean(w19); 
std19=std(w19);

w20=load('ts_w_20.dat');
m20=mean(w20); 
std20=std(w20);
 
w=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]';
avg=[m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19 m20]';
standarddev=[std1 std2 std3 std4 std5 std6 std7 std8 std9 std10...
    std11 std12 std13  std14  std15  std16  std17  std18  std19  std20]';

figure(1)
subplot(2,1,1)
plot(w,avg,'-*')
ylabel('mean','Fontsize',20);
xlabel('w','Fontsize',20);

subplot(2,1,2)
plot(w,standarddev,'-or')
ylabel('Standard deviation','Fontsize',20);
xlabel('w','Fontsize',20);
print -depsc -painters 'U1.0.eps'

