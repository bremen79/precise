% demo.m

close all
clear

n=10000;
delta=0.05;

% i.i.d. uniform samples in [0,1]
x=rand(n,1);

[l_co96,u_co96] = precise_co96(x,delta);
[l_a_co96,u_a_co96] = precise_a_co96(x,delta);
[l_r70,u_r70] = precise_r70(x,delta);

figure
h1=semilogx(l_co96,'r'); hold on; semilogx(u_co96,'r');
h2=semilogx(l_a_co96,'g'); hold on; semilogx(u_a_co96,'g');
h3=semilogx(l_r70,'b'); hold on; semilogx(u_r70,'b');
legend([h1,h2,h3],'PRECiSE-CO96', 'PRECiSE-A-CO96', 'PRECiSE-R70')