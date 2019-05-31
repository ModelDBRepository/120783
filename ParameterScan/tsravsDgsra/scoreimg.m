load scores.txt
v1=scores(1:10,1:10);
v2=scores(11:20,1:10);
v3=scores(21:30,1:10);
v4=scores(31:40,1:10);
v5=scores(41:50,1:10);

v=v1+v2+v3+v4+v5;
v=v/5;
div=(v>=20)*(-0.2)
nondiv=(v<20).*v
% burst=((v>=2)-(v>=20))*1;
% irreg=(v<2)*0.4;
parimag=div+nondiv;
clims=[-0.2 1.0];
figure(1)
subplot(2,3,5);
y=14:14:140;
x=0.0:0.4:3.6;
imagesc(x,y,parimag,clims);
hold on;
plot(1.1, 60, 'k*', 'MarkerSize',8)
xlabel('\Delta g_{sra, Ipc}/g_{Ipc}','FontSize', 14)
ylabel('\tau_{sra, Ipc} (ms)','FontSize', 14)
set(gca,'YDir','normal');
set(gca,'XDir','normal');
set(gca, 'ytick', [0 70 140], 'FontSize', 10);
set(gca, 'xtick', [0.0 1.8 3.6], 'FontSize', 10);


