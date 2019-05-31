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
subplot(2,3,4);
y=10:10:100;
x=0.0:0.02:0.18
imagesc(x,y,parimag,clims);
xlabel('\sigma_{L10} (nA)','FontSize', 14)
ylabel('\Delta_{L10\rightarrow Ipc}','FontSize', 14)
set(gca,'YDir','normal');
set(gca,'XDir','normal');
set(gca, 'ytick', [0 50 100], 'FontSize', 10);
set(gca, 'xtick', [0.0 0.09 0.18], 'FontSize', 10);
hold on;
plot(0.1, 50, 'k*', 'MarkerSize',8)
% axis square
% subplot(2,3,6)
%  colorbar;

