load Ipcspikes.txt
spikes=Ipcspikes;
cellno=length(spikes(1,:))-1;
trialgap=1.5;
xlim=[0 0.1*length(spikes)]
for i=2:(cellno+1)
    times=find(spikes(:,i));
    numspikes=length(times);
    xx=ones(3*numspikes,1)*nan;
    yy=ones(3*numspikes,1)*nan;
    yy(1:3:3*numspikes)=(i-1)*ones(numspikes,1)*trialgap;
    yy(2:3:3*numspikes)=yy(1:3:3*numspikes)+1;
    xx(1:3:3*numspikes)=spikes(times,1);
    xx(2:3:3*numspikes)=spikes(times,1);
    figure(1);
    subplot(3,1,2); 
    h=plot(xx,yy,'r');
    box on;
    hold on;    
end
axis ([0 0.1*length(spikes),0,(cellno)*1.5]);
xlabel('time(ms)','FontSize', 18)
ylabel('Ipc Neuron#','FontSize', 18)
set(gca,'YTick',0:0.75*cellno:1.5*cellno)
set(gca,'YTickLabel', {'0',num2str(cellno/2),num2str(cellno)})


load L10spikes.txt
spikes=L10spikes;
cellno=length(spikes(1,:))-1;
trialgap=1.5;
xlim=[0 0.1*length(spikes)]
for i=2:(cellno+1)
    times=find(spikes(:,i));
    numspikes=length(times);
    xx=ones(3*numspikes,1)*nan;
    yy=ones(3*numspikes,1)*nan;
    yy(1:3:3*numspikes)=(i-1)*ones(numspikes,1)*trialgap;
    yy(2:3:3*numspikes)=yy(1:3:3*numspikes)+1;
    xx(1:3:3*numspikes)=spikes(times,1);
    xx(2:3:3*numspikes)=spikes(times,1);
    subplot(3,1,1); 
    plot(xx,yy, 'r');
    box on;
    hold on;    
end
axis ([0 0.1*length(spikes),0,(cellno)*1.5]);
xlabel('time(ms)','FontSize', 18)
ylabel('L10 Neuron#','FontSize', 18)
set(gca,'YTick',0:0.75*cellno:1.5*cellno)
set(gca,'YTickLabel', {'0',num2str(cellno/2),num2str(cellno)})
