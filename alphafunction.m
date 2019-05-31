x=0;
y=0;
dt=0.1;
alpha=1;
taudac=10;
t=20;
X=[];
Y=[];

for i=0:1000        
	k1=y*dt;			l1=(-alpha*alpha*x-2*alpha*y)*dt+alpha^2*(i==t/dt);
	k2=(y+l1/2)*dt;		l2=(-alpha*alpha*(x+k1/2)-2*alpha*(y+l1/2))*dt+alpha^2*(i==t/dt);
	k3=(y+l2/2)*dt;		l3=(-alpha*alpha*(x+k2/2)-2*alpha*(y+l2/2))*dt+alpha^2*(i==t/dt);
	k4=(y+l3)*dt;		l4=(-alpha*alpha*(x+k3)-2*alpha*(y+l3))*dt+alpha^2*(i==t/dt);
	x=x+k1/6+k2/3+k3/3+k4/6;	
	y=y+l1/6+l2/3+l3/3+l4/6;
    X=[X,x];
    Y=[Y,y];
end
plot((0:1000)*dt,X);