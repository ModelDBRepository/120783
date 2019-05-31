#include "Neuron.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
//Neurons Definitions
Neuron::~Neuron(){}

Neuron::Neuron(double _Er, double _taum, double _Rm, double _Vr, double _Vth, double _Esra, double _Dgsra, double _tausra,
			   double _V, double _gsra, double _sponmean, double _sponvar, double _spongate, int _spkflag, double _Ibias, double _Inet, double _gnet, double _Ispon)
:Er(_Er), taum(_taum), Rm(_Rm), Vr(_Vr), Vth(_Vth), Esra(_Esra), Dgsra(_Dgsra), tausra(_tausra), V(_V), gsra(_gsra), sponmean(_sponmean), sponvar(_sponvar), spongate(_spongate), spkflag(_spkflag), Ibias(_Ibias),
Inet(_Inet), gnet(_gnet), Ispon(_Ispon)
{}                                                                     
//Leak-Integrate-Fire model inti-list
//Neuron(Er, taum, Rm, Vr, Vth, Esra, Dgsra, tausra, V, gsra, spkflag);



Neuron::Neuron(double _Er, double _taum, double _Rm, double _Vr, double _Vth, double _Esra, double _Dgsra, double _tausra, double _tauDAC, double _alpha, double _Adac,double _x, double _y,  double _V, double _gsra, double _sponmean, double _sponvar, double _spongate, int _spkflag, double _Ibias, double _Inet, double _gnet, double _Ispon)
:Er(_Er), taum(_taum), Rm(_Rm), Vr(_Vr), Vth(_Vth), Esra(_Esra), Dgsra(_Dgsra), tausra(_tausra), tauDAC(_tauDAC), alpha(_alpha), Adac(_Adac), x(_x), y(_y), V(_V), gsra(_gsra), sponmean(_sponmean), sponvar(_sponvar), spongate(_spongate), spkflag(_spkflag), Ibias(_Ibias),
Inet(_Inet), gnet(_gnet), Ispon(_Ispon)
{}  


void Neuron::Vout(ofstream &filename){
     filename<<V<<" " ;
     //cout<<V<<" "<<spkflag<<endl; 
   }


void Neuron::Spkout(ofstream &filename){filename<<spkflag<<" ";}


void Neuron::toflowin(double _Inet, double _gnet)
{	
	Inet=_Inet;
	gnet=_gnet;
}

int Neuron::NeuronSignal()
{return spkflag;}
   


void Neuron::Updat(double dt, RNG &randspon)
{

     double Vinf,taumeff;
	 if(spkflag==1){
		 V=Vr;
		 spkflag=0;
	 }

	 else{
		 Vinf=(Er+Rm*gsra*Esra+Rm*Inet+Rm*Ibias+spongate*Rm*randspon.normal(sponmean, sponvar))/(1+Rm*gsra+gnet*Rm);
		 taumeff=taum/(1+Rm*gsra+gnet*Rm);
		 V=Vinf+(V-Vinf)*exp(-dt/taumeff);
		 if (V>Vth){
			 V=1;
			 spkflag=1;
			 gsra+=Dgsra;
		 }
		 else{
			 spkflag=0;
			 gsra=gsra*exp(-dt/tausra);
		 }
	 }
}

void Neuron::UpdatRK(double dt, RNG &randspon)
{

     double Vinf,taumeff;
	 double k1, k2, k3, k4;
	 if(spkflag==1){
		 V=Vr;
		 spkflag=0;
	 }

	 else{
		 Vinf=(Er+Rm*gsra*Esra+Rm*Inet+Rm*Ibias+spongate*Rm*randspon.normal(sponmean, sponvar))/(1+Rm*gsra+gnet*Rm);
		 taumeff=taum/(1+Rm*gsra+gnet*Rm);
		 k1=(-V+Vinf)*dt/taumeff;
		 k2=(-(V+k1/2)+Vinf)*dt/taumeff;
		 k3=(-(V+k2/2)+Vinf)*dt/taumeff;
		 k4=(-(V+k3)+Vinf)*dt/taumeff;
		 V=V+k1/6+k2/3+k3/3+k4/6;		//RungeKuta iteration


		 if (V>Vth){
			 V=1;
			 spkflag=1;
			 gsra+=Dgsra;
		 }
		 else{
			 spkflag=0;
			 gsra=gsra*exp(-dt/tausra);
		 }
	 }
}

void Neuron::UpdatRKwtDAC(double dt, RNG &randspon)
{

     double Vinf,taumeff;
	 double k1, k2, k3, k4;
	 if(spkflag==1){
		 V=Vr;
		 spkflag=0;
	 }

	 else{
		 Vinf=(Er+Rm*gsra*Esra+Rm*Inet+Rm*x*Adac+Rm*Ibias+spongate*Rm*randspon.normal(sponmean, sponvar))/(1+Rm*gsra+gnet*Rm);
		 taumeff=taum/(1+Rm*gsra+gnet*Rm);
		 k1=(-V+Vinf)*dt/taumeff;
		 k2=(-(V+k1/2)+Vinf)*dt/taumeff;
		 k3=(-(V+k2/2)+Vinf)*dt/taumeff;
		 k4=(-(V+k3)+Vinf)*dt/taumeff;
		 V=V+k1/6+k2/3+k3/3+k4/6;		//RungeKuta iteration

		 if (V>Vth){
			 V=1;
			 spkflag=1;
			 gsra+=Dgsra;
		 }
		 else{
			 spkflag=0;
			 gsra=gsra*exp(-dt/tausra);
		 }
	 }
	 //spktrain.push_back(spkflag);
	 //if(spktrain.size()>int(tauDAC/dt))
		// spktrain.erase (spktrain.begin());
	 //DACUpdat(dt, spktrain.at(0));
	 spktrain.push_back(spkflag);
	 if(spktrain.size()<=int(tauDAC/dt))
		 DACUpdat(dt, 0);
	 else{
		 spktrain.erase (spktrain.begin());
		 DACUpdat(dt, spktrain.at(0));
	 }

}

void Neuron::DACUpdat(double dt, double delayspk)
{
	double k1, k2, k3, k4;
	double l1, l2, l3, l4;
	k1=y*dt;			l1=(-alpha*alpha*x-2*alpha*y)*dt+alpha*alpha*delayspk;
	k2=(y+l1/2)*dt;		l2=(-alpha*alpha*(x+k1/2)-2*alpha*(y+l1/2))*dt+alpha*alpha*delayspk;
	k3=(y+l2/2)*dt;		l3=(-alpha*alpha*(x+k2/2)-2*alpha*(y+l2/2))*dt+alpha*alpha*delayspk;
	k4=(y+l3)*dt;		l4=(-alpha*alpha*(x+k3)-2*alpha*(y+l3))*dt+alpha*alpha*delayspk;
	x=x+k1/6+k2/3+k3/3+k4/6;	
	y=y+l1/6+l2/3+l3/3+l4/6;

}



//AccNeuron
//Neurons Definitions
NeuronAcc::~NeuronAcc(){}
NeuronAcc::NeuronAcc(double _Er, double _taum0, double _Rm0, double _Vr, double _Vth0, double _Esra, double _Dgsra, double _tausra,
			   double _V, double _gsra, int _spkflag, double _Ibias, double _Inet, double _gnet)
:Er(_Er), taum0(_taum0), Rm0(_Rm0), Vr(_Vr), Vth0(_Vth0), Esra(_Esra), Dgsra(_Dgsra), tausra(_tausra), V(_V), gsra(_gsra), spkflag(_spkflag), Ibias(_Ibias),
Inet(_Inet), gnet(_gnet)
{
	Vth=Vth0;
	Rm=Rm0;
	DVth=0;
	gAcc=0;
	DRm=0;
	taum=taum0;
}                                                                     
//Leak-Integrate-Fire model inti-list
//Neuron(Er, taum, Rm, Vr, Vth, Esra, Dgsra, tausra, V, gsra, spkflag);

void NeuronAcc::Vout(ofstream &filename){
     filename<<V<<" " <<Vth<<" "<<DVth<<" ";
     //cout<<V<<" "<<spkflag<<endl; 
   }


void NeuronAcc::Spkout(ofstream &filename){
	filename<<spkflag<<" ";
}

void NeuronAcc::toflowin(double _Inet, double _gnet)
{
	Inet=_Inet;
	gnet=_gnet;
}

int NeuronAcc::NeuronSignal()
{
	return spkflag;
}
   
void NeuronAcc::gAccSet(double GACC1, double GACC2)
{
	
	gAcc=0;
}


void NeuronAcc::Updat(double dt)
{

     double Vinf,taumeff;

	// DVth=DVth+(-DVth+(1-exp(-3.5*Rm0*gAcc))*(V-Vr))/(0.5*taum0)*dt;
	 DVth=DVth+(-DVth+3.5*Rm0*gAcc*(V-Vr))/(0.5*taum)*dt;
    //cout<<DVth<<endl;
	 
	 Vth=Vth0+DVth;
	 //if (DVth<0) cout<<Vth<<endl;

	 DRm=DRm+(-DRm+2.0/3.0*Rm0*(1-exp(-(V-Vr)/15.0)))/(taum)*dt;
	 //Rm=Rm0-DRm;
	 Rm=Rm0;
	 taum=taum0*Rm/Rm0;
	 //cout<<DRm<<" "<<Rm<<" "<<taum<<endl;

	 if(spkflag==1){
		 V=Vr;
		 spkflag=0;
	 }

	 else{
		 Vinf=(Er+Rm*gsra*Esra+Rm*Inet+Rm*Ibias)/(1+Rm*gsra+gnet*Rm);
		 taumeff=taum/(1+Rm*gsra+gnet*Rm);
		 V=Vinf+(V-Vinf)*exp(-dt/taumeff);
		 if (V>Vth){
			 V=50;
			 spkflag=1;
			 gsra+=Dgsra;
		 }
		 else{
			 spkflag=0;
			 gsra=gsra*exp(-dt/tausra);
		 }
	 }
}

