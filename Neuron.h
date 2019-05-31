#ifndef _NEURON_H
#define _NEURON_H
//Neuron class
#include <iostream>
#include <fstream>
#include "rng.h"
#include <vector>
using namespace std;
class Neuron
{public:
        Neuron(){}				//constructor
        
		~Neuron();				//Destructor

        Neuron(double _Er, double _taum, double _Rm, double _Vr, double _Vth, double _Esra, double _Dgsra, double _tausra,
			  double _V, double _gsra, double _sponmean, double _sponvar, double _spongate, int _spkflag=0, double _Ibias=0, double _Inet=0, double _gnet=0, double _Ispon=0);                             //constructor of Leaky-integrate-fire model

        Neuron(double _Er, double _taum, double _Rm, double _Vr, double _Vth, double _Esra, double _Dgsra, double _tausra, double _tauDAC, double _alpha, double _Adac, double _x, double _y,  double _V, double _gsra, double _sponmean, double _sponvar, double _spongate, int _spkflag=0, double _Ibias=0, double _Inet=0, double _gnet=0, double _Ispon=0);											//LIF-DAC neuron
        
        void toflowin(double _Inet, double _gnet );			//synapse input (gnet*Enet synaptic source, gnet conductance)
        
        void Updat(double dt, RNG &randspon);				//update neuron with time step dt and random number generator using method as descibed by Dyan Abbott

		void UpdatRK(double dt, RNG &randspon);			//update neuron using Runge-Kuta methods

		void UpdatRKwtDAC(double dt, RNG &randspon);	//update using afterdepolarization current

        void DACUpdat(double dt, double delayspk);		//update after depolarization current

        int NeuronSignal();		//send the spike no of neuron, 1 if the neuron spikes, 0 if the neuron does not.

		void I_inject(double current){Ibias=current;} //injecting biased current, by setting the Ibias varible to certain value.
        
		void Vout(ofstream &filename); //output the membrane potential of neuron to a file in the form of "V<<" ""
		
		void Spkout(ofstream &filename);//output the status of neuron spike to a file 1-spike, 0-nospike.
private:
         int spkflag;														//spike flag
         double Er, taum, Rm, Vr, Vth, tausra, Esra, Dgsra;					//spike flag variable                  
         double V;                                                          //Voltage and Gating Variables                           
         double gsra,gnet;                                                  //spike rate adaptation
         double Ibias, Inet, sponmean, sponvar, Ispon;
		 double spongate;													//spontaneous activity gate, 1-spontaneous activity considered, 0-no spontaneous activity.
		 double x, y;														//DAC current variable;
		 double alpha, tauDAC, Adac;										//DAC current;
		 vector<double> spktrain;											//recording the history of spike of the neuron, with in time period of tauDAC;
       };

class NeuronAcc//define a Neuron class with accommodatioin
{public:
        NeuronAcc(){}
        ~NeuronAcc();
        NeuronAcc(double _Er, double _taum0, double _Rm0, double _Vr, double _Vth0, double _Esra, double _Dgsra, double _tausra,
			  double _V, double _gsra, int _spkflag=0, double _Ibias=0, double _Inet=0, double _gnet=0);                                  //constructor of Leaky-integrate-fire model
        
        void toflowin(double _Inet, double _gnet );
		
		void tauthset(double tauvalue);
        
        void Updat(double dt);
        
       // void LIFSraUpdat(double dt);
        
        int NeuronSignal();		
		void I_inject(double current){Ibias=current;}
        void Vout(ofstream &filename);
		void gAccSet(double GACC1, double GACC2);
		void Spkout(ofstream &filename);
private:
         int spkflag;
         double Er, taum0, taum, Rm, Rm0, Vr, Vth, Vth0, tausra, Esra, Dgsra;						//spike flag variable                  
         double V;                                                               //Voltage and Gating Variables              
		 double DVth, gAcc, DRm;															//accomodation
         double gsra,gnet;                                                     //spike rate adaptation
         double Ibias,Inet, Ispon;
       };

#endif