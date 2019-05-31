#ifndef _SYNAPSE_H
#define _SYNAPSE_H
//Synapse class
#include "Neuron.h"
class Synapse
{
public:
      Synapse(){}
	  ~Synapse();//destructor
	  Synapse(double _tau1, double _tau2, double _gsynap, double _Esynap, 
      double _Ps, double _Ps1, double _Ps2);	//constructor

      void PreSpkSet(int PreSynapSpk){PreSpk=PreSynapSpk;}	//detecting the presynapic neuronal activity, PreSynapSpk;

      double Isynapout(){return (gsynap*Esynap*Ps);}         //output current sources;

      double gsynapout(){return (gsynap*Ps);}                //output conductances;

      void SynpProbUpdat(double dt);						//update each time step dt;

      void SynapDisp();										//output the synapse information into a file;
     
     
private:
      double tau1, tau2, gsynap, Esynap, Ps,  Ps1, Ps2;
      int PreSpk;											//parameters and varialbes.
   
      };
#endif
