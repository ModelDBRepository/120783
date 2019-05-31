#include "Synapse.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
//Synapse definitions
Synapse::~Synapse(){}
Synapse::Synapse(double _tau1, double _tau2, double _gsynap, double _Esynap, 
      double _Ps, double _Ps1, double _Ps2)
:tau1(_tau1), tau2(_tau2), gsynap(_gsynap), Esynap(_Esynap), Ps(_Ps), Ps1(_Ps1), Ps2(_Ps2)
{}
//(tau1, tau2, gsynap, Esynap, Ps, Ps1, Ps2)
//void Synapse::PreSpkSet(Neuron& PreSynapNe){}

void Synapse::SynpProbUpdat(double dt){
     double B, taurise;
     taurise=tau1*tau2/(tau1-tau2);
     B=1/(pow(tau2/tau1,taurise/tau1)-pow(tau2/tau1,taurise/tau2));
     if (PreSpk==0){
          Ps1=Ps1*exp(-dt/tau1);
          Ps2=Ps2*exp(-dt/tau2);
          }
     else{
         Ps1+=B;
         Ps2+=B;
         }
     Ps=Ps1-Ps2;
}

void Synapse::SynapDisp(){
     ofstream myfile("Synapse.txt", ios::app);
     myfile<<Ps<<" "<<PreSpk<<" "<<gsynap*Esynap*Ps<<endl;
     myfile.close();  
   }
