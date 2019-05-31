#ifndef _NETWORK_H
#define _NETWORK_H
//Network Class
#include "Neuron.h"
#include "Synapse.h"
#include <iostream>
#include <fstream>
#include "DataType.h"
#include "rng.h"
using namespace std;

class Network
{
private:
    Neuron  *NeGrp3, *NeGrp4, *NeGrp1, *NeGrp2; //define 4 groups of neuron
	
	Synapse *SypGrp1, *SypGrp2, *SypGrp3, *SypGrp4; //the synapse associated with the neurons.
	
	int NeNo;						//neuron number
	
	double **WmatTp21, **WmatTp31, **WmatTp41, **WmatTp12, **WmatTp23, **WmatTp14;
	//weight matrix from ith gourp to the jth group WmatTpij; 
	
	double fr1,fr2;//firing rate of neuron 1 and 2

public:
	Network(){}
	~Network();


	Network(double **_WmatTp21, double  **_WmatTp31, double **_WmatTp41, double **_WmatTp12, double **_WmatTp23, double **_WmatTp14, Neuron *_NeGrp1, Neuron *_NeGrp2, Neuron *_NeGrp3, Neuron *_NeGrp4, Synapse *_SypGrp1, Synapse *_SypGrp2, Synapse *_SypGrp3, Synapse *_SypGrp4, int _NeNo);
	//network constructor;

/*----more groups of neurons-----*/
	
	int NtWk_4GUpdat(double dt, VecDP &I1, double i1sig, int I_loc1, VecDP &I2, double i2sig, int I_loc2, RNG &randn, ofstream & Ipcstream, 
						   ofstream &L10stream, ofstream &Imc1stream, ofstream &Imc2stream, ofstream &voltstream );
	//network update;
	
	int NtWk_4GUpdatDAC(double dt, VecDP &I1, double i1sig, int I_loc1, VecDP &I2, double i2sig, int I_loc2, RNG &randn, ofstream & Ipcstream, 
						   ofstream &L10stream, ofstream &Imc1stream, ofstream &Imc2stream, ofstream &voltstream );
	//update network considering the ADP

	int NtWk_4G_CN_Updat(double dt,  VecDP &I1, double i1sig, int I_loc1,  VecDP &I2, double i2sig, int I_loc2, RNG &randn, ofstream & Ipcstream, 
						   ofstream &L10stream, ofstream &Imc1stream, ofstream &Imc2stream, ofstream &voltstream );

	double f1rate(){return fr1;}

	void fr1reset(){fr1=0;}

	double f2rate(){return fr2;}

	void fr2reset(){fr2=0;}
//function prototype
//double dotprod(int n, double *vec1, double *vec2);
//function definition

};
#endif