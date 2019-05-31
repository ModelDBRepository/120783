#include "Network.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "rng.h"
#include "matrixpool.h"
#include <ctime>

using namespace std;
Network::~Network(){}

Network::Network(double **_WmatTp21, double  **_WmatTp31, double **_WmatTp41, double **_WmatTp12,  double **_WmatTp23, double **_WmatTp14, Neuron *_NeGrp1, Neuron *_NeGrp2, Neuron *_NeGrp3, Neuron *_NeGrp4, Synapse *_SypGrp1,
				 Synapse *_SypGrp2, Synapse *_SypGrp3, Synapse *_SypGrp4, int _NeNo)
:WmatTp21(_WmatTp21), WmatTp31(_WmatTp31), WmatTp41(_WmatTp41), WmatTp12(_WmatTp12), WmatTp23(_WmatTp23), WmatTp14(_WmatTp14), NeGrp1(_NeGrp1), NeGrp2(_NeGrp2), NeGrp3(_NeGrp3), NeGrp4(_NeGrp4), SypGrp1(_SypGrp1), SypGrp2(_SypGrp2),
SypGrp3(_SypGrp3), SypGrp4(_SypGrp4), NeNo(_NeNo)
{
	
	fr1=0;
	fr2=0;
}																				//4 group system constructor

//Network(double **_WmatTp21, double  **_WmatTp31, double **_WmatTp41, double **_WmatTp12, double **_WmatTp23, double **_WmatTp14, Neuron *_NeGrp1, Neuron *_NeGrp2, Neuron *_NeGrp3, Neuron *_NeGrp4, Synapse *_SypGrp1,
// Synapse *_SypGrp2, Synapse *_SypGrp3, Synapse *_SypGrp4, int _NeNo);


int Network::NtWk_4GUpdat(double dt,  VecDP &I1, double i1sig, int I_loc1,  VecDP &I2, double i2sig, int I_loc2, RNG &randn, ofstream & Ipcstream, 
						   ofstream &L10stream, ofstream &Imc1stream, ofstream &Imc2stream, ofstream &voltstream )
{
		int iGp1=0, iGp2=0, iGp3=0, iGp4, iGp;
		int isspk=0;		
		double *IoutGp1=new double[NeNo];
		double *goutGp1=new double[NeNo];
		double *IoutGp2=new double[NeNo];
		double *goutGp2=new double[NeNo];
		double *IoutGp3=new double[NeNo];
		double *goutGp3=new double[NeNo];
		double *IoutGp4=new double[NeNo];
		double *goutGp4=new double[NeNo];


		double Itot_Gp2Gp1=0,gtot_Gp2Gp1=0;
		double Itot_Gp4Gp1=0,gtot_Gp4Gp1=0;
		double Itot_Gp1Gp2=0,gtot_Gp1Gp2=0;
		double Itot_Gp4Gp2=0,gtot_Gp4Gp2=0;
		double Itot_Gp3Gp2=0,gtot_Gp3Gp2=0;
		double Itot_Gp1Gp3=0,gtot_Gp1Gp3=0;
		double Itot_Gp1Gp4=0,gtot_Gp1Gp4=0;

		for(iGp=0;iGp<NeNo;iGp++){
			SypGrp1[iGp].PreSpkSet(NeGrp1[iGp].NeuronSignal());
			SypGrp1[iGp].SynpProbUpdat(dt);			
			goutGp1[iGp]=SypGrp1[iGp].gsynapout();
			IoutGp1[iGp]=SypGrp1[iGp].Isynapout();	//Projection group 1->other groups

		}					//update every synapse in the synapse group 1; 					
											
	    for(iGp=0;iGp<NeNo; iGp++){
			 SypGrp2[iGp].PreSpkSet(NeGrp2[iGp].NeuronSignal());
			 SypGrp2[iGp].SynpProbUpdat(dt);			
			 goutGp2[iGp]=SypGrp2[iGp].gsynapout();
			 IoutGp2[iGp]=SypGrp2[iGp].Isynapout();	//Projection group 2->other groups 

		  }				  //update every synapse in the synapse group 1; 	

	  
	  for(iGp1=0;iGp1<NeNo;iGp1++){
			  Itot_Gp2Gp1=dotprod(NeNo,WmatTp12[iGp1],IoutGp2);
			  gtot_Gp2Gp1=dotprod(NeNo,WmatTp12[iGp1],goutGp2);		//I&g from Group2 to Goup1
			  Itot_Gp4Gp1=dotprod(NeNo,WmatTp14[iGp1],IoutGp4);
			  gtot_Gp4Gp1=dotprod(NeNo,WmatTp14[iGp1],goutGp4);		//I&g from Group4 to Goup1
	
			  //NeGrp1[iGp1].I_inject(Iin1*exp(-pow(double(iGp1-I_loc1),2)/pow(i1sig,2))+Iin2*exp(-pow(double(iGp1-I_loc2),2)/pow(i2sig,2))); //RGC Gaussian input to L10
			  
			  //NeGrp1[iGp1].I_inject(I1(0)*((abs(iGp1-I_loc1)-i1sig/2)<0)+I2(0)*((abs(iGp1-I_loc2)-i2sig/2)<0)); //RGC square input to L10
			  
			  //NeGrp1[iGp1].I_inject(Iin1*fabs(gsdev(seed))*((abs(iGp1-I_loc1)-i1sig/2)<0)+Iin2*fabs(gsdev(seed))*((abs(iGp1-I_loc2)-i2sig/2)<0)); //RGC Random spatial distribution input to L10
			  
			  NeGrp1[iGp1].I_inject((I1(0)+I1(2)*randn.normal(I1(1),I1(3)))*((abs(iGp1-I_loc1)-i1sig/2)<0)+(I2(0)+I2(2)*randn.normal(I2(1),I2(3)))*((abs(iGp1-I_loc2)-i2sig/2)<0)); 
			  
			 /* NeGrp1[iGp1].I_inject((I1(0)+Irand1)*((abs(iGp1-I_loc1)-i1sig/2)<0)+(I2(0))*((abs(iGp1-I_loc2)-i2sig/2)<0)); */


			  //RGC Random spatial distribution input to L10

			  //NeGrp1[iGp1].I_inject(Iin1*(iGp1==I_loc1)+Iin2*(iGp1==I_loc2));
			  
			  NeGrp1[iGp1].toflowin(Itot_Gp2Gp1+Itot_Gp4Gp1,gtot_Gp2Gp1+gtot_Gp4Gp1);  //synapse current from group2
			
			  NeGrp1[iGp1].UpdatRK(dt, randn);			

		   	  NeGrp1[iGp1].Spkout(L10stream);
		  }						    								//Neuron Group1 Update
      
	  for(iGp2=0; iGp2<NeNo; iGp2++){ 
			  Itot_Gp1Gp2=dotprod(NeNo,WmatTp21[iGp2],IoutGp1);

			  gtot_Gp1Gp2=dotprod(NeNo,WmatTp21[iGp2],goutGp1);		//I&g from Group1 to Goup2	

			  Itot_Gp3Gp2=dotprod(NeNo,WmatTp23[iGp2],IoutGp3);

			  gtot_Gp3Gp2=dotprod(NeNo,WmatTp23[iGp2],goutGp3);		//I&g from Group3 to Goup2

			  NeGrp2[iGp2].toflowin(Itot_Gp1Gp2+Itot_Gp3Gp2, gtot_Gp1Gp2+gtot_Gp3Gp2);


																	//NeGrp2[iGp2].gAccSet(gtot_Gp3Gp2, gtot_Gp1Gp2);
			  NeGrp2[iGp2].UpdatRK(dt, randn);
			
			  NeGrp2[iGp2].Spkout(Ipcstream);
     			if(iGp2==I_loc1) 
				{
					NeGrp2[iGp2].Vout(voltstream);
					if(NeGrp2[iGp2].NeuronSignal()==1)
						isspk=1;

				} 


	  }					//Neuron Group2 Update																

  		delete [] IoutGp1;
		delete [] goutGp1;
		delete [] IoutGp2;
		delete [] goutGp2;
		delete [] IoutGp3;
		delete [] goutGp3;
		delete [] IoutGp4;
		delete [] goutGp4;
		return isspk;
	
}



int Network::NtWk_4G_CN_Updat(double dt,  VecDP &I1, double i1sig, int I_loc1,  VecDP &I2, double i2sig, int I_loc2, RNG &randn, ofstream & Ipcstream, 
						   ofstream &L10stream, ofstream &Imc1stream, ofstream &Imc2stream, ofstream &voltstream )
{
		int iGp1=0, iGp2=0, iGp3=0, iGp4, iGp;
		int isspk=0;		
		double *IoutGp1=new double[NeNo];
		double *goutGp1=new double[NeNo];
		double *IoutGp2=new double[NeNo];
		double *goutGp2=new double[NeNo];
		double *IoutGp3=new double[NeNo];
		double *goutGp3=new double[NeNo];
		double *IoutGp4=new double[NeNo];
		double *goutGp4=new double[NeNo];


		double Itot_Gp2Gp1=0,gtot_Gp2Gp1=0;
		double Itot_Gp4Gp1=0,gtot_Gp4Gp1=0;
		double Itot_Gp1Gp2=0,gtot_Gp1Gp2=0;
		double Itot_Gp4Gp2=0,gtot_Gp4Gp2=0;
		double Itot_Gp3Gp2=0,gtot_Gp3Gp2=0;
		double Itot_Gp1Gp3=0,gtot_Gp1Gp3=0;
		double Itot_Gp1Gp4=0,gtot_Gp1Gp4=0;

		double Irand1=I1(2)*randn.normal(I1(1),I1(3));

		for(iGp=0;iGp<NeNo;iGp++){
			SypGrp1[iGp].PreSpkSet(NeGrp1[iGp].NeuronSignal());
			SypGrp1[iGp].SynpProbUpdat(dt);			
			goutGp1[iGp]=SypGrp1[iGp].gsynapout();
			IoutGp1[iGp]=SypGrp1[iGp].Isynapout();	//Projection group 1->other groups

		}					//update every synapse in the synapse group 1; 					
											
	    for(iGp=0;iGp<NeNo; iGp++){
			 SypGrp2[iGp].PreSpkSet(NeGrp2[iGp].NeuronSignal());
			 SypGrp2[iGp].SynpProbUpdat(dt);			
			 goutGp2[iGp]=SypGrp2[iGp].gsynapout();
			 IoutGp2[iGp]=SypGrp2[iGp].Isynapout();	//Projection group 2->other groups 

		  }				  //update every synapse in the synapse group 1; 	

	  
	  for(iGp1=0;iGp1<NeNo;iGp1++){
			  Itot_Gp2Gp1=dotprod(NeNo,WmatTp12[iGp1],IoutGp2);
			  gtot_Gp2Gp1=dotprod(NeNo,WmatTp12[iGp1],goutGp2);		//I&g from Group2 to Goup1
			  Itot_Gp4Gp1=dotprod(NeNo,WmatTp14[iGp1],IoutGp4);
			  gtot_Gp4Gp1=dotprod(NeNo,WmatTp14[iGp1],goutGp4);		//I&g from Group4 to Goup1
	

			  
			  NeGrp1[iGp1].I_inject((I1(0)+Irand1)*((abs(iGp1-I_loc1)-i1sig/2)<0)+(I2(0))*((abs(iGp1-I_loc2)-i2sig/2)<0)); 


			  //RGC Random spatial distribution input to L10

			  //NeGrp1[iGp1].I_inject(Iin1*(iGp1==I_loc1)+Iin2*(iGp1==I_loc2));
			  
			  NeGrp1[iGp1].toflowin(Itot_Gp2Gp1+Itot_Gp4Gp1,gtot_Gp2Gp1+gtot_Gp4Gp1);  //synapse current from group2
			
			  NeGrp1[iGp1].UpdatRK(dt, randn);			

		   	  NeGrp1[iGp1].Spkout(L10stream);
		  }						    								//Neuron Group1 Update
      
	  for(iGp2=0; iGp2<NeNo; iGp2++){ 
			  Itot_Gp1Gp2=dotprod(NeNo,WmatTp21[iGp2],IoutGp1);

			  gtot_Gp1Gp2=dotprod(NeNo,WmatTp21[iGp2],goutGp1);		//I&g from Group1 to Goup2	

			  Itot_Gp3Gp2=dotprod(NeNo,WmatTp23[iGp2],IoutGp3);

			  gtot_Gp3Gp2=dotprod(NeNo,WmatTp23[iGp2],goutGp3);		//I&g from Group3 to Goup2

			  NeGrp2[iGp2].toflowin(Itot_Gp1Gp2+Itot_Gp3Gp2, gtot_Gp1Gp2+gtot_Gp3Gp2);


																	//NeGrp2[iGp2].gAccSet(gtot_Gp3Gp2, gtot_Gp1Gp2);
			  NeGrp2[iGp2].UpdatRK(dt, randn);
			
			  NeGrp2[iGp2].Spkout(Ipcstream);
     			if(iGp2==I_loc1) 
				{
					NeGrp2[iGp2].Vout(voltstream);
					if(NeGrp2[iGp2].NeuronSignal()==1)
						isspk=1;

				} 


	  }					//Neuron Group2 Update																

  		delete [] IoutGp1;
		delete [] goutGp1;
		delete [] IoutGp2;
		delete [] goutGp2;
		delete [] IoutGp3;
		delete [] goutGp3;
		delete [] IoutGp4;
		delete [] goutGp4;
		return isspk;
	
}



//with DAP
int Network::NtWk_4GUpdatDAC(double dt,  VecDP &I1, double i1sig, int I_loc1,  VecDP &I2, double i2sig, int I_loc2, RNG &randn, ofstream & Ipcstream, 
						   ofstream &L10stream, ofstream &Imc1stream, ofstream &Imc2stream, ofstream &voltstream )
{
		int iGp1=0, iGp2=0, iGp3=0, iGp4, iGp;
		int isspk=0;		
		double *IoutGp1=new double[NeNo];
		double *goutGp1=new double[NeNo];
		double *IoutGp2=new double[NeNo];
		double *goutGp2=new double[NeNo];
		double *IoutGp3=new double[NeNo];
		double *goutGp3=new double[NeNo];
		double *IoutGp4=new double[NeNo];
		double *goutGp4=new double[NeNo];


		double Itot_Gp2Gp1=0,gtot_Gp2Gp1=0;
		double Itot_Gp4Gp1=0,gtot_Gp4Gp1=0;
		double Itot_Gp1Gp2=0,gtot_Gp1Gp2=0;
		double Itot_Gp4Gp2=0,gtot_Gp4Gp2=0;
		double Itot_Gp3Gp2=0,gtot_Gp3Gp2=0;
		double Itot_Gp1Gp3=0,gtot_Gp1Gp3=0;
		double Itot_Gp1Gp4=0,gtot_Gp1Gp4=0;

		for(iGp=0;iGp<NeNo;iGp++){
			SypGrp1[iGp].PreSpkSet(NeGrp1[iGp].NeuronSignal());
			SypGrp1[iGp].SynpProbUpdat(dt);			//Synapse Group	1 Update
			goutGp1[iGp]=SypGrp1[iGp].gsynapout();
			IoutGp1[iGp]=SypGrp1[iGp].Isynapout();	//Projection group 1->other groups

		}										
											
	    for(iGp=0;iGp<NeNo; iGp++){
			 SypGrp2[iGp].PreSpkSet(NeGrp2[iGp].NeuronSignal());
			 SypGrp2[iGp].SynpProbUpdat(dt);			//Synapse Group 2 Update
			 goutGp2[iGp]=SypGrp2[iGp].gsynapout();
			 IoutGp2[iGp]=SypGrp2[iGp].Isynapout();	//Projection group 2->other groups 

		  }									  

	  
	  for(iGp1=0;iGp1<NeNo;iGp1++){
			  Itot_Gp2Gp1=dotprod(NeNo,WmatTp12[iGp1],IoutGp2);
			  gtot_Gp2Gp1=dotprod(NeNo,WmatTp12[iGp1],goutGp2);		//I&g from Group2 to Goup1
			  Itot_Gp4Gp1=dotprod(NeNo,WmatTp14[iGp1],IoutGp4);
			  gtot_Gp4Gp1=dotprod(NeNo,WmatTp14[iGp1],goutGp4);		//I&g from Group4 to Goup1
	
			  //NeGrp1[iGp1].I_inject(Iin1*exp(-pow(double(iGp1-I_loc1),2)/pow(i1sig,2))+Iin2*exp(-pow(double(iGp1-I_loc2),2)/pow(i2sig,2))); //RGC Gaussian input to L10
			  
			  //NeGrp1[iGp1].I_inject(I1(0)*((abs(iGp1-I_loc1)-i1sig/2)<0)+I2(0)*((abs(iGp1-I_loc2)-i2sig/2)<0)); //RGC square input to L10
			  
			  //NeGrp1[iGp1].I_inject(Iin1*fabs(gsdev(seed))*((abs(iGp1-I_loc1)-i1sig/2)<0)+Iin2*fabs(gsdev(seed))*((abs(iGp1-I_loc2)-i2sig/2)<0)); //RGC Random spatial distribution input to L10

			  NeGrp1[iGp1].I_inject((I1(0)+I1(2)*randn.normal(I1(1),I1(3)))*((abs(iGp1-I_loc1)-i1sig/2)<0)+(I2(0)+I2(2)*randn.normal(I2(1),I2(3)))*((abs(iGp1-I_loc2)-i2sig/2)<0)); 
			  
			  //RGC Random spatial distribution input to L10
			  
			  //NeGrp1[iGp1].I_inject(Iin1*(iGp1==I_loc1)+Iin2*(iGp1==I_loc2));
			  
			  NeGrp1[iGp1].toflowin(Itot_Gp2Gp1+Itot_Gp4Gp1,gtot_Gp2Gp1+gtot_Gp4Gp1);  //synapse current from group2			

			  NeGrp1[iGp1].UpdatRK(dt, randn);			
		   	  //NeGrp1[iGp1].Spkout(L10stream);
		  }						    								//Neuron Group1 Update
      
	  for(iGp2=0; iGp2<NeNo; iGp2++){ 
			  Itot_Gp1Gp2=dotprod(NeNo,WmatTp21[iGp2],IoutGp1);
			  gtot_Gp1Gp2=dotprod(NeNo,WmatTp21[iGp2],goutGp1);		//I&g from Group1 to Goup2	
			  Itot_Gp3Gp2=dotprod(NeNo,WmatTp23[iGp2],IoutGp3);
			  gtot_Gp3Gp2=dotprod(NeNo,WmatTp23[iGp2],goutGp3);		//I&g from Group3 to Goup2
			  NeGrp2[iGp2].toflowin(Itot_Gp1Gp2+Itot_Gp3Gp2, gtot_Gp1Gp2+gtot_Gp3Gp2);


																	//NeGrp2[iGp2].gAccSet(gtot_Gp3Gp2, gtot_Gp1Gp2);
			  NeGrp2[iGp2].UpdatRKwtDAC(dt, randn);
			
			  //NeGrp2[iGp2].Spkout(Ipcstream);
     			if(iGp2==I_loc1) 
				{
					//NeGrp2[iGp2].Vout(voltstream);
					if(NeGrp2[iGp2].NeuronSignal()==1)
						isspk=1;

				} 


	  }															

  		delete [] IoutGp1;
		delete [] goutGp1;
		delete [] IoutGp2;
		delete [] goutGp2;
		delete [] IoutGp3;
		delete [] goutGp3;
		delete [] IoutGp4;
		delete [] goutGp4;
		return isspk;
	
}

	  
