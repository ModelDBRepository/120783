// The system have different inputs between L10 to Ipc and Imc neurons
// NO accomdation is considered.
#include "Network.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "matrixpool.h"
#include "Sys.h"
#include "rng.h"
#include <ctime>
using namespace std;

Sys::~Sys(){}
Sys::Sys(double _T, int _ISPOT1, int _ISPOT2, double _AMP1, double _AMP2, double _Anos1, double _Anos2,  double _nctrl1, double _nctrl2, double _nsd1, double _nsd2,
		 double _SIG1, double _SIG2)
:T(_T), ISPOT1(_ISPOT1), ISPOT2(_ISPOT2), AMP1(_AMP1), AMP2(_AMP2), Anos1(_Anos1), Anos2(_Anos2), nctrl1(_nctrl1),  nctrl2(_nctrl2), nsd1(_nsd1),
nsd2(_nsd2),SIG1(_SIG1), SIG2(_SIG2)
{}

double Sys::score(VecDP & VarParSet,VecDP & _L10Par, VecDP &_IpcPar, VecDP & _IpcAdapPar, VecDP & _L10AdapPar, RNG & randn) //Weight matrix, L10 spontaneous parameters, Ipc spontaneous parameters/
{
	  static const int NN=400;
      double dt=0.1;
	  int i=0, iters=(int) T/dt;
      double **Mat21=0, **Mat31=0, **Mat41=0, **Mat12=0, **Mat23=0, **Mat14=0;
	  double r1=0, r2=0; 
	  int isspk=0;
	  int spki=0, prespki=0;
	  double preDspk=0, Dspk=0;
	  double bno=0, rno=0, score=0;
	  int isb=0;
	  int divalert=0;
	  VecDP Iamp1(4);
	  VecDP Iamp2(4);
	  Iamp1(0)=AMP1, Iamp1(1)=Anos1, Iamp1(2)=nctrl1, Iamp1(3)=nsd1;   
	  Iamp2(0)=AMP2, Iamp2(1)=Anos2, Iamp2(2)=nctrl2, Iamp2(3)=nsd2;

	  Mat21=new double *[NN];
	  Mat31=new double *[NN];
  	  Mat41=new double *[NN];
	  Mat12=new double *[NN];
	  Mat23=new double *[NN];
 	  Mat14=new double *[NN];
	  
	  Mat21[0]=new double [NN*NN];
  	  Mat31[0]=new double [NN*NN];
  	  Mat41[0]=new double [NN*NN];
	  Mat12[0]=new double [NN*NN];
	  Mat23[0]=new double [NN*NN];
   	  Mat14[0]=new double [NN*NN];
	  

	  
	  ofstream Ipcf( "Ipcspikes.txt", ios::app);
	  ofstream L10f("L10spikes.txt",ios::app);
	  ofstream voltf( "voltage.txt", ios::app);
	  ofstream Imc1f( "Imc1spikes.txt", ios::app);
	  ofstream Imc2f( "Imc2spikes.txt", ios::app);
	  ofstream Input("input2L10.txt", ios::app);
	  ofstream paraout("Parameters.txt", ios::app);
	  ofstream wmatf("wmatrix.txt",ios::app);


	  NormTopMat(Mat21, VarParSet(0), VarParSet(1), NN, NN);
	  NormTopMat(Mat31, 0.0, 50, NN, NN);
	  NormTopMat(Mat41, 0.0, 50, NN, NN);
	  NormTopMat(Mat12, VarParSet(2), VarParSet(3), NN, NN);
	  //RandTopMat(Mat12, VarParSet(2), VarParSet(3), NN, NN);
	  //nRandTopMat(Mat12, VarParSet(2), VarParSet(3), NN, NN);
	  NulMat(Mat23,  NN, NN);
	  NulMat(Mat14,  NN, NN);
	  
	 // MatDP matrix12(NN, NN, Mat12[0]);
	 // matrix12.disp(wmatf);
	  paraout<<"Stimu Amp "<<AMP1<<"\n"<<"Noise Amp "<<Anos1<<endl;
	  paraout<<"Mat21"<<" "<<VarParSet(0)<<"\t"<<VarParSet(1)<<"\n"<<"Mat12"<<" "<<VarParSet(2)<<"\t"<<VarParSet(3)<<"\n";
	  paraout<<"NN "<<NN<<"\n"<<"Stimu Loc "<<ISPOT1<<"\n"<<"Stimu width "<<SIG1<<endl;
	  paraout<<"L10 spon activity"<<"mean "<<_L10Par(0)<<"\t"<<"sd "<<_L10Par(1)<<endl;
	  paraout<<"Ipc spon activity"<<"mean "<<_IpcPar(0)<<"\t"<<"sd "<<_IpcPar(1)<<endl;
	

	  //CosTopMat(MatT2, 225, NN, NN);
	  //NulMat(MatT2,NN, NN);
	  //ShowMat(MatT2, NN,NN);

      Neuron L10[NN];
	  Neuron Ipc[NN];
	  Neuron ImcT1[NN], ImcT2[NN];         //Neuron(Er, taum, Rm, Vr, Vth, Esra, Dgsra, tausra, V, gsra, spkflag, Ibias);
      Synapse AMPA[NN];
	  Synapse NMDA[NN];
	  Synapse GABA1[NN], GABA2[NN];      //(tau1, tau2, gsynap, Esynap, Ps, Ps1, Ps2)
	  
	//Neuron(Er, taum, Rm, Vr, Vth, Esra, Dgsra, tausra,V, gsra, sponmean, sponvar, spongate);      
	  
	  for(i=0;i<NN;i++){
		  L10[i]=Neuron(-55.0, 104.0, 480.0, -50.0, -39.0, -70.0, _L10AdapPar(0)/480.0, _L10AdapPar(1), -55.0, 0.0, _L10Par(0), _L10Par(1), 0.0);
		
      }      


      for(i=0;i<NN;i++){
		  Ipc[i]=Neuron(-61.0, 25.0, 135.0, -50.0, -40.0, -70.0, _IpcAdapPar(0)/135.0, _IpcAdapPar(1), -61.0, 0.0,  _IpcPar(0), _IpcPar(1), 0.0);
		  
       } //Ipc neuron without DAP


	   // for(i=0;i<NN;i++){
		  //Ipc[i]=Neuron(-61.0, 25.0, 135.0, -50.0, -40.0, -70.0, _IpcAdapPar(0)/135.0, _IpcAdapPar(1), 0.5, 4.5, 0.7, 0.0, 0.0,  -61.0, 0.0,  _IpcPar(0), _IpcPar(1), 1.0);
		  //
    //   } //Ipc neuron with DAP
	  
	  
	  //Neuron Group Initialization
/*	  for(i=0;i<NN;i++){
		  Ipc[i]=Neuron(-61.0, 35.0, 135.0, -60.0, -35.0, -70.0, 0/135.0, 60.0, -61.0, 0.0,  0.0); 
       }	*/		

	  for(i=0;i<NN;i++){
		  ImcT1[i]=Neuron(-58, 150, 130.0, -50.0, -40.0, -70.0, 1.0/130.0, 50.0, -58.0, 0.0,  0.0, 1.0, 0.0);
       }
				  //Neuron(Er, taum, Rm, Vr, Vth, Esra, Dgsra, tausra, V, gsra);
	  for(i=0;i<NN;i++){
		  ImcT2[i]=Neuron(-58, 150, 130.0, -50.0, -40.0, -70.0, 1.0/130.0, 50.0, -58.0, 0.0,  0.0, 1.0, 0.0);
       }	


      for(i=0;i<NN;i++){
		  AMPA[i]=Synapse(5.6, 0.3, 0.25/135.0, 0.0, 0.0, 0.0, 0.0);
      }
    
	  for(i=0;i<NN;i++){
		  NMDA[i]=Synapse(10.0, 1.0, 0.15*0.015/480.0, -5.0, 0.0, 0.0, 0.0);
       }				//Synapse Group Initialization
	  
	  for(i=0;i<NN;i++){
		  GABA1[i]=Synapse(5.6, 0.3, 0/180.0, -80.0, 0.0, 0.0, 0.0);
	  }	

	  for(i=0;i<NN;i++){
		  GABA2[i]=Synapse(7.6, 0.3, 0/180.0, -80.0, 0.0, 0.0, 0.0);
	  }	
      
	  Network L10IpcImcSys(Mat21, Mat31, Mat41, Mat12, Mat23, Mat14, L10, Ipc, ImcT1, ImcT2, AMPA, NMDA, GABA1, GABA2, NN);
		//Network(double **_WmatTp21, double  **_WmatTp31, double **_WmatTp41, double **_WmatTp12, double **_WmatTp23, double **_WmatTp14, Neuron *_NeGrp1, Neuron *_NeGrp2, Neuron *_NeGrp3, Neuron *_NeGrp4, Synapse *_SypGrp1,
		//Synapse *_SypGrp2, Synapse *_SypGrp3, Synapse *_SypGrp4, int _NeNo);

	 for(i=0;i<iters;i++){
		  L10f<<i*dt<<"\t";
		  Ipcf<<i*dt<<"\t";
		  voltf<<i*dt<<"\t";
		  //Imc1f<<i*dt<<"\t";
		  //Imc2f<<i*dt<<"\t";
		  //Iamp1=(AMP1+0.0*fabs(gsdev(seed)))*(i*dt/(dt*iters*1/3)*(i<iters/3)+(i>=iters/3&&i<=iters*2/3)+(i>iters*2/3)*(iters-i)*dt/(dt*iters*1/3));
		  
		  isspk=L10IpcImcSys.NtWk_4GUpdat(dt, Iamp1, SIG1, ISPOT1, Iamp2, SIG2, ISPOT2, randn, Ipcf, L10f, Imc1f, Imc2f, voltf);		//single gaussian input with width SIG1

		  //isspk=L10IpcImcSys.NtWk_4GUpdatDAC(dt, Iamp1, SIG1, ISPOT1, Iamp2, SIG2, ISPOT2, randn, Ipcf, L10f, Imc1f, Imc2f, voltf);		//single gaussian input with width SIG1

		  if(isspk==1){
			  prespki=spki;
			  spki=i;
		  	  preDspk=Dspk;
			  Dspk=(spki-prespki)*dt;
			  if(i*dt>=100){
				  if(isb==0){
					  if(preDspk>10&&Dspk<=4)
					  {
						  isb=1;
						  bno++; //increase the burst number
						  //cout<<"bno="<<bno;
					  }
					  else{
						  rno++; //increase the regular spike number
						  //cout<<"rno="<<rno;
					  }
				  }
				  else{
					  if(Dspk>4)
					  {
						  isb=0;
					  }
				  }	

			  }
			  divalert++;
		  }
		  if (divalert>dt*iters){return 100; } //when the system diverge stop the loop


		  //cout<<fabs(AMP1*gsdev(seed));

		  Input<<Iamp1(0)<<"\t"<<Iamp1(1)<<"\n";
		  L10f<<"\n";
		  Ipcf<<"\n";
		  voltf<<"\n";
		  //Imc1f<<"\n";
		  //Imc2f<<"\n";
	  }
      
		delete [] (Mat21[0]);
		delete [] (Mat21);

		delete [] (Mat31[0]);
		delete [] (Mat31);

		delete [] (Mat41[0]);
		delete [] (Mat41);
		
		delete [] (Mat12[0]);
		delete [] (Mat12);

		delete [] (Mat23[0]);
		delete [] (Mat23);

		delete [] (Mat14[0]);
		delete [] (Mat14);
		L10f.close();
		voltf.close();
		Ipcf.close();
		Imc1f.close();
		Imc2f.close();
		paraout.close();
		Input.close();
		wmatf.close();
		return bno/(rno+bno);

	 
   }
       



void Sys::disp(VecDP & VarParSet,VecDP & _L10Par, VecDP &_IpcPar, VecDP & _IpcAdapPar, VecDP & _L10AdapPar, RNG & randn) //Weight matrix, L10 spontaneous parameters, Ipc spontaneous parameters/
{
	  static const int NN=400;
      double dt=0.1;
	  int i=0, iters=(int) T/dt;
      double **Mat21=0, **Mat31=0, **Mat41=0, **Mat12=0, **Mat23=0, **Mat14=0;
	  double r1=0, r2=0; 
	  int isspk=0;
	  int spki=0, prespki=0;
	  double preDspk=0, Dspk=0;
	  double bno=0, rno=0, score=0;
	  int isb=0;
	  int divalert=0;
	  VecDP Iamp0(4);
	  VecDP Iamp1(4);
	  VecDP Iamp2(4);
	  Iamp1(0)=AMP1, Iamp1(1)=Anos1, Iamp1(2)=nctrl1, Iamp1(3)=nsd1;   
	  Iamp2(0)=AMP2, Iamp2(1)=Anos2, Iamp2(2)=nctrl2, Iamp2(3)=nsd2;
	  Iamp0(0)=0.0, Iamp1(1)=Anos1, Iamp1(2)=nctrl1, Iamp1(3)=nsd1;

	  Mat21=new double *[NN];
	  Mat31=new double *[NN];
  	  Mat41=new double *[NN];
	  Mat12=new double *[NN];
	  Mat23=new double *[NN];
 	  Mat14=new double *[NN];
	  
	  Mat21[0]=new double [NN*NN];
  	  Mat31[0]=new double [NN*NN];
  	  Mat41[0]=new double [NN*NN];
	  Mat12[0]=new double [NN*NN];
	  Mat23[0]=new double [NN*NN];
   	  Mat14[0]=new double [NN*NN];
	  

	  
	  ofstream Ipcf( "Ipcspikes.txt", ios::app);
	  ofstream L10f("L10spikes.txt",ios::app);
	  ofstream voltf( "voltage.txt", ios::app);
	  ofstream Imc1f( "Imc1spikes.txt", ios::app);
	  ofstream Imc2f( "Imc2spikes.txt", ios::app);
	  ofstream Input("input2L10.txt", ios::app);
	  ofstream paraout("Parameters.txt", ios::app);
	  ofstream wmatf("wmatrix.txt",ios::app);


	  NormTopMat(Mat21, VarParSet(0), VarParSet(1), NN, NN);
	  NormTopMat(Mat31, 0.0, 50, NN, NN);
	  NormTopMat(Mat41, 0.0, 50, NN, NN);
	  NormTopMat(Mat12, VarParSet(2), VarParSet(3), NN, NN);
	  //RandTopMat(Mat12, VarParSet(2), VarParSet(3), NN, NN);
	  //nRandTopMat(Mat12, VarParSet(2), VarParSet(3), NN, NN);
	  NulMat(Mat23,  NN, NN);
	  NulMat(Mat14,  NN, NN);
	  
	 // MatDP matrix12(NN, NN, Mat12[0]);
	 // matrix12.disp(wmatf);
	  paraout<<"Stimu Amp "<<AMP1<<"\n"<<"Noise Amp "<<Anos1<<endl;
	  paraout<<"Mat21"<<" "<<VarParSet(0)<<"\t"<<VarParSet(1)<<"\n"<<"Mat12"<<" "<<VarParSet(2)<<"\t"<<VarParSet(3)<<"\n";
	  paraout<<"NN "<<NN<<"\n"<<"Stimu Loc "<<ISPOT1<<"\n"<<"Stimu width "<<SIG1<<endl;
	  paraout<<"L10 spon activity"<<"mean "<<_L10Par(0)<<"\t"<<"sd "<<_L10Par(1)<<endl;
	  paraout<<"Ipc spon activity"<<"mean "<<_IpcPar(0)<<"\t"<<"sd "<<_IpcPar(1)<<endl;
	

	  //CosTopMat(MatT2, 225, NN, NN);
	  //NulMat(MatT2,NN, NN);
	  //ShowMat(MatT2, NN,NN);

      Neuron L10[NN];
	  Neuron Ipc[NN];
	  Neuron ImcT1[NN], ImcT2[NN];         //Neuron(Er, taum, Rm, Vr, Vth, Esra, Dgsra, tausra, V, gsra, spkflag, Ibias);
      Synapse AMPA[NN];
	  Synapse NMDA[NN];
	  Synapse GABA1[NN], GABA2[NN];      //(tau1, tau2, gsynap, Esynap, Ps, Ps1, Ps2)
	  
	//Neuron(Er, taum, Rm, Vr, Vth, Esra, Dgsra, tausra,V, gsra, sponmean, sponvar, spongate);      
	  
	  for(i=0;i<NN;i++){
		  L10[i]=Neuron(-55.0, 104.0, 480.0, -50.0, -39.0, -70.0, _L10AdapPar(0)/480.0, _L10AdapPar(1), -55.0, 0.0, _L10Par(0), _L10Par(1), 1.0);
		
      }      


    //  for(i=0;i<NN;i++){
		  //Ipc[i]=Neuron(-61.0, 25.0, 135.0, -50.0, -40.0, -70.0, _IpcAdapPar(0)/135.0, _IpcAdapPar(1), -61.0, 0.0,  _IpcPar(0), _IpcPar(1), 0.0);
		  //
    //   } //Ipc neuron without DAP


	    for(i=0;i<NN;i++){
		  Ipc[i]=Neuron(-61.0, 25.0, 135.0, -50.0, -40.0, -70.0, _IpcAdapPar(0)/135.0, _IpcAdapPar(1), 0.5, 4.5, 0.7, 0.0, 0.0,  -61.0, 0.0,  _IpcPar(0), _IpcPar(1), 1.0);
		  
       } //Ipc neuron with DAP
	  
	  
	  //Neuron Group Initialization
/*	  for(i=0;i<NN;i++){
		  Ipc[i]=Neuron(-61.0, 35.0, 135.0, -60.0, -35.0, -70.0, 0/135.0, 60.0, -61.0, 0.0,  0.0); 
       }	*/		

	  for(i=0;i<NN;i++){
		  ImcT1[i]=Neuron(-58, 150, 130.0, -50.0, -40.0, -70.0, 1.0/130.0, 50.0, -58.0, 0.0,  0.0, 1.0, 0.0);
       }
				  //Neuron(Er, taum, Rm, Vr, Vth, Esra, Dgsra, tausra, V, gsra);
	  for(i=0;i<NN;i++){
		  ImcT2[i]=Neuron(-58, 150, 130.0, -50.0, -40.0, -70.0, 1.0/130.0, 50.0, -58.0, 0.0,  0.0, 1.0, 0.0);
       }	


      for(i=0;i<NN;i++){
		  AMPA[i]=Synapse(5.6, 0.3, 0.25/135.0, 0.0, 0.0, 0.0, 0.0);
      }
    
	  for(i=0;i<NN;i++){
		  NMDA[i]=Synapse(10.0, 1.0, 0.15*0.015/480.0, -5.0, 0.0, 0.0, 0.0);
       }				//Synapse Group Initialization
	  
	  for(i=0;i<NN;i++){
		  GABA1[i]=Synapse(5.6, 0.3, 0/180.0, -80.0, 0.0, 0.0, 0.0);
	  }	

	  for(i=0;i<NN;i++){
		  GABA2[i]=Synapse(7.6, 0.3, 0/180.0, -80.0, 0.0, 0.0, 0.0);
	  }	
      
	  Network L10IpcImcSys(Mat21, Mat31, Mat41, Mat12, Mat23, Mat14, L10, Ipc, ImcT1, ImcT2, AMPA, NMDA, GABA1, GABA2, NN);
		//Network(double **_WmatTp21, double  **_WmatTp31, double **_WmatTp41, double **_WmatTp12, double **_WmatTp23, double **_WmatTp14, Neuron *_NeGrp1, Neuron *_NeGrp2, Neuron *_NeGrp3, Neuron *_NeGrp4, Synapse *_SypGrp1,
		//Synapse *_SypGrp2, Synapse *_SypGrp3, Synapse *_SypGrp4, int _NeNo);

	 
	  for(i=0;i<iters+1000;i++){
		  L10f<<i*dt<<"\t";
		  Ipcf<<i*dt<<"\t";
		  voltf<<i*dt<<"\t";
		  //Imc1f<<i*dt<<"\t";
		  //Imc2f<<i*dt<<"\t";
		  //Iamp1=(AMP1+0.0*fabs(gsdev(seed)))*(i*dt/(dt*iters*1/3)*(i<iters/3)+(i>=iters/3&&i<=iters*2/3)+(i>iters*2/3)*(iters-i)*dt/(dt*iters*1/3));

		  if(i<500)
		  isspk=L10IpcImcSys.NtWk_4G_CN_Updat(dt, Iamp0, SIG1, ISPOT1, Iamp2, SIG2, ISPOT2, randn, Ipcf, L10f, Imc1f, Imc2f, voltf);		//single gaussian input with width SIG1
		  else if(i>=500&&i<=iters+500)

		  isspk=L10IpcImcSys.NtWk_4G_CN_Updat(dt, Iamp1, SIG1, ISPOT1, Iamp2, SIG2, ISPOT2, randn, Ipcf, L10f, Imc1f, Imc2f, voltf);		//single gaussian input with width SIG1
		  else	  
		  isspk=L10IpcImcSys.NtWk_4G_CN_Updat(dt, Iamp0, SIG1, ISPOT1, Iamp2, SIG2, ISPOT2, randn, Ipcf, L10f, Imc1f, Imc2f, voltf);		//single gaussian input with width SIG1		  

		  //if(i<500)
		  //isspk=L10IpcImcSys.NtWk_4GUpdatDAC(dt, Iamp0, SIG1, ISPOT1, Iamp2, SIG2, ISPOT2, randn, Ipcf, L10f, Imc1f, Imc2f, voltf);		//single gaussian input with width SIG1
		  //else if(i>=500&&i<=iters+500)

		  //isspk=L10IpcImcSys.NtWk_4GUpdatDAC(dt, Iamp1, SIG1, ISPOT1, Iamp2, SIG2, ISPOT2, randn, Ipcf, L10f, Imc1f, Imc2f, voltf);		//single gaussian input with width SIG1
		  //else	  
		  //isspk=L10IpcImcSys.NtWk_4GUpdatDAC(dt, Iamp0, SIG1, ISPOT1, Iamp2, SIG2, ISPOT2, randn, Ipcf, L10f, Imc1f, Imc2f, voltf);		//single gaussian input with width SIG1
//update neurons with ADP
		  Input<<Iamp1(0)<<"\t"<<Iamp1(1)<<"\n";
		  L10f<<"\n";
		  Ipcf<<"\n";
		  voltf<<"\n";
		  //Imc1f<<"\n";
		  //Imc2f<<"\n";
	  }
      
		delete [] (Mat21[0]);
		delete [] (Mat21);

		delete [] (Mat31[0]);
		delete [] (Mat31);

		delete [] (Mat41[0]);
		delete [] (Mat41);
		
		delete [] (Mat12[0]);
		delete [] (Mat12);

		delete [] (Mat23[0]);
		delete [] (Mat23);

		delete [] (Mat14[0]);
		delete [] (Mat14);
		L10f.close();
		voltf.close();
		Ipcf.close();
		Imc1f.close();
		Imc2f.close();
		paraout.close();
		Input.close();
		wmatf.close();

	 
   }