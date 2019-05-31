#include "Sys.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include "DataType.h"
#include "rng.h"

using namespace std;
int main()
{
	int i=0,j=0,k=0;
	ofstream scout("scores.txt", ios::app);
	RNG randgen;
	

	double r=0;	
	double Noisd=0.18;
	VecDP wMatPar(4);	wMatPar(0)=1.0;	wMatPar(1)=50;	wMatPar(2)=1.0;	wMatPar(3)=50;	
	VecDP IpcSponPar(2); IpcSponPar(0)=0.0; IpcSponPar(1)=1.5;
	VecDP L10SponPar(2); L10SponPar(0)=0.0; L10SponPar(1)=0.1;				//default parameters
	VecDP IpcAdapPar(2); IpcAdapPar(0)=1.1; IpcAdapPar(1)=60;
	VecDP L10AdapPar(2); L10AdapPar(0)=0.6; L10AdapPar(1)=50;

	double Ipcsd=0.0, L10sd=0.0;				//Ipc spontaneous activity
	double IpcTsra=14, IpcDgsra=0.0;
	double L10Tsra=15, L10Dgsra=0.0;
	int WidthL10Ipc=0, WidthIpcL10=0;			//width 
	
	double smat[10][10];

	for(i=0;i<10;i++)
		for(j=0;j<10;j++)
			smat[i][j]=0;

	//for (k=0;k<5;k++){					//Ipc sra time const vs Ipc Dgsra
	//	IpcTsra=0;
	//	for (i=0;i<10;i++){
	//		IpcTsra+=14;				//varying the Ipc Sra time constant
	//		IpcAdapPar(1)=IpcTsra;										
	//		IpcDgsra=-0.4;	
	//		for(j=0;j<10; j++){
	//			
	//			IpcDgsra+=0.4;
	//			IpcAdapPar(0)=IpcDgsra;			//varying the Ipc Dgsra

	//			Sys IshmTec(250.0, 199, 1, 0.18, 0.0, 0.0, 0.0, 1.0, 0.0, Noisd, 1.0, 80, 100);
	//			//Sys(_T, _ISPOT1, _ISPOT2, _AMP1, _AMP2, _Anos1, _Anos2, _nctrl1,_nctrl2, _nsd1, _nsd2, _SIG1, _SIG2)
	//			r=IshmTec.score(wMatPar, L10SponPar, IpcSponPar, IpcAdapPar, L10AdapPar, randgen);
	//			smat[i][j]+=r;
	//			scout<<r<<"\t";
	//			cout<<IpcTsra<<"\t"<<IpcDgsra<<"\t"<<r<<"\t"<<smat[i][j]<<endl;
	//		}
	//		scout<<"\n";
	//	}
	//	scout<<"\n";
	//}
	

	//for(i=0;i<10;i++){
	//	for(j=0;j<10;j++){
	//		scout<<smat[i][j]<<"\t";
	//	}
	//	scout<<"\n";
	//}
	//scout.close();

	


	//for (k=0;k<5;k++){					//L10 to Ipc width varying and input noise sd varying
	//	WidthL10Ipc=0;
	//	for (i=0;i<10;i++){
	//		WidthL10Ipc+=10;				//varying the width weight matrix L10->Ipc
	//		wMatPar(1)=WidthL10Ipc;			//changing the passing parameters width of L10->Ipc
	//		Noisd=-0.03;	
	//		for(j=0;j<10; j++){
	//			
	//			Noisd+=0.03;

	//			Sys IshmTec(250.0, 199, 1, 0.18, 0.0, 0.0, 0.0, 1.0, 0.0, Noisd, 1.0, 80, 100);
	//			//Sys(_T, _ISPOT1, _ISPOT2, _AMP1, _AMP2, _Anos1, _Anos2, _nctrl1,_nctrl2, _nsd1, _nsd2, _SIG1, _SIG2)
	//			r=IshmTec.score(wMatPar, L10SponPar, IpcSponPar, IpcAdapPar, L10AdapPar, randgen);
	//			smat[i][j]+=r;
	//			scout<<r<<"\t";
	//			cout<<wMatPar(1)<<"\t"<<Noisd<<"\t"<<r<<"\t"<<smat[i][j]<<endl;
	//		}
	//		scout<<"\n";
	//	}
	//	scout<<"\n";
	//}

	//for(i=0;i<10;i++){
	//	for(j=0;j<10;j++){
	//		scout<<smat[i][j]<<"\t";
	//	}
	//	scout<<"\n";
	//}
	//scout.close();

	
	//for (k=0;k<5;k++){						//L10 to Ipc width varying and Ipc spontaneous activity varying
	//	WidthL10Ipc=0;
	//	for (i=0;i<10;i++){
	//		WidthL10Ipc+=10;				//varying the width weight matrix L10->Ipc
	//		wMatPar(1)=WidthL10Ipc;			//changing the passing parameters width of L10->Ipc
	//		Ipcsd=-0.3;	
	//		for(j=0;j<10; j++){				
	//			Ipcsd+=0.3;
	//			IpcSponPar(1)=Ipcsd;

	//			Sys IshmTec(250.0, 199, 1, 0.18, 0.0, 0.0, 0.0, 1.0, 0.0, Noisd, 1.0, 80, 100);
	//			//Sys(_T, _ISPOT1, _ISPOT2, _AMP1, _AMP2, _Anos1, _Anos2, _nctrl1,_nctrl2, _nsd1, _nsd2, _SIG1, _SIG2)
	//			r=IshmTec.score(wMatPar, L10SponPar, IpcSponPar,IpcAdapPar,L10AdapPar, randgen);
	//			smat[i][j]+=r;
	//			scout<<r<<"\t";
	//			cout<<wMatPar(1)<<"\t"<<IpcSponPar(1)<<"\t"<<r<<"\t"<<smat[i][j]<<endl;
	//		}
	//		scout<<"\n";
	//	}
	//	scout<<"\n";
	//}

	//for(i=0;i<10;i++){
	//	for(j=0;j<10;j++){
	//		scout<<smat[i][j]<<"\t";
	//	}
	//	scout<<"\n";
	//}
	//scout.close();

	//for (k=0;k<5;k++){						//L10 to Ipc width varying and Ipc to L10 Width varying
	//	WidthL10Ipc=0;
	//	for (i=0;i<10;i++){
	//		WidthL10Ipc+=10;				//varying the width weight matrix L10->Ipc
	//		wMatPar(1)=WidthL10Ipc;			//changing the passing parameters width of L10->Ipc
	//		WidthIpcL10=0;	
	//		for(j=0;j<10; j++){				
	//			WidthIpcL10+=10;
	//			wMatPar(3)=WidthIpcL10;	
	//			
	//			Sys IshmTec(250.0, 199, 1, 0.18, 0.0, 0.0, 0.0, 1.0, 0.0, Noisd, 1.0, 80, 100);
	//			//Sys(_T, _ISPOT1, _ISPOT2, _AMP1, _AMP2, _Anos1, _Anos2, _nctrl1,_nctrl2, _nsd1, _nsd2, _SIG1, _SIG2)
	//			r=IshmTec.score(wMatPar, L10SponPar, IpcSponPar, IpcAdapPar, L10AdapPar, randgen);
	//			smat[i][j]+=r;
	//			scout<<r<<"\t";
	//			cout<<wMatPar(1)<<"\t"<<wMatPar(3)<<"\t"<<r<<"\t"<<smat[i][j]<<endl;
	//		}
	//		scout<<"\n";
	//	}
	//	scout<<"\n";
	//}
	//for(i=0;i<10;i++){
	//	for(j=0;j<10;j++){
	//		scout<<smat[i][j]<<"\t";
	//	}
	//	scout<<"\n";
	//}
	//scout.close();


	//for (k=0;k<5;k++){						//L10 to Ipc width varying and L10 spontaneous activity varying
	//	WidthL10Ipc=0;
	//	for (i=0;i<10;i++){
	//		WidthL10Ipc+=10;				//varying the width weight matrix L10->Ipc
	//		wMatPar(1)=WidthL10Ipc;			//changing the passing parameters width of L10->Ipc
	//		L10sd=0.0;	
	//		for(j=0;j<10; j++){				
	//			L10sd+=0.25;
	//			L10SponPar(1)=L10sd;

	//			Sys IshmTec(250.0, 199, 1, 0.18, 0.0, 0.0, 0.0, 1.0, 0.0, Noisd, 1.0, 80, 100);
	//			//Sys(_T, _ISPOT1, _ISPOT2, _AMP1, _AMP2, _Anos1, _Anos2, _nctrl1,_nctrl2, _nsd1, _nsd2, _SIG1, _SIG2)
	//			r=IshmTec.score(wMatPar, L10SponPar, IpcSponPar, IpcAdapPar,L10AdapPar, randgen);
	//			smat[i][j]+=r;
	//			scout<<r<<"\t";
	//			cout<<wMatPar(1)<<"\t"<<L10SponPar(1)<<"\t"<<r<<"\t"<<smat[i][j]<<endl;
	//		}
	//		scout<<"\n";
	//	}
	//	scout<<"\n";
	//}
	//for(i=0;i<10;i++){
	//	for(j=0;j<10;j++){
	//		scout<<smat[i][j]<<"\t";
	//	}
	//	scout<<"\n";
	//}
	//scout.close();



		Sys IshmTec(250.0, 199, 1, 0.18, 0.0, 0.0, 0.0, 1.0, 0.0, Noisd, 1.0, 80, 100);
		//r=IshmTec.score(wMatPar,L10SponPar, IpcSponPar, IpcAdapPar,L10AdapPar, randgen);
		IshmTec.disp(wMatPar,L10SponPar, IpcSponPar, IpcAdapPar,L10AdapPar, randgen);
		



}