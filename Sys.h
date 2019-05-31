// The system have different inputs between L10 to Ipc and Imc neurons
// NO accomdation is considered.
#ifndef _SYS_H
#define _SYS_H
#include "MatVec.h"
#include "DataType.h"
#include "rng.h"
class Sys
{
public:
	Sys(){}
	~Sys();
	Sys(double _T, int _ISPOT1, int _ISPOT2, double _AMP1, double _AMP2, double _Anos1, double _Anos2, double _nctrl1, double _nctrl2,double _nsd1, double _nsd2, double _SIG1, double _SIG2);	
	
	double score(VecDP & VarParSet,VecDP & _L10Par, VecDP &_IpcPar, VecDP & _IpcAdapPar, VecDP & _L10AdapPar, RNG & randn);

	void disp(VecDP & VarParSet,VecDP & _L10Par, VecDP &_IpcPar, VecDP & _IpcAdapPar, VecDP & _L10AdapPar, RNG & randn);

private:
	int ISPOT1, ISPOT2;
	double AMP1, AMP2, SIG1, SIG2, Anos1, Anos2, nctrl1, nctrl2, nsd1, nsd2;
	double T;
};

#endif