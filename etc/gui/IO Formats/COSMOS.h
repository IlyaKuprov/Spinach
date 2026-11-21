#pragma once
#include "IOSpinSystem.h"


using namespace SpinachGUI;


/** \class  COSMOS
*  COSMOS class, inherited from  IOSpinSystem for reading COSMOS format
*  Atoms data, Chemical shields, J-couplings, non framed
*/
public ref class COSMOS:IOSpinSystem 
{
public:
	COSMOS(String^, Model^%  );
	virtual bool LoadFile(void) override;
};

