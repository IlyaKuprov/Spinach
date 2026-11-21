#pragma once
#include "IOSpinSystem.h"

using namespace SpinachGUI;


/** \class  ADF
*  ADF class, inherited from  IOSpinSystem for reading ADF format
*  Atoms data, Chemical shields, J-couplings, framed
*/
public ref class ADF:IOSpinSystem 
{
public:
	ADF(String^, Model^% );
	virtual bool LoadFile(void) override;

private:
	bool read_onceOR;   ///<Flag for checking if multiple frames is being read for coordinates.
    bool read_onceCS;	///<Flag for checking if multiple frames is being read for Chemical shifts.
    bool read_onceJC;	///<Flag for checking if multiple frames is being read for J-couplings.
	bool read_CARTCOORD();
	bool read_J_COUPLINGS();
	bool read_CHEMICAL_SHIELDING();
};
