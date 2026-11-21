#pragma once
#include "IOSpinSystem.h"

using namespace SpinachGUI;

/** \class  GAMESS
*  GAMESS class, inherited from  IOSpinSystem for reading GAMESS format
*  Just atoms data, framed
*/
public ref class GAMESS:IOSpinSystem 
{
public:
	GAMESS(String^, Model^% );
	virtual bool LoadFile(void) override;

private:
	bool read_CARTCOORD(bool bohr);
	bool read_onceOR;

};
