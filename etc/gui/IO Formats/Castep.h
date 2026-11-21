#pragma once
#include "IOSpinSystem.h"

using namespace SpinachGUI;
using namespace System::Collections::Generic;

/** \class  Castep
*  Castep class, inherited from  IOSpinSystem for reading Castep format
*  Atoms data, Chemical shields, Quadrupolar couplings, J-couplings, non framed
*/

public ref class Castep:IOSpinSystem 
{
public:
	Castep(String^, Model^%  );
	virtual bool LoadFile(void) override;

private:
	bool OldFormat(String^ );
	bool NewFormat(String^ );
	bool CheckImportAtom(cli::array<String^>^ , int & );
};

