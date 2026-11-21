#pragma once
#include "IOSpinSystem.h"

using namespace SpinachGUI;

/** \class  MOL
*  MOL class, inherited from  IOSpinSystem for reading MOL format.
*  Just atoms data, no framed.
*/
public ref class MOL:IOSpinSystem 
{
public:
	MOL(String^, Model^% );
    virtual bool LoadFile(void) override;
};