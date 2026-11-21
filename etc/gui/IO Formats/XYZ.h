#pragma once
#include "IOSpinSystem.h"

namespace SpinachGUI {

	/** \class  XYZ
	*  XYZ class, inherited from  IOSpinSystem for reading XYZ format
	*  Just atoms data, framed
	*/
	public ref class XYZ:IOSpinSystem 
	{
	public:
		XYZ(String ^ , Model ^% );
		virtual bool LoadFile(void) override;
	};

}