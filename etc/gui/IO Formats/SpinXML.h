#pragma once
#include "IOSpinSystem.h"
#include "spinxml_schema.h"
#include "spinxml_schema_reference_frame.h"

using namespace SpinachGUI;
using namespace System::Xml;
using namespace System::Xml::Serialization;


/** \class  SpinXML
*  SpinXML class, inherited from  IOSpinSystem for reading and writing SpinXML format
*  Atoms data, Chemical shifts, J-couplings, Hyperfine couplings, exchange couplings,
*  quadropolar couplings, g-tensors, zero-field couplings, non-framed. Difference can 
*  occur to an imported and exported xml file when the interaction matrix is not 
*  symmetrical. It is converted automatically to symmetrical.
*/
ref class SpinXML:IOSpinSystem 
{
public:
	SpinXML(String^, FilterType, Model^% );
	virtual bool LoadFile(void) override;
	void WriteFile(Model^%);

private:

};




