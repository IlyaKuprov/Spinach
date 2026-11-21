#pragma once
#include "IOSpinSystem.h"

using namespace SpinachGUI;


/** \class  Orca
*  Orca class, inherited from  IOSpinSystem for reading Orca format
*  Atoms data, Chemical shifts, J-couplings, Hyperfine couplings, 
*  g-tensors, zero-field couplings, framed
*/
public ref class Orca:IOSpinSystem 
{
public:
	Orca(String^, Model^% );
	virtual bool LoadFile(void) override;

private:
	bool read_onceOR;   ///<Flag for checking if multiple frames is being read for coordinates.
    bool read_onceCS;	///<Flag for checking if multiple frames is being read for Chemical shifts.
    bool read_onceHFC_QUAD;///<Flag for checking if multiple frames is being read for hyperfine couplings.
	bool read_onceGT;   ///<Flag for checking if multiple frames is being read for g-tesnors.
	InteractionsDictionary^ pendingProperties; ///< Tensors from the most recent hyperfine/EFG section.
	bool read_onceZF;   ///<Flag for checking if multiple frames is being read for zero-field couplings.
	bool read_CARTCOORD();
	bool read_HFC_QUAD_COUPLINGS();
	bool read_GTENSOR();
	bool read_CHEMICAL_SHIFTS();
	bool read_ZERO_FIELD();

};

