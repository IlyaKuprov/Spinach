#pragma once
#include "IOSpinSystem.h"

using namespace SpinachGUI;


/** \class  Gauss
*  Gauss class, inherited from  IOSpinSystem for reading GaussView3 format
*  Atoms data, Chemical shifts, J-couplings, Hyperfine couplings, 
*  g-tensors, quadrupolar couplings, framed.
*/
public ref class Gauss:IOSpinSystem 
{
public:
	Gauss(String^, Model^% );
	virtual bool LoadFile(void) override;

private:
	bool read_COORDINATES();
	bool read_ISOTOPES1();
	bool read_JCOUPLINGS();
	bool read_CSTENSOR();
	bool read_ISOCOUPLINGS();
	bool read_ANISCOUPLINGS();
	bool read_GTENSOR();
	bool read_QUADCOUPLINGS();
	bool read_ISOTOPES2();
	bool read_SUSCEPTIBILITY();

	int AtomCount;		///<Number of Atoms.
	int startIndexANIS; ///<Start index for anisotropic couplings.
	bool read_onceJC;	///<Flag for checking if multiple frames is being read for j-couplings.
	bool read_onceISO;	///<Flag for checking if multiple frames is being read for isotropic factors.
	bool read_onceANI;	///<Flag for checking if multiple frames is being read for anisotropic factors.
	bool read_onceSUS; ///<Flag for checking if the magnetic susceptibility has been read.
	bool read_onceCO;   ///<Flag for checking if multiple frames is being read for coordinates.
	bool read_onceCS;	///<Flag for checking if multiple frames is being read for Chemical shifts.   
	bool read_onceGT;   ///<Flag for checking if multiple frames is being read for g-tesnors.
	bool read_onceQC;   ///<Flag for checking if multiple frames is being read for quadropolar couplings.

	ArrayList^ IsotropicCollection;	///<Array that saves the isotropic factors.
	ArrayList ^Elements, ^Mass;		///<Arrays that save elements and mass numbers.
};

