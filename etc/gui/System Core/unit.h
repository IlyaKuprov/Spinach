#pragma once

using namespace System;

/** \class  Unit
*  Unit class for Units API. Original by Gareth Charnock. 
*  The general approach to units of a quantity is stored in a
*  double, then that quantity is in SI. Quantities stored in a double or
*  dimension double can be stored as a user defined unit type. Assigning
*  or implicitly converting a double to a double should result in an SI
*  quantity being assigned to the double (and visa versa with the custom
*  unit type being assigned to the double)
*/
public ref class Unit {
public:
	/**
	* @brief Constructor for unit.
	* @param name Name of unit.
	* @param toSIf Multiply factor to convert to SI units.
	*/
	Unit(String^ name, double toSIf) : mName(name), mToSI(toSIf),mFromSI(1.0f/mToSI) {};

	/**
	* @brief Getter for name.
	* @result Returns the name as string.
	*/
	String^  get_name(){return mName;}

	/**
	* @brief Operator for when converted to string, returns the name.
	* @result Returns the name as string.
	*/
	inline operator const String ^ () {return mName;}

	/**
	* @brief Function for convertin from non SI units(this one) to SI units.
	* @param nonSI Value in non SI unit(this one).
	* @result Returns the value in SI unit.
	*/
	inline double toSI(double nonSI)  {return mToSI*nonSI;}

	/**
	* @brief Function for convertin from SI units to non SI units(this one).
	* @param SI Value in  SI unit.
	* @result Returns the value in non SI unit(this one).
	*/
	inline double fromSI(double SI)  {return mFromSI*SI;}

	/**
	* @brief Getter for factor to SI units.
	* @result Returns the factor to SI units.
	*/
	inline double get_to_si()  {return mToSI;}

	/**
	* @brief Getter for factor from SI units.
	* @result Returns the factorfrom SI units.
	*/
	inline double get_from_si()  {return mFromSI;}

private:
	String^ mName;   ///< String name of the unit
	double mToSI;	 ///< Factor to convert to SI units
	double mFromSI;  ///< Factor to convert from SI units
};



/** \class  Units
* Collection of static units classes where there defined spesicific units
*/
public ref class Units  
{
public:

	//Unknown units
	static Unit^ Unitless = gcnew Unit(L"Unitless",1.0f);
	static Unit^ Unknown= gcnew Unit(L"Unknown",1.0f);

	//Magnetic fluc
	static Unit^ Gauss = gcnew Unit(L"Gauss",1.0f);
	static Unit^ mTesla = gcnew Unit(L"mTesla",10.0f);
	static Unit^ Mi_b = gcnew Unit(L"Bohr magneton",1.0f);
	static Unit^ ppm = gcnew Unit(L"ppm",1.0f);

	//Energy
	static Unit^ Joule = gcnew Unit("Joule",1.0);
	static Unit^ Hz = gcnew Unit("Hz",6.626068e-34);
	static Unit^ rev_cm = gcnew Unit("cm^-1",1.9863e-23);
	static Unit^ KHz = gcnew Unit("KHz",6.626068e-31);
	static Unit^ MHz = gcnew Unit("MHz",6.626068e-28);
	static Unit^ GHz = gcnew Unit("GHz",6.626068e-25);
	static Unit^ eV = gcnew Unit("eV",1.60217646e-19);

	//Time
	static Unit^ seconds = gcnew Unit("Seconds",1.0);
	static Unit^ nsec = gcnew Unit("nsec",1.0e-9);
	static Unit^ mi_sec = gcnew Unit("mi_sec",1.0e-6);

	//Length
	static Unit^ metres = gcnew Unit("Metres",1.0);
	static Unit^ milimetre = gcnew Unit("Milimetre",1e-3);
	static Unit^ micrometre = gcnew Unit("Micrometre",1e-6);
	static Unit^ nanometre = gcnew Unit("Nanometre",1e-9);
	static Unit^ Angstroms = gcnew Unit("Angstroms",1e-10);
	static Unit^ BohrRadius = gcnew Unit("Bohr Radius",52.9177e-12);
};

