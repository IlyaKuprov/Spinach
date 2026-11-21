#pragma once

using namespace System::Collections;

/** \class  Isotope
*  Class for Isotopes using common properties.
*/
public ref class Isotope
{
public:
	/**
	* @brief Constructor for isotope.
	*/
	Isotope(System::String^ element, int atomNumber, int mass, 	int radius, float red, float green,
		float blue, float spin,	float quadrupole_moment, double magnetogyric): 
	Element(element), AtomNumber(atomNumber), Mass(mass), Radius(radius), Red(red), Green(green), 
		Blue(blue), Spin(spin), QuadrupoleMoment(quadrupole_moment), Magnetogyric(magnetogyric){}; 


	System::String^ Element;   ///< Name of element
	int Mass;				   ///< Mass numberof element
	int AtomNumber;			   ///< Atom number of element
	float Red, Green,Blue;	   ///< RBG color of element
	int Radius;				   ///< Atomic Radius of element
	float Spin;				   ///< Spin of element
	float QuadrupoleMoment;	   ///< Quadrupole moment of element
	double Magnetogyric;	   ///< Magnetogyric ratio of element
}; 

/** \class  Isotopes
*  Collection Class with static classes of Isotopes.
*/
public ref class Isotopes
{
public:
	/**
	* @brief Constructor for Isotopes.
	*/
	Isotopes(){};

	/**
	* @brief Getter of isotopes collection.
	* @result Returns the array with the collection of isotopes.
	*/
	static ArrayList^ getAllIsotopes(){Init_Isotopes();return Isotopes_Collection;};

	static Isotope^ FindIsotope(System::String^ element, int atomNumber, int mass);
	static cli::array<System::String^>^ getAllIsotopeMasses(System::String ^element);
	static Isotope^ FindIsotopeWithHigherSpin(Isotope^ isotope);

private:
	static void Init_Isotopes();
	static ArrayList^ Isotopes_Collection; ///< Collection of isotopes
}; 

/** \class  CompareIsotopes
*  Comparator class used by Binary Search of isotopes collection based on atom number and mass number.
*/
ref class CompareIsotopes: public IComparer
{
private:

	// Calls CaseInsensitiveComparer.Compare with the parameters reversed.
	virtual int Compare( Object^ x, Object^ y ) sealed = IComparer::Compare
	{
		//When the atom number is the same, compare the mass number 
		if (((Isotope^)x)->AtomNumber==((Isotope^)y)->AtomNumber)
		{
			if (((Isotope^)x)->Mass==((Isotope^)y)->Mass || ((Isotope^)y)->Mass==0)  return 0;
			else if (((Isotope^)x)->Mass<((Isotope^)y)->Mass)  return -1;
			else if (((Isotope^)x)->Mass>((Isotope^)y)->Mass)  return 1;
		}
		else if (((Isotope^)x)->AtomNumber<((Isotope^)y)->AtomNumber) return -1;
		else if (((Isotope^)x)->AtomNumber>((Isotope^)y)->AtomNumber) return 1;

		return 0; //For avoiding warnings
	}
};

/** \class  CompareIsotopesMass
*  Comparator class used by Binary Search of isotopes collection based on mass number.
*/
ref class CompareIsotopesMass: public IComparer
{
private:

	// Calls CaseInsensitiveComparer.Compare with the parameters reversed.
	virtual int Compare( Object^ x, Object^ y ) sealed = IComparer::Compare
	{
		if (((Isotope^)x)->Mass==((Isotope^)y)->Mass)  return 0;
		else if (((Isotope^)x)->Mass<((Isotope^)y)->Mass)  return -1;
		else if (((Isotope^)x)->Mass>((Isotope^)y)->Mass)  return 1;

		return 0; //For avoiding warnings
	}
};






