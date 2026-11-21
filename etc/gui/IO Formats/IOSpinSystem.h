#pragma once
#include "../System Core/Model.h"
#include "../System Core/isotopes.h"
#include "../System Core/orientation.h"


using namespace System;
using namespace System::IO;
using namespace System::Text;

#define PI           3.141592653589793238462643

namespace SpinachGUI {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	//Numbers in quantum chemistry output use a decimal point whatever the Windows locale
	static inline bool TryParseNumber(System::String^ s, double% x)
	{
		return System::Double::TryParse(s, System::Globalization::NumberStyles::Float|System::Globalization::NumberStyles::AllowThousands,
			System::Globalization::CultureInfo::InvariantCulture, x);
	}

	static String^ ConvertToString(double aNumber, int decDigits, String^ format )
	{
		System::Globalization::NumberFormatInfo^ provider;
		System::Globalization::CultureInfo^ MyCI = gcnew System::Globalization::CultureInfo( "en-US",false );
		provider = MyCI->NumberFormat;

		//number of Decimal digits
		provider->NumberDecimalDigits = decDigits;  

		//seperator for more than a 3 digits
		provider->NumberGroupSeparator = ""; 

		return 	((Double)(aNumber)).ToString(format, provider);
	}

	/**
	* @brief Function for converting DCM to euler angles.
	* @param eigenvec The matrix with the DCM(eigen vectors).
	* @result The resulted Euler angles.
	*/
	static EulerAngles^ EulerAnglesfromEigenVectors(ObservableCollection<double>^ eigenvec)
	{
		SpinachGUI::Matrix3x3^ AA=gcnew SpinachGUI::Matrix3x3();

		//Calculate the determinant of DCM
		double detA= eigenvec[0]*eigenvec[4]*eigenvec[8]
		-eigenvec[0]*eigenvec[5]*eigenvec[7]
		-eigenvec[3]*eigenvec[1]*eigenvec[8]
		+eigenvec[3]*eigenvec[2]*eigenvec[7]
		+eigenvec[6]*eigenvec[1]*eigenvec[5]
		-eigenvec[6]*eigenvec[4]*eigenvec[2];

		//Translate DCM multiplied by its determinant to a Matrix3x3 matrix
		AA->xx=detA*eigenvec[0];AA->xy=detA*eigenvec[1];AA->xz=detA*eigenvec[2];
		AA->yx=detA*eigenvec[3];AA->yy=detA*eigenvec[4];AA->yz=detA*eigenvec[5];
		AA->zx=detA*eigenvec[6];AA->zy=detA*eigenvec[7];AA->zz=detA*eigenvec[8];

		//Setup the initial Orientation object
		SpinachGUI::Orientation^ nOrient= gcnew SpinachGUI::Orientation(AA);

		//Convert to Euler and return
		return nOrient->GetAsEuler();
	}

	/**
	* @brief Function for converting tensor to the lab frame.
	* @param inter Interaction we like to convert.
	* @result The resulted interaction matrix.
	*/
	static ::Matrix3x3^ TransformToLab(Tensor ^ inter)
	{
		Matrix3d A;
		A<<inter->matrix3x3[0], inter->matrix3x3[1], inter->matrix3x3[2],
			inter->matrix3x3[3], inter->matrix3x3[4], inter->matrix3x3[5],
			inter->matrix3x3[6], inter->matrix3x3[7], inter->matrix3x3[8];

		Matrix3d refFrame;

		//Multiplied the reference frames matrices with the interaction matrix
		//Actually Rotation matrix*Interaction Matrix*Transpose rotation matrix
		for(ReferenceFrame^ next=inter->Frame;next!=nullptr;next=next->ParentRefFrame)
		{
			//Convert to Eigen Matrix
			refFrame<<next->matrix3x3[0],next->matrix3x3[1],next->matrix3x3[2],
				next->matrix3x3[3],next->matrix3x3[4],next->matrix3x3[5],
				next->matrix3x3[6],next->matrix3x3[7],next->matrix3x3[8];
			A=refFrame*A*refFrame.transpose();
		}

		//Convert to Matrix3x3 matrix
		::Matrix3x3^ result=gcnew ::Matrix3x3();
		for(int i=0;i<9;i++) result[i]=A(i);

		return result;
	}


	/** \class IOSpinSystem
	* IOSpinSystem abstract class, inherited from  windows form for creating different
	* input/output formats throu inheritance
	*/

	public ref class IOSpinSystem: public System::Windows::Forms::Form
	{
	public:
		IOSpinSystem(): ToShift(false), LineCount(0){};   ///< Constructor

		
		enum class FilterType {
			IMPORT,
			EXPORT,
			IMPORTEXPORT
		};///<Filter type

		//In this verison these parameters would not used. For the moment only for importing
		//virtual FilterType GetFilterType() const=0;
		//virtual const char* GetFilter() const=0;
		//virtual const char* GetFormatName() const=0;

		virtual void SaveFile(String^ filename){};
		virtual bool LoadFile(void){return true;};     ///< Main function for opening and reading file
		Model^ SystemModel;                ///< The  model where imported data will be saved
		bool ToShift;					   ///< Flag to show if the convertion from chemical shielding to shift will be used
		String^ FormatName;				   ///< The name of the format
		String^ FileName;				   ///< The name of the File to read or write
		StreamReader^ ReadingFile;		   ///< The stream that has the file which will be read
		StreamWriter^ WritingFile;		   ///< The stream that has the file which will be written
		int getLineStop() {return LineCount;};   ///< Function that returns the last line that was read


	protected:
		int LineCount;          ///< Line Counter
		String^ ReadLine();
		void ExtractIsotopeAndCoordinates(String ^ , String ^ , bool , int & , int , int , 
			int , int, int , int , int , double, AtomsDictionary^% );
		/**
		* @brief Overload function of ExtractIsotopeAndCoordinates without the AtomNumberIndex parameter
		* when atom number of isotope is unknown.
		*/
		void ExtractIsotopeAndCoordinates(String ^ line, String ^ Format, bool read_once, int & AtomID,
			int FieldsNumber, int IsotopeIndex, int XIndex, int YIndex, int ZIndex, double factor,
			int MassIndex, AtomsDictionary^% AtomCollection)
		{ExtractIsotopeAndCoordinates(line, Format, read_once, AtomID, FieldsNumber, 
		IsotopeIndex, -1, MassIndex, XIndex, YIndex, ZIndex,factor, AtomCollection);};
		/**
		* @brief Overload function of ExtractIsotopeAndCoordinates without the MassIndex parameter
		* when the  mass  of isotope is unknown.
		*/
		void ExtractIsotopeAndCoordinates(String ^ line, String ^ Format, bool read_once, int & AtomID,
			int FieldsNumber, int IsotopeIndex, int XIndex, int YIndex, int ZIndex, double factor,
			AtomsDictionary^% AtomCollection, int AtomNumberIndex)
		{ExtractIsotopeAndCoordinates(line, Format, read_once, AtomID, FieldsNumber, 
		IsotopeIndex, AtomNumberIndex, -1, XIndex, YIndex, ZIndex,factor, AtomCollection);};
		/**
		* @brief Overload function of ExtractIsotopeAndCoordinates without the MassIndex and AtomNumberIndex parameter
		* when the mass and Atom Number  of isotope is unknown.
		*/
		void ExtractIsotopeAndCoordinates(String ^ line, String ^ Format, bool read_once, int & AtomID,
			int FieldsNumber, int IsotopeIndex, int XIndex, int YIndex, int ZIndex, double factor,
			AtomsDictionary^% AtomCollection)
		{ExtractIsotopeAndCoordinates(line, Format, read_once, AtomID, FieldsNumber, 
		IsotopeIndex, -1, -1, XIndex, YIndex, ZIndex,factor, AtomCollection);};

		void ExtractInteractionMatrix3x3(String ^ , String ^ , Atom^ , Atom^ , double , bool , int , int & , 
			int , int , int , int , Unit ^ , InteractionKind , ReferenceFrame ^ , InteractionsDictionary^% );

	};	
	
}