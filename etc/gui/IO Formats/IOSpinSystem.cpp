#include "../StdAfx.h"
#include "IOSpinSystem.h"


/**
*\file IOSpinSystem.cpp
*/

/**
* @brief Reads a line from the stream  ReadingFile keeping track which line reads. 
* @return The read string. Throws string exception.
*/
String^ SpinachGUI::IOSpinSystem::ReadLine()
{
	String ^ result;
	LineCount++;
	
	//Try to read the line and throw exception if something wrong
	try { result=ReadingFile->ReadLine();}
	catch(Exception ^ e) {throw e;};
	return result;
}

/**
* @brief Reads a line, split them to a number of fields which must be equal to FieldsNumber
* and extracts the isotope and coordinates adding to the AtomCollection.
* Can used for multiple frames, but the last frame of data is saved.
* @param line String wih a line of the file.
* @param Format Name of format in order to have it in the exception-error dialogs.
* @param read_once Flag showing that atom data are read more than once.
* @param AtomID Atom ID used when replacing an atom data, not during adding new atoms.
* @param FieldsNumber Number of fields that supposinlgy this format have for the atoms
* data. If not , an exception-errror will be thrown.
* @param IsotopeIndex Zero-Index number which shows the field of the name of isotope.
* @param AtomNumberIndex  Zero-Index number which shows the field of the atom number of isotope.
* @param MassIndex  Zero-Index number which shows the field of the mass of isotope.
* @param XIndex Zero-Index number which shows the field of X coordinate.
* @param YIndex Zero-Index number which shows the field of Y coordinate.
* @param ZIndex Zero-Index number which shows the field of Z coordinate.
* @param factor A multiplying factor for the coordinates.
* @param AtomCollection Collection where the data of the atom will be saved. 
*/
void SpinachGUI::IOSpinSystem::ExtractIsotopeAndCoordinates(String ^ line, String ^ Format, bool read_once, int & AtomID,
								  int FieldsNumber, int IsotopeIndex, int AtomNumberIndex, int MassIndex, int XIndex, int YIndex, int ZIndex,
								  double factor, AtomsDictionary^% AtomCollection)
{
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	

	//Split the line to a number of fields
	numbers = nullptr;
	numbers=line->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

	if(FieldsNumber==numbers->Length) 
	{
		try {
			Atom^ nAtom=gcnew Atom();
			double x,y,z;

			//Try to read coordinates
			if (!SpinachGUI::TryParseNumber(numbers[XIndex], x) || !SpinachGUI::TryParseNumber(numbers[YIndex], y) ||
				!SpinachGUI::TryParseNumber(numbers[ZIndex], z))
				throw gcnew Exception("Error in Atom coordinates");
			else 
			{
				//Save the new atom and coordinates in the model
				nAtom->X=factor*x;nAtom->Y=factor*y;nAtom->Z=factor*z;

				//Can throw exception if does not find the isotope
				if(MassIndex==-1 && AtomNumberIndex==-1) 
					nAtom->isotope=Isotopes::FindIsotope(numbers[IsotopeIndex],0, 0);
				else if(MassIndex==-1) 
					nAtom->isotope=Isotopes::FindIsotope(numbers[IsotopeIndex],
					 Convert::ToInt32(Convert::ToDouble(numbers[AtomNumberIndex], System::Globalization::CultureInfo::InvariantCulture)), 0);
				else if(AtomNumberIndex==-1) 
					nAtom->isotope=Isotopes::FindIsotope(numbers[IsotopeIndex],0,
					 Convert::ToInt32(Convert::ToDouble(numbers[MassIndex], System::Globalization::CultureInfo::InvariantCulture)));
				else
					nAtom->isotope=Isotopes::FindIsotope(numbers[IsotopeIndex],
					 Convert::ToInt32(Convert::ToDouble(numbers[AtomNumberIndex], System::Globalization::CultureInfo::InvariantCulture)), 
					 Convert::ToInt32(Convert::ToDouble(numbers[MassIndex], System::Globalization::CultureInfo::InvariantCulture)));


				//If this is the first time reading this section(read_once==false)
				// just add it in the system else replace the previous frame value(if exists)
				if(read_once) 
					try{AtomCollection[AtomID]=nAtom; AtomID++;}
				catch(String ^ e)  {throw gcnew Exception("Problem in Elements: "+e);}
				else AtomCollection->Add(nAtom);
			};
		}
		catch(Exception ^ e)
		{	
			//Wrong format of the field, where coordinates not giving number and where 
			//isotope not giving valid isotope name
			throw gcnew Exception("Problem in Reading "+Format+ " format in Line "+LineCount+":"+e->Message) ;
		};
	}
	else
	{
		//Wrong number of fields read
		throw gcnew Exception("Problem in Reading Atom and Coordinates in Line "+LineCount+".") ;
	};


}

/**
* @brief Reads three lines, split them to a number of fields which must be equal to FieldsNumber
* (for each line) and extracts the interaction matrix creatting the interaction and 
* adding to the InteractionCollection. Can used for multiple frames, but the last frame of data is saved.
* @param line String of the first line(of three) that contains the matrix of the interaction.
* @param InterName Name of interaction in order to have it in the exception-error dialogs.
* @param A First atom of the interaction.
* @param B Second atom of the interaction.
* @param factor Multiplying factor for the elements of the matrix.
* @param read_once Flag showing that atom data are read more than once.
* @param startIndex Staring index of this kind of interaction. used when multiple frames.
* @param InterID Interaction ID used when replacing an interaction data, not during adding new interaction.
* @param FieldsNumber Number of fields that supposinlgy this format have for the interactions
* data. If not , an exception-errror will be thrown.
* @param M1Index Zero-Index number which shows the field of first element of the matrix.
* @param M2Index Zero-Index number which shows the field of second element of the matrix.
* @param M3Index Zero-Index number which shows the field of third element of the matrix.
* @param unit Unit of the imported interaction matrix data.
* @param IntKind Interaction kind of the imported interaction matrix data.
* @param RefFrame Reference frame of the imported interaction matrix data.
* @param InteractionCollection Collection where the data of the interaction will be saved. 
*/
void SpinachGUI::IOSpinSystem::ExtractInteractionMatrix3x3(String ^ line, String ^ InterName, Atom^ A, Atom^ B, double factor,
			bool read_once, int startIndex, int & InterID, int FieldsNumber, int M1Index, int M2Index, int M3Index, 
			Unit ^ unit, InteractionKind IntKind, ReferenceFrame ^ RefFrame, InteractionsDictionary^% InteractionCollection)
{
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,=";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	
	::Matrix3x3 ^matrix=gcnew ::Matrix3x3(); //3x3 matrix with the interaction
	double x,y,z;                       // For reading the tree element of the matrix
	double test=0.0;


	//               2.0010305    0.0002797    0.0000199
	//			    0.0001454    2.0020587    0.0000201
	//				-0.0001453    0.0004632    2.0045716
	for(int j=0; j<3; j++)
	{
		//The first line is already read, so do it only for hte next two
		if(j>0) line=ReadLine() ;

		//Sometimes the umber have "D" for double precision, so change it to E
	    line=line->Replace(L'D',L'E');

		//Split the line to a number of fields
		numbers = nullptr;
		numbers=line->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

		//Try to read elements and save in the matrix
		if (FieldsNumber!=numbers->Length ||
			!SpinachGUI::TryParseNumber(numbers[M1Index], x) ||
			!SpinachGUI::TryParseNumber(numbers[M2Index], y) ||
			!SpinachGUI::TryParseNumber(numbers[M3Index], z))
			throw gcnew Exception("Problem in Reading "+InterName+ " matrix elements in Line "+LineCount+". " +
			"Either the number of fields in the line is incorrect or "+
			"the format of fields(double) is incorrect.");
		else {matrix[j*3+0]=factor*x;matrix[j*3+1]=factor*y;matrix[j*3+2]=factor*z; test+=abs(x)+abs(y)+abs(z);};
	};

	//If Quadrupolar interactions smaller or equal than zero do not add;
	if(IntKind==InteractionKind::Quadrupolar && test<=0) return;

	//Try to create new Interaction
	Tensor^ nInteraction;
	try{
		nInteraction=gcnew Tensor(A, B,	IntKind, matrix);
		nInteraction->unit=unit;
		nInteraction->Frame=RefFrame;
	}//Throws when an invalid interaction is tried to create.
	catch(Exception ^ e)  
	{ 
		throw gcnew Exception("Problem in Reading Interaction: "+e->Message
			+". Problem in Reading File in Line:"+getLineStop()+" .");
	}
	//If this is the first time reading this section(read_once==false)
	// just add it in the system else replace the previous frame value(if exists)
	if(read_once) 
		try{InteractionCollection[InterID+startIndex]=nInteraction; InterID++;}
	catch(String ^ e)  {throw gcnew Exception("Problem in Elements: "+e);}
	else InteractionCollection->Add(nInteraction);

};