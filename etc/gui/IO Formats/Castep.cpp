#include "../StdAfx.h"
#include "Castep.h"

using namespace System::Collections::Generic;

//Definitions
#define InteractionCollection  SystemModel->InteractionCollection
#define AtomCollection  SystemModel->AtomCollection


/**
*\file Castep.cpp
*/

/**
* @brief Constructor of Castep format
* @param filename Name of the file to be read.
* @param model The model where the imported data will be saved.
*/
Castep::Castep(System::String^ filename, Model^% model)
{
	SystemModel=model;
	FileName=filename;
	ReadingFile = gcnew StreamReader(filename);
	FormatName="CASTEP";
}

/**
* @brief Main function for opening and reading file
* @return Returns false if the format does not seem to be correct.
*/
bool Castep::LoadFile(void)
{
	String^ fileline=gcnew String("");
	bool  OK = true;

	//Turn off Refreshing for faster importing
	AtomCollection->Refreshing=false;


    //read the first line to understand if it is the new format or the old one
	fileline=ReadLine();
	try
	{
		if (fileline->Contains("#$magres-abinitio-v1.0")) 
			OK=NewFormat( fileline); //New CASTEP format
		else OK=OldFormat(fileline); //Old CASTEP Format
	}
	catch(String^ e) 
	{
		throw gcnew Exception("Problem in Reading File in Line "+LineCount+". " +e);
	}


	//Turn on again Refreshing which causes calculation of bonds
	AtomCollection->Refreshing=true;
	ReadingFile->Close();
	return true;
}

/**
* @brief Main function for reading file of old format
* @return Returns false if the format does not seem to be correct.
*/
bool Castep::OldFormat(String^ firstline)
{
	int number; 						//number of fields in each line
	String^ fileline=firstline;         //line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array
	//Dictionary for saving quadropolar constants for specific elements
	Dictionary <String^, double>^ quadrupolecon=gcnew Dictionary<String^, double>();

	//Read until end of file
	while (!ReadingFile->EndOfStream )
	{
		//Bug: case when 2 line have =====
		if(fileline->Contains("============")) 

		{  
			fileline=ReadLine();

			//Reading of quadropolarconstants
			if(fileline->Contains("electric quadrupole moment"))
			{
				while(!(fileline=ReadLine())->Contains("============"))
				{
					//Read the quardopolar constant 
					number=7; 
					numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
					if(number!=numbers->Length) throw gcnew Exception("Problem in quadrupole moment");

					//Try to add quadropolar constant in the dictionary else throw exception
					try 
					{
						quadrupolecon->Add(numbers[0], Convert::ToDouble(numbers[5], System::Globalization::CultureInfo::InvariantCulture));
					}
					catch (ArgumentException^)
					{
						throw gcnew Exception("Quadrupole constant declared twice");
					}

				}
			}
			else if(fileline->Contains("Atom"))
			{
				//Check if an Atom is presented
				number=3; //Number of Strings in the Line describing an Atom
				numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
				if(number!=numbers->Length) throw gcnew Exception("Problem in Reading Atom");
				

				fileline=ReadLine();//============

				//Read Atom and coordinates
				number=7; //Number of Strings in the Line describing an Atom and coordinates
				//H        1 Coordinates      3.499    0.756   -0.603   A
				fileline=ReadLine();
				numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
				if(number!=numbers->Length) throw gcnew Exception("Problem in Reading Atom Coordinates");

				//Check if the atom is already imported
				int AtomID;
				if(!CheckImportAtom(numbers, AtomID)) throw gcnew Exception("Problem in Interaction Atom");
	
				fileline=ReadLine();//Empty line

				//Read Interaction 
				fileline=ReadLine(); //TOTAL Shielding Tensor
				number=2; //Number of Strings in the Line describing kind of interaction
				numbers = nullptr;
				numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

				//Check if number of fields are coorect and if the keyworf TOTAL is  is present
				if(number>numbers->Length || !fileline->Contains("TOTAL"))
				throw gcnew Exception("Problem in Reading Interaction Kind");

				//Check if it is hemical shielding or Quadropolar
				InteractionKind Intkind;
				if(numbers->Length>number && numbers[1]=="Shielding" && numbers[2]=="Tensor") 
				{Intkind=InteractionKind::CShielding;ToShift=true;}
				else if (numbers[1]=="tensor") Intkind=InteractionKind::Quadrupolar;
				fileline=ReadLine();  //Empty line

				//Read quadrupolar interactions
				if(Intkind==InteractionKind::Quadrupolar)
				{
					//Quadropolar factors
					double factor=2.3496e+02;//(efg_atomic*e_charge*1e-28/h_planck)MHz
					double qconst;
					
					//Find nuclear spin quantum number and in case smaller 
					//or equal than 0.5 find other isotope 
					if (((Atom^)AtomCollection[AtomID])->isotope->Spin<=0.5)	
						((Atom^)AtomCollection[AtomID])->isotope=
						Isotopes::FindIsotopeWithHigherSpin(((Atom^)AtomCollection[AtomID])->isotope);
					double NSQN=((Atom^)AtomCollection[AtomID])->isotope->Spin; //nuclear spin quantum number
					
					//Try to find the quadrupolar constant for this element
					try
					{
						qconst=quadrupolecon[((Atom^)AtomCollection[AtomID])->isotope->Element];
					}
					catch (KeyNotFoundException^)
					{
						throw gcnew Exception("Quadrupole constant for "+
							((Atom^)AtomCollection[AtomID])->isotope->Element+" did not declared") ;
					}

					//Reading interaction matrix and creating quadrupolar interaction
					fileline=ReadLine();
					int j;
					//If isotope with no 0.5>spin just ignore
					try
					{
						//Number of Fields: 3, XYZ index: 0,1,2
						ExtractInteractionMatrix3x3(fileline, "Quadrupolar Coupling", 
							(Atom^)AtomCollection[AtomID], (Atom^)AtomCollection[AtomID], 
							factor*qconst/(2*NSQN*(2*NSQN-1)),
							false, 0, j, 3, 0, 1, 2, 
							Units::MHz, Intkind, 
							SystemModel->RefFrameCollection[1], 
							InteractionCollection);
					}
					catch (String^ e) {/*Continue*/};

				}//Read Chemical shield interactions
				else
				{
					//Reading interaction matrix and creating Chemical shield interaction
					fileline=ReadLine();
					int j;
					//Number of Fields: 3, XYZ index: 0,1,2, no factor
					ExtractInteractionMatrix3x3(fileline, "Chemical shielding", 
						(Atom^)AtomCollection[AtomID], (Atom^)AtomCollection[AtomID], 
						1.0,
						false, 0, j, 3, 0, 1, 2, 
						Units::ppm, Intkind, 
						SystemModel->RefFrameCollection[1], 
						InteractionCollection);
				}

				//8 Empty lines
				for (int i=0;i<8;i++) fileline=ReadLine();
			};
		};
		fileline=ReadLine(); //Last thing is read the next line
	};
	return true;
}

/**
* @brief Main function for reading file of new format
* @return Returns false if the format does not seem to be correct.
*/
bool Castep::NewFormat(String^ firstline)
{
	int number;							//number of fields in each line
	String^ fileline=firstline;         //line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,	";          //seperator of fields, includes "tab"
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array

	//Read the file until the end of it 
	while (!ReadingFile->EndOfStream)
	{
		if(fileline->Contains("<atoms>")) //Atom Reading
		{  
			fileline=ReadLine(); //units lattice Angstrom
			fileline=ReadLine(); //lattice    6.0000000000000009E+00   0.0000000
			fileline=ReadLine(); //units atom Angstrom

			while(!(fileline=ReadLine())->Contains("</atoms>"))
			{
				//atom H       H       1        -8.6634114E-01       -1.722893E-01      -6.46057010E-01
				number=7; //Number of Strings in the Line describing an Atom
				numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
				if(number!=numbers->Length)  throw gcnew Exception("Problem in Reading Atom") ;

				//Add atom in the collection
				Atom^ nAtom=gcnew Atom();
				nAtom->X=Convert::ToDouble(numbers[4], System::Globalization::CultureInfo::InvariantCulture);
				nAtom->Y=Convert::ToDouble(numbers[5], System::Globalization::CultureInfo::InvariantCulture);
				nAtom->Z=Convert::ToDouble(numbers[6], System::Globalization::CultureInfo::InvariantCulture);
				nAtom->isotope=Isotopes::FindIsotope(numbers[1],0, 0);

				//Used for recognition of the interaction
				nAtom->Label=numbers[2]+numbers[3];

				AtomCollection->Add(nAtom);
			}
		}
		else if(fileline->Contains("<magres>")) //Interactions Reading
		{
			while(!(fileline=ReadLine())->Contains("</magres>"))	
			{
				if(fileline->Contains("units"))  //Reading Interaction Units
				{
					//Do something with this bloody units
				}
				else if (fileline->Contains("ms") || fileline->Contains("efg") || fileline->Contains("isc "))
				{
					//Read Interaction 
					// ms H   1  3.2E+01  1.2E+00  3.8E+00  1.9E+00   2.7E+01  2.49E+00   4.1E+00  2.2E+00  3.0E+01

					InteractionKind Intkind;
					if(fileline->Contains("ms"))
					{
						Intkind=InteractionKind::CShielding;
						number=12; //Number of Strings in the Line describing ms interaction
						ToShift=true;
					}
					else if(fileline->Contains("efg"))
					{
						Intkind=InteractionKind::Quadrupolar;
						number=12; //Number of Strings in the Line describing efg interaction
					}
					else if(fileline->Contains("isc "))
					{
						Intkind=InteractionKind::Jcoupling;
						number=14; //Number of Strings in the Line describing isc interaction
					}
					else
					{
						throw gcnew Exception("Problem in Reading Interaction Kind:Unknown Interaction") ;
					};


					//Reads the line with the interaction matrix
					numbers = nullptr;
					numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
					if(number!=numbers->Length)
						throw gcnew Exception("Problem in Reading Interaction Kind: Wrong line Format") ;

					//Put the data in the matrix
					::Matrix3x3^ matrix=gcnew ::Matrix3x3();
					int jj=(numbers->Length==12?0:2);
					matrix[0]=Convert::ToDouble(numbers[jj+3], System::Globalization::CultureInfo::InvariantCulture);matrix[1]=Convert::ToDouble(numbers[jj+4], System::Globalization::CultureInfo::InvariantCulture);
					matrix[2]=Convert::ToDouble(numbers[jj+5], System::Globalization::CultureInfo::InvariantCulture);matrix[3]=Convert::ToDouble(numbers[jj+6], System::Globalization::CultureInfo::InvariantCulture);
					matrix[4]=Convert::ToDouble(numbers[jj+7], System::Globalization::CultureInfo::InvariantCulture);matrix[5]=Convert::ToDouble(numbers[jj+8], System::Globalization::CultureInfo::InvariantCulture);
					matrix[6]=Convert::ToDouble(numbers[jj+9], System::Globalization::CultureInfo::InvariantCulture);matrix[7]=Convert::ToDouble(numbers[jj+10], System::Globalization::CultureInfo::InvariantCulture);
					matrix[8]=Convert::ToDouble(numbers[jj+11], System::Globalization::CultureInfo::InvariantCulture);

					//If J-coupling just use the 1/3 trace of matrix
					if(Intkind==InteractionKind::Jcoupling)
					{
						matrix[0]=(matrix[0]+matrix[4]+matrix[8])/3.0f;
						matrix[4]=matrix[0];
						matrix[8]=matrix[0];
						matrix[1]=0.0f;
						matrix[2]=0.0f;
						matrix[3]=0.0f;
						matrix[5]=0.0f;
						matrix[6]=0.0f;
						matrix[7]=0.0f;
					};

					Tensor^ nInteraction;
					int AtomID=0,AtomID2=0;

					//Find the Atom ID from the label
					for each(int i in AtomCollection->Keys)
					{
						if(((Atom^)AtomCollection[i])->Label==numbers[1]+numbers[2]) 
						{
							AtomID=((Atom^)AtomCollection[i])->getID();
							break;
						};
					}

					if(AtomID==0) throw gcnew Exception("Unknown atom in interaction.");
					if(Intkind==InteractionKind::Quadrupolar)
					{
						//Find nuclear spin quantum number and in case smaller 
						//or equal than 0.5 find other isotope 
						if (((Atom^)AtomCollection[AtomID])->isotope->Spin<=0.5f)
							((Atom^)AtomCollection[AtomID])->isotope=
							Isotopes::FindIsotopeWithHigherSpin(((Atom^)AtomCollection[AtomID])->isotope);		

						//Creating Quadrupolar interaction
						//If isotope with no 0.5>spin just ignore
						try
						{
							nInteraction=gcnew Tensor((Atom^)AtomCollection[AtomID], (Atom^)AtomCollection[AtomID],
								Intkind, matrix);
							nInteraction->unit=Units::MHz;
							nInteraction->Frame=SystemModel->RefFrameCollection[1];
							InteractionCollection->Add(nInteraction);
						}
						catch (String^ e) {/*Continue*/};
					}//Creating Chemical shield interaction
					else if(Intkind==InteractionKind::CShielding)
					{
						nInteraction=gcnew Tensor((Atom^)AtomCollection[AtomID], (Atom^)AtomCollection[AtomID],
							Intkind, matrix);
						nInteraction->unit=Units::ppm;
						nInteraction->Frame=SystemModel->RefFrameCollection[1];
						InteractionCollection->Add(nInteraction);

					}//Creating J-coupling interaction
					else if(Intkind==InteractionKind::Jcoupling)
					{
						//Find the second Atom ID from the label
						for each (int i in AtomCollection->Keys)
						{
							if(((Atom^)AtomCollection[i])->Label==numbers[3]+numbers[4]) 
							{
								AtomID2=((Atom^)AtomCollection[i])->getID();
								break;
							};
						}

						//Creating J-coupling interaction when the atoms are different
						if(AtomID2==0) throw gcnew Exception("Unknown atom in coupling.");
						if(AtomID!=AtomID2)
						{
							nInteraction=gcnew Tensor((Atom^)AtomCollection[AtomID],
								(Atom^)AtomCollection[AtomID2],Intkind, matrix);
							nInteraction->unit=Units::ppm;
							nInteraction->Frame=SystemModel->RefFrameCollection[1];
							InteractionCollection->Add(nInteraction);
						};
					}
				}
			}
		};

		fileline=ReadLine(); //Last thing is read the next line
	};

	return true;
}

/**
* @brief Main function for adding an atom to the collection checking if this happen before
* @return Returns always true.
*/
bool Castep::CheckImportAtom(cli::array<String^>^ atomline, int & AtomID)
{

	//check first if the Atom already exists
	for each(int i in AtomCollection->Keys)
	{
		if(((Atom^)AtomCollection[i])->Label==atomline[0]+atomline[1]) 
		{
			AtomID=((Atom^)AtomCollection[i])->getID();
			return true;
		};
	}

	Atom^ nAtom=gcnew Atom();

	nAtom->X=Convert::ToDouble(atomline[3], System::Globalization::CultureInfo::InvariantCulture);
	nAtom->Y=Convert::ToDouble(atomline[4], System::Globalization::CultureInfo::InvariantCulture);
	nAtom->Z=Convert::ToDouble(atomline[5], System::Globalization::CultureInfo::InvariantCulture);

	nAtom->isotope=Isotopes::FindIsotope(atomline[0],0, 0);


	nAtom->Label=atomline[0]+atomline[1];
	AtomCollection->Add(nAtom);
	AtomID=nAtom->getID();
	return true;
}
