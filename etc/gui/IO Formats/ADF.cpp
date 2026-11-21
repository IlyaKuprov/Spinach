#include "../StdAfx.h"
#include "ADF.h"


//Definitions
#define InteractionCollection  SystemModel->InteractionCollection
#define AtomCollection  SystemModel->AtomCollection
//Definitions of keywords-phrases found in this format
#define CMD_CARTCOORD             "G E O M E T R Y  ***  3D  Molecule  ***"
#define CMD_CHEMICAL_SHIFTS       "Stephen K. Wolff & Georg Schreckenbach"
#define CMD_J_COUPLINGS           "*  NMR spin spin coupling tensor calculation  *"

/**
*\file ADF.cpp
*/

/**
* @brief Constructor of ADF format
* @param filename Name of the file to be read.
* @param model The model where the imported data will be saved.
*/
ADF::ADF(System::String^ filename, Model^% model )
{
	SystemModel=model;
	FileName=filename;
	ReadingFile = gcnew StreamReader(filename);
    FormatName="ADF";
}

/**
* @brief Main function for opening and reading file
* @return Returns false if the format does not seem to be correct.
*/
bool ADF::LoadFile(void)
{
	String^ Command=gcnew String("");
	bool  OK = true;

	//Initializing flags to false, meaning that no any frame was yet read
	read_onceOR=false;
	read_onceCS=false;
	read_onceJC=false;

	//Turn off Refreshing for faster importing
	AtomCollection->Refreshing=false;


	//Read the file until the end of it and as long flag OK is true
	while (!ReadingFile->EndOfStream  && OK)
	{
		Command= ReadLine() ;

		//Reading different sections of the file according to keywords-phrases
		if(OK)
		{
			//Read coordinates
			if(Command->Contains(CMD_CARTCOORD))
			{  
				OK=read_CARTCOORD();
				read_onceOR=true; //read at least once
			}//Read chemical shields	
			else if(Command->Contains( CMD_CHEMICAL_SHIFTS)) 
			{ 
				OK=read_CHEMICAL_SHIELDING();
				read_onceCS=true;//read at least once
			} //Read J-couplings
			else if(Command->Contains(CMD_J_COUPLINGS)) 
			{
				OK=read_J_COUPLINGS();
				read_onceJC=true;//read at least once
			}
		};

		//Wrong format inside senctions
		if(!OK) throw gcnew Exception("Problem in Reading File in Line "+LineCount+".");

	};

	// Check if file is a valid ADF format, having at least one of the sections
	if(!(read_onceOR ||  read_onceCS || read_onceJC )) return false;

	//Turn on again Refreshing which causes calculation of bonds
	AtomCollection->Refreshing=true;
	ReadingFile->Close();
	return true;
}

/**
* @brief Main function for reading coordinates from the file
* @return Returns false if the format does not seem to be correct.
*/
bool ADF::read_CARTCOORD()
{
	String^ fileline=gcnew String("");  //line
	int AtomID=1;						//Saving the atom ID

	//Readin initial dummy lines
	fileline=ReadLine() ; 
	// ===============
	fileline=ReadLine() ; 
	//
	fileline=ReadLine() ; 
	//
	fileline=ReadLine() ; 
	//ATOMS
	fileline=ReadLine() ; 
	//=====                            X Y Z                    CHARGE
	fileline=ReadLine() ; 
	//                             (Angstrom)             Nucl     +Core       At.Mass
	fileline=ReadLine() ; 
	//                    --------------------------    ----------------       -------

	fileline="";
	//Read the section until the end of it which is an line with less than 2 chars
	while ((fileline=ReadLine())->Length>2 )
	{
	//Number of Fields: 8, Isotope index: 1, Mass index: 7, XYZ index: 2,3,4
	 ExtractIsotopeAndCoordinates(fileline, FormatName, read_onceOR, AtomID,
								   8, 1, 2, 3, 4, 1.0, 7, AtomCollection);
	};

	return true;

}

/**
* @brief Main function for reading Chemical shieldings from the file
* @return Returns false if the format does not seem to be correct.
*/
bool ADF::read_CHEMICAL_SHIELDING()
{
	int number; 						//number of fields in each line
	String^ fileline=gcnew String("");  //line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array
	int AtomID;						  //Saving the atom ID

	//Saving the starting index of this interaction for case of multiple frames
	int j=1;
	static int startIndexCS;
	if(!read_onceCS) startIndexCS=InteractionCollection->Count;

	//It is chemical shielding, so conversion to shift is needed
	ToShift= true;

	//Read until end of senction
	while (!(fileline=ReadLine())->Contains("NORMAL TERMINATION"))
	{
		
		//****  N U C L E U S : C(2)
		if (fileline->Contains("****  N U C L E U S :"))
		{ 
			//(blank)
			fileline=ReadLine() ;
			//ADF-GUI atom: C(2)
			fileline=ReadLine() ;

			number=3;
			numbers = nullptr;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
			if(number>numbers->Length) throw gcnew Exception("Problem in Reading CHEMICAL SHIELDING");

			//Naking of "C(2)"
			for(int ii=0;numbers[2][ii]!=Convert::ToChar("(");) 
			{
				numbers[2]=numbers[2]->Remove(ii,1);
			}
			numbers[2]=numbers[2]->Remove(0,1);
			numbers[2]=numbers[2]->Remove(numbers[2]->Length-1);
			AtomID=Convert::ToInt32(numbers[2]);

			while (!(fileline=ReadLine())->Contains(
				"********************************************************************************"))
			{
				if(fileline->Contains("===  TOTAL NMR SHIELDING TENSOR (ppm)"))
				{
					fileline=ReadLine() ; 
					//
					fileline=ReadLine() ; 
					//        ***********************************
					fileline=ReadLine() ;
					//        CARTESIAN AXIS REPRESENTATION
					fileline=ReadLine() ;
					//	
					fileline=ReadLine() ;
					//        ==== total shielding tensor
					fileline=ReadLine() ;
					//         -----------------------------------

					//Reading interaction matrix and creating interaction
					fileline=ReadLine();
					//Number of Fields: 3, XYZ index: 0,1,2
					ExtractInteractionMatrix3x3(fileline, "Chemical Shielding", 
						(Atom^)AtomCollection[AtomID], (Atom^)AtomCollection[AtomID], 1.0,
						read_onceCS, startIndexCS, j, 3, 0, 1, 2, 
						Units::ppm, InteractionKind::CShielding, 
						SystemModel->RefFrameCollection[1], 
						InteractionCollection);
						
				};// "total shielding tensor"
			};//Searching lines  for "total shielding tensor"

		}; //Nucleus

	}; //CHEMICAL SHIELDING

	return true;
}

/**
* @brief Main function for reading J-Couplings from the file
* @return Returns false if the format does not seem to be correct.
*/
bool ADF::read_J_COUPLINGS()
{

	int number; 						//number of fields in each line
	String^ fileline=gcnew String("");  //line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array
	int AtomID1, AtomID2; //Saving ID of the first and the second atom
	double value;		  //Scalar value of the j-coupling
	int existingID; //Flag for if an j-coupling is duplicated

	//Saving the starting index of this interaction for case of multiple frames

	//Read until end of senction
	while (!(fileline=ReadLine())->Contains("NORMAL TERMINATION"))
	{
		//  Nucleus 1 (C) perturbing nucleus 2 (C)
		if (fileline->Contains("Nucleus") && fileline->Contains("perturbing nucleus") )
		{ 
			//Read the atom IDs which interact through j-coupling
			number=7;
			numbers = nullptr;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
			if(number>numbers->Length) throw gcnew Exception("Problem in Reading J-Couplings");

			AtomID1=Convert::ToInt32(numbers[1]);
			AtomID2=Convert::ToInt32(numbers[5]);

			//Check for doubles, especially the ones with reverse order
			//quite slow, checking one by one
			existingID=0;
			for each (int i in InteractionCollection->Keys)
			{
				if(((Interaction^ )InteractionCollection[i])->IntKind==InteractionKind::Jcoupling)
					if((((Tensor^ )InteractionCollection[i])->A->ID==AtomID1 &&
				    ((Tensor^ )InteractionCollection[i])->B->ID==AtomID2) ||
				   (((Tensor^ )InteractionCollection[i])->A->ID==AtomID2 &&
					((Tensor^ )InteractionCollection[i])->B->ID==AtomID1)) 
					 existingID=i;
			};

			while (!(fileline=ReadLine())->Contains("Principal Axis Representation"))
			{
				//read only if it is not a repeated data
				if(fileline->Contains("total calculated spin-spin coupling:") )
				{
					::Matrix3x3 ^ J_Coupling=gcnew ::Matrix3x3();
					
					//Try read scalar quantity
					//isotropic    k=     180.984                              j=     137.514      
					number=5;
					fileline=ReadLine() ;
					numbers = nullptr;
					numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
					if(number!=numbers->Length)  throw gcnew Exception("Problem in Reading J-Coupling matrix");

					if (!Double::TryParse(numbers[4], value)) return false;

					//Set the scalar value to the diagonal of the matrix
					J_Coupling[0]=value; J_Coupling[4]=value; J_Coupling[8]=value;
					J_Coupling[1]=0.0; J_Coupling[2]=0.0; J_Coupling[3]=0.0; 
					J_Coupling[5]=0.0; J_Coupling[6]=0.0; J_Coupling[7]=0.0;

					//Create j-coupling interaction, unit and reference frame
					Tensor^ nInteraction=gcnew Tensor((Atom^)AtomCollection[AtomID1],
						(Atom^)AtomCollection[AtomID2],
						InteractionKind::Jcoupling, J_Coupling);
					nInteraction->unit=Units::Hz;
					nInteraction->Frame=SystemModel->RefFrameCollection[1];

					//If this is the first time reading this section(read_onceJC==false)
					// just add it in the system else replace the previous frame value(if exists)
					if(existingID!=0)
						try{InteractionCollection[existingID]=nInteraction;}
					catch(String ^ e)  {throw gcnew Exception("Problem in Elements: "+e);}
					else InteractionCollection->Add(nInteraction);

				};// "total calculated spin-spin coupling"
			};//Searching lines  for "total calculated spin-spin coupling"

		}; //Nucleus

	}; //J-Couplings

	return true;
}



