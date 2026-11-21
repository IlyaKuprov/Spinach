#include "../StdAfx.h"
#include "Orca.h"

//Definitions
#define InteractionCollection  SystemModel->InteractionCollection
#define AtomCollection  SystemModel->AtomCollection

//Definitions of keywords-phrases found in this format
#define CMD_CARTCOORD             "CARTESIAN COORDINATES (A.U.)"
#define CMD_GTENSOR               "The g-matrix:"
#define CMD_ZERO_FIELD            "ZERO-FIELD-SPLITTING TENSOR"
#define CMD_CHEMICAL_SHIFTS        "CHEMICAL SHIFTS"
#define CMD_HFC_QUAD_COUPLINGS    "ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE"
//Definition of factor from convertion from atomic units to amstrong
#define  AUtoAmstrong             0.529177249

/**
*\file Orca.cpp
*/

/**
* @brief Constructor of Orca format
* @param filename Name of the file to be read.
* @param model The model where the imported data will be saved.
*/
Orca::Orca(System::String^ filename, Model^% model)
{
	SystemModel=model;
	FileName=filename;
	ReadingFile = gcnew StreamReader(filename);
	FormatName="Orca";
}

/**
* @brief Main function for opening and reading file
* @return Returns false if the format does not seem to be correct.
*/
bool Orca::LoadFile(void)
{
 String^ Command;
 
 bool  OK = true;

 //Initializing flags to false, meaning that no any frame was yet read
 read_onceOR=false;
 read_onceCS=false;
 read_onceHFC_QUAD=false;
 read_onceGT=false;
 pendingProperties=gcnew InteractionsDictionary();
 read_onceZF=false;

 //Turn off Refreshing for faster importing
 AtomCollection->Refreshing=false;




 //Read the file
 while (!ReadingFile->EndOfStream  && OK )
 {
	 //if(Command==nullptr)		
	 //  {
	 Command=gcnew String("");
	 Command= ReadLine() ;
	 //  };


	 if(OK)
	 {
		 if ( String::Compare(Command->Trim(), CMD_CARTCOORD  ) == 0)
		 {  OK=read_CARTCOORD();read_onceOR=true; }	
		 else if ( String::Compare(Command->Trim(), CMD_GTENSOR) == 0 ) 
		 {OK=read_GTENSOR();read_onceGT=true; } 
		 else if ( String::Compare(Command->Trim(), CMD_ZERO_FIELD ) == 0 ) 
		 {OK=read_ZERO_FIELD();read_onceZF=true; } 
		 else if (String::Compare(Command->Trim(), CMD_CHEMICAL_SHIFTS ) == 0 ) 
		 {
			 OK=read_CHEMICAL_SHIFTS();read_onceCS=true;
		 } 
		 else if (String::Compare(Command->Trim(), CMD_HFC_QUAD_COUPLINGS)== 0) 
		 {OK=read_HFC_QUAD_COUPLINGS();read_onceHFC_QUAD=true;}
	 };

	 if(!OK) throw gcnew Exception("Problem in Reading File in Line "+LineCount+",");


 };

  // Check if Orca format
 if(!(read_onceOR || read_onceGT || read_onceZF || read_onceCS || read_onceHFC_QUAD)) return false;

//  sww->Close();
  ReadingFile->Close();
  if(read_onceHFC_QUAD)
  {

	int electronindex=0;

	if(AtomCollection->get_electron_number()<1)
	{

	  //Create electron
	  Atom^ nAtom=gcnew Atom();
	  nAtom->isotope=Isotopes::FindIsotope("e",0, 0); 
	  nAtom->X=0.0f;
	  nAtom->Y=0.0f;
	  nAtom->Z=0.0f;
	  AtomCollection->Add(nAtom);
	  electronindex=nAtom->getID();

	  //Create gtensor
	  ::Matrix3x3 ^g_tensor=gcnew ::Matrix3x3();
	  g_tensor[0]=2.0023f; g_tensor[4]=2.0023f;g_tensor[8]=2.0023f;
	  g_tensor[1]=0.0f; g_tensor[2]=0.0f; g_tensor[3]=0.0f; 
	  g_tensor[5]=0.0f; g_tensor[6]=0.0f; g_tensor[7]=0.0f;
	  Tensor^ nInteraction=gcnew Tensor(nAtom, nAtom,
										  InteractionKind::GTensor, g_tensor);
	 nInteraction->unit=Units::Mi_b;
	 nInteraction->Frame=SystemModel->RefFrameCollection[1];
	 InteractionCollection->Add(nInteraction);

	}
	else
		for each(int i in AtomCollection->Keys) 
			if(((Atom^)AtomCollection[i])->isotope->Element=="e") electronindex=i;

	// Commit only the final property section. Hyperfine and quadrupolar rows
	// may interleave, so neither kind has a contiguous range of model IDs.
	for each(int id in pendingProperties->Keys)
	{
		Tensor^ tensor=(Tensor^)pendingProperties[id];
		if(tensor->IntKind==InteractionKind::spinrotation)
		{
			tensor=gcnew Tensor(tensor->A,AtomCollection[electronindex],
				InteractionKind::HFC,tensor->matrix3x3);
			tensor->unit=Units::MHz;
			tensor->Frame=SystemModel->RefFrameCollection[1];
		}
		InteractionCollection->Add(tensor);
	}
	pendingProperties->Clear();

  };

   AtomCollection->Refreshing=true; 
  return true;
	
}

/**
* @brief Main function for reading coordinates from the file
* @return Returns false if the format does not seem to be correct.
*/
bool Orca::read_CARTCOORD()
{
	String^ fileline=gcnew String("");  //line
	int AtomID=1;						//Saving the atom ID

	//dummy lines
	//----------------------------
	fileline=ReadLine() ; 
	// NO LB      ZA    FRAG    MASS        X           Y           Z
	fileline=ReadLine() ;

	//Read the section until the end of it
	while (!(fileline=ReadLine())->Contains("--------") &&
		fileline!="")
	{
		//Number of Fields: 8, Isotope index: 1, AtomNumber Index: 2, Mass index: 4, XYZ index: 5,6,7
		ExtractIsotopeAndCoordinates(fileline, FormatName, read_onceOR, AtomID,
			8, 1, 2, 4, 5, 6, 7, AUtoAmstrong, AtomCollection);
	};

	return true;
}

/**
* @brief Main function for reading Chemical shifts from the file
* @return Returns false if the format does not seem to be correct.
*/
bool Orca::read_CHEMICAL_SHIFTS()
{
	int number;							//Number of fields in a line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	
	String^ fileline=gcnew String("");  //line
	int AtomID=1;						//Saving the atom ID

	//Saving the starting index of this interaction for case of multiple frames
	int j=1;
	static int startIndexCS;
	if(!read_onceCS) startIndexCS=InteractionCollection->Count;

	//Dummy lines
	//----------------
	//
	// ---------------
	fileline=ReadLine() ;
	fileline=ReadLine() ;
	fileline=ReadLine() ;

	fileline="";

	//Read until end of senction or end of file
	while (!(fileline=ReadLine())->Contains("-----------------------------------------") &&
		   !fileline->Contains("****ORCA TERMINATED NORMALLY****"))
	{
		
		//Nucleus   0C : A
		if (fileline->Contains("Nucleus"))
		{ 
			//Split the line to a number of fields
			number=2;
			numbers = nullptr;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
			if(number>numbers->Length) throw gcnew Exception("Problem in Reading CHEMICAL SHIFTS");

			//Clean string leaving only atom ID
			for(int i=numbers[1]->Length-1;i>-1;i--) 
				if(Char::IsLetter(numbers[1][i])) numbers[1]=numbers[1]->Remove(i);
			AtomID=Convert::ToInt32(numbers[1]);

			// ----------------------------------------------------------- 
			fileline=ReadLine();

			if((fileline=ReadLine())->Contains("Raw-matrix :"))
			{
				//Reading interaction matrix and creating interaction
				fileline=ReadLine();
				//Number of Fields: 3, XYZ index: 0,1,2
				ExtractInteractionMatrix3x3(fileline, "Chemical Shielding", 
					(Atom^)AtomCollection[AtomID+1], (Atom^)AtomCollection[AtomID+1], 1.0,
					read_onceCS, startIndexCS, j, 3, 0, 1, 2, 
					Units::MHz, InteractionKind::Shift, 
					SystemModel->RefFrameCollection[1], 
					InteractionCollection);
			}
		}; //Nucleus

	}; //CHEMICAL SHIFTS

	return true;

}

/**
* @brief Main function for reading Hyperfine couplings and quadropolar couplings from the file
* @return Returns false if the format does not seem to be correct.
*/
bool Orca::read_HFC_QUAD_COUPLINGS()
{

	int number;							//Number of fields in a line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	
	String^ fileline=gcnew String("");  //line
	int AtomID=1;						//Saving the atom ID
	double factor=2.3496e+02;           //(efg_atomic*e_charge*1e-28/h_planck)MHz
	double Qbarn;						//Quadropolar factor
	pendingProperties->Clear();
	int j=1;

	//Dummy lines
	//-----------------------------------------
	//
	//-----------------------------------------------------------
	fileline=ReadLine() ;
	fileline=ReadLine() ;
	fileline=ReadLine() ;

	//Read until end of senction
	while (!(fileline=ReadLine())->Contains("****ORCA TERMINATED NORMALLY****") )
	{
		
		//Nucleus   0C : A:ISTP=   13 I=  0.5 P=134.1900 MHz/au**3
		//fileline= ReadLine() ;
		if (fileline->Contains("Nucleus"))
		{ 
			//Split the line to a number of fields
			number=5;
			numbers = nullptr;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
			if(number>numbers->Length) 
				throw gcnew Exception("Problem in Reading ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE");

			//Clean string leaving only atom ID
			for(int i=numbers[1]->Length-1;i>-1;i--) 
				if(Char::IsLetter(numbers[1][i])) numbers[1]=numbers[1]->Remove(i);
			AtomID=Convert::ToInt32(numbers[1]);

			//Reading Qbarn factor
			//                Q:ISTP=   13 I=  0.5 Q=  0.0000 barn  
			fileline=ReadLine();
			number=7;
			numbers = nullptr;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
			if(number!=numbers->Length) 
				throw gcnew Exception("Problem in Reading ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE");
			if(!TryParseNumber(numbers[5],Qbarn))
				throw gcnew FormatException("Invalid quadrupole moment in ORCA output.");

			// ----------------------------------------------------------- 
			fileline=ReadLine();

			//Read until end of senction of end of file
			while (!(fileline=ReadLine())->
				Contains(" -----------------------------------------------------------") &&
				!fileline->Contains("****ORCA TERMINATED NORMALLY****"))
			{
				//Read Hyperfine interactions
				if(fileline->Contains("Raw HFC matrix (all values in MHz):"))
				{
					//  ------------------------------
					fileline=ReadLine();

					//Reading interaction matrix and creating interaction
					fileline=ReadLine();
					//Number of Fields: 3, XYZ index: 0,1,2
					ExtractInteractionMatrix3x3(fileline, "Hyeperfine Couplings", 
						(Atom^)AtomCollection[AtomID+1], (Atom^)AtomCollection[AtomID+1], 1.0,
						false, 0, j, 3, 0, 1, 2, 
						Units::MHz, InteractionKind::spinrotation, 
						SystemModel->RefFrameCollection[1], 
						pendingProperties);

				}//Read quadropolar interactions
				else if(fileline->Contains("Raw EFG matrix (all values in a.u.**-3):"))
				{
					//Find nuclear spin quantum number and in case smaller 
					//or equal than 0.5 find other isotope 
					if (((Atom^)AtomCollection[AtomID+1])->isotope->Spin<=0.5)	
						((Atom^)AtomCollection[AtomID+1])->isotope=
						Isotopes::FindIsotopeWithHigherSpin(((Atom^)AtomCollection[AtomID+1])->isotope);
					double NSQN=((Atom^)AtomCollection[AtomID+1])->isotope->Spin; //nuclear spin quantum number


					// Keep any nonzero EFG, including a traceless diagonal tensor.
					::Matrix3x3^ Quad_Coupling=gcnew ::Matrix3x3();
					const double scale=NSQN>0.5 ? factor*Qbarn/(2*NSQN*(2*NSQN-1)) : 0.0;
					bool nonzero=false;
					for(int row=0; row<3; row++)
					{
						numbers=ReadLine()->Split(delimiter,StringSplitOptions::RemoveEmptyEntries);
						if(numbers->Length!=3) throw gcnew FormatException("Invalid ORCA EFG row.");
						for(int col=0; col<3; col++)
						{
							double value;
							if(!TryParseNumber(numbers[col],value)) throw gcnew FormatException("Invalid ORCA EFG element.");
							Quad_Coupling[3*row+col]=scale*value;
							if(scale*value!=0.0) nonzero=true;
						}
					}
					if(nonzero)
					{
						//If isotope with no 0.5>spin just ignore
						try
						{
							Tensor^ nInteraction=gcnew Tensor((Atom^)AtomCollection[AtomID+1], 
								(Atom^)AtomCollection[AtomID+1],
								InteractionKind::Quadrupolar, Quad_Coupling);
							nInteraction->unit=Units::MHz;
							nInteraction->Frame=SystemModel->RefFrameCollection[1];
							pendingProperties->Add(nInteraction);
						}
						catch (String^ e) {/*Continue*/};
					};

				};


			};

			//in order to stop while
			if (fileline->Contains("****ORCA TERMINATED NORMALLY****")) return true;
		}; //Nucleus

	}; //ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE

	return true;

}

/**
* @brief Main function for reading g-tensor couplins from the file
* @return Returns false if the format does not seem to be correct.
*/
bool Orca::read_GTENSOR()
{
	String^ fileline=gcnew String(""); //line
	int j=1;
	//Saving the starting index of this interaction for case of multiple frames
	static int startIndexGT;
	if(!read_onceGT) startIndexGT=InteractionCollection->Count;

	//Creating new atom for electron
	Atom^ nAtom=gcnew Atom();
	nAtom->X=0.0f;
	nAtom->Y=0.0f;
	nAtom->Z=0.0f;
	nAtom->isotope=Isotopes::FindIsotope("e",0, 0); //Electron mass=0

	//add it to the collection
	if(read_onceGT) AtomCollection[AtomCollection->Count]=nAtom;
	else AtomCollection->Add(nAtom);

	//Reading interaction matrix and creating interaction
	fileline=ReadLine();
	//Number of Fields: 3, XYZ index: 0,1,2
	ExtractInteractionMatrix3x3(fileline, "G-Tensor", 
		nAtom, nAtom, 1.0,
		read_onceGT, startIndexGT, j, 3, 0, 1, 2,
		Units::Mi_b, InteractionKind::GTensor, 
		SystemModel->RefFrameCollection[1], 
		InteractionCollection);

	return true;
}

/**
* @brief Main function for reading zero-field couplings from the file
* @return Returns false if the format does not seem to be correct.
*/
bool Orca::read_ZERO_FIELD()
{
	String^ fileline=gcnew String(""); //line
	int j=1;
	//Saving the starting index of this interaction for case of multiple frames
	static int startIndexZF;
	if(!read_onceZF) startIndexZF=InteractionCollection->Count;

	//Dummy lines
	//---------------------------
	//
	//raw-matrix : 
	fileline=ReadLine() ;
	fileline=ReadLine() ;
	fileline=ReadLine() ;

	//Check if electrons were created before
	Atom^ nAtom=gcnew Atom();
	if(AtomCollection->get_electron_number()>0 || read_onceGT || read_onceZF)
	{
		for each(int i in AtomCollection->Keys)
			if(AtomCollection[i]->isotope->Element=="e") nAtom=AtomCollection[i];
	}
	else
	{
		nAtom->X=0.0f;
		nAtom->Y=0.0f;
		nAtom->Z=0.0f;
		nAtom->isotope=Isotopes::FindIsotope("e",0, 0); //Electron mass=0
		AtomCollection->Add(nAtom);
	}

	//Reading interaction matrix and creating interaction
	fileline=ReadLine();
	//Number of Fields: 3, XYZ index: 0,1,2
	ExtractInteractionMatrix3x3(fileline, "Zero-Field coupling", 
		nAtom, nAtom, 1.0,
		read_onceZF, startIndexZF, j, 3, 0, 1, 2,
		Units::rev_cm, InteractionKind::ZFS, 
		SystemModel->RefFrameCollection[1], 
		InteractionCollection);

	return true;
}
