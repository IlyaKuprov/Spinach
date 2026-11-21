#include "../StdAfx.h"
#include "Gauss.h"

using namespace System;

//Definitions
#define InteractionCollection  SystemModel->InteractionCollection
#define AtomCollection  SystemModel->AtomCollection

//Definitions of keywords-phrases found in this format
#define CMD_ISOTOPES1   		  "Multiplicity ="//
#define CMD_ISOTOPES2		      "Distance matrix (angstroms):"//
#define CMD_COORDINATES1          "Standard orientation:"//
#define CMD_COORDINATES2          "Input orientation:"//
#define CMD_GTENSOR               "g tensor [g = g_e + g_RMC + g_DC + g_OZ/SOC]:"
#define CMD_SUSCEPTIBILITY "Magnetic susceptibility (cgs-ppm):"
#define CMD_CSTENSOR              "SCF GIAO Magnetic shielding tensor (ppm):"//
#define CMD_JCOUPLINGS            "Total nuclear spin-spin coupling J (Hz):"//
#define CMD_ISOCOUPLINGS          "Isotropic Fermi Contact Couplings"
#define CMD_ANISCOUPLINGS         "Anisotropic Spin Dipole Couplings in Principal Axis System"
#define CMD_QUADCOUPLINGS         "Nuclear quadrupole coupling constants [Chi] (MHz):"

/**
*\file Gauss.cpp
*/

/**
* @brief Constructor of Gauss format
* @param filename Name of the file to be read.
* @param model The model where the imported data will be saved.
*/
Gauss::Gauss(System::String^ filename, Model^% model)
{
	SystemModel=model;
	FileName=filename;
	ReadingFile = gcnew StreamReader(filename);
	IsotropicCollection=gcnew ArrayList();
	Elements= gcnew ArrayList();
	Mass= gcnew ArrayList();
	FormatName="Gauss";
}

/**
* @brief Main function for opening and reading file
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::LoadFile(void)
{
 String^ Command=gcnew String("");
 bool  OK = true;

 //Initializing flags to false, meaning that no any frame was yet read
 read_onceCO=false;
 read_onceCS=false;
 read_onceJC=false;
 read_onceISO=false;
 read_onceANI=false;
 read_onceGT=false;
 read_onceQC=false;
 read_onceSUS=false;

 //Turn off Refreshing for faster importing
 AtomCollection->Refreshing=false;


 //Read the file
 while (!ReadingFile->EndOfStream  && OK )
 {
	 Command= ReadLine() ;

	 //Reading different sections of the file according to keywords-phrases
	 if(OK)
	 {
		 //Read isotopes
		 if (Command->Contains(CMD_ISOTOPES1))	
		 { 
			 OK=read_ISOTOPES1();
		 }//Read isotopes - alternative syntax
		 else if (Command->Contains(CMD_ISOTOPES2))	
		 {
			 OK=read_ISOTOPES2();
		 } //Read coordinates
		 else if (Command->Contains(CMD_COORDINATES1) || Command->Contains(CMD_COORDINATES2))	
		 { 
			 OK=read_COORDINATES(); 
			 read_onceCO=true; //read at least once
		 }	//Read g-tensors
		 else if (Command->Contains(CMD_GTENSOR)) 
		 {
			 OK=read_GTENSOR();
			 read_onceGT=true;  //read at least once
		 } //Read the molecular magnetic susceptibility
		 else if (Command->Contains(CMD_SUSCEPTIBILITY)) 
		 {
			 OK=read_SUSCEPTIBILITY();
			 read_onceSUS=true;  //read at least once
		 } //Read Chemical Shields
		 else if (Command->Contains(CMD_CSTENSOR)) 
		 {
			 OK=read_CSTENSOR();
			 read_onceCS=true; //read at least once
		 } //Read j-couplings
		 else if (Command->Contains(CMD_JCOUPLINGS)) 
		 {
			 OK=read_JCOUPLINGS();
			 read_onceJC=true; //read at least once
		 } //Read isotropic factors
		 else if (Command->Contains(CMD_ISOCOUPLINGS)) 
		 {
			 OK=read_ISOCOUPLINGS();
			 read_onceISO=true; //read at least once
		 } //Read Quadrupolar couplings
		 else if (Command->Contains(CMD_QUADCOUPLINGS)) 
		 {
			 OK=read_QUADCOUPLINGS();
			 read_onceQC=true; //read at least once
		 } //Read anisotropic couplings for hfc
		 else if (Command->Contains(CMD_ANISCOUPLINGS)) 
		 {
			 OK=read_ANISCOUPLINGS();
			 read_onceANI=true; //read at least once
		 }

	 };

	 //Wrong format inside senctions
	 if(!OK) throw gcnew Exception("Problem in Reading File in Line "+LineCount+".") ;
 };

 // Check if Gaussian format, at least one of the section to be read
 if(!(read_onceCO || read_onceGT || read_onceCS ||
	 read_onceJC || read_onceISO || read_onceQC || read_onceANI || read_onceSUS)) return false;


 //Compine isotropic factor and anisotropic coupligs to create HFC
 if(read_onceANI && read_onceISO)
 {
	 int electronindex;

	 //If not any electron, create one and a gtensor
	 if(AtomCollection->get_electron_number()<1)
	 {
		 //Create electron
		 Atom^ nAtom=gcnew Atom();
		 nAtom->isotope=Isotopes::FindIsotope("e",0, 0); 
		 nAtom->X=0.0f;
		 nAtom->Y=0.0f;
		 nAtom->Z=0.0f;

		 //Add electron to collection
		 AtomCollection->Add(nAtom);

		 //Save the ID of the electron
		 electronindex=nAtom->getID();

		 //Create gtensor
		 ::Matrix3x3 ^g_tensor=gcnew ::Matrix3x3();
		 g_tensor[0]=2.0023f; g_tensor[4]=2.0023f;g_tensor[8]=2.0023f;
		 g_tensor[1]=0.0f; g_tensor[2]=0.0f; g_tensor[3]=0.0f; 
		 g_tensor[5]=0.0f; g_tensor[6]=0.0f; g_tensor[7]=0.0f;


		 //Try to create new Interaction
		 Tensor^ nInteraction;
		 try{
			 nInteraction=gcnew Tensor(nAtom, nAtom,
				 InteractionKind::GTensor, g_tensor);
			 nInteraction->unit=Units::Mi_b;
			 nInteraction->Frame=SystemModel->RefFrameCollection[1];
		 }//Throws when an invalid interaction is tried to create.
		 catch(Exception ^ e)  
		 { 
			 throw gcnew Exception("Problem in creating gtensor"+e->Message
				 +". Problem in Reading File in Line:"+getLineStop()+" .");
		 };

		 //Add interaction to collection
		 InteractionCollection->Add(nInteraction);

	 }
	 else //find alectron and save its ID
		 for each(int i in AtomCollection->Keys) 
			 if(((Atom^)AtomCollection[i])->isotope->Element=="e") electronindex=i;


	 int jj=0;

	 //Loop over all interactions
	 for each(int i in gcnew List<int>(InteractionCollection->Keys)) 
	 {
		 //Find the spin rotations in order to change them to HFC
		 if(((Interaction ^)InteractionCollection[i])->IntKind==InteractionKind::spinrotation)
		 {
			 ::Matrix3x3 ^ HF_Couplings=gcnew ::Matrix3x3();

			 //The interaction matrix of HFC is the anistropix part 
			 //plus the istropic factor in the diagonal elements
			 HF_Couplings=((Tensor ^)InteractionCollection[i])->matrix3x3;
			 HF_Couplings[0]+=(double)IsotropicCollection[jj];
			 HF_Couplings[4]+=(double)IsotropicCollection[jj];
			 HF_Couplings[8]+=(double)IsotropicCollection[jj];

			 //Try to create new Interaction
			 Tensor^ nInteraction;
			 try{
				 nInteraction=gcnew Tensor(((Interaction ^)InteractionCollection[i])->A , 
					 (Atom^)AtomCollection[electronindex],
					 InteractionKind::HFC, 
					 HF_Couplings);
				 nInteraction->unit=Units::Gauss;
				 nInteraction->Frame=SystemModel->RefFrameCollection[1];
			 }//Throws when an invalid interaction is tried to create.
			 catch(Exception ^ e)  
			 { 
				 throw gcnew Exception("Problem in creating gtensor"+e->Message
					 +". Problem in Reading File in Line:"+getLineStop()+" .");
			 };

			 //Try to replace existing interaction with new HFC
			 try{InteractionCollection[i]=nInteraction; jj++;}
			 catch(String ^ e)  {throw gcnew Exception("Problem in Elements: "+e);}
		 };
	 };
 };

 //Turn on again Refreshing which causes calculation of bonds
 AtomCollection->Refreshing=true;
 ReadingFile->Close();
 return true;
}

/**
* @brief Main function for opening  isotope section and save the elements and mass numbers if exist.
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::read_ISOTOPES1()
{

	int number=1;   //number of fields in the a line, here can be anything
	String^ fileline=gcnew String(""); //line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,";            //seperator of fields =
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	

	//Sometimes next line contains the word symbolic, so ignore and read next line
	if((fileline=ReadLine())->Contains("Symbolic")) fileline=ReadLine();

	//Read as long as they are more than zero fields in the line
	while (!(fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries)->Length<1))
	{
		// C                     0.4865   -0.1848    1.3191   //
		//Split the line to a number of fields
		numbers = nullptr;
		numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

		//Check if the number fields is correct
		if(number>numbers->Length) 
		{
			throw gcnew Exception("Problem in Reading Isotope."
				+" Problem in Reading File in Line: "+getLineStop()+".") ;
		};

		//Change to capital, for sometims the first letter is not capital
		Char^ cap=Char::ToUpper(numbers[0][0]);
		numbers[0]=cap+numbers[0]->Remove(0,1);

		//Add element to the element array
		Elements->Add(numbers[0]);

		//Read next line
		fileline=ReadLine();
	};

	//We create an array with masses. For the case of " H(Iso=2) " take this number a mass
	//else mass of element zero.
	for(int i=0;i<Elements->Count;i++)
	{
		//Put mass of isotope to  zero if eveything OK
		Mass->Add((Int32 ^)0);

		//Check every character in the string of the element
		for(int jj=0;jj<((String^)Elements[i])->Length;jj++)
		{
			//In case of something is not a letter  or it is a number
			if(!Char::IsLetter(((String^)Elements[i])[jj]) ||
				Char::IsDigit(((String^)Elements[i])[jj])) 
			{
				//Find if it contains the word Iso
				if(((String^)Elements[i])->Substring(jj)->Contains("Iso"))
				{
					//Find index of mass which is where "Iso=" ends
					int num=((String^)Elements[i])->IndexOf("Iso")+4;

					//Create a substring from this index until the end
					String^ temp=((String^)Elements[i])->Substring(num);

					//Find the first no number character and cut it until the end
					for(int tt=0; tt<temp->Length;tt++)
						if(!Char::IsDigit(temp[tt])) temp=temp->Substring(0,tt);

					//Convert to number
					Mass[i]=Int32::Parse(temp);
				}

				//At the end cut the rest string and keep the beginning only
				Elements[i]=((String^)Elements[i])->Substring(0,jj);
			}
		}
	};

	return true;
}


/**
* @brief Alternative function for opening  isotope section and save the elements
* and mass numbers if exist. At the end not used
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::read_ISOTOPES2()
{
	int j=0;
	int number=3;   //number of fields in the a line
	String^ fileline=gcnew String(""); //line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,";            //seperator of fields=
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	
	
	//                    1          2          3          4          5
	fileline=ReadLine();

	//Read as long as they are thre or more fields in the line
	while (!(fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries)->Length>2))
	{
		//1  C    0.000000
		//Split the line to a number of fields
		numbers = nullptr;
		numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

		//Check if the numebr fields is correct
		if(number!=numbers->Length) 
		{
			throw gcnew Exception("Problem in Reading Isotope."
				+" Problem in Reading File in Line: "+getLineStop()+".") ;
		};

		//If elements and masses was not added in the arrays before,  add them now else overwrite them.
		if(++j==Elements->Count) 
		{
			Elements->Add(numbers[1]);
			Mass->Add((Int32 ^)0);
		}
		else 
		{
			Elements[j]=numbers[1];
			Mass[j]=(Int32 ^)0;
		};

		//Read next line
		fileline=ReadLine();
	};

	//We create an array with masses. For the case of " H(Iso=2) " take this number a mass
	//else mass of element zero.
	for(int i=0;i<Elements->Count;i++)
	{
		//Check every character in the string of the element
		for(int jj=0;jj<((String^)Elements[i])->Length;jj++)
		{
			//In case of something is not a letter  or it is a number
			if(!Char::IsLetter(((String^)Elements[i])[jj]) ||
				Char::IsDigit(((String^)Elements[i])[jj])) 
			{
				//Find if it contains the word Iso
				if(((String^)Elements[i])->Substring(jj)->Contains("Iso"))
				{
					//Find index of mass which is where "Iso=" ends
					int num=((String^)Elements[i])->IndexOf("Iso")+4;

					//Create a substring from this index until the end
					String^ temp=((String^)Elements[i])->Substring(num);

					//Find the first no number character and cut it until the end
					for(int tt=0; tt<temp->Length;tt++)
						if(!Char::IsDigit(temp[tt])) temp=temp->Substring(0,tt);

					//Convert to number
					Mass[i]=Int32::Parse(temp);
				}

				//At the end cut the rest string and keep the beginning only
				Elements[i]=((String^)Elements[i])->Substring(0,jj);
			}
		}
	};

	return true;
}

/**
* @brief Main function for reading coordinates and add atoms in the collection
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::read_COORDINATES()
{
 	int number=6;   //number of fileds in the a line
	double value;  //readedd value
	String^ fileline=gcnew String(""); //line
	int AtomID=1;   
	int index=0;    //Index for the elements array
	double x,y,z;   //Coordinates

	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,=";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	

	// ---------------------------------------------------------------------
	// Center     Atomic      Atomic             Coordinates (Angstroms)
	// Number     Number       Type             X           Y           Z
	// ---------------------------------------------------------------------
	fileline=ReadLine() ;
	fileline=ReadLine() ;
	fileline=ReadLine() ;
	fileline=ReadLine() ;

	//Read until end of section
	while (!(fileline=ReadLine())->Contains("--------"))
	{
		//Split the line to a number of fields
		numbers = nullptr;
		numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

		//Check if the numebr fields is correct
		if(number!=numbers->Length) 
		{
			throw gcnew Exception("Problem in Reading Coordinates."
				+" Problem in Reading File in Line: "+getLineStop()+".") ;
		};

		//Try to read coordinates
		if (!SpinachGUI::TryParseNumber(numbers[3], x) || !SpinachGUI::TryParseNumber(numbers[4], y) ||
			!SpinachGUI::TryParseNumber(numbers[5], z))
			throw gcnew Exception("Error in Atom coordinates");

		Atom^ nAtom=gcnew Atom();

		//Save the new atom and coordinates in the model
		nAtom->X=x; nAtom->Y=y; nAtom->Z=z;

		//Find isotope from the elements matrix
		try{
			nAtom->isotope=Isotopes::FindIsotope(((String^) Elements[index]),
				Convert::ToInt32(numbers[1]), 
				(int)((Int32^) Mass[index]));
			index++;
		}
		catch(Exception ^ e)
		{
			throw gcnew Exception("Problem in Elements:"+e->Message
				+". Problem in Reading File in Line:"+getLineStop()+" .");
		}

		//If this is the first time reading this section(read_once==false)
		// just add it in the system else replace the previous frame value(if exists)
		if(read_onceCO) 
			try{AtomCollection[AtomID]=nAtom; AtomID++;}
		catch(Exception ^ e)  {throw gcnew Exception("Problem in Elements: "+e->Message);}
		else AtomCollection->Add(nAtom);
	};

	//Save nuber of atoms
	AtomCount=AtomCollection->Count;

	return true;
}

/**
* @brief Function for reading isotropic factors  for each atom to used for hyperfine couplings later
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::read_ISOCOUPLINGS()
{

	int number=6;   //number of fields in the a line
	double value; //readed value
	String^ fileline=gcnew String(""); //line
	int jj=0;
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,=";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	

	//Atom                 a.u.       MegaHertz       Gauss      10(-4) cm-1
	fileline=ReadLine() ;

	//Read until end of section
	while (!(fileline=ReadLine())->Contains("--------"))
	{
		//Split the line to a number of fields
		numbers = nullptr;
		numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

		//Check if the numebr fields is correct
		if(number!=numbers->Length)
		{
			throw gcnew Exception("Problem in Reading Isotropic couplings"
				+". Problem in Reading File in Line: "+getLineStop()+".") ;
		};

		//Try to read elements and save in the matrix
		if (!SpinachGUI::TryParseNumber(numbers[4], value)) 
			throw gcnew Exception("Problem in Reading Isotropic factor in Line "+
			LineCount+". The format of field(double) is incorrect.");

		//If this is the first time reading this section(read_once==false)
		// just add it in the system else replace the previous frame value(if exists)
		if(read_onceISO) 
			try{IsotropicCollection[jj]=value; jj++;}
		catch(String ^ e)  {throw gcnew Exception("Problem in Elements: "+e);}
		else  IsotropicCollection->Add(value);
	};

	return true;
}

/**
* @brief Main function for reading J-Couplings from the file
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::read_JCOUPLINGS()
{
	int number;   //number of fields in the a line
	double value; //readed value
	String^ fileline=gcnew String(""); //line
	int jj=0;

	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,=";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	

	//Saving the starting index of this interaction for case of multiple frames
	static int startIndex;
	if(!read_onceJC) startIndex=InteractionCollection->Count;


	//Loop over the parts of matrix  
	for(int k=0;k<AtomCount/5+(int)(AtomCount%5>0); k++)
	{
		//        1             2             3             4             5 
		fileline=ReadLine() ;

		//Read line with up to 6 columns
		for(int i=0; i<AtomCount-5*k; i++)
		{
			fileline=ReadLine();

			//Sometimes the umber have "D" for double precision, so change it to E
			fileline=fileline->Replace(L'D',L'E');

			//If line number >4th line the fields are 6 else line number +2
			if((i+2) < 6) number=i+2; else number=6;

			//Split the line to a number of fields
			numbers = nullptr;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

			//Check if the numebr fields is correct
			if(number!=numbers->Length) 
			{
				throw gcnew Exception("Problem in Reading Total Couplings"
					+". Problem in Reading File in Line: "+getLineStop()+".") ;
			};

			//      4  0.601661D+00  0.182273D+01  0.147429D+03  0.385007D+12
			//Loop over the fields in line
			for(int j=1; j<number-((int)(i<5)); j++)
			{

				//Try to read elements and save in the matrix
				if (!SpinachGUI::TryParseNumber(numbers[j], value)) 
					throw gcnew Exception("Problem in Reading J-Couplings matrix elements in Line "+
					LineCount+". The format of fields(double) is incorrect.");

				//Set the matrix
				::Matrix3x3 ^Matrix=gcnew ::Matrix3x3();
				Matrix[0]=value; Matrix[4]=value; Matrix[8]=value;
				Matrix[1]=0.0; Matrix[2]=0.0; Matrix[3]=0.0; 
				Matrix[5]=0.0; Matrix[6]=0.0; Matrix[7]=0.0;

				//Try to create new Interaction
				Tensor^ nInteraction;
				try{
					nInteraction=gcnew Tensor((Atom^)AtomCollection[i+5*k+1],
						(Atom^)AtomCollection[j+5*k],
						InteractionKind::Jcoupling, 
						Matrix);
					nInteraction->unit=Units::Hz;
					nInteraction->Frame=SystemModel->RefFrameCollection[1];
				}//Throws when an invalid interaction is tried to create.
				catch(Exception ^ e)  
				{ 
					throw gcnew Exception("Problem in Reading J-Couplings:"+e->Message
						+". Problem in Reading File in Line:"+getLineStop()+" .");
				};

				//If this is the first time reading this section(read_once==false)
				// just add it in the system else replace the previous frame value(if exists)
				if(read_onceJC) 
					try{InteractionCollection[startIndex+jj]=nInteraction; jj++;}
				catch(String ^ e)  {throw gcnew Exception("Problem in Elements: "+e);}
				else InteractionCollection->Add(nInteraction);
			};

		};

	};

	return true;
}

/**
* @brief Main function for reading Chemical shieldings from the file
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::read_CSTENSOR()
{
	String^ fileline=gcnew String(""); //line
	int ii=1;

	//Saving the starting index of this interaction for case of multiple frames
	static int startIndex;
	if(!read_onceCS) startIndex=InteractionCollection->Count;

	//Convert to Shifts
	ToShift=true;

	//Loop over all atoms
	for each(int i in AtomCollection->Keys)
	{
		// 1  C    Isotropic =    57.9607   Anisotropy =   128.1526
		fileline=ReadLine();

		//Reading interaction matrix and creating interaction
		//  XX=   124.2135   YX=    17.8802   ZX=   -52.5149
		//  XY=    16.8895   YY=    49.5584   ZY=    50.7000
		//  XZ=   -52.1651   YZ=    56.3286   ZZ=     0.1103
		fileline=ReadLine();
		//Number of Fields: 6, XYZ index: 1,3,5
		ExtractInteractionMatrix3x3(fileline, "Chemical Shielding", 
			(Atom^)AtomCollection[i], (Atom^)AtomCollection[i], 1.0,
			read_onceCS, startIndex, ii, 6, 1, 3, 5, 
			Units::ppm, InteractionKind::CShielding, 
			SystemModel->RefFrameCollection[1], 
			InteractionCollection);

		//Dummy
		//Eigenvalues:   -50.5935    81.0799   143.3958
		fileline=ReadLine() ;
	};   
	return true;
}

/**
* @brief Main function for reading Anisotropic eigenvalues and DCMs and save them as spin 
* rotation couplings in order to create hyperfine interactions. Actually we create the 
* interaction matrix from them.
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::read_ANISCOUPLINGS()
{
	int number=8;   //number of fields in the a line
	String^ fileline=gcnew String(""); //line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,=";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	

	//Saving the starting index of this interaction for case of multiple frames
	if(!read_onceANI) startIndexANIS=InteractionCollection->Count;

	// ---------------------------------------------------------------------------------
	//
	//    Atom             a.u.   MegaHertz   Gauss  10(-4) cm-1        Axes
	//
	fileline=ReadLine() ;
	fileline=ReadLine() ;
	fileline=ReadLine() ;

	//Loop over all atoms excepts electrons
	for each(int i in AtomCollection->Keys) 
		if(((Atom^)AtomCollection[i])->isotope->Element!="e")
		{
			double *eigen_Anis_Couplings=new double[9]; //eigen values
			::Matrix3x3 ^ eigenvec_Anis_Couplings=gcnew ::Matrix3x3(); //dcm
			double x,y,z,k;

			// <Empty line>
			fileline=ReadLine();

			//          Baa    -0.2281   -30.604   -10.920   -10.209  0.0000 -0.4863  0.8738
			fileline=ReadLine();

			//Split the line to a number of fields
			numbers = nullptr;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

			//Try to read elements and save in the matrix
			// Reject wrong field counts before indexing the fields.
			if (number!=numbers->Length || 
				!SpinachGUI::TryParseNumber(numbers[3], x) ||
				!SpinachGUI::TryParseNumber(numbers[5], y) ||
				!SpinachGUI::TryParseNumber(numbers[6], z) ||
				!SpinachGUI::TryParseNumber(numbers[7], k))
				throw gcnew Exception("Problem in Reading Anisotropic Couplings in Line "+LineCount+". " +
				"Either the number of fields in the line is incorrect or "+
				"the format of fields(double) is incorrect.");
			else {
				eigen_Anis_Couplings[1]=eigen_Anis_Couplings[2]=0.0f;
				eigen_Anis_Couplings[0]=x;
				eigenvec_Anis_Couplings[0]=y;
				eigenvec_Anis_Couplings[1]=z;
				eigenvec_Anis_Couplings[2]=k;
			};

			// 1 C(13)  Bbb    -0.2215   -29.718   -10.604    -9.913  0.0000  0.8738  0.4863
			fileline=ReadLine();

			//Split the line to a number of fields
			numbers = nullptr;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

			//Try to read elements and save in the matrix
			// Reject wrong field counts before indexing the fields.
			if ((number+2)!=numbers->Length ||
				!SpinachGUI::TryParseNumber(numbers[5], x) ||
				!SpinachGUI::TryParseNumber(numbers[7], y) ||
				!SpinachGUI::TryParseNumber(numbers[8], z) ||
				!SpinachGUI::TryParseNumber(numbers[9], k))
				throw gcnew Exception("Problem in Reading Anisotropic Couplings in Line "+LineCount+". " +
				"Either the number of fields in the line is incorrect or "+
				"the format of fields(double) is incorrect.");
			else {
				eigen_Anis_Couplings[3]=eigen_Anis_Couplings[5]=0.0f;
				eigen_Anis_Couplings[4]=x;
				eigenvec_Anis_Couplings[3]=y;
				eigenvec_Anis_Couplings[4]=z;
				eigenvec_Anis_Couplings[5]=k;
			};

			//          Bcc     0.4495    60.322    21.524    20.121  1.0000  0.0000  0.0000
			fileline=ReadLine();

			//Split the line to a number of fields
			numbers = nullptr;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );

			//Try to read elements and save in the matrix
			// Reject wrong field counts before indexing the fields.
			if (number!=numbers->Length || 
				!SpinachGUI::TryParseNumber(numbers[3], x) ||
				!SpinachGUI::TryParseNumber(numbers[5], y) ||
				!SpinachGUI::TryParseNumber(numbers[6], z) ||
				!SpinachGUI::TryParseNumber(numbers[7], k))
				throw gcnew Exception("Problem in Reading Anisotropic Couplings in Line "+LineCount+". " +
				"Either the number of fields in the line is incorrect or "+
				"the format of fields(double) is incorrect.");
			else {
				eigen_Anis_Couplings[6]=eigen_Anis_Couplings[7]=0.0f;
				eigen_Anis_Couplings[8]=x;
				eigenvec_Anis_Couplings[6]=y;
				eigenvec_Anis_Couplings[7]=z;
				eigenvec_Anis_Couplings[8]=k;
			};

			//Generate Tensor from eigenvectors and eigen values
			Matrix3d result,eigenVectors,eigenValues;
			for(int j=0;j<9;j++) 
			{
				eigenVectors(j)=eigenvec_Anis_Couplings[j];
				eigenValues(j)=eigen_Anis_Couplings[j];
			};

			//Make the calculation
			result= (eigenVectors*eigenValues)*eigenVectors.transpose();

			//Transfer to the matrix
			for(int j=0;j<9;j++) eigenvec_Anis_Couplings[j]=result(j);


			//Try tocreate temporary spirotation then convert to HFC
			Tensor^ nInteraction;
			try{
				nInteraction=gcnew Tensor((Atom^)AtomCollection[i], 
					(Atom^)AtomCollection[i], 
					InteractionKind::spinrotation, 
					eigenvec_Anis_Couplings);
			}
			catch(Exception ^ e) 
			{
				throw gcnew Exception("Problem in Reading Hyperfine: "+e->Message
					+". Problem in Reading File in Line:"+getLineStop()+" .");
			};

			//If this is the first time reading this section(read_once==false)
			// just add it in the system else replace the previous frame value(if exists)
			if(read_onceANI)
				try{InteractionCollection[i+startIndexANIS]=nInteraction;}
			catch(String ^ e)  {throw gcnew Exception("Problem in Elements: "+e);}
			else InteractionCollection->Add(nInteraction);		
		};

		return true;
}

/**
* @brief Main function for reading G-tesnor from the file
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::read_GTENSOR()
{
   	String^ fileline=gcnew String(""); //line
	int ii=1;

	//Saving the starting index of this interaction and electron for case of multiple frames
	static int startIndex;
	static int ElectronstartIndex;
	if(!read_onceGT) 
	{
		startIndex=InteractionCollection->Count;
		ElectronstartIndex=AtomCollection->Count;
	}

	//Create electron
	Atom^ nAtom=gcnew Atom();
	nAtom->isotope=Isotopes::FindIsotope("e",0, 0); //Electron mass=0
	nAtom->X=0.0f;
	nAtom->Y=0.0f;
	nAtom->Z=0.0f;

	//If this is the first time reading this section(read_once==false)
	// just add it in the system else replace the previous frame value(if exists)
	if(read_onceGT) AtomCollection[ElectronstartIndex+1]=nAtom;
	else AtomCollection->Add(nAtom);

	//nAtom->ID=AtomCollection->Count;

	//Reading interaction matrix and creating interaction
	//  XX= 124.2135  YX=    17.8802   ZX=-52.5149
	//  XY= 16.8895   YY=    49.5584   ZY= 50.7000
	//  XZ=-52.1651   YZ=    56.3286   ZZ=  0.1103
	fileline=ReadLine();
	
	//Number of Fields: 6, XYZ index: 1,3,5
	ExtractInteractionMatrix3x3(fileline, "G-tesnor", 
		nAtom, nAtom, 1.0,
		read_onceGT, startIndex, ii, 6, 1, 3, 5, 
		Units::Mi_b, InteractionKind::GTensor, 
		SystemModel->RefFrameCollection[1], 
		InteractionCollection);

	return true;
}

/**
* @brief Main function for reading Quadrupolar interactions from the file
* @return Returns false if the format does not seem to be correct.
*/
bool Gauss::read_QUADCOUPLINGS()
{
	int number;   //number of fields in the a line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,=";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array	

	String^ fileline=gcnew String(""); //line
	int ii=1;

	//Saving the starting index of this interaction for case of multiple frames
	static int startIndex;
	if(!read_onceQC) startIndex=InteractionCollection->Count;

	//Loop over all atoms
	for each(int i in AtomCollection->Keys)
	{
		//    1  C(12) 
		fileline=ReadLine();

		//Split the line
	   numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
	   
	   //check if all the atoms have quadropolar and save atom ID
	   int AtomID;
	   if(numbers->Length!=2 || !Int32::TryParse(numbers[0], AtomID)) return true;

		//For Quadrupolar interactions, try to find proper isotope
		try{
			if(((Atom^)AtomCollection[AtomID])->isotope->Spin<=0.5f) 
				((Atom^)AtomCollection[AtomID])->isotope=
				Isotopes::FindIsotopeWithHigherSpin(((Atom^)AtomCollection[i])->isotope);

			//Reading interaction matrix and creating interaction
			//aa=     0.0000   ba=     0.0000   ca=     0.0000
			//ab=     0.0000   bb=     0.0000   cb=     0.0000
			//ac=     0.0000   bc=     0.0000   cc=     0.0000
			fileline=ReadLine();
			//Number of Fields: 6, XYZ index: 1,3,5
			ExtractInteractionMatrix3x3(fileline, "Quadrupolar Interaction", 
				(Atom^)AtomCollection[AtomID], (Atom^)AtomCollection[AtomID], 1.0,
				read_onceQC, startIndex, ii, 6, 1, 3, 5, 
				Units::MHz, InteractionKind::Quadrupolar, 
				SystemModel->RefFrameCollection[1], 
				InteractionCollection);
		}//In case we can not find isotope with  spin>0.5
		catch (String^ e) {/*Continue*/};
	};   
	return true;
}

/**
* @brief Reads the molecular magnetic susceptibility tensor printed by Gaussian NMR jobs.
*
* Gaussian 16 prints
*    Magnetic susceptibility (cgs-ppm):
*      Isotropic =   -49.0654   Anisotropy=    36.5994
*      XX=   -24.6658   YX=    -0.0000   ZX=     0.0000
*      XY=     0.0000   YY=   -24.6658   ZY=     0.0000
*      XZ=    -0.0000   YZ=    -0.0000   ZZ=   -97.8646
* Elements are placed by their labels, so the order of the rows does not matter. The tensor
* is a property of the molecule (cgs-ppm, in the frame and about the origin of the coordinates
* in the file), not a spin interaction; a later block replaces an earlier one.
* @return false unless all nine elements are found, once each, right after the header.
*/
bool Gauss::read_SUSCEPTIBILITY()
{
	::Matrix3x3 ^chi=gcnew ::Matrix3x3();
	cli::array<bool> ^found=gcnew cli::array<bool>(9);
	int nfound=0;
	System::Text::RegularExpressions::Regex ^element=
		gcnew System::Text::RegularExpressions::Regex("([XYZ])([XYZ])=\\s*([-+]?\\d+\\.\\d+)");

	//The isotropy/anisotropy line and the three rows follow the header
	for (int line=0; line<6 && nfound<9; line++)
	{
		String ^fileline=ReadLine();
		if (fileline==nullptr) return false;
		for each (System::Text::RegularExpressions::Match ^m in element->Matches(fileline))
		{
			int k=3*(m->Groups[1]->Value[0]-L'X')+(m->Groups[2]->Value[0]-L'X');
			double value;
			if (found[k] || !Double::TryParse(m->Groups[3]->Value, System::Globalization::NumberStyles::Float,
				System::Globalization::CultureInfo::InvariantCulture, value)) return false;
			chi[k]=value;
			found[k]=true;
			nfound++;
		}
	}
	if (nfound!=9) return false;

	//The carrier atom is not added to the model, so the tensor never enters the atom or interaction collections
	Atom ^carrier=gcnew Atom();
	SystemModel->Susceptibility=gcnew Tensor(carrier, carrier, InteractionKind::Susceptibility, chi);
	return true;
}
