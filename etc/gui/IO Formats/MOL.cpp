#include "../StdAfx.h"
#include "MOL.h"


//Definitions
#define AtomCollection  SystemModel->AtomCollection


/**
*\file MOL.cpp
*/

/**
* @brief Constructor of MOL format
* @param filename Name of the file to be read.
* @param model The model where the imported data will be saved.
*/
MOL::MOL(System::String^ filename, Model^% model)
{
	SystemModel=model;
	FileName=filename;
	ReadingFile = gcnew StreamReader(filename);
	FormatName="MOL";
}

/**
* @brief Main function for opening and reading file
* @return Returns false if the format does not seem to be correct.
*/
bool MOL::LoadFile(void)
{
	bool  OK = true;

	//Turn off Refreshing for faster importing
	AtomCollection->Refreshing=false;

	
	//Read the file until the end of it and as long flag OK is true
	while (!ReadingFile->EndOfStream  && OK )
	{
		int totalAtoms=0;					//total number of atoms in the system
		int totalBonds=0;					//total number of bonds in the system
		String^ fileline=gcnew String("");  //line

		try {
			//Read comments 3 lines
			fileline=ReadLine();
			fileline=ReadLine();
			fileline=ReadLine();

			//Read number of atoms and Bonds
			fileline=ReadLine();
            if(fileline==nullptr || fileline->Length<6 || fileline->Contains("V3000") ||
                !Int32::TryParse(fileline->Substring(0,3),totalAtoms) ||
                !Int32::TryParse(fileline->Substring(3,3),totalBonds) || totalAtoms<0 || totalBonds<0)
                OK=false;
		}
		catch(Exception ^ e)
		{
			throw gcnew Exception("Unexpected end of file in Line "+LineCount+":"+e->Message) ;
		};

		
        if(!OK) break;
		//Read the atom format: 0.4865   -0.1848    1.3191 C 0 0 0 0 0 0 0 0 0 0 0 0              //
		for(int i=0;i<totalAtoms;i++)
		{
			//Try to read a line from the file
			try {fileline=ReadLine();}
			catch(Exception ^ e)
			{
				throw gcnew Exception("Unexpected end of file in Line "+LineCount+":"+e->Message);
			};

            if(fileline==nullptr || fileline->Length<34)
                throw gcnew Exception("Problem in Reading Atom and Coordinates in Line "+LineCount+".");
            Atom^ atom=gcnew Atom();
            atom->X=Double::Parse(fileline->Substring(0,10),System::Globalization::CultureInfo::InvariantCulture);
            atom->Y=Double::Parse(fileline->Substring(10,10),System::Globalization::CultureInfo::InvariantCulture);
            atom->Z=Double::Parse(fileline->Substring(20,10),System::Globalization::CultureInfo::InvariantCulture);
            atom->isotope=Isotopes::FindIsotope(fileline->Substring(31,3)->Trim(),0,0);
            AtomCollection->Add(atom);
		};

	   //Read the bond format: 2  1  1  0  0  0  0          //
		for(int i=0;i<totalBonds;i++)
		{
			fileline=ReadLine();
            int first,second;
            if(fileline==nullptr || fileline->Length<9 ||
                !Int32::TryParse(fileline->Substring(0,3),first) ||
                !Int32::TryParse(fileline->Substring(3,3),second) ||
                first<1 || first>totalAtoms || second<1 || second>totalAtoms)
                throw gcnew Exception("Problem in Reading Bonds in Line "+LineCount+".");
		};
		// Skip optional properties before the mandatory record terminator.
        while((fileline=ReadLine())!=nullptr && fileline->Trim()!="M  END") {}
        OK=(fileline!=nullptr);
        break;

	}

	//Turn on again Refreshing which causes calculation of bonds
	AtomCollection->Refreshing=true;
	ReadingFile->Close();
	return OK;

}

