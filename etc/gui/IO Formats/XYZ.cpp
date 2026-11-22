#include "../StdAfx.h"
#include "XYZ.h"

//Definitions
#define AtomCollection  SystemModel->AtomCollection


/**
*\file XYZ.cpp
*/

/**
* Constructor of XYZ format which can use frames
* @param filename Name of the file to be read.
* @param model The model where the imported data will be saved.
*/
SpinachGUI::XYZ::XYZ(String^ filename, Model^% model)
{
	SystemModel=model;
	FileName=filename;
	ReadingFile = gcnew StreamReader(filename);
    FormatName="XYZ";
}

/**
* Main function for opening and reading file
* @return Returns false if the format does not seem to be correct.
*/
bool SpinachGUI::XYZ::LoadFile(void)
{
	bool  OK = true;
	bool read_once=false;	//Initializing flags to false, meaning that no any frame was yet read
	
	//Turn off Refreshing for faster importing
	AtomCollection->Refreshing=false;

	//Read the file until the end of it and as long flag OK is true
	while (!ReadingFile->EndOfStream  && OK )
	{
		int totalAtoms=0;					//total number of atoms in the system
		String^ fileline=gcnew String("");  //line
		int AtomID=1;						//Saving the atom ID
		

		//Read number of atoms, checking first if it is an empty line to ignore
		if((fileline=ReadLine())->Length==0) continue;
		if (!Int32::TryParse(fileline, totalAtoms ) || totalAtoms<0) { OK=false; break; } 

		AtomsDictionary::CheckAtomCount(totalAtoms);
		// The second line is a comment, even when it resembles an atom.
        fileline=ReadLine();
        for(int i=0;i<totalAtoms;i++)
		{
			//Try to read a line from the file
			try {fileline=ReadLine();}
			catch(Exception ^ e)
			{
				throw gcnew Exception("Unexpected end of file in Line "+LineCount+":"+e->Message);
			};

			//Number of Fields: 4, Isotope index: 0, XYZ index: 1,2,3
			ExtractIsotopeAndCoordinates(fileline, FormatName, read_once, AtomID, 4, 0, 1, 2, 3, 1.0, AtomCollection);
		};

		read_once=true;//read at least once
	};

	//Turn on again Refreshing which causes calculation of bonds
	AtomCollection->Refreshing=true;
	ReadingFile->Close();
	return OK;

}

