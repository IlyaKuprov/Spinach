#include "../StdAfx.h"
#include "GAMESS.h"

//Definitions
#define AtomCollection  SystemModel->AtomCollection
//Definitions of keywords-phrases found in this format
#define CMD_CARTCOORD             "***** EQUILIBRIUM GEOMETRY LOCATED *****"


/**
*\file GAMESS.cpp
*/

/**
* @brief Constructor of GAMESS format 
* @param filename Name of the file to be read.
* @param model The model where the imported data will be saved.
*/
GAMESS::GAMESS(System::String^ filename, Model^% model)
{
	SystemModel=model;
	FileName=filename;
	ReadingFile = gcnew StreamReader(filename);
    FormatName="GAMESS";

}

/**
* @brief Main function for opening and reading file
* @return Returns false if the format does not seem to be correct.
*/
bool GAMESS::LoadFile(void)
{
	String^ Command=gcnew String("");
	bool  OK = true;

	//Turn off Refreshing for faster importing
	AtomCollection->Refreshing=false;


	//Initializing flags to false, meaning that no any frame was yet read
	read_onceOR=false;


	//Read the file until the end of it and as long flag OK is true
	while (!ReadingFile->EndOfStream  && OK )
	{

		Command= ReadLine() ;

		bool initial=!read_onceOR && Command->TrimStart()->StartsWith("ATOM") &&
			Command->Contains("ATOMIC") && Command->TrimEnd()->EndsWith("COORDINATES (BOHR)");
		if(initial || Command->Contains(CMD_CARTCOORD))
		{
			OK=read_CARTCOORD(initial);
			read_onceOR=OK;
		}

		//Wrong format inside senctions
		if(!OK)  throw gcnew Exception("Problem in Reading File in Line "+LineCount+".") ;
	};

	AtomCollection->Refreshing=true;
	ReadingFile->Close();
	return read_onceOR;
}

/**
* @brief Main function for reading coordinates from the file
* @return Returns false if the format does not seem to be correct.
*/
bool GAMESS::read_CARTCOORD(bool bohr)
{
	String^ fileline=gcnew String("");  //line
	int AtomID=1;						//Saving the atom ID

	// Initial coordinates have one more header line. Equilibrium coordinates
	// have a heading, column labels and a rule, and are already in angstroms.
	for(int i=0; i<(bohr ? 1 : 3); i++) ReadLine();
	bool found=false;
	while(!ReadingFile->EndOfStream)
	{
		fileline=ReadLine();
		if(String::IsNullOrWhiteSpace(fileline)) break;
		ExtractIsotopeAndCoordinates(fileline,FormatName,read_onceOR,AtomID,
			5,0,2,3,4,bohr ? 0.529177249 : 1.0,AtomCollection);
		found=true;
	}

	return found;

}


