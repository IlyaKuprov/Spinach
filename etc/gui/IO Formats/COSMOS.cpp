#include "../stdafx.h"
#include "COSMOS.h"


//Definitions
#define InteractionCollection  SystemModel->InteractionCollection
#define AtomCollection  SystemModel->AtomCollection


/**
*\file COSMOS.cpp
*/

/**
* Constructor of COSMOS format
* @param filename Name of the file to be read.
* @param model The model where the imported data will be saved.
*/
COSMOS::COSMOS(System::String^ filename, Model^% model)
{
	SystemModel=model;
	FileName=filename;
	ReadingFile = gcnew StreamReader(filename);
    FormatName="COSMOS";
}

/**
* Main function for opening and reading file
* @return Returns false if the format does not seem to be correct.
*/
bool COSMOS::LoadFile(void)
{
	bool  OK = true;

	//Turn off Refreshing for faster importing
	AtomCollection->Refreshing=false;


	int number=2;						//number of fields in each line
	int totalAtoms=0;					//total number of atoms in the system
	int totalBonds=0;					//total number of bonds in the system
	String^ fileline=gcnew String("");  //line
	cli::array<String^>^numbers = nullptr;   // Matrix of the fields after splitting of a line
	String^ delimStr = " ,";            //seperator of fields
	cli::array<Char>^delimiter = delimStr->ToCharArray(); //Conversion to char array


	//Read the file until the end of it and as long flag OK is true
	while (!ReadingFile->EndOfStream  && OK )
	{
		fileline=ReadLine();

		//Read atom section
		//ATOMS 72
		if(fileline->Contains("ATOMS"))
		{
			//Read number of atoms and Bonds
			number=2;
			numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
			if(number==numbers->Length) 
			{
				if (!Int32::TryParse(numbers[1], totalAtoms ) || totalAtoms<0) OK=false;
			}
			else OK=false;
			if(!OK) { AtomCollection->Refreshing=true; ReadingFile->Close(); return false; }



			number=18;
			//Read the atom format: H12_c1 1  0.82629  7.03135  4.14038 csmb  0.24052 9 2 0 0 0 0 0 0 0 0 0  //
			for(int i=0;i<totalAtoms;i++)
			{
				fileline=ReadLine();

				//Split the line to a number of fields
				numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
				if(number==numbers->Length) 
				{
					try {
						Atom^ nAtom=gcnew Atom();
						int AtomicN;
						double x,y,z;

						//Try to read coordinates and atomic number
						if (!Int32::TryParse(numbers[1], AtomicN) || !Double::TryParse(numbers[2], x) || 
							!Double::TryParse(numbers[3], y ) || !Double::TryParse(numbers[4], z )) 
							throw gcnew Exception("Error in Atom coordinates and atomic number");
						else 
						{
							//Save  coordinates in the model
							nAtom->X=x;nAtom->Y=y;nAtom->Z=z;

							//Check if the second char is lower than n it means that  it 
							//belongs to the name of the element  else also remove it
							if(Char::IsLower(numbers[0][1])) numbers[0]=numbers[0]->Remove(2);
							else numbers[0]=numbers[0]->Remove(1);

							//Can throw exception if does not find the isotope
							try{
								nAtom->isotope=Isotopes::FindIsotope(numbers[0],AtomicN, 0);
							}
							catch(String ^ e)
							{
								nAtom->isotope=Isotopes::FindIsotope(nullptr,AtomicN, 0);
							};
							AtomCollection->Add(nAtom);
						};
					}
					catch(String ^ e)
					{
						//Wrong format of the field, where coordinates not giving number and where 
						//isotope not giving valid isotope name
						throw gcnew Exception("Problem in Reading COSMO format:"+e) ;
					}
				} //Wrong number of fields read
				else throw gcnew Exception("Problem in Reading Atom and Coordinates in Line "+LineCount+".") ;
			};
		}//Reading interactions-couplings section
		else if((fileline->Contains("CALCULATED_PROPERTIES")))
		{
			fileline=ReadLine();

			//Loop over all atoms
			for(int i=0;i<totalAtoms;i++)
			{
				number=2;
				//Split the line to a number of fields
				numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
				if(number<=numbers->Length) 
				{   
					//Set the Atom which the interaction is read
					int atomID;
					Int32::TryParse(numbers[0], atomID);
					Atom ^ A=AtomCollection[atomID];

					//Read until end of senction
					OK=true;
					while (fileline!="END_CALCULATED_PROPERTIES"   && OK )
					{
						//Read line which is shielding or j-couplings or others
						fileline=ReadLine();


						//Read chemical shielding
						if(fileline->Contains("CS-TENSOR"))
						{
							//It is chemical shielding, so conversion to shift is needed
							ToShift= true;

							//Split the line to a number of fields
							number=7;
							numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
							try {

								//Try to read matrix
								::Matrix3x3 ^a_tensor=gcnew ::Matrix3x3();
								a_tensor[5]=0.0f; a_tensor[6]=0.0f; a_tensor[7]=0.0f; 
								double xx,yy,zz,xy,xz,yz;
								if (number!=numbers->Length ||
									!Double::TryParse(numbers[1], xx) || 
									!Double::TryParse(numbers[2], yy) || 
									!Double::TryParse(numbers[3], zz) ||
									!Double::TryParse(numbers[4], xy) ||
									!Double::TryParse(numbers[5], xz) ||
									!Double::TryParse(numbers[6], yz) ) 
									throw gcnew String("Error in CS-Tensor");
								else {
									a_tensor[0]=xx;a_tensor[1]=xy;a_tensor[2]=xz;
									a_tensor[3]=xz;a_tensor[4]=yy;a_tensor[8]=zz;

									//Add chemical shielding interaction, unit and reference frame
									Tensor^ nInteraction=gcnew Tensor(A, A, InteractionKind::Shift, a_tensor);
									nInteraction->unit=Units::ppm;
									nInteraction->Frame=SystemModel->RefFrameCollection[1];
									InteractionCollection->Add(nInteraction);
								};
							}
							catch(String ^ e)
							{	
								//Wrong number of fields read or Wrong format of fields
								throw gcnew Exception("Problem in Reading COSMO format:"+e) ;
							};
						}//Read J-COUPLINGS
						else if(fileline->Contains("J-COUPLINGS"))
						{
							//Split the line to a number of fields
							number=3;
							numbers=fileline->Split(delimiter,StringSplitOptions::RemoveEmptyEntries );
							try {

								int coupled;
								double value;
								Atom ^ B;

								//Try to read matrices for each combination
								::Matrix3x3 ^Matrix=gcnew ::Matrix3x3();
								Matrix[1]=0.0; Matrix[2]=0.0; Matrix[3]=0.0; 
								Matrix[5]=0.0; Matrix[6]=0.0; Matrix[7]=0.0;
								for(int j=0;j<((numbers->Length-1)/2);j++)
								{
									if (!Int32::TryParse(numbers[j*2+1], coupled) || 
										!Double::TryParse(numbers[j*2+2], value) ) 
										throw gcnew Exception("Error in J-COUPLINGS");
									else {	
										Matrix[0]=value; Matrix[4]=value; Matrix[8]=value;
										B=AtomCollection[coupled];

										//Add chemical shielding interaction, unit and reference frame
										Tensor^ nInteraction=gcnew Tensor(A, B, InteractionKind::Jcoupling, Matrix);
										nInteraction->unit=Units::Hz;
										nInteraction->Frame=SystemModel->RefFrameCollection[1];
										InteractionCollection->Add(nInteraction);
									};
								};

							}
							catch(String ^ e)
							{	
								//Wrong format of fields
								throw gcnew Exception("Problem in Reading COSMO format:"+e) ;
							};
						}
						else if(fileline->Contains("RDC-TENSOR") ||
							fileline->Contains("ORDER-TENSOR") ||
							fileline->Contains("QC-TENSOR") ||
							fileline->Contains("DD-TENSOR") ||
							fileline->Contains("DISTANCES"))
						{
							//Do Nothing
						}
						else OK=false; //Do Nothing
					};

				};

			};

		};

	};

	//Turn on again Refreshing which causes calculation of bonds
	AtomCollection->Refreshing=true;
	ReadingFile->Close();
	return true;
}


