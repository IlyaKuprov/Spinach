#include "../stdafx.h"
#include "SpinEvolution.h"

//Definitions
#define To_String(aNumber) SpinachGUI::ConvertToString(aNumber, 4, "N")

/**
*\file SpinEvolution.cpp
*/

/**
* @brief Initializator of Spinevolution class.
* @param filename Name of the file to write into.
* @param model The model which the data will be extracted.
*/
void SpinEvolution::Initialize_SpinEvolution(Model^ model, String^ filename)
{
	SystemModel=model;
	ChanneldataGridView->Rows->Clear();
	TitleFile=filename;

	String^ temp=gcnew String("");			//channel strings
	String^ ChannelList=gcnew String("");   //Total channels string
	int j=0;

	//Creating channel list and dtagrid with channels
	for each(int i in SystemModel->AtomCollection->Keys) 
	{
		//Create the channel string from the elements and mass 
		temp=" "+((Atom^)SystemModel->AtomCollection[i])->isotope->Element+
			((Atom^)SystemModel->AtomCollection[i])->isotope->Mass;

		//Check if does not exist already
		if(!ChannelList->Contains(temp)) 
		{
			//Add line in the datagridview and fill the line with the channel string
			ChanneldataGridView->Rows->Add();
			ChanneldataGridView[1,j]->Value=temp;
			ChanneldataGridView[0,j]->Value=true;
			j++;

			//Add the channel string in the total channels string 
			ChannelList+=temp;
		};
	};
}

/**
* @brief Function for writing in the file all the data.
* @param filename Name of the file to write into.
*/
void SpinEvolution::WriteFile(System::String^ filename)
{
	//Create and open the file
	WritingFile=gcnew FileStream(filename,FileMode::Create);

	//Euler angles object, values and isotropic, anisotropic etc
	EulerAngles^ euleran;
	double alpha, beta, gamma;
	double iso, aniso, hta;

	String ^ FileText=gcnew String(""); //The Text written in the file

	//Some signature staff
	FileText="* Produced by Spinach GUI "+ SpinachGUIver+" \n";
	FileText+="*--------------------------------------------------------------------\n";
	FileText+="*System Name:"+TitleFile+" \n";
	FileText+="**A basic SpinEvolution simulation route is appended after \n";
	FileText+="**'The System' section with random values.\n";
	FileText+="**Please edit as appropriate for your specific case.\n";
	FileText+="*--------------------------------------------------------------------\n";
	FileText+="****** The System *******\n";
	FileText+="spectrometer(MHz)  500\n";
	FileText+="spinning_freq(kHz) 10.0\n";

	//Print Atoms
	String^ AtomList=gcnew String("");
	String^ IndexList=gcnew String("");
	String^ SpinachGUIIndex=gcnew String("");
	String^ ChannelList=gcnew String("");
	String^ temp=gcnew String("");
	String^ temp2=gcnew String("");
	Dictionary <int, int>^ AtomIDs=gcnew Dictionary<int, int>();
	ArrayList^ SortedAtomList=gcnew ArrayList();
	AtomsDictionary^ AtomCollection=gcnew AtomsDictionary(SystemModel->InteractionCollection);
	int j=1;

	//First add atoms in a new array
	for each(int i in SystemModel->AtomCollection->Keys) 
	{
		SortedAtomList->Add(((Atom^)SystemModel->AtomCollection[i]));
	};

	//Sorting according to element name
	AtomCollection=SystemModel->AtomCollection;

	//Loop over all channels
	for(int i=0;i<ChanneldataGridView->Rows->Count;i++) 
	{ 
		//Write in the channel list anly the ones that are on
		if(Convert::ToInt32(ChanneldataGridView[0,i]->Value)) ChannelList+=ChanneldataGridView[1,i]->Value;
	}

	//Loop over all atoms except electrons
	for each(int i in SystemModel->AtomCollection->Keys) 
	{
		if(((Atom^)AtomCollection[i])->isotope->Element!="e" &&
			ChannelList->Contains(((Atom^)AtomCollection[i])->isotope->Element+((Atom^)AtomCollection[i])->isotope->Mass))
		{
			//Create a string with mass number and element, and add in the atomList
			temp=" "+((Atom^)AtomCollection[i])->isotope->Element+
				((Atom^)AtomCollection[i])->isotope->Mass;
			AtomList+=temp;


			//Create spaced String with the same length as the Atom
			temp2 = gcnew String(" ",0, temp->Length);
			temp2=temp2->Insert(1, Convert::ToString(j));
			IndexList+=temp2->Remove(temp->Length); //Reduced it to original size

			//Similarly for writing SpinachGUI Indexes
			temp2 = gcnew String(" ",0, temp->Length);
			temp2=temp2->Insert(1, Convert::ToString(((Atom^)AtomCollection[i])->ID));
			SpinachGUIIndex+=temp2->Remove(temp->Length); //Reduced it to original size

			//Try to save a list of all unique IDs
			try 
			{
				AtomIDs->Add(((Atom^)AtomCollection[i])->ID, j++);
			}
			catch (ArgumentException^)
			{
				MessageBox::Show("AtomID declared twice", "Error", 
					MessageBoxButtons::OK, MessageBoxIcon::Error) ;
			};
		};
	};

	//Write all the above in the stream in the correct format
	FileText+="channels          "+ChannelList+"\n";
	FileText+="nuclei            "+AtomList+"\n";
	FileText+="atomic_coords      ";

	//If more than one atom loop over all atoms and write them
	if(AtomIDs->Count<=1) FileText+="*";
	else 
		for each(int i in SystemModel->AtomCollection->Keys)
			if(AtomIDs->ContainsKey(((Atom^)AtomCollection[i])->ID))
				FileText+=((Atom^)AtomCollection[i])->X.ToString("G17", System::Globalization::CultureInfo::InvariantCulture)+" "+
					((Atom^)AtomCollection[i])->Y.ToString("G17", System::Globalization::CultureInfo::InvariantCulture)+" "+((Atom^)AtomCollection[i])->Z.ToString("G17", System::Globalization::CultureInfo::InvariantCulture)+" ";

	FileText+="\n";

	//Print Interactions
	String ^ CSIsotropic=gcnew String("");
	cli::array<double>^ isotropicShifts=gcnew cli::array<double>(AtomIDs->Count);
	String ^ ChemicalShift=gcnew String("");
	String ^ JCoupling=gcnew String("");
	String ^ Quadrupolar=gcnew String("");
	String ^ A=gcnew String("");
	String ^ B=gcnew String("");

	//Loop over all interactions
	for each(int i in SystemModel->InteractionCollection->Keys) 
	{
		//Short name for the A Atom string
		A=((Interaction^)SystemModel->InteractionCollection[i])->A->isotope->Element+
			((Interaction^)SystemModel->InteractionCollection[i])->A->isotope->Mass;

		//Short name for the B Atom string
		if(((Interaction^)SystemModel->InteractionCollection[i])->B!=nullptr)
			B=((Interaction^)SystemModel->InteractionCollection[i])->B->isotope->Element+
			((Interaction^)SystemModel->InteractionCollection[i])->B->isotope->Mass;
		else B="dummy";

		//For Chemical shifts, jcouplings and quadropolar interactions and when channel is activated
		if((((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Shift ||
			((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Jcoupling||
			((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Quadrupolar) &&
			ChannelList->Contains(A))
		{
			//Calculate euler angles object form this function
			euleran=EulerAnglesfromEigenVectors(((Interaction^)SystemModel->InteractionCollection[i])->dcm);

			//Convert angles to degrees 
			alpha=float(180.0f*euleran->alpha/PI);
			beta=float(180.0f*euleran->beta/PI);
			gamma=float(180.0f*euleran->gamma/PI);

			//Isotorpic Factor
			iso=(((Interaction^)SystemModel->InteractionCollection[i])->eigen[0]+
				((Interaction^)SystemModel->InteractionCollection[i])->eigen[1]+
				((Interaction^)SystemModel->InteractionCollection[i])->eigen[2])/3.0f;

			double yy=1000000.0f,xx,zz=iso;
			//Find from abs(s_zz-s_iso)>=abs(s_xx-s_iso)>=abs(s_yy-s_iso)
			for(int jj=0;jj<3;jj++)
			{
				if(abs(((Interaction^)SystemModel->InteractionCollection[i])->eigen[jj]-iso)>=abs(zz-iso))
					zz=((Interaction^)SystemModel->InteractionCollection[i])->eigen[jj];
				if(abs(((Interaction^)SystemModel->InteractionCollection[i])->eigen[jj]-iso)<=abs(yy-iso))
					yy=((Interaction^)SystemModel->InteractionCollection[i])->eigen[jj];
			};
			for(int jj=0;jj<3;jj++)
			{
				if(((Interaction^)SystemModel->InteractionCollection[i])->eigen[jj]!=zz &&
					((Interaction^)SystemModel->InteractionCollection[i])->eigen[jj]!=yy)
					xx=((Interaction^)SystemModel->InteractionCollection[i])->eigen[jj];
			};

			//Anisotorpic Factor etc
			aniso=zz-iso;
			if(abs(aniso)<1.0e-10) aniso=hta=0.0f;
			else hta=(yy-xx)/aniso;

			//Write Chemical shift in appropriate format
			if(((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Shift)
			{
				isotropicShifts[AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->A->ID]-1]=iso;
				ChemicalShift+="csa_parameters     "+AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->A->ID]+
					" "+To_String(aniso)+" "+To_String(hta)+" "+
					To_String(alpha)+" "+ To_String(beta)+" "+ To_String(gamma)+" ppm \n";
			};

			//Write J-coupling in appropriate format
			if(((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Jcoupling
				&& ChannelList->Contains(B))
			{
				JCoupling+="j_coupling         "+
					AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->A->ID]+" "+
					AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->B->ID]+" "+
					To_String(iso)+" "+ To_String(aniso)+" "+To_String(hta)+" "+
					To_String(alpha)+" "+ To_String(beta)+" "+ To_String(gamma)+" \n";
			};

			//Write Quadrupolar coupling in appropriate format
			if(((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Quadrupolar)
			{
				Quadrupolar+="quadrupole         "+
					AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->A->ID]+" "+
					To_String(aniso)+" "+
					To_String(hta)+" "+
					To_String(alpha)+" "+ To_String(beta)+" "+ To_String(gamma)+" \n";
			};
		} 	
	};

	//Write shifts in exported nucleus order, including zero for absent shifts.
	for each(double shift in isotropicShifts) CSIsotropic+=To_String(shift)+" ";
	//Write all the above in the stream in the correct format
	FileText+="cs_isotropic       "+(CSIsotropic->Length>1 ? CSIsotropic+" ppm" : "*")+" \n";
	FileText+=(ChemicalShift->Length>1 ? ChemicalShift : "csa_parameters     * \n");
	FileText+=(JCoupling->Length>1 ? JCoupling : "j_coupling         * \n");
	FileText+=(Quadrupolar->Length>1 ? Quadrupolar : "quadrupole         * \n");

	//This is a basic Spinevolution parameters
	FileText+="dip_switchboard    *\n";
	FileText+="csa_switchboard    *\n";
	FileText+="exchange_nuclei    *\n";
	FileText+="bond_len_nuclei    *\n";
	FileText+="bond_ang_nuclei    *\n";
	FileText+="tors_ang_nuclei    *\n";
	FileText+="groups_nuclei      *\n";

	//Printing of the pulse sequence according to the channels we have activated
	FileText+="******* Pulse Sequence ******************************\n";
	int jj=1;
	for(int i=0;i<ChanneldataGridView->Rows->Count;i++) 
	{
		if(ChannelList->Contains(Convert::ToString(ChanneldataGridView[1,i]->Value)))
		{
			FileText+="CHN "+jj+"\n";
			FileText+="timing(usec)       (4   4)\n";
			FileText+="power(kHz)         125 125\n";  
			FileText+="phase(deg)         0   15\n"; 
			FileText+="freq_offs(kHz)     0   0\n";
			jj++;
		};
	};

	//This is a basic Spinevolution parameters which append at the end of our data in order to run in the SpinEvolution  
	FileText+="******* Variables ***********************************\n";
	FileText+="******* Options *************************************\n";
	FileText+="rho0               F1x\n";
	FileText+="observables        F1p\n";
	FileText+="EulerAngles        *\n";
	FileText+="n_gamma            *\n";
	FileText+="line_broaden(Hz)   *\n";
	FileText+="zerofill           *\n";
	FileText+="FFT_dimensions     *\n";
	FileText+="options            *\n";
	FileText+="******************************************************\n";

	//Convert and write to the file      
	cli::array<Byte>^ info = (gcnew UTF8Encoding( true ))->GetBytes(FileText);
	WritingFile->Write( info, 0, info->Length  );
	WritingFile->Close();
}

/**
* @brief EventHandler for opening the save dialoguebox and  writing the stream to file.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
Void SpinEvolution::Exportbutton_Click(System::Object^  sender, System::EventArgs^  e) 
{
	//Setting the format and initila filename for the save dialoguebox 
	saveFileDialog->Filter="SpinEvolution File|*.in";
	saveFileDialog->FileName=Path::GetFileNameWithoutExtension(TitleFile)+".in";

	//Opens the save dialogue box and check if the OK button clicked then save to file
	if (saveFileDialog->ShowDialog() == System::Windows::Forms::DialogResult::OK )
		WriteFile(saveFileDialog->FileName);
}