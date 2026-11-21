#include "../StdAfx.h"
#include "Simpson.h"

//Definitions
#define To_String(aNumber) SpinachGUI::ConvertToString(aNumber, 4, "N")

/**
*\file Simpson.cpp
*/

/**
* @brief Initializator of Simpson class.
* @param filename Name of the file to write into.
* @param model The model which the data will be extracted.
*/
void Simpson::Initialize_Simpson(Model^ model, String^ filename)
{
	SystemModel=model;
	ChanneldataGridView->Rows->Clear();
	FileName=filename;

	String^ temp=gcnew String("");			//channel strings
	String^ ChannelList=gcnew String("");   //Total channels string
	int j=0;

	//Creating channel list and dtagrid with channels
	for each(int i in SystemModel->AtomCollection->Keys) 
	{
		//Create the channel string from the mass and elements
		temp=" "+((Atom^)SystemModel->AtomCollection[i])->isotope->Mass+
			((Atom^)SystemModel->AtomCollection[i])->isotope->Element;

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
void Simpson::WriteFile(String^ filename)
{
	//Create and open the file
	WritingFile=gcnew FileStream(filename,FileMode::Create);

	//Euler angles object, values and isotropic, anisotropic etc
	EulerAngles^ euleran;
	double alpha, beta, gamma;
	double iso, aniso, hta;

	//Quadropolar order read by the combo box option
	int QuadOrder=(OrdercomboBox->Text=="1st"? 1: 2);

	String ^ FileText=gcnew String("");  //The Text written in the file

	//Some signature staff
	FileText="# Produced by Spinach GUI "+ SpinachGUIver+" \r\n \r\n";
	FileText+="#--------------------------------------------------------------------\r\n";
	FileText+="#System Name:"+FileName+" \r\n";
	FileText+="#--------------------------------------------------------------------\r\n \r\n";
	FileText+="spinsys { \r\n";

	//Print Atoms
	String^ AtomList=gcnew String("");
	String^ IndexList=gcnew String("");
	String^ SpinachGUIIndex=gcnew String("");
	String^ ChannelList=gcnew String("");
	String^ temp=gcnew String("");
	String^ temp2=gcnew String("");
	Dictionary <int, int>^ AtomIDs=gcnew Dictionary<int, int>();
	ArrayList^ SortedAtomList=gcnew ArrayList();
	int j=1;

	//First add atoms in a new array
	for each(int i in SystemModel->AtomCollection->Keys) 
	{
		SortedAtomList->Add(((Atom^)SystemModel->AtomCollection[i]));
	};

	//Sort a view only; do not renumber the source atoms or add source bonds.
	if(SortcheckBox->Checked) SortedAtomList->Sort(gcnew ComparerElement());

	//Loop over all atoms except electrons
	for each(Atom^ atom in SortedAtomList)
	{
		if(atom->isotope->Element!="e")
		{
			//Create a string with mass number and element, and add in the atomList
			temp=" "+atom->isotope->Mass+
				atom->isotope->Element;
			AtomList+=temp;


			//Create spaced String with the same length as the Atom
			temp2 = gcnew String(" ",0, temp->Length);
			temp2=temp2->Insert(1, Convert::ToString(j));
			IndexList+=temp2->Remove(temp->Length); //Reduced it to original size

			//Similarly for writing SpinachGUI Indexes
			temp2 = gcnew String(" ",0, temp->Length);
			temp2=temp2->Insert(1, Convert::ToString(atom->getID()));
			SpinachGUIIndex+=temp2->Remove(temp->Length); //Reduced it to original size

			//Try to save a list of all unique IDs
			try 
			{
				AtomIDs->Add(atom->ID, j++);
			}
			catch (ArgumentException^)
			{
				MessageBox::Show("AtomID declared twice", "Error", 
					MessageBoxButtons::OK, MessageBoxIcon::Error) ;
			};
		};
	};

	//Loop over all channels
	for(int i=0;i<ChanneldataGridView->Rows->Count;i++) 
	{ 
		//Write in the channel list anly the ones that are on
		if(Convert::ToInt32(ChanneldataGridView[0,i]->Value)) ChannelList+=ChanneldataGridView[1,i]->Value;
	}

	//Write all the above in the stream in the correct format
	FileText+="	channels "+ChannelList+"\r\n";
	FileText+="	nuclei "+AtomList+"\r\n";
	FileText+="# Atom Index "+ IndexList+"\r\n";
	FileText+="# SpinachGUI "+ SpinachGUIIndex+"\r\n";

	//For printing Interactions
	String ^ ChemicalShift=gcnew String("");
	String ^ JCoupling=gcnew String("");
	String ^ Dipolar=gcnew String("");
	String ^ Quadrupolar=gcnew String("");

	//Loop over all interactions
	for each(int i in SystemModel->InteractionCollection->Keys)
	{
		//For Chemical shifts, jcouplings and quadropolar interactions
		if(((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Shift ||
			((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Jcoupling||
			((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Quadrupolar)
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
				ChemicalShift+="	shift "+AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->A->ID]+" "+
					To_String(iso)+"p "+ To_String(aniso)+"p "+To_String(hta)+" "+
					To_String(alpha)+" "+ To_String(beta)+" "+ To_String(gamma)+" \r\n";
			};

			//Write J-coupling in appropriate format
			if(((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Jcoupling)
			{
				JCoupling+="	jcoupling "+
					AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->A->ID]+" "+
					AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->B->ID]+" "+
					To_String(iso)+" "+ To_String(aniso)+" "+To_String(hta)+" "+
					To_String(alpha)+" "+ To_String(beta)+" "+ To_String(gamma)+" \r\n";
			};

			//Write Quadrupolar coupling in appropriate format
			if(((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::Quadrupolar)
			{
				//nuclear spin quantum number
				double NSQN=((Interaction^)SystemModel->InteractionCollection[i])->A->isotope->Spin; 
				double C_q;

				//Calculation of C_q factor =2I*(2I-1)*1.0e-6 Hz
				C_q=((Tensor^)SystemModel->InteractionCollection[i])->matrix3x3[8]*
					(2*NSQN*(2*NSQN-1))*1.0e-6;


				Quadrupolar+="	quadrupole "+
					AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->A->ID]+" "+
					Convert::ToString(QuadOrder)+" "+ 
					To_String(C_q)+" "+
					To_String(hta)+" "+
					To_String(alpha)+" "+ To_String(beta)+" "+ To_String(gamma)+" \r\n";
			};

		} //For writing dipolar couplings that are higher than the threshold
		else if(((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::CBond
			&& ((Tensor^)SystemModel->InteractionCollection[i])->matrix3x3[0]<=Convert::ToDouble(DipolarThreshold->Value))
		{
			//Calculate euler angles object form this function
			euleran=EulerAnglesfromEigenVectors(((Interaction^)SystemModel->InteractionCollection[i])->dcm);

			alpha=0.0f;
			//Calculate the distance in each direction
			double x=((Tensor^)SystemModel->InteractionCollection[i])->A->X-((Tensor^)SystemModel->InteractionCollection[i])->B->X;
			double y=((Tensor^)SystemModel->InteractionCollection[i])->A->Y-((Tensor^)SystemModel->InteractionCollection[i])->B->Y;
			double z=((Tensor^)SystemModel->InteractionCollection[i])->A->Z-((Tensor^)SystemModel->InteractionCollection[i])->B->Z;

			//For bonds-distnces larger than zero
			if(((Tensor^)SystemModel->InteractionCollection[i])->matrix3x3[0] > std::numeric_limits<double>::epsilon())
				beta=acos(z/((Tensor^)SystemModel->InteractionCollection[i])->matrix3x3[0])*(180/PI); //Theta and convert to degrees
			else beta=0.0;

			gamma=atan2(y,x)*(180/PI);//Fai and convert to degrees
			double bij=-1.054571726e-11* //h_bar*e-7/(amstrong*amstrong*amstrong)
				((Tensor^)SystemModel->InteractionCollection[i])->A->isotope->Magnetogyric*
				((Tensor^)SystemModel->InteractionCollection[i])->B->isotope->Magnetogyric
				/(2*PI*((Tensor^)SystemModel->InteractionCollection[i])->matrix3x3[0]*
				((Tensor^)SystemModel->InteractionCollection[i])->matrix3x3[0]*
				((Tensor^)SystemModel->InteractionCollection[i])->matrix3x3[0]);

			//Write dipolar coupling in appropriate format
			Dipolar+="	dipole "+
				AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->A->ID]+" "+
				AtomIDs[((Tensor^)SystemModel->InteractionCollection[i])->B->ID]+" "+
				To_String(bij)+" "+ 
				To_String(alpha)+" "+ To_String(beta)+" "+ To_String(gamma)+" \r\n";
		};

	};

	//Write all the above in the stream in the correct format
	FileText+=ChemicalShift;
	FileText+=JCoupling;
	FileText+=Quadrupolar;
	FileText+=Dipolar;
	FileText+="\r\n  \r\n } \r\n";

	//This is a basic simpson script which append at the end of our data in order to run in the simpson
	FileText+="#--------------------------------------------------------------------\r\n";
	FileText+="# A basic SIMPSON simulation route is appended below,\r\n";
	FileText+="# please edit as appropriate for your specific case.\r\n";
	FileText+="#--------------------------------------------------------------------\r\n";

	FileText+="\r\n";
	FileText+="par {\r\n";
	FileText+="proton_frequency 600e6 \r\n";
	FileText+="spin_rate        0 \r\n";
	FileText+="sw               20000 \r\n";
	FileText+="start_operator   I1x \r\n";
	FileText+="detect_operator  I1p \r\n";
	FileText+="gamma_angles     1 \r\n";
	FileText+="verbose          1101 \r\n";
	FileText+="np               256 \r\n";
	FileText+="crystal_file     zcw4180 \r\n";
	FileText+="} \r\n";

	FileText+="\r\n";
	FileText+="proc pulseq {} {\r\n";
	FileText+="     acq_block { \r\n";
	FileText+="         delay 20 \r\n";
	FileText+="     } \r\n";
	FileText+="} \r\n";


	FileText+="\r\n";
	FileText+="proc main {} {\r\n";
	FileText+="global par \r\n";
	FileText+="set f [fsimpson] \r\n";
	FileText+="fsave $f $par(name).fid \r\n";
	FileText+="faddlb $f 50 0 \r\n";
	FileText+="fzerofill $f 4096 \r\n";
	FileText+="fft $f \r\n";
	FileText+="fsave $f $par(name).spe \r\n";
	FileText+="funload $f \r\n";
	FileText+="} \r\n";

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
Void Simpson::Exportbutton_Click(System::Object^  sender, System::EventArgs^  e)
{
	//Setting the format and initila filename for the save dialoguebox 
	saveFileDialog->Filter="Simpson File|*.in";
	saveFileDialog->FileName=Path::GetFileNameWithoutExtension(FileName)+".in";

	//Opens the save dialogue box and check if the OK button clicked then save to file
	if (saveFileDialog->ShowDialog() == System::Windows::Forms::DialogResult::OK )
		WriteFile(saveFileDialog->FileName);
}



