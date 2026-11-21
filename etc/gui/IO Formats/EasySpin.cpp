#include "../StdAfx.h"
#include "EasySpin.h"

//Definitions
#define To_String(aNumber) SpinachGUI::ConvertToString(aNumber, 6, "N")
#define logtext(text)  EasySpinNumber(pow(10.0,EasySpinDouble(text)))
#define InteractionCollection SystemModel->InteractionCollection
#define AtomCollection SystemModel->AtomCollection

//Text fields and generated MATLAB code use decimal points, independent of locale.
static double EasySpinDouble(String^ value)
{
    return Double::Parse(value, System::Globalization::CultureInfo::InvariantCulture);
}

static String^ EasySpinNumber(double value)
{
    return value.ToString("G17", System::Globalization::CultureInfo::InvariantCulture);
}

//A red field reports malformed pasted/typed input without throwing from a UI event.
static bool EasySpinValidNumber(Object^ sender)
{
    System::Windows::Forms::TextBox^ box=dynamic_cast<System::Windows::Forms::TextBox^>(sender);
    double value;
    if(box!=nullptr && box->Text!="" &&
       (!Double::TryParse(box->Text, System::Globalization::NumberStyles::Float,
                         System::Globalization::CultureInfo::InvariantCulture, value) ||
        Double::IsNaN(value) || Double::IsInfinity(value)))
    {
        box->BackColor=System::Drawing::Color::MistyRose;
        return false;
    }
    return true;
}

/**
*\file EasySpin.cpp
*/

/**
* @brief Initializator of EasySpin class.
* @param aCurrentTitleFile Name of the file to write into.
* @param model The model which the data will be extracted.
*/
System::Void SpinachGUI::EasySpin::EasySpin_Initialization(Model^% model, String^ aCurrentTitleFile)
{
	//Set initial valuables
	SystemModel=model;
	CurrentTitleFile=aCurrentTitleFile;
	PreviewtextBox->AutoSize=true;

	//Set the position of panels
	DynamicsgroupBox->Size=System::Drawing::Size(396,220);
	RotationallogcheckBox->Location=System::Drawing::Point(12, 185);
	OrderinggroupBox->Location=System::Drawing::Point(3, 225);
	Heisenbergpanel->Location=System::Drawing::Point(6, 335);

	//Initial state garlic etc
	radioButton_CheckedChanged(gcnew RadioButton (), gcnew System::EventArgs());
	radioButton_diffusion_CheckedChanged(gcnew RadioButton (), gcnew System::EventArgs());
	ExperimenttextBoxKeyPress(gcnew TextBox(),gcnew KeyPressEventArgs((char)13));
	Exp_comboBox_SelectedIndexChanged(gcnew ComboBox(),gcnew System::EventArgs());
	Expnumeric_ValueChanged(gcnew NumericUpDown(),gcnew System::EventArgs());
	CrystalSymmetryBoxKeyPress(gcnew TextBox(),gcnew KeyPressEventArgs((char)13));

	//Initial state of options
	Optioncombo_SelectedIndexChanged(gcnew ComboBox(), gcnew System::EventArgs()) ;
	Optnumeric_ValueChanged(gcnew NumericUpDown(), gcnew System::EventArgs());
	Opt_textBox_KeyPress(gcnew TextBox(), gcnew KeyPressEventArgs((char)13));

	//For conversion of convolutions
	TextBox^ textbox=gcnew TextBox();
	textbox->Text="0";

	//Initial state for convolutions in broadening tab
	textbox->Tag=Convert::ToString(1);
	convolutionBoxKeyPress(textbox, gcnew KeyPressEventArgs((char)13));
	textbox->Tag=Convert::ToString(3);
	convolutionBoxKeyPress(textbox, gcnew KeyPressEventArgs((char)13));
	
	//Initial dynamics tabs parameters
    Anisotropic_text_KeyPress(gcnew TextBox(), gcnew KeyPressEventArgs((char)13));
	textBoxDynamics_KeyPress(gcnew TextBox(), gcnew KeyPressEventArgs((char)13));
    orderheisentext_KeyPress(gcnew TextBox(), gcnew KeyPressEventArgs((char)13));

	//Update the final strings
	WriteHeader();
	WriteFooter();
}

/**
* @brief Write begin of the file(Atoms and interactions) with the Easyspin format.
*/
System::Void SpinachGUI::EasySpin::WriteHeader()
{
	String ^ FileText=gcnew String("");
	Dictionary<int, String^>^ A=gcnew Dictionary<int, String^>();
	Dictionary<int, String^>^ Q=gcnew Dictionary<int, String^>();

	//Initial Name etc
	FileText="clear Sys Exp Opt; \r\n";
	FileText+="Sys.S=1/2; \r\n";

	//Loop over all interactions
	for each(int i in InteractionCollection->Keys)
	{
		//Write nucleus and HFC
		if(((Interaction^)InteractionCollection[i])->IntKind==InteractionKind::HFC)
		{
			String ^ astring=gcnew String("");
			::Matrix3x3^ matrix3x3=TransformToLab(((Tensor^)InteractionCollection[i]));

			Unit^ unit=InteractionCollection[i]->unit;
			bool field=(unit==nullptr || unit==Units::Unknown || unit==Units::Unitless ||
				unit==Units::Gauss || unit==Units::mTesla);
			double scale=field ? (unit==Units::mTesla ? 1.0 : 0.1) : Units::MHz->fromSI(unit->toSI(1.0));
			astring="Sys = nucspinadd(Sys, '"+InteractionCollection[i]->A->isotope->Mass+
				InteractionCollection[i]->A->isotope->Element+"', ";
			if(field) astring+="mt2mhz(";
			astring+="[";
			for(int j=0;j<9;j++)
			{
				if(j>0) astring+=(j%3==0 ? "; " : ", ");
				astring+=To_String(scale*matrix3x3[j]);
			}
			astring+="]";
			if(field) astring+=")";

			//Add to the final string
			A[((Interaction^)InteractionCollection[i])->A->ID]=astring;
		};

		//Write quadropolar
		if(((Interaction^)InteractionCollection[i])->IntKind==InteractionKind::Quadrupolar)
		{
			String ^ qstring=gcnew String("");
			::Matrix3x3^ matrix3x3=TransformToLab(((Tensor^)InteractionCollection[i]));

			Unit^ unit=InteractionCollection[i]->unit;
			double scale=(unit==nullptr || unit==Units::Unitless || unit==Units::Unknown) ?
				1.0 : Units::MHz->fromSI(unit->toSI(1.0));
			qstring=", [], [";
			for(int j=0;j<9;j++)
			{
				if(j>0) qstring+=(j%3==0 ? "; " : ", ");
				qstring+=To_String(scale*matrix3x3[j]);
			}
			qstring+="]";

			//Add to the final string
			Q[((Interaction^)InteractionCollection[i])->A->ID]=qstring;
		};

		//Write G-Tensor
		if(((Interaction^)InteractionCollection[i])->IntKind==InteractionKind::GTensor)
		{
			::Matrix3x3^ matrix3x3=TransformToLab(((Tensor^)InteractionCollection[i]));
			FileText+="Sys.g = ["+
				To_String(matrix3x3[0])+", "+
				To_String(matrix3x3[1])+", "+
				To_String(matrix3x3[2])+"; "+
				To_String(matrix3x3[3])+", "+
				To_String(matrix3x3[4])+", "+
				To_String(matrix3x3[5])+"; "+
				To_String(matrix3x3[6])+", "+
				To_String(matrix3x3[7])+", "+
				To_String(matrix3x3[8])+" "+
				"];\r\n";
		};

		//Write Zero-Splitting
		if(((Interaction^)InteractionCollection[i])->IntKind==InteractionKind::ZFS)
		{
			//Transform to Mhz
			Unit^ funit=Units::MHz;
			::Matrix3x3^ matrix3x3=TransformToLab(((Tensor^)InteractionCollection[i]));
			double con=funit->fromSI(((Tensor^)InteractionCollection[i])->unit->toSI(1.0f));

			//Write matrix
			FileText+="Sys.D = ["+
				To_String(con*matrix3x3[0])+", "+
				To_String(con*matrix3x3[1])+", "+
				To_String(con*matrix3x3[2])+"; "+
				To_String(con*matrix3x3[3])+", "+
				To_String(con*matrix3x3[4])+", "+
				To_String(con*matrix3x3[5])+"; "+
				To_String(con*matrix3x3[6])+", "+
				To_String(con*matrix3x3[7])+", "+
				To_String(con*matrix3x3[8])+" "+
				"];\r\n";
		};
	};

	//Emit every nucleus once, matching tensors by atom identity rather than order.
	for each(int i in AtomCollection->Keys)
	{
		Atom^ atom=(Atom^)AtomCollection[i];
		if(atom->isotope->Element=="e") continue;
		String^ hyperfine;
		if(!A->TryGetValue(atom->ID, hyperfine))
			hyperfine="Sys = nucspinadd(Sys, '"+atom->isotope->Mass+
				atom->isotope->Element+"', [0 0 0]";
		String^ quadrupole;
		FileText+=hyperfine+(Q->TryGetValue(atom->ID,quadrupole) ? quadrupole : "")+");\r\n";
	}

	//Write in the final string
	PreviewtextBox->Text=FileText+"\r\n";
}

/**
* @brief Write end of the file(Parameters) with the  Easyspin format.
*/
System::Void SpinachGUI::EasySpin::WriteFooter()
{
	String ^ FileText=gcnew String("");

	//Write the Header for the parameters
	if (easyspinfunction==0) //garlic
	{
		FileText+="\r\n garlic(Sys,Exp,Opt);";
	}
	else if (easyspinfunction==1) //chili
	{
		FileText+="\r\n chili(Sys,Exp,Opt);";
	}
	else if (easyspinfunction==2) //pepper
	{
		FileText+="\r\n pepper(Sys,Exp,Opt);";
	};

	//write all the parameters strings
	PreviewtextBox->Text+=experimentstring;
	PreviewtextBox->Text+=experiment1string;
	PreviewtextBox->Text+=experiment2string;
	if(easyspinfunction==2) PreviewtextBox->Text+=symmetrystring;
	PreviewtextBox->Text+=optionstring;
	PreviewtextBox->Text+=option1string;
	PreviewtextBox->Text+=option2string;
	PreviewtextBox->Text+=broadeningstring;
	if(easyspinfunction==2) PreviewtextBox->Text+=broadening1string;
	if(easyspinfunction!=2) PreviewtextBox->Text+=dynamicstring;
	if(easyspinfunction==1) PreviewtextBox->Text+=dynamic1string;
	PreviewtextBox->Text+=FileText;


	int maxWidth=0;
	int fontsize=11;

	//Find the maximum line width
	for(int i=0;i<PreviewtextBox->Lines->Length;i++) 
		if(PreviewtextBox->Lines[i]->Length>maxWidth) maxWidth=PreviewtextBox->Lines[i]->Length;

	//Change the size of the preview panel
	PreviewtextBox->Width=fontsize*maxWidth;
	PreviewtextBox->Height=PreviewtextBox->Lines->Length*(fontsize+7);
}

#pragma region Others

/**
* @brief Convert units of a quantity.
* @param value The number to convert. String form to double.
* @param result The unit we want to convert to.
* @param selected The initially unit in string form.
* result The converted number in string form.
*/
String^ SpinachGUI::EasySpin::unitconv(String^ value, Unit^ result, String^ selected)
{
	//Set the initial unit
	Unit^ unit;
	if (selected=="ns" ) unit=Units::nsec;
	else if (selected=="mTesla" ) unit=Units::mTesla;
	else if (selected=="Gauss" ) unit=Units::Gauss;
	else if (selected=="MHz" ) unit=Units::MHz;
	else if (selected=="GHz" ) unit=Units::GHz;
	else if (selected=="Hz" ) unit=Units::Hz;
	else if (selected=="Seconds" ) unit=Units::seconds;
	else  unit=Units::mi_sec;/*if (selected=="?s" )*/

	
	double nvalue;
	//Convert number to SI units
	if(value=="") return ""; 
	else nvalue=unit->toSI(EasySpinDouble(value));

	//Convert to final unit
	return EasySpinNumber(result->fromSI(nvalue));
}

/**
* @brief EventHandler for radio button of EasySpin Function
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::radioButton_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{
	if(!((RadioButton^)sender)->Checked) return;
	//Choose between garlic, chili, pepper and activate  proper  elements
	switch(Convert::ToInt32((String^)(((RadioButton^)sender)->Tag))) {
	case 1:  //garlic
		pepperBroadeningspanel->Visible=false;
		Optionsgarlicpanel->Visible=true;
		Optionsgarlicpanel->Location=System::Drawing::Point(6, 3);
		Optionschilipanel->Visible=false;
		Optionspepperpanel->Visible=false;
		Isotropicradio->Enabled=true;
		Isotropicradio->Checked=true;
		Axialradio->Enabled=false;
		Rhombicradio->Enabled=false;
		IsotropicDiffussionpanel->Visible=false;
		axialpanel->Enabled=true;
		rhombicpanel->Enabled=true;
		isotropicpanel->Enabled=true;
		RotationallogcheckBox->Enabled=true;
		OrderinggroupBox->Visible=false;
		Heisenbergpanel->Visible=false;
		PepperExppanel->Visible=false;
		easyspinfunction=0; //garlic
		break;
	case 2:  //chili
		pepperBroadeningspanel->Visible=false;
		Optionsgarlicpanel->Visible=false;
		Optionschilipanel->Location=System::Drawing::Point(6, 3);
		Optionschilipanel->Visible=true;
		Optionspepperpanel->Visible=false;
		Isotropicradio->Enabled=true;
		Axialradio->Enabled=true;
		Rhombicradio->Enabled=true;
		IsotropicDiffussionpanel->Visible=true;
		axialpanel->Enabled=true;
		rhombicpanel->Enabled=true;
		isotropicpanel->Enabled=true;
		RotationallogcheckBox->Enabled=true;
		OrderinggroupBox->Visible=true;
		Heisenbergpanel->Visible=true;
		PepperExppanel->Visible=false;
		easyspinfunction=1; //chili
		break;
	case 3:  //pepper
		pepperBroadeningspanel->Visible=true;
		Optionsgarlicpanel->Visible=false;
		Optionschilipanel->Visible=false;
		Optionspepperpanel->Location=System::Drawing::Point(6, 3);
		Optionspepperpanel->Visible=true;
		Isotropicradio->Enabled=false;
		Axialradio->Enabled=false;
		Rhombicradio->Enabled=false;
		IsotropicDiffussionpanel->Visible=false;
		axialpanel->Enabled=false;
		rhombicpanel->Enabled=false;
		isotropicpanel->Enabled=false;
		RotationallogcheckBox->Enabled=false;
		OrderinggroupBox->Visible=false;
		Heisenbergpanel->Visible=false;
		PepperExppanel->Visible=true;
		easyspinfunction=2; //pepper
		break;

	};

	//Rebuild all mode-dependent fields rather than keeping the previous mode.
	Optioncombo_SelectedIndexChanged(sender, e);
	Optnumeric_ValueChanged(sender, e);
	Opt_textBox_KeyPress(gcnew TextBox(), gcnew KeyPressEventArgs((char)13));
	Exp_comboBox_SelectedIndexChanged(sender, e);
	CrystalSymmetryBoxKeyPress(gcnew TextBox(), gcnew KeyPressEventArgs((char)13));
	Anisotropic_text_KeyPress(gcnew TextBox(), gcnew KeyPressEventArgs((char)13));

	//Update the final strings
	WriteHeader();
	WriteFooter();
}

/**
* @brief EventHandler for Export button of EasySpin 
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Exportbutton_Click(System::Object^  sender, System::EventArgs^  e) 
{		
	//Setting the format and initialize filename for the save dialoguebox 
	saveFileDialog->Filter="EasySpin File|*.m";
	saveFileDialog->FileName=Path::GetFileNameWithoutExtension(CurrentTitleFile)+".m";

	//Opens the save dialogue box and check if the OK button clicked then save to file
	if (saveFileDialog->ShowDialog() == System::Windows::Forms::DialogResult::OK )
	{	 
		//Create export file
		FileStream^ fs=gcnew FileStream(saveFileDialog->FileName,FileMode::Create );

		//Convert and write to file
		cli::array<Byte>^ info = (gcnew UTF8Encoding( true ))->GetBytes( PreviewtextBox->Text);
		fs->Write( info, 0, info->Length  );

		fs->Close();
	};
}

#pragma endregion

#pragma region Experiment tab

/**
* @brief EventHandler for top text boxes in Experiment tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::ExperimenttextBoxKeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	if(e->KeyChar==(char)13 && !EasySpinValidNumber(sender)) return;
	if(e->KeyChar==(char)13)
		for each(TextBox^ box in gcnew cli::array<TextBox^>{Centertext,Sweeptext,Minimumtext,Maximumtext,Frequencytext,Modetext})
			if(!EasySpinValidNumber(box)) return;
	//Clear the string
	experimentstring="";

	//If key pressed is not number, or . or backspace or enter do nothing
	if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '-' && e->KeyChar != '+' && e->KeyChar != 'e' && e->KeyChar != 'E' && e->KeyChar != '\b' && e->KeyChar !=(char)13)
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted
		if (e->KeyChar ==(char)13)
		{
			//Returned to original color of the textbox
			((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

			//If one the Center and Sweep text fields is not empty, update text field
			if(Centertext->Text!="" && Sweeptext->Text!="")
			{
				//Write the center-sweep parameters after convertion to proper units
				experimentstring="Exp.CenterSweep = ["+unitconv(Centertext->Text,Units::mTesla,ExpUnitscombo->Text)+
					" "+ unitconv(Sweeptext->Text,Units::mTesla,ExpUnitscombo->Text)+"];\r\n";
			}
			else //If one the min-max text fields is not empty, update text field
			{
				if(Minimumtext->Text!="" && Maximumtext->Text!="")

					//Write the min-max parameters after convertion to proper units
						experimentstring="Exp.Range = ["+unitconv(Minimumtext->Text,Units::mTesla,ExpUnitscombo->Text)+
						" "+ unitconv(Maximumtext->Text,Units::mTesla,ExpUnitscombo->Text)+"];\r\n"; 
			}

			//Write the microwave frequence value in the aproppriate format
			experimentstring+="Exp.mwFreq = "+ Frequencytext->Text+";\r\n";

			//Write the Phase(degrees) value in the aproppriate format
			experimentstring+="Exp.mwPhase = "+ Phasetext->Text+";\r\n";

			//Write the Mod amp value in the aproppriate format
			experimentstring+=(ModcheckBox->Checked? "Exp.ModAmp = "+
				unitconv(Modetext->Text,Units::mTesla,ModUnitscombo->Text)+";\r\n" : "");

			//Update the final strings
			WriteHeader();
			WriteFooter();

		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for buttons for X-, Q-, W- Band of microwave freq  in Experiment tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::buttonFreqClick(System::Object^  sender, System::EventArgs^  e) 
{
	//Set the below values in Microwave freq. for X-, Q-, W- band respectively
	switch(Convert::ToInt32((String^)(((Button^)sender)->Tag))) {
	case 1:
		Frequencytext->Text="9.5";
		break;
	case 2:
		Frequencytext->Text="35";
		break;
	case 3:
		Frequencytext->Text="95";
		break;
	};

	//Update the value
	ExperimenttextBoxKeyPress(gcnew TextBox(),gcnew KeyPressEventArgs((char)13));
}

/**
* @brief EventHandler for combo boxes in Experiment tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Exp_comboBox_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e) 
{
	//Clear the string
	experiment1string="";

	//For pepper, write the crystal symmetry  mode
	if(easyspinfunction==2) experiment1string+="Exp.Mode = '"+ Modecombo->SelectedItem+"' ;\r\n";

	//Write the harmonic and numbero of point ofr the detection setting panel
	experiment1string+="Exp.Harmonic = "+ Harmoniccombo->SelectedItem+";\r\n";
	experiment1string+="Exp.nPoints = "+ Pointscombo->SelectedItem+";\r\n";

	//Update the final strings
	WriteHeader();
	WriteFooter();
}

/**
* @brief EventHandler for Temperature numeric box in Experiment tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Expnumeric_ValueChanged(System::Object^  sender, System::EventArgs^  e) 
{
	//Write the new value of the temperature in the string if checked
	experiment2string=(TempcheckBox->Checked? "Exp.Temperature = "+ Decimal::ToDouble(Tempnumeric->Value).ToString(System::Globalization::CultureInfo::InvariantCulture)+";\r\n" : "");

	//Update the final strings
	WriteHeader();
	WriteFooter();
}

/**
* @brief EventHandler for check boxes in Experiment tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Exp_checkBox_CheckedChanged(System::Object^  sender, System::EventArgs^  e)
{
	//Enable Temperature text box
	Tempnumeric->Enabled=TempcheckBox->Checked;
   
	//Enable mod. amp and units
	Modetext->Enabled=ModcheckBox->Checked;
	ModUnitscombo->Enabled=ModcheckBox->Checked;

	//Update the values
	Expnumeric_ValueChanged(gcnew NumericUpDown(), gcnew System::EventArgs());
	ExperimenttextBoxKeyPress(gcnew TextBox(), gcnew KeyPressEventArgs((char)13));
	
	//Update the final strings
	WriteHeader();
	WriteFooter();
}

/**
* @brief EventHandler for crystal symmetry parameter in Experiment tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::CrystalSymmetryBoxKeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	//Clear the string
	symmetrystring="";
	Int32 test; 

	//If key pressed is not number, or . or backspace or enter do nothing
	if (!Char::IsLetterOrDigit(e->KeyChar) && e->KeyChar != '/' && e->KeyChar != '\b' && e->KeyChar !=(char)13)
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted and pepper is enabled
		if (e->KeyChar ==(char)13 && easyspinfunction==2)
		{
			//Write crystal symmetry parameter parameters in correct format
			symmetrystring="Exp.CrystalSymmetry = "+ 
				(Int32::TryParse(CrystalSymmetrytext->Text, test)?
				CrystalSymmetrytext->Text : "'"+CrystalSymmetrytext->Text+"'")+";\r\n";

			//Update the final strings
			WriteHeader();
			WriteFooter();
		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

#pragma endregion

#pragma region Option tab

/**
* @brief EventHandler for combo boxes in Options tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Optioncombo_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e) 
{	
	//Clear the string
	optionstring="";
	
	//Garlic function , print the Method comobox (1st panel)
	if (easyspinfunction==0)
		optionstring="Opt.Method = '"+ GarlicMethodcombo->SelectedItem+"' ;\r\n";

	//Chilli function , print the Output and Verbosity comobox (2nd panel)
	if (easyspinfunction==1) 
	{
		optionstring="Opt.Output = '"+ chiliOutputcombo->SelectedItem+"' ;\r\n";
		optionstring+="Opt.Verbosity = "+ chiliVerbcombo->SelectedIndex+" ;\r\n";
	};

	//Pepper function , print the Method, Output, Symmetry and Verbosity comobox (3rd panel)
	if (easyspinfunction==2) 
	{
		optionstring="Opt.Method = '"+ PepperMethodcombo->SelectedItem+"' ;\r\n";
		optionstring+="Opt.Output = '"+ PepperOutputcombo->SelectedItem+"' ;\r\n";
		optionstring+="Opt.Symmetry = '"+ Symmetrycombo->SelectedItem+"' ;\r\n";
		optionstring+="Opt.Verbosity = "+ PepperVerbcombo->SelectedIndex+" ;\r\n";
	};

	//Update the final strings
	WriteHeader();
	WriteFooter();
}

/**
* @brief EventHandler for Numerix boxes in Options tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Optnumeric_ValueChanged(System::Object^  sender, System::EventArgs^  e)
{
	//Clear the string
	option1string="";

	//Chilli function , print the MOMD, knots and Director tilt values (2nd panel)
	if (easyspinfunction==1)
	{
		option1string="Exp.MOMD = "+(MOMDcheckBox->Checked? "1": "0")+";\r\n"+
			(MOMDcheckBox->Checked? "Opt.nKnots = "+Decimal::ToDouble(ChillKnotsnumeric->Value).ToString(System::Globalization::CultureInfo::InvariantCulture)+";\r\n":"")+
			(MOMDcheckBox->Checked? "":"Exp.psi = "+psitextBox->Text+";\r\n");
	}
	
	//Pepper function , print the threshold, knots, interpolate, cutoff and intensity values (3nd panel)
	if (easyspinfunction==2) //pepper
	{
		option1string="Opt.nKnots = ["+ pepKnotsnumeric->Value+
			(IntercheckBox->Checked? " "+ Interpnumeric->Value+"]" : " 0]")+";\r\n"+
			"Opt.Threshold = "+EasySpinNumber(pow(10.0,Decimal::ToDouble(Thresnumeric->Value)))+";\r\n"+
			"Opt.IsoCutoff = "+EasySpinNumber(pow(10.0,Decimal::ToDouble(Isonumeric->Value)))+";\r\n"+
			"Opt.Intensity = "+(IntenscheckBox->Checked? "'on'": "'off'")+";\r\n";
	};

	//Update the final strings
	WriteHeader();
	WriteFooter();
}

/**
* @brief EventHandler for check boxes in Options tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Opt_checkBox_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{
	//Chili function, depending the MOMD check box, enable or disable different field (2nd panel)
	if (easyspinfunction==1) 
	{
		ChillKnotsnumeric->Enabled=MOMDcheckBox->Checked;
		Knotslabel->Enabled=MOMDcheckBox->Checked;
		psilabel->Enabled=!MOMDcheckBox->Checked;
		psitextBox->Enabled=!MOMDcheckBox->Checked;
	}

	//Pepper function, depending the interpolate check box , enable numeric field next to it (3nd panel)
	else if (easyspinfunction==2) 
	{
		Interpnumeric->Enabled=IntercheckBox->Checked;
	};

	//Update hte numerics values
	Optnumeric_ValueChanged(gcnew NumericUpDown(), gcnew System::EventArgs());

	//Update the final strings
	WriteHeader();
	WriteFooter();
}

/**
* @brief EventHandler for text boxes in Options tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Opt_textBox_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	if(e->KeyChar==(char)13 && !EasySpinValidNumber(sender)) return;
	//Clear the string
	option2string="";

	//If key pressed is not number, or . or backspace or enter do nothing
	if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '-' && e->KeyChar != '+' && e->KeyChar != 'e' && e->KeyChar != 'E' && e->KeyChar != '\b' && e->KeyChar !=(char)13)
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted and chilli function
		if (e->KeyChar ==(char)13 && easyspinfunction==1)
		{
			//Returned to original color of the textbox
			((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

			//If one the below text fields is not empty, update text field
			if(EvenLMaxtext->Text!="" || OddLMaxtext->Text!="" || EvenKMaxtext->Text!=""||
				OddMMaxtext->Text!="" )
			{
				//For the empty texboxes set value zero
				if(EvenLMaxtext->Text=="") EvenLMaxtext->Text="0";
				if(OddLMaxtext->Text=="") OddLMaxtext->Text="0";
				if(EvenKMaxtext->Text=="") EvenKMaxtext->Text="0";
				if(OddMMaxtext->Text=="") OddMMaxtext->Text="0";

				//Write LKLM  parameters in correct format
				option2string="Opt.LLKM = ["+EvenLMaxtext->Text+" "+
					OddLMaxtext->Text+" "+ EvenKMaxtext->Text+" "+OddMMaxtext->Text+"];\r\n";
			}

			//Write Director tilt parameter in correct format if MOMD checked
			option2string+=(MOMDcheckBox->Checked? "":"Exp.psi = "+psitextBox->Text+";\r\n");

			//Update the final strings
			WriteHeader();
			WriteFooter();

		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

#pragma endregion

#pragma region Broadeninigs tab

/**
* @brief EventHandler for convolutions parameters in Broadenings tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::convolutionBoxKeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	if(e->KeyChar==(char)13 && !EasySpinValidNumber(sender)) return;
	if(e->KeyChar==(char)13)
		for each(TextBox^ box in gcnew cli::array<TextBox^>{GaussianFWHMtext,Gaussianpptext,LorentzianFWHMtext,Lorentzianpptext})
			if(!EasySpinValidNumber(box)) return;
	//Clear the string
	broadeningstring="";

	//If key pressed is not number, or . or backspace or enter do nothing
	if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '-' && e->KeyChar != '+' && e->KeyChar != 'e' && e->KeyChar != 'E' && e->KeyChar != '\b' && e->KeyChar !=(char)13 )
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted
		if (e->KeyChar ==(char)13 )
		{
			//Returned to original color of the textbox
			((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

			//Choose between Gaussian and Lorentzian FWHM or pp
			switch(Convert::ToInt32((String^)(((TextBox^)sender)->Tag))) {
			case 1://Gaussian FWHM to pp

				//For not empty text box
				if(((TextBox^)sender)->Text!="")
				{
					//Convert from FWHM to pp
					Gaussianpptext->Text=To_String(EasySpinDouble(
						GaussianFWHMtext->Text)/sqrt(2.0*log(2.0)));
				}
				else ///if FWHM text empty set also pp empty
					Gaussianpptext->Text="";

				break;

			case 2://Gaussian pp to FWHM

				//For not empty text box
				if(((TextBox^)sender)->Text!="")
				{
					//Convert from pp to FWHM
					GaussianFWHMtext->Text=To_String(EasySpinDouble(
						Gaussianpptext->Text)*sqrt(2.0*log(2.0)));
				}
				else ///if pp text empty set also FWHM empty
					GaussianFWHMtext->Text="";

				break;

			case 3://Lorentzian FWHM to pp

				//For not empty text box
				if(((TextBox^)sender)->Text!="")
				{
					//Convert from FWHM to pp
					Lorentzianpptext->Text=To_String(EasySpinDouble(
						LorentzianFWHMtext->Text)/sqrt(3.0));
				}
				else ///if FWHM text empty set also pp empty
					Lorentzianpptext->Text="";

				break;

			case 4://Lorentzian pp to FWHM

				//For not empty text box
				if(((TextBox^)sender)->Text!="")
				{
					//Convert from pp to FWHM
					LorentzianFWHMtext->Text=To_String(EasySpinDouble(
						Lorentzianpptext->Text)*sqrt(3.0));
				}
				else ///if pp text empty set also FWHM empty
					LorentzianFWHMtext->Text!="";

				break;

			};

			//If Lorentzian is empty, and Gaussian ont write Gaussian
			if(Lorentzianpptext->Text=="" && Gaussianpptext->Text!="")
			{
				//Write Gaussian in the appropriate format
				broadeningstring="Sys.lwpp = " +unitconv(Gaussianpptext->Text,Units::Gauss,BroadUnitscombo->Text)+";\r\n";
			}//If Lorentzian is not empty, write both if Gaussian not empty
			else if(Lorentzianpptext->Text!="")
			{
				//Write Lorentzian and Gaussian in the appropriate format
				broadeningstring="Sys.lwpp = ["+(Gaussianpptext->Text==""?"0":
					unitconv(Gaussianpptext->Text,Units::Gauss,BroadUnitscombo->Text))+
					" "+ unitconv(Lorentzianpptext->Text,Units::Gauss,BroadUnitscombo->Text)+"];\r\n";
			}

			//Update the final strings
			WriteHeader();
			WriteFooter();

		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for anisotropic parameters in Broadenings tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Anisotropic_text_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	if(e->KeyChar==(char)13 && !EasySpinValidNumber(sender)) return;
	//Clear the string
	broadening1string="";

	//If key pressed is not number, or . or backspace or enter do nothing
	if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '-' && e->KeyChar != '+' && e->KeyChar != 'e' && e->KeyChar != 'E' && e->KeyChar != '\b' && e->KeyChar !=(char)13)
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted and pepper is choosen
		if (e->KeyChar ==(char)13 && easyspinfunction==2)
		{
			//Returned to original color of the textbox
			((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

			//if any of Hstrain parameters are not empty
			if(HStrainxtext->Text!="" || HStrainytext->Text!="" ||	HStrainztext->Text!="" )
			{
				//Set the other parameters which are empty to zero
				if(HStrainxtext->Text=="") HStrainxtext->Text="0";
				if(HStrainytext->Text=="") HStrainytext->Text="0";
				if(HStrainztext->Text=="") HStrainztext->Text="0";

				//Write the final string in the appropriate format
				broadening1string+="Sys.HStrain = ["+ HStrainxtext->Text+" "+
					HStrainytext->Text+" "+
					HStrainztext->Text+"];\r\n";
			}; 

			//if any of gstrain parameters are not empty
			if(gStrainxtext->Text!="" || gStrainytext->Text!="" || gStrainztext->Text!="" )
			{
				//Set the other parameters which are empty to zero
				if(gStrainxtext->Text=="") gStrainxtext->Text="0";
				if(gStrainytext->Text=="") gStrainytext->Text="0";
				if(gStrainztext->Text=="") gStrainztext->Text="0";

				//Write the final string in the appropriate format
				broadening1string+="Sys.gStrain = ["+ gStrainxtext->Text+" "+
					gStrainytext->Text+" "+
					gStrainztext->Text+"];\r\n";
			};

			//if any of Astrain parameters are not empty
			if(AStrainxtext->Text!="" || AStrainytext->Text!="" || AStrainztext->Text!="" )
			{
				//Set the other parameters which are empty to zero
				if(AStrainxtext->Text=="") AStrainxtext->Text="0";
				if(AStrainytext->Text=="") AStrainytext->Text="0";
				if(AStrainztext->Text=="") AStrainztext->Text="0";

				//Write the final string in the appropriate format
				broadening1string+="Sys.AStrain = ["+ AStrainxtext->Text+" "+
					AStrainytext->Text+" "+
					AStrainztext->Text+"];\r\n";
			};

			//if any of Dstrain parameters are not empty
			if(DStrainDtext->Text!="" || DStrainEtext->Text!="" )
			{
				//Set the other parameters which are empty to zero
				if(DStrainDtext->Text=="") DStrainDtext->Text="0";
				if(DStrainEtext->Text=="") DStrainEtext->Text="0";

				//Write the final string in the appropriate format if corr is empty else write evrything
				if(DStraincorrtext->Text=="") broadening1string+="Sys.DStrain = ["+ DStrainDtext->Text+" "+ DStrainEtext->Text+"];\r\n";
				else broadening1string+="Sys.DStrain = ["+ DStrainDtext->Text+" "+
					DStrainEtext->Text+" "+
					DStraincorrtext->Text+"];\r\n";
			};

			//Update the final strings
			WriteHeader();
			WriteFooter();

		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

#pragma endregion

#pragma region Dynamics tab

/**
* @brief EventHandler for rotational radio button of Dynamics tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::radioButton_diffusion_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{
	if(!((RadioButton^)sender)->Checked) return;
	//Choose between isotorpic or axial or rhombic diffusion to enable respective panels
	switch(Convert::ToInt32((String^)(((RadioButton^)sender)->Tag))) {
	case 1:  //isotropic
		isotropicpanel->Visible=true;
		axialpanel->Visible=false;
		rhombicpanel->Visible=false;
		isotropicpanel->Location=System::Drawing::Point(3, 83);
		break;
	case 2:  //axial
		isotropicpanel->Visible=false;
		axialpanel->Location=System::Drawing::Point(3, 83);
		axialpanel->Visible=true;
		rhombicpanel->Visible=false;
		break;
	case 3:  //rhombic
		isotropicpanel->Visible=false;
		axialpanel->Visible=false;
		rhombicpanel->Visible=true;
		rhombicpanel->Location=System::Drawing::Point(3, 83);
		break;

	};

	//Set the postion of the paneles
	DynamicsgroupBox->Size=System::Drawing::Size(396,220);
	RotationallogcheckBox->Location=System::Drawing::Point(12, 185);
	OrderinggroupBox->Location=System::Drawing::Point(3, 225);
	Heisenbergpanel->Location=System::Drawing::Point(6, 335);
	if(SystemModel!=nullptr) textBoxDynamics_KeyPress(isocorrtext, gcnew KeyPressEventArgs((char)13));
}

/**
* @brief EventHandler for diffusion parameters in Dynamics tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::textBoxDynamics_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	if(e->KeyChar==(char)13 && !EasySpinValidNumber(sender)) return;
	if(e->KeyChar==(char)13)
	{
		cli::array<TextBox^>^ inputs=Convert::ToInt32(((TextBox^)sender)->Tag)==1 ?
			gcnew cli::array<TextBox^>{isocorrtext,axcorrxytext,axcorrztext,corrxtext,corrytext,corrztext} :
			gcnew cli::array<TextBox^>{isodifftext,axdiffxytext,axdiffztext,diffxtext,diffytext,diffztext};
		for each(TextBox^ box in inputs) if(!EasySpinValidNumber(box)) return;
	}
	dynamicstring="";	//Clear the string
	String^ temp;
	Unit^ unit;

	//If key pressed is not number, or . or backspace or enter do nothing
	if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '-' && e->KeyChar != '+' && e->KeyChar != 'e' && e->KeyChar != 'E' && e->KeyChar != '\b' && e->KeyChar !=(char)13)
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted
		if (e->KeyChar ==(char)13)
		{
			//Returned to original color of the textbox
			((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

			//for choosing between Corelation time and Diffusion tensor
			switch(Convert::ToInt32((String^)(((TextBox^)sender)->Tag))) {
			case 1: //When writing correlation time, we change also diffusion auto

				//Isotropic diffusion and when the text is not empty
				if(isocorrtext->Text!="" )
				{
					//Set the units for diffusion tensor
					if(DiffIsoUnitscombo->Text=="GHz") unit=Units::GHz;
					else  unit=Units::MHz;

					//Save the text for time after conversion to chosen SI units 
					temp=unitconv(RotationallogcheckBox->Checked? logtext(isocorrtext->Text) : isocorrtext->Text,Units::seconds,IsoUnitscombo->Text);

					//Calulate the diffusion tensor from the correlation 
					//time and convert to selected units
					isodifftext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::Hz->get_name());
				} //if time text empty set also tensor empty
				else
					isodifftext->Text="";

				//axial anisotropic diffusion xy and when the text is not empty
				if(axcorrxytext->Text!="")
				{
					//Set the units for diffusion tensor
					if(DiffAxUnitscombo->Text=="GHz") unit=Units::GHz;
					else  unit=Units::MHz;

					//Save the text for time after conversion to chosen SI units 
					temp=unitconv(RotationallogcheckBox->Checked? logtext(axcorrxytext->Text) : axcorrxytext->Text,Units::seconds,AxUnitscombo->Text);

					//Calulate the diffusion tensor from the correlation 
					//time and convert to selected units
					axdiffxytext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::Hz->get_name());
				} //if time text empty set also tensor empty
				else 
					axdiffxytext->Text="";

				//axial anisotropic diffusion z and when the text is not empty
				if(axcorrztext->Text!="")
				{
					//Set the units for diffusion tensor
					if(DiffAxUnitscombo->Text=="GHz") unit=Units::GHz;
					else  unit=Units::MHz;

					//Save the text for time after conversion to chosen SI units 
					temp=unitconv(RotationallogcheckBox->Checked? logtext(axcorrztext->Text) : axcorrztext->Text,Units::seconds,AxUnitscombo->Text);

					//Calulate the diffusion tensor from the correlation 
					//time and convert to selected units
					axdiffztext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::Hz->get_name());
				} //if time text empty set also tensor empty
				else
					axdiffztext->Text="";

				//rhombic anisotropic diffusion x and when the text is not empty
				if(corrxtext->Text!="")
				{
					//Set the units for diffusion tensor
					if(DiffRhUnitscombo->Text=="GHz") unit=Units::GHz;
					else  unit=Units::MHz;

					//Save the text for time after conversion to chosen SI units 
					temp=unitconv(RotationallogcheckBox->Checked? logtext(corrxtext->Text) : corrxtext->Text,Units::seconds,RhUnitscombo->Text);

					//Calulate the diffusion tensor from the correlation 
					//time and convert to selected units
					diffxtext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::Hz->get_name());
				} //if time text empty set also tensor empty
				else
					diffxtext->Text="";

				//rhombic anisotropic diffusion y and when the text is not empty
				if(corrytext->Text!="")
				{
					//Set the units for diffusion tensor
					if(DiffRhUnitscombo->Text=="GHz") unit=Units::GHz;
					else  unit=Units::MHz;

					//Save the text for time after conversion to chosen SI units 
					temp=unitconv(RotationallogcheckBox->Checked? logtext(corrytext->Text) : corrytext->Text,Units::seconds,RhUnitscombo->Text);

					//Calulate the diffusion tensor from the correlation 
					//time and convert to selected units
					diffytext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::Hz->get_name());
				} //if time text empty set also tensor empty
				else
					diffytext->Text="";

				//rhombic anisotropic diffusion z and when the text is not empty
				if(corrztext->Text!="")
				{
					//Set the units for diffusion tensor
					if(DiffRhUnitscombo->Text=="GHz") unit=Units::GHz;
					else  unit=Units::MHz;

					//Save the text for time after conversion to chosen SI units 
					temp=unitconv(RotationallogcheckBox->Checked? logtext(corrztext->Text) : corrztext->Text,Units::seconds,RhUnitscombo->Text);

					//Calulate the diffusion tensor from the correlation 
					//time and convert to selected units
					diffztext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::Hz->get_name());
				}
				else//if time text empty set also tensor empty
					diffztext->Text="";

				break;	
				//////////////////////////////////////////////////////////

			case 2: //When writing diffusion tensor, we change also time auto

				///Isotropic time and when the text is not empty
				if(isodifftext->Text!="")
				{
					//Set the units for correlation time
					if(IsoUnitscombo->Text=="ns") unit=Units::nsec;
					else  unit=Units::mi_sec;

					//Save the text for tensor after conversion to chosen SI units 
					temp=unitconv(RotationallogcheckBox->Checked? logtext(isodifftext->Text) : isodifftext->Text,Units::Hz,DiffIsoUnitscombo->Text);

					//Calulate the correlation time from the diffusion 
					//tensor and convert to selected units
					isocorrtext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::seconds->get_name());
				} //if tensor text empty set also time empty
				else
					isocorrtext->Text="";

				//Axial anisotropic time xy and when the text is not empty
				if(axdiffxytext->Text!="")
				{
					//Set the units for correlation time
					if(AxUnitscombo->Text=="ns") unit=Units::nsec;
					else  unit=Units::mi_sec;

					//Save the text for tensor after conversion to chosen SI units
					temp=unitconv(RotationallogcheckBox->Checked? logtext(axdiffxytext->Text) : axdiffxytext->Text,Units::Hz,DiffAxUnitscombo->Text);

					//Calulate the correlation time from the diffusion 
					//tensor and convert to selected units
					axcorrxytext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::seconds->get_name());
				}//if tensor text empty set also time empty
				else 
					axcorrxytext->Text="";

				//Axial anisotropic time z and when the text is not empty
				if(axdiffztext->Text!="")
				{
					//Set the units for correlation time
					if(AxUnitscombo->Text=="ns") unit=Units::nsec;
					else  unit=Units::mi_sec;

					//Save the text for tensor after conversion to chosen SI units
					temp=unitconv(RotationallogcheckBox->Checked? logtext(axdiffztext->Text) : axdiffztext->Text,Units::Hz,DiffAxUnitscombo->Text);

					//Calulate the correlation time from the diffusion 
					//tensor and convert to selected units
					axcorrztext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::seconds->get_name());
				}//if tensor text empty set also time empty
				else
					axcorrztext->Text="";

				//rhombic anisotropic time x and when the text is not empty
				if(diffxtext->Text!="")
				{
					//Set the units for correlation time
					if(RhUnitscombo->Text=="ns") unit=Units::nsec;
					else  unit=Units::mi_sec;

					//Save the text for tensor after conversion to chosen SI units
					temp=unitconv(RotationallogcheckBox->Checked? logtext(diffxtext->Text) : diffxtext->Text,Units::Hz,DiffRhUnitscombo->Text);

					//Calulate the correlation time from the diffusion 
					//tensor and convert to selected units
					corrxtext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::seconds->get_name());
				}//if tensor text empty set also time empty
				else
					corrxtext->Text="";

				//rhombic anisotropic time y and when the text is not empty
				if(diffytext->Text!="")
				{
					//Set the units for correlation time
					if(RhUnitscombo->Text=="ns") unit=Units::nsec;
					else  unit=Units::mi_sec;

					//Save the text for tensor after conversion to chosen SI units
					temp=unitconv(RotationallogcheckBox->Checked? logtext(diffytext->Text) : diffytext->Text,Units::Hz,DiffRhUnitscombo->Text);

					//Calulate the correlation time from the diffusion 
					//tensor and convert to selected units
					corrytext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::seconds->get_name());
				}//if tensor text empty set also time empty
				else
					corrytext->Text="";

				//rhombic anisotropic time z and when the text is not empty
				if(diffztext->Text!="")
				{
					//Set the units for correlation time
					if(RhUnitscombo->Text=="ns") unit=Units::nsec;
					else  unit=Units::mi_sec;

					//Save the text for tensor after conversion to chosen SI units
					temp=unitconv(RotationallogcheckBox->Checked? logtext(diffztext->Text) : diffztext->Text,Units::Hz,DiffRhUnitscombo->Text);

					//Calulate the correlation time from the diffusion 
					//tensor and convert to selected units
					corrztext->Text=unitconv(EasySpinNumber(1.0/(6.0*EasySpinDouble(temp))),
						unit,Units::seconds->get_name());
				}//if tensor text empty set also time empty
				else
					corrztext->Text="";

				break;
			};

			//The checkbox specifies logarithms for both correlation and diffusion fields.
			if(RotationallogcheckBox->Checked)
			{
				cli::array<TextBox^>^ converted=Convert::ToInt32(((TextBox^)sender)->Tag)==1 ?
					gcnew cli::array<TextBox^>{isodifftext,axdiffxytext,axdiffztext,diffxtext,diffytext,diffztext} :
					gcnew cli::array<TextBox^>{isocorrtext,axcorrxytext,axcorrztext,corrxtext,corrytext,corrztext};
				for each(TextBox^ box in converted)
					if(box->Text!="" && EasySpinDouble(box->Text)>0.0) box->Text=EasySpinNumber(log10(EasySpinDouble(box->Text)));
			}

			//When the isotropic diffusion is selected, write parameters to the 
			//appropriate format is SI unit. If log is enable write them in log form.
			if(Isotropicradio->Checked && isocorrtext->Text!="") dynamicstring="Sys.tcorr = " +
				(RotationallogcheckBox->Checked?
				unitconv(logtext(isocorrtext->Text),Units::seconds,IsoUnitscombo->Text):
				unitconv(isocorrtext->Text,Units::seconds,IsoUnitscombo->Text))+";\r\n";

			//When the axial anisotropic diffusion is selected, write parameters to the 
			//appropriate format is SI unit. If log is enable write them in log form.
			if(Axialradio->Checked && axcorrxytext->Text!="" && axcorrztext->Text!="") dynamicstring="Sys.tcorr = [" +
				(RotationallogcheckBox->Checked?
				unitconv(logtext(axcorrxytext->Text),Units::seconds,AxUnitscombo->Text):
				unitconv(axcorrxytext->Text,Units::seconds,AxUnitscombo->Text))+" "+
				(RotationallogcheckBox->Checked?
				unitconv(logtext(axcorrztext->Text),Units::seconds,AxUnitscombo->Text):
				unitconv(axcorrztext->Text,Units::seconds,AxUnitscombo->Text))+"];\r\n";

			//When the rhombic aniisotropic diffusion is selected, write parameters to the 
			//appropriate format is SI unit. If log is enable write them in log form.
			if(Rhombicradio->Checked && corrxtext->Text!="" && corrytext->Text!="" && corrztext->Text!="") dynamicstring="Sys.tcorr =[ " +
				(RotationallogcheckBox->Checked?
				unitconv(logtext(corrxtext->Text),Units::seconds,RhUnitscombo->Text):
				unitconv(corrxtext->Text,Units::seconds,RhUnitscombo->Text))+" "+
				(RotationallogcheckBox->Checked?
				unitconv(logtext(corrytext->Text),Units::seconds,RhUnitscombo->Text):
				unitconv(corrytext->Text,Units::seconds,RhUnitscombo->Text))+" "+
				(RotationallogcheckBox->Checked?
				unitconv(logtext(corrztext->Text),Units::seconds,RhUnitscombo->Text):
				unitconv(corrztext->Text,Units::seconds,RhUnitscombo->Text))+" "+"];\r\n";

			//Update the final strings
			WriteHeader();
			WriteFooter();

		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for ordering potential parameters in Dynamics tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::orderheisentext_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	if(e->KeyChar==(char)13 && !EasySpinValidNumber(sender)) return;
	//Clear the string
	dynamic1string="";

	//If key pressed is not number, or . or backspace or enter do nothing
	if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '-' && e->KeyChar != '+' && e->KeyChar != 'e' && e->KeyChar != 'E' && e->KeyChar != '\b' && e->KeyChar !=(char)13)
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted and chili is enabled
		if (e->KeyChar ==(char)13 && easyspinfunction==1)
		{
			//Returned to original color of the textbox
			((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

			//If one the below text fields is not empty, update text field
			if(l20text->Text!="" || l22text->Text!="" || l40text->Text!=""||
				l42text->Text!="" || l44text->Text!="" )
			{
				//For the empty texboxes set value zero
				if(l20text->Text=="") l20text->Text="0";
				if(l22text->Text=="") l22text->Text="0";
				if(l40text->Text=="") l40text->Text="0";
				if(l42text->Text=="") l42text->Text="0";
				if(l44text->Text=="") l44text->Text="0";

				//Write lambda  parameters in correct format
				dynamic1string="Sys.labda = ["+ l20text->Text+" "+l22text->Text+" "+l40text->Text+
					" "+l42text->Text+" "+l44text->Text+" "+"];\r\n";
			} 

			//If heisenberg textbox is not empty write/update in the final string
			if(heisentext->Text!="") dynamic1string+="Sys.Exchange = "+heisentext->Text+";\r\n";

			//Update the final strings
			WriteHeader();
			WriteFooter();

		} //When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}
  
#pragma endregion

#pragma region Units

/**
* @brief EventHandler for when a unit is selected using units combobox
* in Broadenings tab and Experiment tab.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::Unitscombo_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e)
{
	switch(Convert::ToInt32((String^)(((ComboBox^)sender)->Tag))) 
	{
	case 1: //For updating convolutions textboxes in Broadenings tab, hitting enter (char)13
		convolutionBoxKeyPress(gcnew TextBox(),gcnew KeyPressEventArgs((char)13));
		break;
	case 2: //For updating top textboxes in Experiment tab, hitting enter (char)13
		ExperimenttextBoxKeyPress(gcnew TextBox(),gcnew KeyPressEventArgs((char)13));
		break;
	};
}

/**
* @brief EventHandler for when a unit is selected using units combobox in Dynamics Tab
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
 System::Void SpinachGUI::EasySpin::DynamicsUnitscombo_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e) 
 {
	 //Choose tag 
	 TextBox^ textbox=gcnew TextBox();
	 textbox->Tag=(String^)(((ComboBox^)sender)->Tag);

	 //Update textboxes in Dynamics tab, hitting enter (char)13
	 textBoxDynamics_KeyPress(textbox,gcnew KeyPressEventArgs((char)13));
 }

#pragma endregion

#pragma region Preview Box

 /**
* @brief EventHandler for when we move the mouse with mouse button down, hand effect.
* Not perfect from times to times crash.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
 System::Void SpinachGUI::EasySpin::PreviewtextBox_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) 
 {
	 //Only for left button
	 if( e->Button == System::Windows::Forms::MouseButtons::Left)
	 {
		 //for the X coordinate
		 int tempX=panel2->HorizontalScroll->Value+2*(InitialPoint->X-e->Location.X);

		 //The value of position should be inside the range of horizontal scroll
		 if(tempX>panel2->HorizontalScroll->Minimum && tempX<panel2->HorizontalScroll->Maximum)
		 {
			 //Set the value of the Horizontal scroll
			 panel2->HorizontalScroll->Value+=2*(InitialPoint->X-e->Location.X);

			 //Reset the initial point
			 InitialPoint->X=e->Location.X;
		 }

		 //for the Y coordinate
		 int tempY=panel2->VerticalScroll->Value+2*(InitialPoint->Y-e->Location.Y);

		 //The value of position should be inside the range of Vertical scroll
		 if(tempY>panel2->VerticalScroll->Minimum && tempY<panel2->VerticalScroll->Maximum)
		 {
			 //Set the value of the Vertical scroll
			 panel2->VerticalScroll->Value+=2*(InitialPoint->Y-e->Location.Y);

			 //Reset the initial point
			 InitialPoint->Y=e->Location.Y;
		 }
	 }
 }

 /**
* @brief EventHandler for when we press mouse button down to save initial position.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::EasySpin::PreviewtextBox_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) 
{
	//Save initial position and then fire mouseMove
	InitialPoint=e->Location;
	PreviewtextBox_MouseMove(sender,e);
}

#pragma endregion