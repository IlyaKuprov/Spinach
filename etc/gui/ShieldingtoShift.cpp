#include "stdafx.h"
#include "ShieldingtoShift.h"

//Definitions
#define MainICollection  SystemModel->InteractionCollection
#define MainCollection  SystemModel->AtomCollection

/**
*\file ShieldingtoShift.cpp
*/


/**
* @brief Function that initializes the default list of reference substances.
* Also makes channel list and initialize the datagridview.
*/
void SpinachGUI::ShieldingtoShift::Initialization()
{
	//strings for Channel list and other things
	String^ temp=gcnew String("");
	List<String^>^ ChannelList=gcnew List<String^>();
	int j=0;

	//Containers for default substances
	reference_shielding=gcnew Dictionary<String^, double>();
	reference_substance=gcnew Dictionary<String^, String^>();


	//Default values for reference shifts substances
	//1H
	reference_substance->Add("1H","TMS");
	reference_shielding->Add("1H",33.440);
	//2H
	reference_substance->Add("2H","TMS");
	reference_shielding->Add("2H",33.440);
	//13C
	reference_substance->Add("13C","TMS");
	reference_shielding->Add("13C",186.38);
	//14N
	reference_substance->Add("14N","CH3NO2");
	reference_shielding->Add("14N",-132.83);
	//15N
	reference_substance->Add("15N","CH3NO2");
	reference_shielding->Add("15N",-134.42);
	//17O
	reference_substance->Add("17O","D2O");
	reference_shielding->Add("17O",290.66);
	//19F
	reference_substance->Add("19F","CFCl3");
	reference_shielding->Add("19F",192.97);
	//29Si
	reference_substance->Add("29Si","TMS");
	reference_shielding->Add("29Si",378.94);
	//31P
	reference_substance->Add("31P","H3PO4");
	reference_shielding->Add("31P",332.07);
	//33S
	reference_substance->Add("33S","(NH4)2SO4");
	reference_shielding->Add("33S",205.17);


	//Loop over all elements to find shieldings and create channel list
	for each(int i in SystemModel->InteractionCollection->Keys)
	{
		//Only for chemical shiefts
		if(((Interaction^)SystemModel->InteractionCollection[i])->IntKind==InteractionKind::CShielding)
		{
			//Create a temporary channel
			temp=((Interaction^)SystemModel->InteractionCollection[i])->A->isotope->Mass+
				((Interaction^)SystemModel->InteractionCollection[i])->A->isotope->Element;

			//If this channel is not already inside the list
			if(!ChannelList->Contains(temp)) 
			{
				//Add a new line and name to the datagridview and to the channel list
				ChanneldataGridView->Rows->Add();
				ChanneldataGridView[0,j]->Value=temp;
				ChannelList->Add(temp);


				//If reference is given
				if(((Interaction^)SystemModel->InteractionCollection[i])->reference!=nullptr)
				{

					try //Try to find channel-element in the default list of substances
					{
						//If channel has the same reference with the default list
						if(reference_substance[temp]==
							((Interaction^)SystemModel->InteractionCollection[i])->reference)
						{
							//set the value and the name
							ChanneldataGridView[2,j]->Value=reference_substance[temp];
							ChanneldataGridView[1,j]->Value=reference_shielding[temp];
						}
						else //If channel has different reference substance from default list
						{
							//change the reference in the default list
							reference_substance[temp]=((Interaction^)SystemModel->InteractionCollection[i])->reference;

							//refresh value in the datagridview
							ChanneldataGridView[2,j]->Value=reference_substance[temp];

							//set the value zero in order to insert it
							reference_shielding[temp]=NULL;
						};
					}
					catch (KeyNotFoundException^) //if channel-element does not exist
					{
						//Set it in the datagrid view
						ChanneldataGridView[2,j]->Value=((Interaction^)SystemModel->InteractionCollection[i])->reference;

						//Add it  hte new channel in the dafault list
						reference_substance->Add(temp,
							((Interaction^)SystemModel->InteractionCollection[i])->reference);
					}
				} //if reference substance not given
				else 
				{
					//Try to find if channel-element is existing in defaults
					try 
					{
						//Update datagridview with substance and value
						ChanneldataGridView[2,j]->Value=reference_substance[temp];
						ChanneldataGridView[1,j]->Value=reference_shielding[temp];

					}//channel-element not in list
					catch (KeyNotFoundException^) 
					{
						//Add the channel in the default list with no name and keep column 1 and 2 empty
						reference_substance->Add(temp,"");
					};
				};

				j++;
			};
		};
	};

	//Check if all the data are filled, enable OK button
	CheckShiftCompletness();
}

/**
* @brief EventHandler for when finishing editing the first or second column.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::ShieldingtoShift::ChanneldataGridView_CellEndEdit(System::Object^  sender, System::Windows::Forms::DataGridViewCellEventArgs^  e) 
{
	//Change reference substance
	if (e->ColumnIndex==2) 
	{
		//Try to add channel-element in the default list of substances
		try 
		{
			reference_substance->Add(((String^)ChanneldataGridView[0,e->RowIndex]->Value),
				((String^)ChanneldataGridView[2,e->RowIndex]->Value));
		}//If already exist, just change the value
		catch (ArgumentException^)
		{
			reference_substance[((String^)ChanneldataGridView[0,e->RowIndex]->Value)]=
				((String^)ChanneldataGridView[2,e->RowIndex]->Value);
		};

		//Check if all the data are filled, enable OK button
		CheckShiftCompletness();

	}//Change reference shielding
	else if(e->ColumnIndex==1)  
	{
		double result;

		//Check if the value is a valid number
		if(Double::TryParse((String^)ChanneldataGridView[e->ColumnIndex,e->RowIndex]->Value, result))
		{
			//In order to convert it to double
			ChanneldataGridView[e->ColumnIndex,e->RowIndex]->Value=result;

			//Try to add channel-element in the default list of shielding value
			try 
			{
				reference_shielding->Add(((String^)ChanneldataGridView[0,e->RowIndex]->Value),result);
			}//If already exist, just change the value
			catch (ArgumentException^)
			{
				reference_shielding[((String^)ChanneldataGridView[0,e->RowIndex]->Value)]=result;

			};

			//Check if all the data are filled, enable OK button
			CheckShiftCompletness();

		} //if not a valid number
		else
		{
			//throw error message
			MessageBox::Show("Wrong Number Format for Coordinates", "Error", 
				MessageBoxButtons::OK, MessageBoxIcon::Error) ;

			//After hiiting OK in the error message,  edit is reinvoked for this cell
			SetColumnIndex^ method= gcnew SetColumnIndex(this, &ShieldingtoShift::DelegateMethod);
			ChanneldataGridView->BeginInvoke(method,e->ColumnIndex,e->RowIndex);
		};
	};
}

/**
* @brief Function that checks if all the cells in the first and second
* column are filled in order to enabled OK button
*/
System::Void SpinachGUI::ShieldingtoShift::CheckShiftCompletness()
{
	bool test=true;

	//Check each line and first and second column are filled with the appropriate data
	for(int i=0;i<ChanneldataGridView->Rows->Count;i++)
	{
		test=test && (Double^)ChanneldataGridView[1,i]->Value!=nullptr && 
			         (String^)ChanneldataGridView[2,i]->Value!=nullptr;
	};
    
	OKbutton->Enabled=test;
}

/**
* @brief EventHandler for setting up the keypress handler to a cell in the first or second column.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::ShieldingtoShift::ChanneldataGridView_EditingControlShowing(System::Object^  sender, System::Windows::Forms::DataGridViewEditingControlShowingEventArgs^  e) 
{
	//For only first and second column
	if(ChanneldataGridView->CurrentCellAddress.X==1 || ChanneldataGridView->CurrentCellAddress.X==2)
	{
		//Setup the text for selected cell and set the keypress event handler 
		TextBox^ txtbox = (TextBox^)(e->Control);
		if (txtbox != nullptr) txtbox->KeyPress += gcnew KeyPressEventHandler(this,&ShieldingtoShift::ChanneldataGridView_KeyPress);
	};
}

/**
* @brief EventHandler for text in cells.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::ShieldingtoShift::ChanneldataGridView_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	//For second column(reference shielding) only, check if it is number
	if (ChanneldataGridView->CurrentCellAddress.X==1)
	{
		//If key pressed is not number, or . or - or backspace or enter do nothing
		if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '\b'  && e->KeyChar != '-' && e->KeyChar !=(char)13)
		{
			e->Handled = true;
		};
	};

	//When enter is hitted call the endedit handler
	if(e->KeyChar==(char)13) ChanneldataGridView->EndEdit(DataGridViewDataErrorContexts::Commit);
}

/**
* @brief EventHandler for when OK button click with finishing the converstion from shielding to shifting tensors
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::ShieldingtoShift::OKbutton_Click(System::Object^  sender, System::EventArgs^  e)
{
	String ^element=gcnew String("");

	//Loop over all interactions
	for each(int i in MainICollection->Keys)
	{
		//Find all chemical shieldings
		if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::CShielding)
		{
			//Change of the the shielding to shift kind
			((Interaction^)MainICollection[i])->IntKind=InteractionKind::Shift;

			//Create channel string for easyness
			element=((Interaction^)SystemModel->InteractionCollection[i])->A->isotope->Mass+
				((Interaction^)SystemModel->InteractionCollection[i])->A->isotope->Element;

			//Set the substance for the interaction in case  it doesnot have
			((Interaction^)MainICollection[i])->reference=reference_substance[element];

			//Calculate shifting (ref_Shield_parameter*I)-Matrix
			((Tensor^)MainICollection[i])->matrix3x3[0]=reference_shielding[element]-
				((Tensor^)MainICollection[i])->matrix3x3[0];
			((Tensor^)MainICollection[i])->matrix3x3[4]=reference_shielding[element]-
				((Tensor^)MainICollection[i])->matrix3x3[4];
			((Tensor^)MainICollection[i])->matrix3x3[8]=reference_shielding[element]-
				((Tensor^)MainICollection[i])->matrix3x3[8];

			((Tensor^)MainICollection[i])->matrix3x3[1]=-((Tensor^)MainICollection[i])->matrix3x3[1];
			((Tensor^)MainICollection[i])->matrix3x3[2]=-((Tensor^)MainICollection[i])->matrix3x3[2];
			((Tensor^)MainICollection[i])->matrix3x3[3]=-((Tensor^)MainICollection[i])->matrix3x3[3];
			((Tensor^)MainICollection[i])->matrix3x3[5]=-((Tensor^)MainICollection[i])->matrix3x3[5];
			((Tensor^)MainICollection[i])->matrix3x3[6]=-((Tensor^)MainICollection[i])->matrix3x3[6];
			((Tensor^)MainICollection[i])->matrix3x3[7]=-((Tensor^)MainICollection[i])->matrix3x3[7];
		};
	};
}

/**
* @brief Function which used with the delegate in order to redit a cell with wrong format
* @param columnIndex Index of the column of the edited cell.
* @param rowIndex Index of the row of the edited cell.
*/
System::Void SpinachGUI::ShieldingtoShift::DelegateMethod(int columnIndex, int rowIndex)
{
	//Set the current cell as the one clicked
	ChanneldataGridView->CurrentCell = ChanneldataGridView[columnIndex,rowIndex];

	//Begin edit of the cell
	ChanneldataGridView->BeginEdit(true);

}