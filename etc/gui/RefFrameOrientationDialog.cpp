#include "stdafx.h"
#include "RefFrameOrientationDialog.h"

//Definitions
#define To_String(aNumber) SpinachGUI::ConvertToString(aNumber, 6, "N")

/**
*\file RefFrameOrientationDialog.cpp
*/

/**
* @brief Function for setting up orientation object from rotation matrix.
*/
void SpinachGUI::RefFrameOrientationDialog::CreateOrient()
{
	Matrix3x3^ AA=gcnew Matrix3x3();

	//Calculate the determinant of DCM
	double detA= rotation[0]*rotation[4]*rotation[8]
	-rotation[0]*rotation[5]*rotation[7]
	-rotation[3]*rotation[1]*rotation[8]
	+rotation[3]*rotation[2]*rotation[7]
	+rotation[6]*rotation[1]*rotation[5]
	-rotation[6]*rotation[4]*rotation[2];

	detA=1.0f;

	//Translate DCM multiplied by its determinant to a Matrix3x3 matrix
	AA->xx=detA*rotation[0];AA->xy=detA*rotation[1];AA->xz=detA*rotation[2];
	AA->yx=detA*rotation[3];AA->yy=detA*rotation[4];AA->yz=detA*rotation[5];
	AA->zx=detA*rotation[6];AA->zy=detA*rotation[7];AA->zz=detA*rotation[8];

	//Setup the initial Orientation object
	nOrient= gcnew SpinachGUI::Orientation(AA);
}

/**
* @brief Function for updating Euler angle textboxes after calculating them from orientation.
*/
void SpinachGUI::RefFrameOrientationDialog::UpdateEulerAngles()
{
	//Get euler angles from orientation object
	EulerAngles^ euleran=nOrient->GetAsEuler();

	//Convert agles from grad to degrees
	alpha=float(180.0f*euleran->alpha/PI);
	beta=float(180.0f*euleran->beta/PI);
	gamma=float(180.0f*euleran->gamma/PI);

	//set the euler textboxes
	textBoxalpha->Text= To_String(alpha);
	textBoxbeta->Text= To_String(beta);
	textBoxgamma->Text= To_String(gamma);
}

/**
* brief Function for updating quaternion textboxes after calculating them from orientation.
*/
void SpinachGUI::RefFrameOrientationDialog::UpdateQuaternion()
{
	//Get quaternion parameters from orientation object
	Quaternion^ quater=nOrient->GetAsQuaternion();

	//set the quaternion textboxes
	textBoxquat1->Text= To_String(quater->w);
	textBoxquat2->Text= To_String(quater->x);
	textBoxquat3->Text= To_String(quater->y);
	textBoxquat4->Text= To_String(quater->z);
}

/**
* brief Function for updating angle axis textboxes after calculating them from orientation.
*/
void SpinachGUI::RefFrameOrientationDialog::UpdateAngleAxis()
{
	//Get angle axis parameters from orientation object
	AngleAxis^ angleax=nOrient->GetAsAngleAxis();

	//set the angle axis textboxes
	textBoxangle1->Text= To_String(angleax->x);
	textBoxangle2->Text= To_String(angleax->y);
	textBoxangle3->Text= To_String(angleax->z);
	textBoxangle4->Text= To_String(180.0f*angleax->angle/PI);
}

/**
* @brief Function for updating rotation matrix textboxes.
*/
void SpinachGUI::RefFrameOrientationDialog::UpdateRotation()
{
	//set the rotation matrix textboxes
	XXtextBox->Text= To_String(rotation[0]);
	XYtextBox->Text= To_String(rotation[1]);
	XZtextBox->Text= To_String(rotation[2]);
	YXtextBox->Text= To_String(rotation[3]);
	YYtextBox->Text= To_String(rotation[4]);
	YZtextBox->Text= To_String(rotation[5]);
	ZXtextBox->Text= To_String(rotation[6]);
	ZYtextBox->Text= To_String(rotation[7]);
	ZZtextBox->Text= To_String(rotation[8]);
}

/**
* @brief EventHandler for text inside angle axis , angle textbox.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::RefFrameOrientationDialog::textBoxangle_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	//If key pressed is not number, or .  or - or backspace or enter do nothing
	if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '-' && e->KeyChar != '\b' && e->KeyChar !=(char)13)
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted 
		if (e->KeyChar ==(char)13)
		{
			double result;

			//if valid number is entered
			if(Double::TryParse(((TextBox^)sender)->Text, System::Globalization::NumberStyles::Float | System::Globalization::NumberStyles::AllowThousands, System::Globalization::CultureInfo::InvariantCulture, result))
			{
				//Returned to original color of the textbox
				((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

				//Collect angle axis parameters from textboxes
				nangleaxis=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture)*PI/180.0f;
				Vector3d naxis= Vector3d(Convert::ToDouble(textBoxangle1->Text, System::Globalization::CultureInfo::InvariantCulture),
					Convert::ToDouble(textBoxangle2->Text, System::Globalization::CultureInfo::InvariantCulture),
					Convert::ToDouble(textBoxangle3->Text, System::Globalization::CultureInfo::InvariantCulture));

				//Set up orientation object from angle axis parameters
				nOrient= gcnew SpinachGUI::Orientation(gcnew AngleAxis(naxis.x(),
					naxis.y(),
					naxis.z(),
					nangleaxis));

				//Convert to rotation matrix
				Matrix3x3^ dc=nOrient->GetAsDCM();
				rotation[0]=dc->xx;rotation[1]=dc->xy;rotation[2]=dc->xz;
				rotation[3]=dc->yx;rotation[4]=dc->yy;rotation[5]=dc->yz;
				rotation[6]=dc->zx;rotation[7]=dc->zy;rotation[8]=dc->zz;

				//Update rotation matrix, quaternion, euler-angles texttboxes
				CreateOrient();
				UpdateEulerAngles();
				UpdateQuaternion();
				UpdateRotation();

				//Update 3D visualization
				OpenGL->UpdateDisplay(UpdateType::ELLIPSOIDS);
			}//if not valid number do nothing
			else e->Handled = true;
		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for text inside euler angles textboxes.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::RefFrameOrientationDialog::textBoxEuler_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	//If key pressed is not number, or .  or backspace or enter do nothing
	if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '-' && e->KeyChar != '\b' && e->KeyChar !=(char)13)
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted 
		if (e->KeyChar ==(char)13)
		{
			double textfield;

			//if valid number is entered
			if(Double::TryParse(((TextBox^)sender)->Text, System::Globalization::NumberStyles::Float | System::Globalization::NumberStyles::AllowThousands, System::Globalization::CultureInfo::InvariantCulture, textfield))
			{
				//Returned to original color of the textbox
				((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

				//Choose between alpha , beta , gamma
				switch(Convert::ToInt32((String^)(((TextBox^)sender)->Tag))) {
				case 1:
					alpha=textfield;
					break;
				case 2:
					beta=textfield;
					break;
				case 3:
					gamma=textfield;
					break;

				};

				//Set up orientation object from euler angles
				nOrient= gcnew SpinachGUI::Orientation(gcnew EulerAngles(alpha*PI/180.0f,beta*PI/180.0f, gamma*PI/180.0f));

				//Convert to rotation matrix
				Matrix3x3^ dc=nOrient->GetAsDCM();
				rotation[0]=dc->xx;rotation[1]=dc->xy;rotation[2]=dc->xz;
				rotation[3]=dc->yx;rotation[4]=dc->yy;rotation[5]=dc->yz;
				rotation[6]=dc->zx;rotation[7]=dc->zy;rotation[8]=dc->zz;

				//Update rotation matrix, quaternion, angle axis texttboxes
				UpdateRotation();
				UpdateQuaternion();
				UpdateAngleAxis();

				//Update 3D visualization
				OpenGL->UpdateDisplay(UpdateType::ELLIPSOIDS);
			}//if not valid number do nothing
			else e->Handled = true;
		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}