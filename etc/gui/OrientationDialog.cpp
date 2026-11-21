#include "StdAfx.h"
#include "OrientationDialog.h"

//Definitions
#define To_String(aNumber) SpinachGUI::ConvertToString(aNumber, 6, "F")
#define To_StringE(aNumber) SpinachGUI::ConvertToString(aNumber, 6, "e3")

/**
*\file OrientationDialog.cpp
*/

/**
* @brief Initialization of orientation dialog.
* @param aTensor The tensor that will change.
* @param aOpenGL The OpenGL panel which will be redrawn in every update.
* @param aCheck Special fuction that control the update
*/
void SpinachGUI::OrientationDialog::Initialization(Tensor^% aTensor, OpenGLForm::COpenGL^% aOpenGL, Check^ aCheck)
{
	//set the parameters of the interaction, OpenGl and update function
	tensor=aTensor;
	OpenGL=aOpenGL;
	nAtom=tensor->A;
	ID=tensor->ID;
	unit=tensor->unit;
	refID=tensor->Frame->getID();
	InterKind=tensor->IntKind;
	Update=aCheck;

	//Construct spherical  matrices of different forms
	Re_Spherical_Tens=new double[9];
	Im_Spherical_Tens=new double[9];
	Interact_Mat=gcnew ::Matrix3x3();
	eigenvec=gcnew ObservableCollection<double>();
	eigenval=gcnew ObservableCollection<double>();

	//Set theXY,XZ,YY,YZ,ZZ elemtns textboxes disable for j-couplings and exchange interactions
	if(InterKind==InteractionKind::Jcoupling || InterKind==InteractionKind::Exchange)
	{
		textBoxXYTensor->Enabled=false;
		textBoxXZTensor->Enabled=false;
		textBoxYYTensor->Enabled=false;
		textBoxYZTensor->Enabled=false;
		textBoxZZTensor->Enabled=false;
	} 
	else
	{
		textBoxXYTensor->Enabled=true;
		textBoxXZTensor->Enabled=true;
		textBoxYYTensor->Enabled=true;
		textBoxYZTensor->Enabled=true;
		textBoxZZTensor->Enabled=true;
	} 

	//Initialize tensor, eigen values and dcm
	//In order to not use the automatic update
	for (int i=0;i<9;i++)
	{
		Interact_Mat[i]=tensor->matrix3x3[i];
		eigenvec->Add(tensor->dcm[i]);
	}
	eigenval->Add(tensor->eigen[0]);
	eigenval->Add(tensor->eigen[1]);
	eigenval->Add(tensor->eigen[2]);

	//if interaction is chemical shift appear at specific position and set the textbox
	if(InterKind==InteractionKind::Shift)
	{
		groupBox1->Width=groupBox1->Width+270;
		Referencelabel->Visible=true;
		ReferencetextBox->Visible=true;
		ReferencetextBox->Text=tensor->reference;
	}

	//Calculate and update all textboxes
	UpdateAll();
}

#pragma region KeyPress Handlers

/**
* @brief EventHandler for text inside interaction tensor matrix.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationDialog::InteractionMatrix_KeyPress(System::Object^  sender, KeyPressEventArgs^  e) 
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

				//Depending the tag of the text box, update the interaction matrix
				switch(Convert::ToInt32((String^)(((TextBox^)sender)->Tag))) {
				case 1: //for XX

					//Convert to double
					Interact_Mat[0]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);

					//If j-coupling or exchange interaction the YY and ZZ textboes 
					//will be locked so update them automaticaly
					if(InterKind==InteractionKind::Jcoupling || InterKind==InteractionKind::Exchange)
					{
						Interact_Mat[4]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
						Interact_Mat[8]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					}
					break;
				case 2: //For XY, convert to double and update automatically YX
					Interact_Mat[1]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					Interact_Mat[3]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					break;
				case 3://For XZ, convert to double and update automatically ZX
					Interact_Mat[2]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					Interact_Mat[6]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					break;
				case 4: //For YY, convert to double
					Interact_Mat[4]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					break;
				case 5: //For YZ, convert to double and update automatically ZY
					Interact_Mat[5]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					Interact_Mat[7]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					break;
				case 6: //For ZZ, convert to double
					Interact_Mat[8]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					break;
				};

				//Update the text boxes, if jcoupling or exchange interaction only YY and ZZ
				if(InterKind==InteractionKind::Jcoupling || InterKind==InteractionKind::Exchange)
				{
					textBoxYYTensor->Text=Convert::ToString(Interact_Mat[0]);
					textBoxZZTensor->Text=Convert::ToString(Interact_Mat[0]);
				} //Update the text boxes for XY,ZX,ZY
				else
				{
					textBoxYXTensor->Text=Convert::ToString(Interact_Mat[3]);
					textBoxZXTensor->Text=Convert::ToString(Interact_Mat[6]);
					textBoxZYTensor->Text=Convert::ToString(Interact_Mat[7]);
				}

				//Calculate all the forms and update the text boxes, the model and 3d model
				Calculate_Eigen();
				UpdateInteraction();
				UpdateEigenState();
				UpdateIsotropic();
				UpdateSphericaTensor();
				CreateOrient();
				UpdateEulerAngles();
				UpdateWignerMatrix();
				UpdateQuaternion();
				UpdateAngleAxis();
				UpdateParameters();
				UpdateModel();
			}//if not valid number do nothing
			else e->Handled = true;
		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for text inside Euler angles boxes.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationDialog::textBoxeuler_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
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
			double textfield;

			//if valid number is entered
			if(Double::TryParse(((TextBox^)sender)->Text, System::Globalization::NumberStyles::Float | System::Globalization::NumberStyles::AllowThousands, System::Globalization::CultureInfo::InvariantCulture, textfield))
			{
				//Returned to original color of the textbox 
				((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

				//Depending the tag of the text box, update the euler angles
				switch(Convert::ToInt32((String^)(((TextBox^)sender)->Tag))) {
				case 1://alpha
					alpha=textfield;
					break;
				case 2://beta
					beta=textfield;
					break;
				case 3://gamma
					gamma=textfield;
					break;
				};

				//In order to work correcly, according to Spinach, convert to gradiants
				nOrient= gcnew SpinachGUI::Orientation(gcnew EulerAngles(gamma*PI/180.0f,beta*PI/180.0f,alpha*PI/180.0f));

				//Convert to dcm
				Matrix3x3^ dc=nOrient->GetAsDCM();

				//Update dcm matrix
				eigenvec[0]=dc->xx;
				eigenvec[1]=dc->xy;
				eigenvec[2]=dc->xz;
				eigenvec[3]=dc->yx;
				eigenvec[4]=dc->yy;
				eigenvec[5]=dc->yz;
				eigenvec[6]=dc->zx;
				eigenvec[7]=dc->zy;
				eigenvec[8]=dc->zz;

				//Calculate all the forms and update the text boxes, the model and 3d model
				Calculate_TensorfromEigen();
				UpdateInteraction();
				UpdateEigenState();
				UpdateIsotropic();
				UpdateSphericaTensor();
				UpdateWignerMatrix();
				UpdateQuaternion();
				UpdateAngleAxis();
				UpdateModel();

			}//if not valid number do nothing
			else e->Handled = true;
		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for text inside Spherical tensor boxes.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationDialog::textBoxSpherical_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	//save tag of texbox
	int texttag=Convert::ToInt32((String^)(((TextBox^)sender)->Tag));

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
			double textfield;

			//if valid number is entered
			if(Double::TryParse(((TextBox^)sender)->Text, System::Globalization::NumberStyles::Float | System::Globalization::NumberStyles::AllowThousands, System::Globalization::CultureInfo::InvariantCulture, textfield))
			{
				//Returned to original color of the textbox
				((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Window;

				//if the tag is even number the real part will saved else the imaginary
				if(texttag%2==0) Re_Spherical_Tens[texttag/2]= textfield;
				else Im_Spherical_Tens[texttag/2]= textfield;

				//Update Symmetrical terms
				textBoxSpherical2n1R->Text= To_String(Re_Spherical_Tens[7]=-Re_Spherical_Tens[5]);
				textBoxSpherical2n1I->Text= To_String(Im_Spherical_Tens[7]=Im_Spherical_Tens[5]);
				textBoxSpherical2n2R->Text= To_String(Re_Spherical_Tens[8]=Re_Spherical_Tens[4]);
				textBoxSpherical2n2I->Text= To_String(Im_Spherical_Tens[8]=-Im_Spherical_Tens[4]);

				//if all of this expression are different from 0 then the tensor matrix will imaginary, not good
				if( (3.0f* Im_Spherical_Tens[4] - sqrt(6.0f)* Im_Spherical_Tens[6] + 3.0f* Im_Spherical_Tens[8])/6.0f!=0 &&
					(sqrt(2.0f)*Re_Spherical_Tens[2] + Re_Spherical_Tens[4] - Re_Spherical_Tens[8])/2.0f!=0 &&
					(Im_Spherical_Tens[1]+Im_Spherical_Tens[3] -Im_Spherical_Tens[5]+ Im_Spherical_Tens[7])/2.0f!=0 &&
					(-sqrt(2.0f)* Re_Spherical_Tens[2] + Re_Spherical_Tens[4] - Re_Spherical_Tens[8])/2.0f!=0 &&
					(- 3.0f* Im_Spherical_Tens[4] -sqrt(6.0f)*Im_Spherical_Tens[6]- 3.0f* Im_Spherical_Tens[8])/6.0f!=0 &&
					(Re_Spherical_Tens[1]-Re_Spherical_Tens[3] -Re_Spherical_Tens[5]- Re_Spherical_Tens[7])/2.0f!=0 &&
					(-Im_Spherical_Tens[1] - Im_Spherical_Tens[3] -Im_Spherical_Tens[5] + Im_Spherical_Tens[7])/2.0f!=0 &&
					(-Re_Spherical_Tens[1] + Re_Spherical_Tens[3] -Re_Spherical_Tens[5] + Re_Spherical_Tens[7])/2.0f!=0 &&
					-sqrt(2.0f/3.0f)*Im_Spherical_Tens[6]!=0) 
				{
					System::Windows::Forms::MessageBox::Show("Complex Interaction Matrix Element", "Error", 
						System::Windows::Forms::MessageBoxButtons::OK,
						System::Windows::Forms::MessageBoxIcon::Error);}
				else //if the above conditio is false, we can have real tensor
				{
					//Calculation of XX from the spherical tensor
					Interact_Mat[0]=(-2.0f*sqrt(3.0f)*Re_Spherical_Tens[0] + 3.0f*Re_Spherical_Tens[4] 
					-sqrt(6.0f)*Re_Spherical_Tens[6] + 3.0f*Re_Spherical_Tens[8])/6.0f;

					//Calculation of XY from the spherical tensor
					Interact_Mat[1]=(-sqrt(2.0f)* Im_Spherical_Tens[2]- Im_Spherical_Tens[4]+ 
						Im_Spherical_Tens[8] )/2.0f;

					//Calculation of XZ from the spherical tensor
					Interact_Mat[2]=(Re_Spherical_Tens[1] + Re_Spherical_Tens[3] 
					- Re_Spherical_Tens[5] + Re_Spherical_Tens[7])/2.0f;

					//Calculation of YX from the spherical tensor
					Interact_Mat[3]=(sqrt(2.0f)* Im_Spherical_Tens[2] - Im_Spherical_Tens[4] + 
						Im_Spherical_Tens[8])/2.0f;

					//Calculation of YY from the spherical tensor
					Interact_Mat[4]= (-2.0f*sqrt(3.0f)*Re_Spherical_Tens[0]- 3.0f* Re_Spherical_Tens[4] - 
						sqrt(6.0f)*Re_Spherical_Tens[6]- 3.0f* Re_Spherical_Tens[8])/6.0f;

					//Calculation of YZ from the spherical tensor
					Interact_Mat[5]=(-Im_Spherical_Tens[1] + Im_Spherical_Tens[3]
					+Im_Spherical_Tens[5] + Im_Spherical_Tens[7])/2.0f;

					//Calculation of ZX from the spherical tensor
					Interact_Mat[6]=(-Re_Spherical_Tens[1] - Re_Spherical_Tens[3] 
					-Re_Spherical_Tens[5] + Re_Spherical_Tens[7])/2.0f;

					//Calculation of ZY from the spherical tensor
					Interact_Mat[7]=(Im_Spherical_Tens[1] - Im_Spherical_Tens[3]
					+Im_Spherical_Tens[5] + Im_Spherical_Tens[7])/2.0f;

					//Calculation of ZZ from the spherical tensor
					Interact_Mat[8]=(-Re_Spherical_Tens[0]+sqrt(2.0f)*Re_Spherical_Tens[6])/sqrt(3.0);

					//Calculate all the forms and update the text boxes, the model and 3d model
					Calculate_Eigen();
					UpdateInteraction();
					UpdateEigenState();
					UpdateIsotropic();
					CreateOrient();
					UpdateEulerAngles();
					UpdateWignerMatrix();
					UpdateQuaternion();
					UpdateAngleAxis();
					UpdateParameters();
					UpdateModel();
					UpdateSphericaTensor();
				};
			}//if not valid number do nothing
			else e->Handled = true;
		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for text inside Eigen values boxes.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationDialog::textBoxEigen_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
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

				//Depending the tag of the text box, update the eigenvalues
				switch(Convert::ToInt32((String^)(((TextBox^)sender)->Tag))) {
				case 1://1st
					eigenval[0]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					break;
				case 2://2nd
					eigenval[1]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					break;
				case 3://3rd
					eigenval[2]=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture);
					break;
				};

				//Calculate all the forms and update the text boxes, the model and 3d model
				Calculate_TensorfromEigen();
				UpdateInteraction();
				UpdateEigenState();
				UpdateIsotropic();
				UpdateSphericaTensor();
				CreateOrient();
				UpdateWignerMatrix();
				UpdateEulerAngles();
				UpdateQuaternion();
				UpdateAngleAxis();
				UpdateParameters();
				UpdateModel();
			}//if not valid number do nothing
			else e->Handled = true;
		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for text inside angle axis, angle box.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationDialog::textBoxAngle_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
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
				nOrient= gcnew SpinachGUI::Orientation(gcnew AngleAxis(naxis.x(), naxis.y(), naxis.z(), nangleaxis));

				//Convert to dcm matrix
				Matrix3x3^ dc=nOrient->GetAsDCM();
				eigenvec[0]=dc->xx;eigenvec[1]=dc->xy;eigenvec[2]=dc->xz;
				eigenvec[3]=dc->yx;eigenvec[4]=dc->yy;eigenvec[5]=dc->yz;
				eigenvec[6]=dc->zx;eigenvec[7]=dc->zy;eigenvec[8]=dc->zz;

				//Calculate all the forms and update the text boxes, the model and 3d model
				CreateOrient();
				UpdateEulerAngles();
				UpdateQuaternion();
				UpdateWignerMatrix();
				UpdateEigenState();
				UpdateIsotropic();
				Calculate_TensorfromEigen();
				UpdateInteraction();
				UpdateSphericaTensor();
				UpdateModel();
			}//if not valid number do nothing
			else e->Handled = true;
		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for text insideisotropic conventions, isotropic box.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationDialog::textBoxIso_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
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

				//Calculate trace and then subtrack from the diagonal elements
				double trace=(Interact_Mat[0]+Interact_Mat[4]+Interact_Mat[8])/3.0f;
				Interact_Mat[0]+=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture)-trace;
				Interact_Mat[4]+=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture)-trace;
				Interact_Mat[8]+=Convert::ToDouble(((TextBox^)sender)->Text, System::Globalization::CultureInfo::InvariantCulture)-trace;

				//Calculate all the forms and update the text boxes, the model and 3d model
				Calculate_Eigen();
				UpdateInteraction();
				UpdateEigenState();
				UpdateIsotropic();
				UpdateSphericaTensor();
				CreateOrient();
				UpdateEulerAngles();
				UpdateWignerMatrix();
				UpdateQuaternion();
				UpdateAngleAxis();
				UpdateParameters();
				UpdateModel();
			}//if not valid number do nothing
			else e->Handled = true;
		}//When normal writing is happenning the Textbox color change to blue
		else ((TextBox^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}
#pragma endregion

#pragma region Updates of Textboxes

/**
* @brief Function for updating interaction tensor matrix textboxes.
*/
void SpinachGUI::OrientationDialog::UpdateInteraction()
{
	InteractionText=gcnew String("interactionmatrix=[ ");

	//set the interaction tensor matrix textboxes
	textBoxXXTensor->Text= To_String(Interact_Mat[0]);
	textBoxXYTensor->Text= To_String(Interact_Mat[1]);
	textBoxXZTensor->Text= To_String(Interact_Mat[2]);
	textBoxYXTensor->Text= To_String(Interact_Mat[3]);
	textBoxYYTensor->Text= To_String(Interact_Mat[4]);
	textBoxYZTensor->Text= To_String(Interact_Mat[5]);
	textBoxZXTensor->Text= To_String(Interact_Mat[6]);
	textBoxZYTensor->Text= To_String(Interact_Mat[7]);
	textBoxZZTensor->Text= To_String(Interact_Mat[8]);


	//Writing to matlab format
	InteractionText+=To_StringE(Interact_Mat[0])+" "+To_StringE(Interact_Mat[1])+" "+To_StringE(Interact_Mat[2])+"\r\n"+
		"                    "+To_StringE(Interact_Mat[3])+" "+To_StringE(Interact_Mat[4])+" "+To_StringE(Interact_Mat[5])+"\r\n"+
		"                    "+To_StringE(Interact_Mat[6])+" "+To_StringE(Interact_Mat[7])+" "+To_StringE(Interact_Mat[8])+" ];\r\n";	
}

/**
* @brief Function for updating eigenvalues and eigenvector matrix textboxes.
*/
void SpinachGUI::OrientationDialog::UpdateEigenState()
{
	EigenValuesText=gcnew String("eigenvalues=[ ");
	EigenVectorsText=gcnew String("eigenvectors=[ ");

	//EigenValues
	textBoxEigenVal1->Text= To_String(eigenval[0]);
	textBoxEigenVal2->Text= To_String(eigenval[1]);
	textBoxEigenVal3->Text= To_String(eigenval[2]);

	//Writing to matlab format
	EigenValuesText+=To_StringE(eigenval[0])+" "+To_StringE(eigenval[1])+" "+To_StringE(eigenval[2])+" ];\r\n";

	//EigenVector
	textBoxXXEigenVec->Text= To_String(eigenvec[0]);
	textBoxXYEigenVec->Text= To_String(eigenvec[1]);
	textBoxXZEigenVec->Text= To_String(eigenvec[2]);
	textBoxYXEigenVec->Text= To_String(eigenvec[3]);
	textBoxYYEigenVec->Text= To_String(eigenvec[4]);
	textBoxYZEigenVec->Text= To_String(eigenvec[5]);
	textBoxZXEigenVec->Text= To_String(eigenvec[6]);
	textBoxZYEigenVec->Text= To_String(eigenvec[7]);
	textBoxZZEigenVec->Text= To_String(eigenvec[8]);

	//Writing to matlab format
	EigenVectorsText+=To_StringE(eigenvec[0])+" "+To_StringE(eigenvec[1])+" "+To_StringE(eigenvec[2])+" \r\n"+
		"              "+To_StringE(eigenvec[3])+" "+To_StringE(eigenvec[4])+" "+To_StringE(eigenvec[5])+" \r\n"+
		"              "+To_StringE(eigenvec[6])+" "+To_StringE(eigenvec[7])+" "+To_StringE(eigenvec[8])+" ];\r\n";
}

/**
* @brief Function for calculating and updateing isotropic parameters textboxes. 
* Call after calculating eigenvalues
*/
void SpinachGUI::OrientationDialog::UpdateIsotropic()
{
	double temp=0.0;
	ConventionsText=gcnew String("");

	//Isotropic: calculation, updating textbox, Write to matlab format
	temp=(eigenval[0]+eigenval[1]+eigenval[2])/3.0;
	textBoxIso->Text= To_String(temp);
	ConventionsText+="isotropic= "+To_StringE(temp)+" ; \r\n";

	//Axiality: calculation, updating textbox, Write to matlab format
	temp=2.0*eigenval[2]-(eigenval[0]+eigenval[1]);
	textBoxAx->Text= To_String(temp);
	ConventionsText+="axiality= "+To_StringE(temp)+" ; \r\n";

	//Rhombicity: calculation, updating textbox, Write to matlab format
	temp=eigenval[0]-eigenval[1];
	textBoxRh->Text= To_String(temp);
	ConventionsText+="rhombicity= "+To_StringE(temp)+" ; \r\n";

	//Span: calculation, updating textbox, Write to matlab format
	temp=eigenval[2]-eigenval[0];
	textBoxSp->Text= To_String(temp);
	ConventionsText+="span= "+To_StringE(temp)+" ; \r\n";

	//Skew: calculation, updating textbox, Write to matlab format
	temp=(eigenval[0]-2.0*eigenval[1]+eigenval[2])/(2*temp);
	textBoxSk->Text= To_String(temp);
	ConventionsText+="skew= "+To_StringE(temp)+" ; \r\n";
}

/**
* @brief Function for calculating and updateing of spherical tensor textboxes. 
* Call after calculating tensor
*/
void SpinachGUI::OrientationDialog::UpdateSphericaTensor()
{
	//Strings for the matlab export
	String^ sign, ^rank0, ^rank1,^rank2;
	SphericalTensorText=gcnew String("sphericaltensor");
	rank0=gcnew String(".rank0= ");
	rank1=gcnew String(".rank1=[ ");
	rank2=gcnew String(".rank2=[ ");

	//T2,2 : calculation, updating textbox
	Re_Spherical_Tens[4]=(Interact_Mat[0]-Interact_Mat[4])/2.0;
	Im_Spherical_Tens[4]=(-Interact_Mat[1]-Interact_Mat[3])/2.0;
	textBoxSpherical22R->Text= To_String(Re_Spherical_Tens[4]);
	textBoxSpherical22I->Text= To_String(Im_Spherical_Tens[4]);

	//Write to matlab format, checking if needed the + sign first
	sign=(Im_Spherical_Tens[4]<0?"":"+");
	rank2+=To_StringE(Re_Spherical_Tens[4])+sign+To_StringE(Im_Spherical_Tens[4])+"i \r\n";

	//T2,1 : calculation, updating textbox
	Re_Spherical_Tens[5]=-(Interact_Mat[2]+Interact_Mat[6])/2.0;
	Im_Spherical_Tens[5]=-(-Interact_Mat[5]-Interact_Mat[7])/2.0;
	textBoxSpherical21R->Text= To_String(Re_Spherical_Tens[5]);
	textBoxSpherical21I->Text= To_String(Im_Spherical_Tens[5]);

	//Write to matlab format, checking if needed the + sign first
	sign=(Im_Spherical_Tens[5]<0?"":"+");
	rank2+="                        "+To_StringE(Re_Spherical_Tens[5])+sign+To_StringE(Im_Spherical_Tens[5])+"i \r\n";

	//T1,1 : calculation, updating textbox
	Re_Spherical_Tens[1]=-(Interact_Mat[6]-Interact_Mat[2])/2.0;
	Im_Spherical_Tens[1]=-(-Interact_Mat[7]+Interact_Mat[5])/2.0;
	textBoxSpherical11R->Text= To_String(Re_Spherical_Tens[1]);
	textBoxSpherical11I->Text= To_String(Im_Spherical_Tens[1]);

	//Write to matlab format, checking if needed the + sign first
	sign=(Im_Spherical_Tens[1]<0?"":"+");
	rank1+=To_StringE(Re_Spherical_Tens[1])+sign+To_StringE(Im_Spherical_Tens[1])+"i \r\n";

	//T2,0 : calculation, updating textbox
	Re_Spherical_Tens[6]=(2.0*Interact_Mat[8]-Interact_Mat[0]-Interact_Mat[4])/sqrt(6.0);
	Im_Spherical_Tens[6]=0.0;
	textBoxSpherical20R->Text= To_String(Re_Spherical_Tens[6]);
	textBoxSpherical20I->Text= To_String(Im_Spherical_Tens[6]);

	//Write to matlab format, checking if needed the + sign first
	sign=(Im_Spherical_Tens[6]<0?"":"+");
	rank2+="                        "+To_StringE(Re_Spherical_Tens[6])+sign+To_StringE(Im_Spherical_Tens[6])+"i \r\n";

	//T1,0 : calculation, updating textbox
	Re_Spherical_Tens[2]=0.0;
	Im_Spherical_Tens[2]=-(Interact_Mat[1]-Interact_Mat[3])/sqrt(2.0);
	textBoxSpherical10R->Text= To_String(Re_Spherical_Tens[2]);
	textBoxSpherical10I->Text= To_String(Im_Spherical_Tens[2]);

	//Write to matlab format, checking if needed the + sign first
	sign=(Im_Spherical_Tens[2]<0?"":"+");
	rank1+="                        "+To_StringE(Re_Spherical_Tens[2])+sign+To_StringE(Im_Spherical_Tens[2])+"i  \r\n";

	//T0,0 : calculation, updating textbox
	Re_Spherical_Tens[0]=-(Interact_Mat[0]+Interact_Mat[4]+Interact_Mat[8])/sqrt(3.0);
	Im_Spherical_Tens[0]=0.0;
	textBoxSpherical00R->Text= To_String(Re_Spherical_Tens[0]);
	textBoxSpherical00I->Text= To_String(Im_Spherical_Tens[0]);

	//Write to matlab format, checking if needed the + sign first
	sign=(Im_Spherical_Tens[0]<0?"":"+");
	rank0+=To_StringE(Re_Spherical_Tens[0])+sign+To_StringE(Im_Spherical_Tens[0])+"i ;\r\n";


	//T2,-1 : calculation, updating textbox
	Re_Spherical_Tens[7]=(Interact_Mat[2]+Interact_Mat[6])/2.0;
	Im_Spherical_Tens[7]=(Interact_Mat[5]+Interact_Mat[7])/2.0;
	textBoxSpherical2n1R->Text= To_String(Re_Spherical_Tens[7]);
	textBoxSpherical2n1I->Text= To_String(Im_Spherical_Tens[7]);

	//Write to matlab format, checking if needed the + sign first
	sign=(Im_Spherical_Tens[7]<0?"":"+");
	rank2+="                        "+To_StringE(Re_Spherical_Tens[7])+sign+To_StringE(Im_Spherical_Tens[7])+"i  \r\n";

	//T1,-1 : calculation, updating textbox
	Re_Spherical_Tens[3]=-(Interact_Mat[6]-Interact_Mat[2])/2.0;
	Im_Spherical_Tens[3]=-(Interact_Mat[7]-Interact_Mat[5])/2.0;
	textBoxSpherical1n1R->Text= To_String(Re_Spherical_Tens[3]);
	textBoxSpherical1n1I->Text= To_String(Im_Spherical_Tens[3]);

	//Write to matlab format, checking if needed the + sign first
	sign=(Im_Spherical_Tens[3]<0?"":"+");
	rank1+="                        "+To_StringE(Re_Spherical_Tens[3])+sign+To_StringE(Im_Spherical_Tens[3])+"i ];\r\n";

	//T2,-2 : calculation, updating textbox
	Re_Spherical_Tens[8]=(Interact_Mat[0]-Interact_Mat[4])/2.0;
	Im_Spherical_Tens[8]=(Interact_Mat[1]+Interact_Mat[3])/2.0;
	textBoxSpherical2n2R->Text= To_String(Re_Spherical_Tens[8]);
	textBoxSpherical2n2I->Text= To_String(Im_Spherical_Tens[8]);

	//Write to matlab format, checking if needed the + sign first
	sign=(Im_Spherical_Tens[8]<0?"":"+");
	rank2+="                        "+To_StringE(Re_Spherical_Tens[8])+sign+To_StringE(Im_Spherical_Tens[8])+"i ]; \r\n";

	//Write the final matlab format
	SphericalTensorText=SphericalTensorText+rank2+
		SphericalTensorText+rank1+
		SphericalTensorText+rank0;
}

/**
* @brief Function for updating Euler angle textboxes after calculating them from orientation.
*/
void SpinachGUI::OrientationDialog::UpdateEulerAngles()
{
	EulerAnglesText=gcnew String("eulerangles=[ ");

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

	//Writing to matlab format
	EulerAnglesText+=To_StringE(alpha)+" "+To_StringE(beta)+" "+ To_StringE(gamma)+" ]; \r\n";
}

/**
* @brief Function for calculating and updateing Wigner Matrix textboxes. 
* Call after calculating eigenvalues
*/
void SpinachGUI::OrientationDialog::UpdateWignerMatrix()
{
	//Strings for the matlab export
	WignerMatrixText=gcnew String("wignermatrix=[ ");
	double Realtemp=0.0;
	double Imaginarytemp=0.0;
	double alpha, beta, gamma;
	String ^sign;

	//Take euler angles and convert them to gradiants.
	EulerAngles^ euleran=gcnew EulerAngles(this->alpha*PI/180.0f,this->beta*PI/180.0f,this->gamma*PI/180.0f);
	
	//Retain the full Euler angles; rotations are not quadrant symmetric.
	alpha=euleran->alpha;
	beta=euleran->beta;
	gamma=euleran->gamma;

//M[1,1] : calculation, updating textboxes
	Realtemp=pow(cos(beta/2.0),4)*cos(2.0*(alpha + gamma));
	Imaginarytemp=-(pow(cos(beta/2.),4)*sin(2*(alpha + gamma)));
	textBox60->Text= To_String(Realtemp);
	textBox59->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[1,2] : calculation, updating textboxes
	Realtemp=(cos(alpha + 2*gamma)*(2*sin(beta) + sin(2*beta)))/4.;
	Imaginarytemp=-((2*sin(beta) + sin(2*beta))*sin(alpha + 2*gamma))/4.;
	textBox65->Text= To_String(Realtemp);
	textBox64->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[1,3] : calculation, updating textboxes
	Realtemp=(sqrt(1.5)*cos(2*gamma)*pow(sin(beta),2))/2.;
	Imaginarytemp=-(sqrt(1.5)*cos(gamma)*pow(sin(beta),2)*sin(gamma));
	textBox112->Text= To_String(Realtemp);
	textBox107->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[1,4] : calculation, updating textboxes
	Realtemp=cos(alpha - 2*gamma)*pow(sin(beta/2.),2)*sin(beta);
	Imaginarytemp=pow(sin(beta/2.),2)*sin(beta)*sin(alpha - 2*gamma);
	textBox111->Text= To_String(Realtemp);
	textBox110->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[1,5] : calculation, updating textboxes
	Realtemp=cos(2*alpha - 2*gamma)*pow(sin(beta/2.),4);
	Imaginarytemp=pow(sin(beta/2.),4)*sin(2*alpha - 2*gamma);
	textBox109->Text= To_String(Realtemp);
	textBox108->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i \r\n";

	////////////////////////////////////////////////////
	//M[2,1] : calculation, updating textboxes
	Realtemp=-(cos(2*alpha + gamma)*(2*sin(beta) + sin(2*beta)))/4.;
	Imaginarytemp=((2*sin(beta) + sin(2*beta))*sin(2*alpha + gamma))/4.;
	textBox47->Text= To_String(Realtemp);
	textBox46->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+="               "+To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[2,2] : calculation, updating textboxes
	Realtemp=((cos(beta) + cos(2*beta))*cos(alpha + gamma))/2.;
	Imaginarytemp=-((cos(beta) + cos(2*beta))*sin(alpha + gamma))/2.;
	textBox58->Text= To_String(Realtemp);
	textBox57->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[2,3] : calculation, updating textboxes
	Realtemp=sqrt(1.5)*cos(beta)*cos(gamma)*sin(beta);
	Imaginarytemp=-(sqrt(1.5)*cos(beta)*sin(beta)*sin(gamma));
	textBox102->Text= To_String(Realtemp);
	textBox97->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[2,4] : calculation, updating textboxes
	Realtemp=((cos(beta) - cos(2*beta))*cos(alpha - gamma))/2.;
	Imaginarytemp=((cos(beta) - cos(2*beta))*sin(alpha - gamma))/2.;
	textBox106->Text= To_String(Realtemp);
	textBox105->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[2,5] : calculation, updating textboxes
	Realtemp=cos(2*alpha - gamma)*pow(sin(beta/2.),2)*sin(beta);
	Imaginarytemp=pow(sin(beta/2.),2)*sin(beta)*sin(2*alpha - gamma);
	textBox104->Text= To_String(Realtemp);
	textBox103->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i \r\n";

	////////////////////////////////////////////////////
	//M[3,1] : calculation, updating textboxes
	Realtemp=(sqrt(1.5)*cos(2*alpha)*pow(sin(beta),2))/2.;
	Imaginarytemp=-(sqrt(1.5)*cos(alpha)*sin(alpha)*pow(sin(beta),2));
	textBox62->Text= To_String(Realtemp);
	textBox61->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+="               "+To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[3,2] : calculation, updating textboxes
	Realtemp=-(sqrt(1.5)*cos(alpha)*cos(beta)*sin(beta));
	Imaginarytemp=sqrt(1.5)*cos(beta)*sin(alpha)*sin(beta);
	textBox50->Text= To_String(Realtemp);
	textBox49->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[3,3] : calculation, updating textboxes
	Realtemp=(1 + 3*cos(2*beta))/4.;
	Imaginarytemp=0.0;
	textBox56->Text= To_String(Realtemp);
	textBox53->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[3,4] : calculation, updating textboxes
	Realtemp=sqrt(1.5)*cos(alpha)*cos(beta)*sin(beta);
	Imaginarytemp=sqrt(1.5)*cos(beta)*sin(alpha)*sin(beta);
	textBox96->Text= To_String(Realtemp);
	textBox95->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[3,5] : calculation, updating textboxes
	Realtemp=(sqrt(1.5)*cos(2*alpha)*pow(sin(beta),2))/2.;
	Imaginarytemp=sqrt(1.5)*cos(alpha)*sin(alpha)*pow(sin(beta),2);
	textBox94->Text= To_String(Realtemp);
	textBox93->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i \r\n";

	////////////////////////////////////////////////////
	//M[4,1] : calculation, updating textboxes
	Realtemp=-(cos(2*alpha - gamma)*pow(sin(beta/2.),2)*sin(beta));
	Imaginarytemp=pow(sin(beta/2.),2)*sin(beta)*sin(2*alpha - gamma);
	textBox67->Text= To_String(Realtemp);
	textBox66->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+="               "+To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[4,2] : calculation, updating textboxes
	Realtemp=((cos(beta) - cos(2*beta))*cos(alpha - gamma))/2.;
	Imaginarytemp=((-cos(beta) + cos(2*beta))*sin(alpha - gamma))/2.;
	textBox45->Text= To_String(Realtemp);
	textBox44->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[4,3] : calculation, updating textboxes
	Realtemp=-(sqrt(1.5)*cos(beta)*cos(gamma)*sin(beta));
	Imaginarytemp=-(sqrt(1.5)*cos(beta)*sin(beta)*sin(gamma));
	textBox48->Text= To_String(Realtemp);
	textBox43->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[4,4] : calculation, updating textboxes
	Realtemp=((cos(beta) + cos(2*beta))*cos(alpha + gamma))/2.;
	Imaginarytemp=((cos(beta) + cos(2*beta))*sin(alpha + gamma))/2.;
	textBox79->Text= To_String(Realtemp);
	textBox76->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[4,5] : calculation, updating textboxes
	Realtemp=(cos(2*alpha + gamma)*(2*sin(beta) + sin(2*beta)))/4.;
	Imaginarytemp=((2*sin(beta) + sin(2*beta))*sin(2*alpha + gamma))/4.;
	textBox73->Text= To_String(Realtemp);
	textBox70->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i \r\n";

	////////////////////////////////////////////////////
	//M[5,1] : calculation, updating textboxes
	Realtemp=cos(2*alpha - 2*gamma)*pow(sin(beta/2.),4);
	Imaginarytemp=-(pow(sin(beta/2.),4)*sin(2*alpha - 2*gamma));
	textBox52->Text= To_String(Realtemp);
	textBox51->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+="               "+To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[5,2] : calculation, updating textboxes
	Realtemp=-(cos(alpha - 2*gamma)*pow(sin(beta/2.),2)*sin(beta));
	Imaginarytemp=pow(sin(beta/2.),2)*sin(beta)*sin(alpha - 2*gamma);
	textBox55->Text= To_String(Realtemp);
	textBox54->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[5,3] : calculation, updating textboxes
	Realtemp=(sqrt(1.5)*cos(2*gamma)*pow(sin(beta),2))/2.;
	Imaginarytemp=sqrt(1.5)*cos(gamma)*pow(sin(beta),2)*sin(gamma);
	textBox92->Text= To_String(Realtemp);
	textBox63->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[5,4] : calculation, updating textboxes
	Realtemp=-(cos(alpha + 2*gamma)*(2*sin(beta) + sin(2*beta)))/4.;
	Imaginarytemp=-((2*sin(beta) + sin(2*beta))*sin(alpha + 2*gamma))/4.;
	textBox101->Text= To_String(Realtemp);
	textBox100->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i  ";

	//M[5,5] : calculation, updating textboxes
	Realtemp=pow(cos(beta/2.),4)*cos(2*(alpha + gamma));
	Imaginarytemp=pow(cos(beta/2.),4)*sin(2*(alpha + gamma));
	textBox99->Text= To_String(Realtemp);
	textBox98->Text= To_String(Imaginarytemp);

	//Write to matlab format, checking if needed the + sign first
	sign=(Imaginarytemp<0?"":"+");
	WignerMatrixText+=To_StringE(Realtemp)+sign+To_StringE(Imaginarytemp)+"i ]; \r\n";
}

/**
* brief Function for updating quaternion textboxes after calculating them from orientation.
*/
void SpinachGUI::OrientationDialog::UpdateQuaternion()
{
	QuaternionText=gcnew String("quaternion=[ ");

	//Get quaternion parameters from orientation object
	Quaternion^ quater=nOrient->GetAsQuaternion();

	//set the quaternion textboxes
	textBoxquat1->Text= To_String(quater->w);
	textBoxquat2->Text= To_String(quater->x);
	textBoxquat3->Text= To_String(quater->y);
	textBoxquat4->Text= To_String(quater->z);

	//Writing to matlab format
	QuaternionText+=To_StringE(quater->w)+" "+To_StringE(quater->x)+" "+
		To_StringE(quater->y)+" "+ To_StringE(quater->z)+" ]; \r\n";
}

/**
* brief Function for updating angle axis textboxes after calculating them from orientation.
*/
void SpinachGUI::OrientationDialog::UpdateAngleAxis()
{
	AngleAxisText=gcnew String("angleaxis=[ ");

	//Get angle axis parameters from orientation object
	AngleAxis^ angleax=nOrient->GetAsAngleAxis();

	//set the angle axis textboxes
	textBoxangle1->Text= To_String(angleax->x);
	textBoxangle2->Text= To_String(angleax->y);
	textBoxangle3->Text= To_String(angleax->z);
	textBoxangle4->Text= To_String(180.0f*angleax->angle/PI);

	//Writing to matlab format
	AngleAxisText+=To_StringE(angleax->x)+" "+To_StringE(angleax->y)+" "+
		To_StringE(angleax->z)+" "+ To_StringE(180.0f*angleax->angle/PI)+" ]; \r\n";
}

/**
* brief Function for updating interaction  parameters: ID, interaction kind, unit, ref frame ID.
*/
void SpinachGUI::OrientationDialog::UpdateParameters()
{
 //Update parameters textboxes
 textBoxID->Text=Convert::ToString(ID);					     //ID
 textBoxInteractType->Text=InterKind.ToString()->ToLower();	 //Interaction Kind
 UnittextBox->Text=unit->get_name();						 //Unit
 RefFrametextBox->Text=Convert::ToString(refID);			 //Ref Frame ID
}

/**
* brief Function for updating all the dialog's textboxes after calculating the appropriate  orientation.
*/
void SpinachGUI::OrientationDialog::UpdateAll()
{
	//Refresh the snapshot after aligning the model to this interaction.
	for (int i=0;i<9;i++) Interact_Mat[i]=tensor->matrix3x3[i];

	//Interaction matrix and eigen
	UpdateInteraction();
	Calculate_Eigen();
	UpdateEigenState();
	UpdateIsotropic();
	UpdateSphericaTensor();

	//Orientation, euler angles, Wingner Matrix, quaternion and angles axis parameters
	CreateOrient();
	UpdateEulerAngles();
	UpdateWignerMatrix();
	UpdateQuaternion();
	UpdateAngleAxis();

	//Parameters
	UpdateParameters();
}

#pragma endregion

#pragma region "Calculations and 3d model update"

/**
* @brief Function for calculating  eigen vector and eigen values from tensor matrix.
*/
void SpinachGUI::OrientationDialog::Calculate_Eigen()
{
	//Migrate matrix to eigen matrix
	Matrix3d A;
	A<<Interact_Mat[0],Interact_Mat[1],Interact_Mat[2],
		Interact_Mat[3],Interact_Mat[4],Interact_Mat[5],
		Interact_Mat[6],Interact_Mat[7],Interact_Mat[8];

	//Set eigen solver ans solve 
	SelfAdjointEigenSolver<Matrix3d> eigensolver(A);
	if (eigensolver.info() != Success) abort();
	//Migrate eigen matrix to matrix, eigen values
	eigenval[0]=eigensolver.eigenvalues()[0];
	eigenval[1]=eigensolver.eigenvalues()[1];
	eigenval[2]=eigensolver.eigenvalues()[2];

	//dcm
	//In order to work correcly, according to Spinach
	eigenvec[0]=eigensolver.eigenvectors()(0);
	eigenvec[1]=eigensolver.eigenvectors()(3);
	eigenvec[2]=-eigensolver.eigenvectors()(6);
	eigenvec[3]=eigensolver.eigenvectors()(1);
	eigenvec[4]=eigensolver.eigenvectors()(4);
	eigenvec[5]=-eigensolver.eigenvectors()(7);
	eigenvec[6]=eigensolver.eigenvectors()(2);
	eigenvec[7]=eigensolver.eigenvectors()(5);
	eigenvec[8]=-eigensolver.eigenvectors()(8);
}

/**
* @brief Function for calculating  tensor matrix from eigen vector(dcm) and eigen values.
*/
void SpinachGUI::OrientationDialog::Calculate_TensorfromEigen()
{
	Matrix3d result,eigenVectors,eigenValues;

	//Migrate matrix to eigen matrix
	//In order to work correcly, according to Spinach
	eigenVectors(0)=eigenvec[0];
	eigenVectors(1)=eigenvec[3];
	eigenVectors(2)=eigenvec[6];
	eigenVectors(3)=eigenvec[1];
	eigenVectors(4)=eigenvec[4];
	eigenVectors(5)=eigenvec[7];
	eigenVectors(6)=eigenvec[2];
	eigenVectors(7)=eigenvec[5];
	eigenVectors(8)=eigenvec[8];

	//Migrate matrix to eigen matrix
	eigenValues(0)=eigenval[0];
	eigenValues(1)=0.0;
	eigenValues(2)=0.0;
	eigenValues(3)=0.0;
	eigenValues(4)=eigenval[1];
	eigenValues(5)=0.0;
	eigenValues(6)=0.0;
	eigenValues(7)=0.0;
	eigenValues(8)=eigenval[2];

	//Multiply D*E*D^T
	result= (eigenVectors*eigenValues)*eigenVectors.transpose();

	//Migrate eigen matrix to matrix
	for(int j=0;j<9;j++) Interact_Mat[j]=result(j);
}

/**
* @brief Function for setting up orientation object from dcm matrix.
*/
void SpinachGUI::OrientationDialog::CreateOrient()
{
	Matrix3x3^ AA=gcnew Matrix3x3();

	//Calculate the determinant of DCM
	double detA= eigenvec[0]*eigenvec[4]*eigenvec[8]
	-eigenvec[0]*eigenvec[5]*eigenvec[7]
	-eigenvec[3]*eigenvec[1]*eigenvec[8]
	+eigenvec[3]*eigenvec[2]*eigenvec[7]
	+eigenvec[6]*eigenvec[1]*eigenvec[5]
	-eigenvec[6]*eigenvec[4]*eigenvec[2];

	//setup sigh of det
	if(detA>0) detA=1.0f; else detA=-1.0f;

	//Translate DCM multiplied by its determinant to a Matrix3x3 matrix
	AA->xx=detA*eigenvec[0];AA->xy=detA*eigenvec[1];AA->xz=detA*eigenvec[2];
	AA->yx=detA*eigenvec[3];AA->yy=detA*eigenvec[4];AA->yz=detA*eigenvec[5];
	AA->zx=detA*eigenvec[6];AA->zy=detA*eigenvec[7];AA->zz=detA*eigenvec[8];

	//Setup the initial Orientation object
	nOrient= gcnew SpinachGUI::Orientation(AA);
}

/**
* @brief Function for updating the tensor matrix in the model and repaint 3d model
*/
void SpinachGUI::OrientationDialog::UpdateModel() 
{
	//update tensor matrix in the model
	for (int i=0;i<9;i++) tensor->matrix3x3[i]=Interact_Mat[i];

	//Repaint only interaction draw lists in OpenGL
	Update(true);
	OpenGL->UpdateDisplay(UpdateType::ELLIPSOIDS);
	OpenGL->UpdateDisplay(UpdateType::BONDS_LIKE);
}

#pragma endregion

#pragma region Button clicks

/**
* @brief EventHandler for Export button of the orientation. 
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationDialog::Exportbutton_Click(System::Object^  sender, System::EventArgs^  e) 
{
	//Create and show the dialog box for exporting the orientation
	SpinachGUI::OrientationExport^ orientExport=gcnew SpinachGUI::OrientationExport(this);
	orientExport->ShowDialog();
}

/**
* @brief EventHandler for button for removing anisotropic part from the spherical tensor coeefficients.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationDialog::RemoveAnisbutton_Click(System::Object^  sender, System::EventArgs^  e)
{ 
	//Set T[2,2] real and imaginary part to zero and update
	textBoxSpherical22I->Text= To_String(0.00000);		 
	textBoxSpherical_KeyPress(textBoxSpherical22I,gcnew KeyPressEventArgs((char)13)) ;
	textBoxSpherical22R->Text= To_String(0.00000);
	textBoxSpherical_KeyPress(textBoxSpherical22R,gcnew KeyPressEventArgs((char)13)) ;

	//Set T[2,1]  real and imaginary part to zero and update
	textBoxSpherical21I->Text= To_String(0.00000);
	textBoxSpherical_KeyPress(textBoxSpherical21I,gcnew KeyPressEventArgs((char)13)) ;
	textBoxSpherical21R->Text= To_String(0.00000);
	textBoxSpherical_KeyPress(textBoxSpherical21R, gcnew KeyPressEventArgs((char)13)) ;

	//Set T[2,0] real part to zero and update
	textBoxSpherical20R->Text= To_String(0.00000);
	textBoxSpherical_KeyPress(textBoxSpherical20R, gcnew KeyPressEventArgs((char)13)) ;
}

/**
* @brief EventHandler for button for making this interaction the lab frame 
* and update all the orientation in this dialog and 3d Model.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationDialog::Alignmentbutton_Click(System::Object^  sender, System::EventArgs^  e)
{
	//Make this interaction lab frame
	OpenGL->SystemModel->MakeInteractionLabFrame(tensor->getID());

	//update the whole dialog
	UpdateAll();
	//Update drawing list in openGL
	OpenGL->UpdateDisplay(UpdateType::ALL);
}

#pragma endregion