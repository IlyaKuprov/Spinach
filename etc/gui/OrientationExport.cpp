#include "stdafx.h"
#include "OrientationExport.h"

/**
*\file OrientationExport.cpp
*/

/**
* @brief Writes-updates the preview box with the text.
*/
System::Void SpinachGUI::OrientationExport::Write()
{
	String ^ FileText=gcnew String("");

	//If checkbox checked, write interaction matrix in  preview box
	FileText+=(checkBoxInteract->Checked?((OrientationDialog^ )orientation)->InteractionText : "")+"\r\n";

	//If checkbox checked, write eigen values in  preview box
	FileText+=(checkBoxEigen->Checked?((OrientationDialog^ )orientation)->EigenValuesText : "");

	//If checkbox checked, write eigen vectors in  preview box
	FileText+=(checkBoxEigenVector->Checked?((OrientationDialog^ )orientation)->EigenVectorsText : "")+"\r\n";

	//If checkbox checked, write axiality , rmombicity, skew etc in  preview box
	FileText+=(checkBoxConvention->Checked?((OrientationDialog^ )orientation)->ConventionsText : "");

	//If checkbox checked, write euler angles in  preview box
	FileText+=(checkBoxEuler->Checked?((OrientationDialog^ )orientation)->EulerAnglesText : "");

	//If checkbox checked, write quaternion parameters in  preview box
	FileText+=(checkBoxQuaternion->Checked?((OrientationDialog^ )orientation)->QuaternionText : "");

	//If checkbox checked, write angle axis parameters in  preview box
	FileText+=(checkBoxAngle->Checked?((OrientationDialog^ )orientation)->AngleAxisText : "")+"\r\n";

	//If checkbox checked, write spherical tensor matrix in  preview box
	FileText+=(checkBoxSpherical->Checked?((OrientationDialog^ )orientation)->SphericalTensorText : "")+"\r\n";

	//If checkbox checked, write Wigner matrix in  preview box
	FileText+=(checkBoxWigner->Checked?((OrientationDialog^ )orientation)->WignerMatrixText : "")+"\r\n";

	//Write final preview box
	PreviewtextBox->Text=FileText+"\r\n";
}

/**
* @brief EventHandler for Export button of the orientation 
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::OrientationExport::Exportbutton_Click(System::Object^  sender, System::EventArgs^  e) 
{
	//Setting the format and initialize filename for the save dialoguebox 
	ExportFileDialog->Filter="MatLab File|*.m";
	ExportFileDialog->FileName="InteractionID"+((OrientationDialog^ )orientation)->ID+".m";

	//Opens the save dialogue box and check if the OK button clicked then save to file
	if (ExportFileDialog->ShowDialog() == System::Windows::Forms::DialogResult::OK )
	{	 
		//Create export file
		FileStream^ fs=gcnew FileStream(ExportFileDialog->FileName,FileMode::Create );

		//Convert and write to file
		cli::array<Byte>^ info = (gcnew  System::Text::UTF8Encoding( true ))->GetBytes( PreviewtextBox->Text);
		fs->Write( info, 0, info->Length  );

		fs->Close();
	};
}

