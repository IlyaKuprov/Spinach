// SpinachGUI.cpp : main project file.

#include "stdafx.h"
#include "MainForm.h"

using namespace SpinachGUI;


[STAThreadAttribute]
int main(cli::array<System::String ^> ^args)
{
	// Enabling visual styles before any controls are created
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false); 

	// Create the main window and run it, without any file
	if(args->Length==0) 
	{
		MainForm ^ mainform=gcnew MainForm(gcnew String(""));
		mainform->Visible=false;
		Application::Run(mainform);
	} //or open with the file
	else Application::Run(gcnew MainForm((String^)args[0]));
	return 0;
}        

