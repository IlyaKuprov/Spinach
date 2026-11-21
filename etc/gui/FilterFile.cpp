#include "StdAfx.h"
#include "FilterFile.h"

//Definitions
#define MainICollection  SystemModel->InteractionCollection
#define MainCollection  SystemModel->AtomCollection

/**
*\file FilterFile.cpp
*/

/**
* @brief Function for initialization of parameters and creating distributions.
*/
void SpinachGUI::FilterFile::FilterFile_Initialization()
{
	//initialization of parameters
	MinJC=MinHFC=MinQC=Double::MaxValue;                 //Minima
	MaxJC=MaxHFC=MaxQC=-10000.0;                //Maxima
	maxhistoJC=maxhistoHFC=maxhistoQC=-10000;   //Maxima valus of histograms
	bool flagJC,flagHFC,flagQC;                 //Flags for showing panels
	flagJC=flagHFC=flagQC=false;
	JCbigpanel->Enabled=false;
	HFCbigpanel->Enabled=false;
	QCbigpanel->Enabled=false;

	//Initialization of histograms matrices
	HistogramJC=new int[20];
	HistogramHFC=new int[20];
	HistogramQC=new int[20];
	for(int i=0;i<20;i++) {	HistogramJC[i]=0;HistogramHFC[i]=0;HistogramQC[i]=0;};

	//Loop over interactions to find minima and maxima
	for each(int i in MainICollection->Keys)
	{
		//Find Maximum and Minimum Value of abs(J-Coupling)
		if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::Jcoupling)
		{   
			flagJC=true;

			//save min and max
			if(abs(((Tensor^)MainICollection[i])->matrix3x3[0]) < MinJC)
				MinJC=abs(((Tensor^)MainICollection[i])->matrix3x3[0]);
			if(abs(((Tensor^)MainICollection[i])->matrix3x3[0]) > MaxJC)
				MaxJC=abs(((Tensor^)MainICollection[i])->matrix3x3[0]);
		};

		//Find Maximum and Minimum Value of z-norm for HFC coupling
		if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::HFC)
		{
			flagHFC=true;
			double z_norm=znorm(((Tensor^)MainICollection[i])->matrix3x3);

			//save min and max
			if(z_norm < MinHFC) 	MinHFC=z_norm;
			if(z_norm > MaxHFC)  MaxHFC=z_norm;
		};

		//Find Maximum and Minimum Value of z-norm for Quadrupolar coupling
		if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::Quadrupolar)
		{
			flagQC=true;
			double z_normq=znorm(((Tensor^)MainICollection[i])->matrix3x3);

			//save min and max
			if(z_normq < MinQC)  MinQC=z_normq;
			if(z_normq > MaxQC)  MaxQC=z_normq;
		};
	};

	//Loop over interactions to create distribution matrices
	for each(int i in MainICollection->Keys)
	{
		//Histogram of abs(J-Coupling)
		if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::Jcoupling)
		{   
			HistogramJC[int(20*(abs(((Tensor^)MainICollection[i])->matrix3x3[0])-MinJC)/
				(MaxJC-MinJC+0.00001))]++;
		};

		//Histogram of z-norm for HFC coupling
		if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::HFC)
		{
			Matrix3d result,eigenVectors;
			double z_norm=znorm(((Tensor^)MainICollection[i])->matrix3x3);
			HistogramHFC[int(20*(z_norm-MinHFC)/(MaxHFC-MinHFC+0.00001))]++;
		};

		//Histogram of z-norm for Quadrupolar coupling
		if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::Quadrupolar)
		{
			Matrix3d resultq,eigenVectorsq;
			double znormq=znorm(((Tensor^)MainICollection[i])->matrix3x3);
			HistogramQC[int(20*(znormq-MinQC)/(MaxQC-MinQC+0.00001))]++;
		};
	};

	//Find histograms maxima
	for(int i=0;i<20;i++) 
	{
		if(HistogramJC[i]>maxhistoJC) maxhistoJC=HistogramJC[i];
		if(HistogramHFC[i]>maxhistoHFC) maxhistoHFC=HistogramHFC[i];
		if(HistogramQC[i]>maxhistoQC) maxhistoQC=HistogramQC[i];
	};

	//For Formating numbers
	System::Globalization::CultureInfo^ MyCI = gcnew System::Globalization::CultureInfo( "en-US",false );
	System::Globalization::NumberFormatInfo^ provider = MyCI->NumberFormat;//gcnew System::Globalization::NumberFormatInfo;
	provider->NumberDecimalDigits = 3;  //number of Decimal digits

	//For appearring the j-coupling panels
	if(flagJC)
	{
		JCbigpanel->Enabled=true;

		//Set the parameters for the numeric control
		JCnumericUpDown->Maximum=System::Decimal(MaxJC+0.00001);
		JCnumericUpDown->Minimum=System::Decimal(0.0);
		JCnumericUpDown->Value=System::Decimal(0.0);

		JCthreshold=0.0;//MinJC;

		//Set the label for the maximum, minimum and maximum height
		JCMinlabel->Text=((Double)(MinJC)).ToString("N", provider);
		JCMaxlabel->Text=((Double)(MaxJC)).ToString("N", provider);
		JCMaxFreqlabel->Text=Convert::ToString(maxhistoJC);
	};

	//For appearring the HFC panels
	if(flagHFC)
	{
		HFCbigpanel->Enabled=true;

		//Set the parameters for the numeric control
		HFCnumericUpDown->Maximum=System::Decimal(MaxHFC+0.00001);
		HFCnumericUpDown->Minimum=System::Decimal(0.0);
		HFCnumericUpDown->Value=System::Decimal(0.0);

		HFCthreshold=0.0;//MinHFC;

		//Set the label for the maximum, minimum and maximum height
		HFCMinlabel->Text=((Double)(MinHFC)).ToString("N", provider);
		HFCMaxlabel->Text=((Double)(MaxHFC)).ToString("N", provider);
		HFCMaxFreqlabel->Text=Convert::ToString(maxhistoHFC);
	};

	//For appearring the quadropolar panels
	if(flagQC)
	{
		QCbigpanel->Enabled=true;

		//Set the parameters for the numeric control
		QCnumericUpDown->Maximum=System::Decimal(MaxQC+0.00001);
		QCnumericUpDown->Minimum=System::Decimal(0.0);
		QCnumericUpDown->Value=System::Decimal(0.0);

		QCthreshold=0.0;//MinQC;

		//Set the label for the maximum, minimum and maximum height
		QCMinlabel->Text=((Double)(MinQC)).ToString("N", provider);
		QCMaxlabel->Text=((Double)(MaxQC)).ToString("N", provider);
		QCMaxFreqlabel->Text=Convert::ToString(maxhistoQC);
	};
}

/**
* @brief Function for making the final filtering.
*/
void SpinachGUI::FilterFile::Filtering()
{
	//Create temporary interaction collection with the undeleted interactions
	InteractionsDictionary ^newInteractionDic=gcnew InteractionsDictionary();

	//Loop over all interactions and keep undelete intereactions for bigger than thresholds
	for each(int i in MainICollection->Keys)
	{
		//If not J-Coupling, true-> continue to the next and add to the colelction
		//If J-coupling and bigger than j-coupling threshold, true -> continue to the next 
		//and add to the colelction
		if(!(((Interaction^)MainICollection[i])->IntKind==InteractionKind::Jcoupling &&
			abs(((Tensor^)MainICollection[i])->matrix3x3[0]) <= JCthreshold))
		{   
			//if not HFC and quadropolar, add them to the list
			if(!(((Interaction^)MainICollection[i])->IntKind==InteractionKind::HFC ||
				((Interaction^)MainICollection[i])->IntKind==InteractionKind::Quadrupolar))
			{
				newInteractionDic->Add((Interaction ^)MainICollection[i]);	
			} //when interaction is HFC and quadropolar
			else
			{
				//Calulating z-norm
				double z_norm=znorm(((Tensor^)MainICollection[i])->matrix3x3);;

				//if interaction HFC and bigger than threshold add to collection
				if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::HFC &&
					(z_norm >= HFCthreshold))
				{
					newInteractionDic->Add((Interaction ^)MainICollection[i]);	
				}//if interaction Quadrupolar and bigger than threshold add to collection
				else if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::Quadrupolar &&
					(z_norm >= QCthreshold))
				{
					newInteractionDic->Add((Interaction ^)MainICollection[i]);	
				};
			};
		};
	};

	//Clear interaction collection and add the remaining interactions form temporary collection
	MainICollection->Clear();
	for each(int i in  newInteractionDic->Keys) MainICollection->Add(newInteractionDic[i]);


	//if Remove Orphaned Atoms checked
	if(RemoveAtomsCheckBox->Checked)
	{
		//Create a boolean matrix which flags which atoms to removed, initialize everything to true
		bool *removeAtom=new bool[MainCollection->Count+1];
		for each(int j in MainCollection->Keys) removeAtom[j]=true;

		//Loop over interaction collection and check one by one ineraction marking which atoms exist
		for each(int i in MainICollection->Keys)
		{
			//Do not check bonds, as being created automatically
			if(((Interaction^)MainICollection[i])->IntKind!=InteractionKind::CBond)
			{
				//mark atom A of interactions as existing so false for removing
				int i1=((Interaction^)MainICollection[i])->A->ID; 
				removeAtom[i1]=false;

				//mark atom B(if valid) of interactions as existing so false for removing
				if(((Interaction^)MainICollection[i])->B!=nullptr)
				{int i2=((Interaction^)MainICollection[i])->B->ID; removeAtom[i2]=false;};
			};
		};

		int jj=0;

		//Loop over all atoms checking which one is for deleting and delete it
		//keeping track of the index
		for(int j=1, count=MainCollection->Count;j<=count;j++)
		{
			if(removeAtom[j]) 
			{
				MainCollection->RemoveAt(j-jj);
				jj++;
			};
		};
		delete[] removeAtom;
	};
}

/**
* @brief EventHandler for text in numerics.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::FilterFile::numericUpDown_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) 
{
	//If key pressed is not number, or .  or backspace or enter do nothing
	if (!Char::IsNumber(e->KeyChar) && e->KeyChar != '.' && e->KeyChar != '\b' && e->KeyChar !=(char)13)
	{
		e->Handled = true;
	}
	else
	{
		//When enter is hitted 
		if (e->KeyChar ==(char)13)
		{
			//Returned to original color of the textbox
			((NumericUpDown^)sender)->BackColor=System::Drawing::SystemColors::Window;

			//Choose between j-couplings, HFC, quadropolar and save the value and update
			switch(Convert::ToInt32((String^)(((NumericUpDown^)sender)->Tag))) {
			case 1:
				JCthreshold=Convert::ToDouble(((NumericUpDown^)sender)->Value);
				JCpanel->Refresh();
				break;
			case 2:
				HFCthreshold=Convert::ToDouble(((NumericUpDown^)sender)->Value);
				HFCpanel->Refresh();
				break;
			case 3:
				QCthreshold=Convert::ToDouble(((NumericUpDown^)sender)->Value);
				QCpanel->Refresh();
				break;
			};

		}//When normal writing is happenning the Textbox color change to blue
		else ((NumericUpDown^)sender)->BackColor=System::Drawing::SystemColors::Highlight;
	};
}

/**
* @brief EventHandler for numerics.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::FilterFile::numericUpDown_ValueChanged(System::Object^  sender, System::EventArgs^  e) 
{
	//Choose between j-couplings, HFC, quadropolar and save the value and update
	switch(Convert::ToInt32((String^)(((NumericUpDown^)sender)->Tag))) {
	case 1:
		JCthreshold=Convert::ToDouble(((NumericUpDown^)sender)->Value);
		JCpanel->Refresh();
		break;
	case 2:
		HFCthreshold=Convert::ToDouble(((NumericUpDown^)sender)->Value);
		HFCpanel->Refresh();
		break;
	case 3:
		QCthreshold=Convert::ToDouble(((NumericUpDown^)sender)->Value);
		QCpanel->Refresh();
		break;
	};

}

/**
* @brief EventHandler for Painting quadropolar interacions distribution bars.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::FilterFile::QCpanel_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) 
{
	//Setting the graphics of the sender
	Graphics^ g = e->Graphics;

	//In order toavoid weird result, set  to 1
	if (maxhistoQC<1) maxhistoQC=1;

	//the initial sizes of the bars
	int x1 = 1;
	int y1 = 51;
	int y2 = 50;
	int x2 = 0;

	//Loop over the 20 bars of the distribution
	for(int i=0; i<20; i++) 
	{
		//Setting the width of the bar, constant=12
		x2  = x1 + 12;

		//setting the height of the bar, normalized by the maximum of the distribution
		// and multiplied by the maximu height =50
		y2=int(50*float(HistogramQC[i])/float(maxhistoQC));

		//setting a special brush for when we have strippes
		System::Drawing::Drawing2D::HatchBrush^ myHatchBrush = 
			gcnew System::Drawing::Drawing2D::HatchBrush(
			System::Drawing::Drawing2D::HatchStyle::WideUpwardDiagonal, Color::Gray,  Color::Red);

		//Setting the bar rectangle size
		System::Drawing::Rectangle myRectangle=	System::Drawing::Rectangle(x1,y1-y2,10,y2);

		//Fill the bar rectagle with normal color
		//or with the stripes when the values is smaller than the threshold
		if(i<int(20*(QCthreshold-MinQC)/(MaxQC-MinQC+0.00001)+0.00001))
			g->FillRectangle( myHatchBrush, myRectangle );
		else g->FillRectangle(Brushes::Red, myRectangle);

		//set the end of yhis bar, the beginning of the next one
		x1=x2;
	};
}

/**
* @brief EventHandler for Painting j-couplings interacions distribution bars.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
System::Void SpinachGUI::FilterFile::JCpanel_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) 
{
	//Setting the graphics of the sender
	Graphics^ g = e->Graphics;

	//In order toavoid weird result, set  to 1
	if (maxhistoJC<1) maxhistoJC=1;

	//the initial sizes of the bars
	int x1 = 1;
	int y1 = 51;
	int y2 = 50;
	int x2 = 0;

	//Loop over the 20 bars of the distribution
	for(int i=0; i<20; i++) 
	{
		//Setting the width of the bar, constant=12
		x2  = x1 + 12;

		//setting the height of the bar, normalized by the maximum of the distribution
		// and multiplied by the maximu height =50
		y2=int(50*float(HistogramJC[i])/float(maxhistoJC));

		//setting a special brush for when we have strippes
		System::Drawing::Drawing2D::HatchBrush^ myHatchBrush = 
			gcnew System::Drawing::Drawing2D::HatchBrush(
			System::Drawing::Drawing2D::HatchStyle::WideUpwardDiagonal, Color::Gray,  Color::Red);

		//Setting the bar rectangle size
		System::Drawing::Rectangle myRectangle=	System::Drawing::Rectangle(x1,y1-y2,10,y2);

		//Fill the bar rectagle with normal color
		//or with the stripes when the values is smaller than the threshold
		if(i<int(20*(JCthreshold-MinJC)/(MaxJC-MinJC+0.00001)+0.00001))
			g->FillRectangle( myHatchBrush, myRectangle );
		else g->FillRectangle(Brushes::Red, myRectangle);

		//set the end of yhis bar, the beginning of the next one
		x1=x2;
	};
}

/**
* @brief EventHandler for Painting hyperfine coupling interacions distribution bars.
* @param sender The object that fire this event.
* @param e The event arguments for this event handler.
*/
 System::Void SpinachGUI::FilterFile::HFCpanel_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) 
 {
 	//Setting the graphics of the sender
	Graphics^ g = e->Graphics;

	//In order toavoid weird result, set  to 1
	if (maxhistoHFC<1) maxhistoHFC=1;

	//the initial sizes of the bars
	int x1 = 1;
	int y1 = 51;
	int y2 = 50;
	int x2 = 0;

	//Loop over the 20 bars of the distribution
	for(int i=0; i<20; i++) 
	{
		//Setting the width of the bar, constant=12
		x2  = x1 + 12;

		//setting the height of the bar, normalized by the maximum of the distribution
		// and multiplied by the maximu height =50
		y2=int(50*float(HistogramHFC[i])/float(maxhistoHFC));

		//setting a special brush for when we have strippes
		System::Drawing::Drawing2D::HatchBrush^ myHatchBrush = 
			gcnew System::Drawing::Drawing2D::HatchBrush(
			System::Drawing::Drawing2D::HatchStyle::WideUpwardDiagonal, Color::Gray,  Color::Red);

		//Setting the bar rectangle size
		System::Drawing::Rectangle myRectangle=	System::Drawing::Rectangle(x1,y1-y2,10,y2);

		//Fill the bar rectagle with normal color
		//or with the stripes when the values is smaller than the threshold
		if(i<int(20*(HFCthreshold-MinHFC)/(MaxHFC-MinHFC+0.00001)+0.00001))
			g->FillRectangle( myHatchBrush, myRectangle );
		else g->FillRectangle(Brushes::Red, myRectangle);

		//set the end of yhis bar, the beginning of the next one
		x1=x2;
	};
 }

 /**
* @brief Function calculating the z-norm from a 3x3 matrix.
* @param matrix3x3 The 3x3 matrix which is the z-norm is calculated.
* @result The z-norm value.
*/
double SpinachGUI::FilterFile::znorm(Matrix3x3^ matrix3x3)
{
	//Calulating z-norm
	Matrix3d result,eigenVectors;
	double z_norm;
	for(int j=0;j<9;j++) 
	{
		//Migrate to eigen matrix
		eigenVectors(j)=matrix3x3[j];
	};

	//Calulate M^T*M
	result= eigenVectors.transpose()*eigenVectors;

	//if  sum of diagonal z-norm bigger or equal than zero, 
	//z-norm=sqrt sum of diagonal elements, else throw error
	if ((result(0)+result(4)+result(8))>=0) z_norm=sqrt(result(0)+result(4)+result(8));
	else MessageBox::Show("Problem with Filtering:Negative znorm", "Error", 
		MessageBoxButtons::OK, MessageBoxIcon::Error);

	return z_norm;
}