#include "../StdAfx.h"
#include "Model.h"

/**
*\file Model.cpp
*/

#pragma region "Referenceframe, Tensor and Interaction"

/**
* @brief Constructor of Reference frame.
*/
ReferenceFrame::ReferenceFrame()	 
{ 
	//Set the rotation matrix as ID
	Matrix3x3 ^Matrix=gcnew Matrix3x3();
	Matrix[0]=1.0; Matrix[4]=1.0; Matrix[8]=1.0;
	Matrix[1]=0.0; Matrix[2]=0.0; Matrix[3]=0.0; 
	Matrix[5]=0.0; Matrix[6]=0.0; Matrix[7]=0.0;
	rotation=Matrix;

	//Fire an change event
	OnPropertyChanged(gcnew PropertyChangedEventArgs("Rotation"));

	//Connect the event with the event handler for the rotation matrix
	rotation->PropertyChanged+=
		gcnew System::ComponentModel::PropertyChangedEventHandler(this,&ReferenceFrame::RotationChange);

};

/**
* @brief Constructor of Reference frame using rotation matrix.
* @param  anotherMatrix3x3 The rotation matrix of the tensor.
*/
ReferenceFrame::ReferenceFrame(Matrix3x3^ anotherMatrix3x3)	
{ 
	//Set the rotation matrix
	rotation=anotherMatrix3x3;

	//Fire an change event
	OnPropertyChanged(gcnew PropertyChangedEventArgs("Rotation"));

	//Connect the event with the event handler for the rotation matrix
	rotation->PropertyChanged+=
		gcnew System::ComponentModel::PropertyChangedEventHandler(this,&ReferenceFrame::RotationChange);
};

/**
* @brief Constructor of Tensor inherited by interaction using basic properties.
* @param A Atom A of interaction.
* @param B Atom B of interaction.
* @param Kind Kind of the interaction.
* @param  anotherMatrix3x3 The rotation matrix of the tensor.
*/
Tensor::Tensor(Atom^ A, Atom^ B, InteractionKind Kind, ::Matrix3x3^ anotherMatrix3x3):Interaction(A,B,Kind)
{ 
	//Always Symmetrize-Make Diagonal
	aMatrix3x3=gcnew ::Matrix3x3();
	aMatrix3x3[0]=anotherMatrix3x3[0];
	aMatrix3x3[4]=anotherMatrix3x3[4];
	aMatrix3x3[8]=anotherMatrix3x3[8];
	aMatrix3x3[1]=(anotherMatrix3x3[1]+anotherMatrix3x3[3])/2.0f;
	aMatrix3x3[3]=aMatrix3x3[1];
	aMatrix3x3[2]=(anotherMatrix3x3[2]+anotherMatrix3x3[6])/2.0f;
	aMatrix3x3[6]=aMatrix3x3[2];
	aMatrix3x3[5]=(anotherMatrix3x3[5]+anotherMatrix3x3[7])/2.0f;
	aMatrix3x3[7]=aMatrix3x3[5];

	//Calculate the eigenvalues and eigenvectors of the interaction
	Calculate_EigenStates();

	//Fire an change event
	OnPropertyChanged(gcnew PropertyChangedEventArgs("Matrix3x3"));

	//Connect the event with the event handler for the matrix
	aMatrix3x3->PropertyChanged+=gcnew System::ComponentModel::PropertyChangedEventHandler(this,&Tensor::TensorChange);
};

/**
* @brief Function for calculating the eigenvalues and 
* eigenvectors of the interaction. Eigenvectors are calculated to the root frame.
*/
System::Void Tensor::Calculate_EigenStates()
{
	//Pass matrix to the eigen library
	Matrix3d A;
	A<<aMatrix3x3[0],aMatrix3x3[1],aMatrix3x3[2],
		aMatrix3x3[3],aMatrix3x3[4],aMatrix3x3[5],
		aMatrix3x3[6],aMatrix3x3[7],aMatrix3x3[8];

	//Solve the eigenvalue problem
	SelfAdjointEigenSolver<Matrix3d> eigensolver(A);
	if (eigensolver.info() != Success) abort();

	//Pass it back to the eigen value matrix
	eigen[0]=eigensolver.eigenvalues()[0];
	eigen[1]=eigensolver.eigenvalues()[1];
	eigen[2]=eigensolver.eigenvalues()[2];

	Matrix3d refFrame,result;
	result=eigensolver.eigenvectors();

	//Loop over all parent ref frames, in order to multiply them with the matrix
	for(ReferenceFrame^ next=frame;next!=nullptr;next=next->ParentRefFrame)
	{
		//Pass ref frame matrix to the eigen library
		refFrame<<next->matrix3x3[0],next->matrix3x3[1],next->matrix3x3[2],
			next->matrix3x3[3],next->matrix3x3[4],next->matrix3x3[5],
			next->matrix3x3[6],next->matrix3x3[7],next->matrix3x3[8];

		//Multiply the matrix with the eigenvector matric (dcm)
		result=refFrame*result;
	}

	//Pass it back to the eigen vector matrix
	for(int i=0;i<9;i++) dcm[i]=result(i);
	//Or this
	/*
	dcm[0]=result(0); dcm[1]=result(3); dcm[2]=result(6);
	dcm[3]=result(1); dcm[4]=result(4); dcm[5]=result(7);
	dcm[6]=result(2); dcm[7]=result(5);	dcm[8]=result(8);
	*/
}

/**
* @brief Constructor of Interaction using basic properties.
* @param a Atom A of interaction.
* @param b Atom B of interaction.
* @param kind Kind of the interaction.
*/
Interaction::Interaction(Atom^% a, Atom^% b, InteractionKind kind)
{
	//Try to create valid interaction else throw message
	try{
		//Checks if the attempted creation of interaction is valid
		ValidCreationChange(a, b, kind, "create ");

		//Make 3x1 Matrix with 0 eigen values
		aeigen= gcnew ObservableCollection<double>();
		aeigen->Add(0.0f);aeigen->Add(0.0f);aeigen->Add(0.0f);

		//Make 9x1 Matrix with eigenvectors to 0
		adcm= gcnew ObservableCollection<double>();
		for (int i=0;i<9;i++) adcm->Add(0.0f);
		//Set atom A, B, kind of interaction
		A=a;B=b;
		IntKind=kind;

		//Connect the events with the event handlers
		aeigen->CollectionChanged+=gcnew System::Collections::Specialized::NotifyCollectionChangedEventHandler(this,&Interaction::InteractionChange);
		adcm->CollectionChanged+=gcnew System::Collections::Specialized::NotifyCollectionChangedEventHandler(this,&Interaction::InteractionChange);	  
	}
	catch(String^ e){throw e;};
};

/**
* @brief Function that checks if the attempted interaction is valid.
* @param tA Atom A of interaction.
* @param tB Atom B of interaction.
* @param aIntKind Kind of the interaction.
* @param function Function that is attempted, a string usually create or change.
*/
System::Void Interaction::ValidCreationChange(Atom^ tA, Atom ^tB, InteractionKind aIntKind, String ^ function)
{
	//Quadropolar and  are valid for isotopes only
	if(aIntKind==InteractionKind::Quadrupolar && !(tA->isotope->Spin!=0.5 && 
		tA->isotope->Element!="e" && tB->isotope->Element!="e" && tA==tB))
		throw gcnew System::String("Cannot "+function+"this interaction: Valid only for spin>0.5 and no electrons.");
	//Chemical shielding and Shift are valid for isotopes only
	else if((aIntKind==InteractionKind::CShielding || aIntKind==InteractionKind::Shift ) && 
		!(tA->isotope->Element!="e" && tB->isotope->Element!="e" && tA==tB))
		throw gcnew System::String("Cannot "+function+"this interaction: Valid only for no electrons.");
	//Zerofield tensor and Gtensor are valid on electrons only
	else if((aIntKind==InteractionKind::ZFS || aIntKind==InteractionKind::GTensor) && 
		!(tA->isotope->Element=="e" && tB->isotope->Element=="e" && tA==tB))
		throw gcnew System::String("Cannot "+function+"this interaction: Valid only for electrons.");
	//Hyperfine is valid when isotope coupled to electron
	else if(aIntKind==InteractionKind::HFC &&  !(tA->isotope->Element!="e" && tB->isotope->Element=="e" && tA!=tB))
		throw gcnew System::String("Cannot "+function+"this interaction: Must couple isotope with electron.");
	//For changing Exhange we need two different electrons
	else if(aIntKind==InteractionKind::Exchange &&  !(tA->isotope->Element=="e" && tB->isotope->Element=="e" && tA!=tB))
		throw gcnew System::String("Cannot "+function+"this interaction: Must couple two different electrons."); 
	//For changing to Jcoupling we need two different elements, no electrons
	else if(aIntKind==InteractionKind::Jcoupling &&  !(tA->isotope->Element!="e" && tB->isotope->Element!="e" && tA!=tB))
		throw gcnew System::String("Cannot "+function+" Jcoupling: Must couple two different elements, no electrons.");
};

#pragma endregion

#pragma region AtomsDictionary

System::Void AtomsDictionary::Clear()
{
	newID=1;
	electron_number=0;
	for each(Atom^ atom in aDictionary->Values)
		atom->PropertyChanged-=gcnew PropertyChangedEventHandler(this,&AtomsDictionary::AtomsDictionaryChange);
	aDictionary->Clear();
	MinX=MinY=MinZ=dMax;
	MaxX=MaxY=MaxZ=-dMax;
};

AtomsDictionary::AtomsDictionary(InteractionsDictionary ^%aInteractions)
{
	aDictionary=gcnew Dictionary <int, Atom^>();
	newID=1;
	electron_number=0;
	Interactions=aInteractions;
	//Synchronize AtomsDictionary_electron_number
	//pin_ptr<int> temp= &electron_number;
	//Interactions->AtomsDictionary_electron_number=temp;
	//Initialization of Minimum and maximum
	MinX=MinY=MinZ=dMax;
	MaxX=MaxY=MaxZ=-dMax;
	refresh=true;
};

System::Void AtomsDictionary::Add(Atom^ atom)
{
	if(atom->isotope->Element!="e") CheckAtomCount(aDictionary->Count-electron_number+1);
	//Add atom in the end of Dictionary
	atom->ID=newID;
	aDictionary->Add(newID,atom); 
	//Check if new atom electron, then increase number of electrons
	if(atom->isotope->Element=="e") electron_number++;
	//Add bonds for the new atom which was added at the end
	Atom^ A=atom;
	for(int i=1;i<aDictionary->Count;i++) 
	{
		Atom^ B=(Atom^)aDictionary[i];
		if (A->isotope->Element!="e" && B->isotope->Element!="e")
		{
			double distance=sqrt((A->X-B->X)*(A->X-B->X)+(A->Y-B->Y)*(A->Y-B->Y)+(A->Z-B->Z)*(A->Z-B->Z));
			::Matrix3x3 ^Matrix=gcnew ::Matrix3x3();
			Matrix[0]=distance; Matrix[4]=distance; Matrix[8]=distance;
			Matrix[1]=0.0; Matrix[2]=0.0; Matrix[3]=0.0; 
			Matrix[5]=0.0; Matrix[6]=0.0; Matrix[7]=0.0;
			Tensor^ nInteraction=gcnew Tensor(A,B,InteractionKind::CBond,Matrix);
			Interactions->Add(nInteraction);
		};
	};

	//Check if Atom is maximum or minimum and save ID
	if(atom->isotope->Element!="e")
	{
		if(atom->X < MinX) {MinX=atom->X; iMinX=atom->ID;}
		if(atom->X > MaxX) {MaxX=atom->X; iMaxX=atom->ID;}
		if(atom->Y < MinY) {MinY=atom->Y; iMinY=atom->ID;}
		if(atom->Y > MaxY) {MaxY=atom->Y; iMaxY=atom->ID;}
		if(atom->Z < MinZ) {MinZ=atom->Z; iMinZ=atom->ID;}
		if(atom->Z > MaxZ) {MaxZ=atom->Z; iMaxZ=atom->ID;}
	};

	//Add event handler to atom	
	//tempAtom=atom;
	atom->PropertyChanged+=gcnew System::ComponentModel::PropertyChangedEventHandler(
		this,&AtomsDictionary::AtomsDictionaryChange);

	//Event Handler
	OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Add, 
		gcnew KeyValuePair<int, Atom^>(newID, atom)));
	//Increase new ID number
	newID++;

}

System::Void AtomsDictionary::RemoveAt(int aID)
{
	//Check existance with a simple test
	try	{int test=((Atom ^)aDictionary[aID])->getID();}
	catch(KeyNotFoundException^){throw gcnew System::String("Atom with ID "+aID+" does not exist.");};
	//Check if electron then reduce number of electrons
	if(((Atom ^)aDictionary[aID])->isotope->Element=="e") electron_number--;
	//Remove interactions with this atom

	InteractionsDictionary ^newInteractionDic=gcnew InteractionsDictionary();
	for each(int i in  Interactions->Keys)
	{
		if(((Interaction^)Interactions[i])->B==nullptr)
		{
			if(((Interaction^)Interactions[i])->A->ID!=aID) 
				newInteractionDic->Add((Interaction ^)Interactions[i]);	
		}
		else if(((Interaction^)Interactions[i])->B->ID!=aID && 
			((Interaction^)Interactions[i])->A->ID!=aID) 
			newInteractionDic->Add((Interaction ^)Interactions[i]);
	};
	Interactions->Clear();
	for each(int i in  newInteractionDic->Keys) Interactions->Add(newInteractionDic[i]);

	//Event Handler
	OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Remove, 
		gcnew KeyValuePair<int, Atom^>(aID, ((Atom ^)this[aID]))));

	//Move eveything one place before ovewriting the desired deleted element
	aDictionary[aID]->PropertyChanged-=gcnew PropertyChangedEventHandler(this,&AtomsDictionary::AtomsDictionaryChange);
	for(int i=aID;i<aDictionary->Count;i++) {aDictionary[i]=aDictionary[i+1];((Atom ^)aDictionary[i])->ID=i;}
	//Remove the last element-duplicate as pre-last 
	aDictionary->Remove(--newID);


	//Recompute bounds and their IDs after compacting the dictionary.
	MinX=MinY=MinZ=dMax;
	MaxX=MaxY=MaxZ=-dMax;
	for each(int i in aDictionary->Keys)
	{
		Atom^ atom=aDictionary[i];
		if(atom->X<MinX) {MinX=atom->X;iMinX=i;}
		if(atom->Y<MinY) {MinY=atom->Y;iMinY=i;}
		if(atom->Z<MinZ) {MinZ=atom->Z;iMinZ=i;}
		if(atom->X>MaxX) {MaxX=atom->X;iMaxX=i;}
		if(atom->Y>MaxY) {MaxY=atom->Y;iMaxY=i;}
		if(atom->Z>MaxZ) {MaxZ=atom->Z;iMaxZ=i;}
	}
};

System::Void AtomsDictionary::AtomsDictionaryChange(System::Object ^ sender,
													System::ComponentModel::PropertyChangedEventArgs^ e) 
{
	OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Remove, 
		gcnew KeyValuePair<int, Atom^>(
		((Atom^) sender)->getID(), 
		((Atom^) sender))));

	//Changes in X,Y,Z
	if(refresh && (e->PropertyName->Contains("X") || e->PropertyName->Contains("Y") || e->PropertyName->Contains("Z")))
	{
		//Recalculate bond distances
		for(int i=-1;i>Interactions->bondID;i--)
		{
			Atom^ A=((Interaction^)Interactions[i])->A;
			Atom^ B=((Interaction^)Interactions[i])->B;
			if(A->getID()==((Atom^) sender)->getID() || B->getID()==((Atom^) sender)->getID())
			{
				double distance=sqrt((A->X-B->X)*(A->X-B->X)+(A->Y-B->Y)*(A->Y-B->Y)+(A->Z-B->Z)*(A->Z-B->Z));
				::Matrix3x3 ^Matrix=gcnew ::Matrix3x3();
				Matrix[0]=distance; Matrix[4]=distance; Matrix[8]=distance;
				Matrix[1]=0.0; Matrix[2]=0.0; Matrix[3]=0.0; 
				Matrix[5]=0.0; Matrix[6]=0.0; Matrix[7]=0.0;
				((Tensor^)Interactions[i])->matrix3x3=Matrix;
			};
		};

		//Recalculate maximums and minimums
		if(e->PropertyName->Contains("X"))
		{
			//If any of elements change to larger than current Max, then this is the new Max
			//else if the element which is Max is changed(to smaller) search all for new Max.
			if(((Atom^) sender)->X > MaxX ) {MaxX=((Atom^) sender)->X;iMaxX=((Atom^) sender)->ID;}
			else if(((Atom^) sender)->ID==iMaxX)
			{
				MaxX=-dMax;
				for each(int i in aDictionary->Keys) if(MaxX<this[i]->X) {MaxX=this[i]->X;iMaxX=i;};
			};			
			//Similarly Min	
			if(((Atom^) sender)->X < MinX ) {MinX=((Atom^) sender)->X;iMinX=((Atom^) sender)->ID;}
			else if(((Atom^) sender)->ID==iMinX)
			{
				MinX=dMax;
				for each(int i in aDictionary->Keys) if(MinX>this[i]->X) {MinX=this[i]->X;iMinX=i;};
			};
		}
		if(e->PropertyName->Contains("Y"))
		{
			//If any of elements change to larger than current Max, then this is the new Max
			//else if the element which is Max is changed(to smaller) search all for new Max.
			if(((Atom^) sender)->Y > MaxY ) {MaxY=((Atom^) sender)->Y;iMaxY=((Atom^) sender)->ID;}
			else if(((Atom^) sender)->ID==iMaxY)
			{
				MaxY=-dMax;
				for each(int i in aDictionary->Keys) if(MaxY<this[i]->Y) {MaxY=this[i]->Y;iMaxY=i;};
			};			
			//Similarly Min	
			if(((Atom^) sender)->Y < MinY ) {MinY=((Atom^) sender)->Y;iMinY=((Atom^) sender)->ID;}
			else if(((Atom^) sender)->ID==iMinY)
			{
				MinY=dMax;
				for each(int i in aDictionary->Keys) if(MinY>this[i]->Y) {MinY=this[i]->Y;iMinY=i;};
			};
		}
		if(e->PropertyName->Contains("Z"))
		{
			//If any of elements change to larger than current Max, then this is the new Max
			//else if the element which is Max is changed(to smaller) search all for new Max.
			if(((Atom^) sender)->Z > MaxZ ) {MaxZ=((Atom^) sender)->Z;iMaxZ=((Atom^) sender)->ID;}
			else if(((Atom^) sender)->ID==iMaxZ)
			{
				MaxZ=-dMax;
				for each(int i in aDictionary->Keys) if(MaxZ<this[i]->Z) {MaxZ=this[i]->Z;iMaxZ=i;};
			};			
			//Similarly Min	
			if(((Atom^) sender)->Z < MinZ ) {MinZ=((Atom^) sender)->Z;iMinZ=((Atom^) sender)->ID;}
			else if(((Atom^) sender)->ID==iMinZ)
			{
				MinZ=dMax;
				for each(int i in aDictionary->Keys) if(MinZ>this[i]->Z) {MinZ=this[i]->Z;iMinZ=i;};
			};
		}; 
	};

	if(e->PropertyName->Contains("Isotope"))
	{
		String ^ isotopeName=e->PropertyName->Substring(e->PropertyName->LastIndexOf("$"))->Replace("$","");
		//New non electron /old electron, decrease number of electrons
		if(isotopeName=="e" && ((Atom ^)sender)->isotope->Element!="e") electron_number--;
		//New electron /old non electron, increase number of electrons
		if(isotopeName!="e" && ((Atom ^)sender)->isotope->Element=="e") electron_number++;

		//If the isotope change to isotope with spin=0.5 remove quadropolar interactions of it
		if(e->PropertyName->Contains("SpinChange") && ((Atom^) sender)->isotope->Spin==0.5 )
		{
			InteractionsDictionary ^newInteractionDic=gcnew InteractionsDictionary();
			for each(int i in Interactions->Keys)
			{
				if((((Atom^) sender)!=((Interaction^)Interactions[i])->A &&
					((Atom^) sender)!=((Interaction^)Interactions[i])->B) || 
					((Interaction^)Interactions[i])->IntKind!=InteractionKind::Quadrupolar)
					newInteractionDic->Add((Interaction ^)Interactions[i]);
			};
			Interactions->Clear();
			for each(int i in  newInteractionDic->Keys) Interactions->Add(newInteractionDic[i]);
		};

		//Change to non electron
		if(isotopeName=="e" && ((Atom^) sender)->isotope->Element!="e")
		{
			InteractionsDictionary ^newInteractionDic=gcnew InteractionsDictionary();
			for each(int i in Interactions->Keys)
			{
				if(!(
					//Zerofield tensor and Gtensor are valid on electrons only
						(((Atom^) sender)==((Interaction^)Interactions[i])->A ||
						((Atom^) sender)==((Interaction^)Interactions[i])->B) && 
						(((Interaction^)Interactions[i])->IntKind==InteractionKind::ZFS ||
						((Interaction^)Interactions[i])->IntKind==InteractionKind::GTensor)
						||
						//Hyperfine is valid when isotope coupled to electron
						((Atom^) sender)==((Interaction^)Interactions[i])->B && 
						((Interaction^)Interactions[i])->IntKind==InteractionKind::HFC
						||
						//Exhange we need two different electrons
						(((Atom^) sender)==((Interaction^)Interactions[i])->A ||
						((Atom^) sender)==((Interaction^)Interactions[i])->B) && 
						((Interaction^)Interactions[i])->IntKind==InteractionKind::Exchange)
						)
						newInteractionDic->Add((Interaction ^)Interactions[i]);
			};
			Interactions->Clear();
			for each(int i in  newInteractionDic->Keys) Interactions->Add(newInteractionDic[i]);

			Atom^ A=((Atom^) sender);
			for(int i=1;i<=aDictionary->Count;i++) 
			{
				Atom^ B=(Atom^)aDictionary[i];
				if (B!=A && B->isotope->Element!="e")
				{
					double distance=sqrt((A->X-B->X)*(A->X-B->X)+(A->Y-B->Y)*(A->Y-B->Y)+(A->Z-B->Z)*(A->Z-B->Z));
					::Matrix3x3 ^Matrix=gcnew ::Matrix3x3();
					Matrix[0]=distance; Matrix[4]=distance; Matrix[8]=distance;
					Matrix[1]=0.0; Matrix[2]=0.0; Matrix[3]=0.0; 
					Matrix[5]=0.0; Matrix[6]=0.0; Matrix[7]=0.0;
					Tensor^ nInteraction=gcnew Tensor(A,B,InteractionKind::CBond,Matrix);
					Interactions->Add(nInteraction);
				};
			};

			//Check if Atom is maximum or minimum and save ID
			if(A->X < MinX) {MinX=A->X; iMinX=A->ID;}
			if(A->X > MaxX) {MaxX=A->X; iMaxX=A->ID;}
			if(A->Y < MinY) {MinY=A->Y; iMinY=A->ID;}
			if(A->Y > MaxY) {MaxY=A->Y; iMaxY=A->ID;}
			if(A->Z < MinZ) {MinZ=A->Z; iMinZ=A->ID;}
			if(A->Z > MaxZ) {MaxZ=A->Z; iMaxZ=A->ID;}

		};

		//Change to electron
		if(isotopeName!="e" && ((Atom^) sender)->isotope->Element=="e")
		{

			InteractionsDictionary ^newInteractionDic=gcnew InteractionsDictionary();
			for each(int i in Interactions->Keys)
			{
				if(!(
					//For adding Jcoupling we need two different elements, no electrons
						(((Atom^) sender)==((Interaction^)Interactions[i])->A ||
						((Atom^) sender)==((Interaction^)Interactions[i])->B) && 
						(((Interaction^)Interactions[i])->IntKind==InteractionKind::Jcoupling ||
						((Interaction^)Interactions[i])->IntKind==InteractionKind::CBond)
						||
						//Hyperfine is valid when isotope coupled to electron
						(((Atom^) sender)==((Interaction^)Interactions[i])->A && 
						(((Interaction^)Interactions[i])->IntKind==InteractionKind::HFC||
						((Interaction^)Interactions[i])->IntKind==InteractionKind::CShielding||
						((Interaction^)Interactions[i])->IntKind==InteractionKind::Shift||
						((Interaction^)Interactions[i])->IntKind==InteractionKind::Quadrupolar)))
						)
						newInteractionDic->Add((Interaction ^)Interactions[i]);
			};
			Interactions->Clear();
			for each(int i in  newInteractionDic->Keys) Interactions->Add(newInteractionDic[i]);
		};


	};
};


#pragma endregion

#pragma region InteractionsDictionary

InteractionsDictionary::InteractionsDictionary()
{
	aDictionary=gcnew Dictionary <int, Interaction^>();
	newID=1;bondID=-1;
	InterKindCount=gcnew Dictionary<InteractionKind, int>();
	InterKindMax=gcnew Dictionary<InteractionKind, double>();
	InterKindMin=gcnew Dictionary<InteractionKind, double>();
	iInterKindMax=gcnew Dictionary<InteractionKind, int>();
	iInterKindMin=gcnew Dictionary<InteractionKind, int>() ;
	//Initiate the number, Max and Min of each interaction kind
	for each (InteractionKind IntKind in Enum::GetValues(InteractionKind::typeid))
	{
		InterKindCount->Add(IntKind,0);
		InterKindMax->Add(IntKind,-1*dMax);
		InterKindMin->Add(IntKind,dMax);
		iInterKindMax->Add(IntKind,0);
		iInterKindMin->Add(IntKind,0);
	};
};

System::Void InteractionsDictionary::InteractionsDictionaryChange(System::Object ^ sender,
																  System::ComponentModel::PropertyChangedEventArgs^ e) 
{
	//Changes in interaction kind
	if(e->PropertyName->Contains("InteractionKind"))
	{ 
		//Increase the number of new kind interaction
		InterKindCount[((Interaction^) sender)->IntKind]++;
		UpdateEigenMinMax((Interaction ^)sender, ((Interaction^) sender)->IntKind);
		//Reduce the number of old kind interaction, stored in the property
		InterKindCount[(InteractionKind) Enum::Parse(InteractionKind::typeid, 
			e->PropertyName->Replace("InteractionKind ",""))]--;
		UpdateEigenMinMax((Interaction ^)sender, (InteractionKind) Enum::Parse(InteractionKind::typeid, 
			e->PropertyName->Replace("InteractionKind ","")));
	};
	//Changes in eigen values
	if(e->PropertyName=="eigen") UpdateEigenMinMax((Interaction ^)sender, ((Interaction^) sender)->IntKind);

	OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Remove, 
		gcnew KeyValuePair<int, Interaction^>(
		((Interaction^) sender)->getID(), 
		((Interaction^) sender))));
};


System::Void InteractionsDictionary::UpdateEigenMinMax(Interaction^ aInteraction, InteractionKind intKind)
{
	//If any of elements change to larger than current Max, then this is the new Max
	//else if the element which is Max is changed(to smaller) search all for new Max.
	//Keep valid cached extrema; rescan only when their owner or value changed.
	int minIndex=iInterKindMin[intKind], maxIndex=iInterKindMax[intKind];
	int minID=(minIndex<0 ? (minIndex-2)/3 : minIndex/3);
	int maxID=(maxIndex<0 ? (maxIndex-2)/3 : maxIndex/3);
	if(aInteraction!=nullptr && aDictionary->ContainsKey(minID) && aDictionary->ContainsKey(maxID) &&
		aDictionary[minID]->IntKind==intKind && aDictionary[maxID]->IntKind==intKind &&
		abs(aDictionary[minID]->eigen[minIndex-3*minID])==InterKindMin[intKind] &&
		abs(aDictionary[maxID]->eigen[maxIndex-3*maxID])==InterKindMax[intKind])
	{
		if(aInteraction->IntKind==intKind)
			for(int j=0;j<3;j++)
			{
				double value=abs(aInteraction->eigen[j]);
				if(value<InterKindMin[intKind]) {InterKindMin[intKind]=value;iInterKindMin[intKind]=3*aInteraction->getID()+j;}
				if(value>InterKindMax[intKind]) {InterKindMax[intKind]=value;iInterKindMax[intKind]=3*aInteraction->getID()+j;}
			}
		return;
	}
	InterKindMin[intKind]=dMax;
	InterKindMax[intKind]=-dMax;
	for each(int i in Keys)
		if(this[i]->IntKind==intKind)
			for(int j=0;j<3;j++)
			{
				double value=abs(this[i]->eigen[j]);
				if(value<InterKindMin[intKind]) {InterKindMin[intKind]=value;iInterKindMin[intKind]=3*i+j;}
				if(value>InterKindMax[intKind]) {InterKindMax[intKind]=value;iInterKindMax[intKind]=3*i+j;}
			}
};

System::Void InteractionsDictionary::RemoveAt(int aID) 
{
	InteractionKind kind=aDictionary[aID]->IntKind;
	InterKindCount[kind]--;
	aDictionary[aID]->PropertyChanged-=gcnew PropertyChangedEventHandler(this,&InteractionsDictionary::InteractionsDictionaryChange);

			OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Remove, 
				gcnew KeyValuePair<int, Interaction^>(aID, ((Interaction ^)aDictionary[aID]))));

			if(aID<0) //Bond
			{
				for(int i=aID;i>bondID+1;i--) {aDictionary[i]=aDictionary[i-1];((Interaction ^)aDictionary[i])->ID=i;}
				try	{aDictionary->Remove(++bondID);}
				catch(KeyNotFoundException^){throw gcnew System::String("Interaction with ID "+aID+" does not exist.");};
			}
			else //Other Kinds
			{
				//Move eveything one place before ovewriting the desired deleted interaction
				for(int i=aID;i<newID-1;i++) {aDictionary[i]=aDictionary[i+1];((Interaction ^)aDictionary[i])->ID=i;}
				//Remove the last interaction-duplicate as pre-last 
				try	{aDictionary->Remove(--newID);}
				catch(KeyNotFoundException^){throw gcnew System::String("Interaction with ID "+aID+" does not exist.");};
			}
	UpdateEigenMinMax(nullptr,kind);
};
System::Void InteractionsDictionary::Clear()
{
	newID=1;bondID=-1;
	for each(Interaction^ inter in aDictionary->Values)
		inter->PropertyChanged-=gcnew PropertyChangedEventHandler(this,&InteractionsDictionary::InteractionsDictionaryChange);
	aDictionary->Clear();
	//Initiate the number, Max and Min of each interaction kind
	for each (InteractionKind IntKind in Enum::GetValues(InteractionKind::typeid))
	{
		InterKindCount[IntKind]=0;
		InterKindMax[IntKind]=-1*dMax;
		InterKindMin[IntKind]=dMax;
		iInterKindMax[IntKind]=0;
		iInterKindMin[IntKind]=0;
	};

};
System::Void InteractionsDictionary::Add(Interaction^ inter)
{
	//Check if bond or other kind
	if(inter->IntKind==InteractionKind::CBond)
	{
		inter->ID=bondID;
		aDictionary->Add(bondID,inter); 
		bondID--;
	}
	else 	//Everything else coupling
	{
		inter->ID=newID;
		aDictionary->Add(newID,inter); 
		newID++;

	};

	//For the new Interaction add event handler
	inter->PropertyChanged+=gcnew System::ComponentModel::PropertyChangedEventHandler(
		this,&InteractionsDictionary::InteractionsDictionaryChange);
	//Increase the number of the specific kind
	InterKindCount[inter->IntKind]++;
	//Check if Interaction is minimum and save ID
	for(int j=0;j<3;j++) if(abs(inter->eigen[j]) < InterKindMin[inter->IntKind]) 
	{
		InterKindMin[inter->IntKind]=abs(inter->eigen[j]); 
		iInterKindMin[inter->IntKind]=3*inter->ID+j;
	};
	//Check if Interaction is maximum and save ID
	for(int j=0;j<3;j++) if(abs(inter->eigen[j]) > InterKindMax[inter->IntKind]) 
	{
		InterKindMax[inter->IntKind]=abs(inter->eigen[j]); 
		iInterKindMax[inter->IntKind]=3*inter->ID+j;
	};

	//Event Handler
	OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Add, 
		gcnew KeyValuePair<int, Interaction^>(inter->getID(), inter)));


};

#pragma endregion

#pragma region ReferenceFramesDictionary
ReferenceFramesDictionary::ReferenceFramesDictionary(InteractionsDictionary ^%aInteractions)
{
	aDictionary=gcnew Dictionary <int, ReferenceFrame^>();
	newID=2;
	Interactions=aInteractions;

	//Create lab-root reference frame
	::Matrix3x3 ^Matrix=gcnew ::Matrix3x3();
	Matrix[0]=1.0; Matrix[4]=1.0; Matrix[8]=1.0;
	Matrix[1]=0.0; Matrix[2]=0.0; Matrix[3]=0.0; 
	Matrix[5]=0.0; Matrix[6]=0.0; Matrix[7]=0.0;
	ReferenceFrame^ frame=gcnew ReferenceFrame(Matrix);
	frame->ID=1;
	frame->ParentRefFrame=nullptr;
	frame->Label="Lab Frame";
	frame->PropertyChanged+=gcnew System::ComponentModel::PropertyChangedEventHandler(
		this,&ReferenceFramesDictionary::ReferenceFramesDictionaryChange);
	aDictionary->Add(1,frame); 


};



static Matrix3d readMatrix(::Matrix3x3^ source)
{
	Matrix3d result;
	for(int row=0;row<3;row++)
		for(int col=0;col<3;col++) result(row,col)=source[3*row+col];
	return result;
}

static ::Matrix3x3^ modelMatrix(const Matrix3d& source)
{
	::Matrix3x3^ result=gcnew ::Matrix3x3();
	for(int row=0;row<3;row++)
		for(int col=0;col<3;col++) result[3*row+col]=source(row,col);
	return result;
}

System::Void ReferenceFramesDictionary::TransformToFrame(ReferenceFrame^ frame,Tensor ^% inter)
{
	Matrix3d matrix=readMatrix(inter->matrix3x3);
	for(ReferenceFrame^ next=inter->Frame;next!=frame->ParentRefFrame;next=next->ParentRefFrame)
	{
		Matrix3d R=readMatrix(next->matrix3x3);
		matrix=R*matrix*R.transpose();
	}
	inter->matrix3x3=modelMatrix(matrix);
}

System::Void ReferenceFramesDictionary::MakeFrameLabFrame(ReferenceFrame^ frame)
{
	if(frame->ParentRefFrame==nullptr) return;

	Matrix3d absolute=Matrix3d::Identity();
	for(ReferenceFrame^ next=frame;next!=nullptr;next=next->ParentRefFrame)
		absolute=readMatrix(next->matrix3x3)*absolute;

	// The laboratory coordinates of atoms and physical tensors stay unchanged.
	// Absorb the selected frame's old orientation into its tensors and children.
	for each(int id in Interactions->Keys)
	{
		Tensor^ tensor=(Tensor^)Interactions[id];
		if(id>0 && tensor->Frame==frame)
			tensor->matrix3x3=modelMatrix(absolute*readMatrix(tensor->matrix3x3)*absolute.transpose());
	}
	for each(ReferenceFrame^ child in aDictionary->Values)
		if(child->ParentRefFrame==frame)
			child->matrix3x3=modelMatrix(absolute*readMatrix(child->matrix3x3));

	// Reverse the path to the old root. Only the first parent needs an
	// absolute rotation; higher reversed edges use the old relative inverse.
	Matrix3d rotation=absolute*readMatrix(frame->matrix3x3).transpose();
	ReferenceFrame^ parent=frame->ParentRefFrame;
	Dictionary<ReferenceFrame^,Matrix3x3^>^ rotations=gcnew Dictionary<ReferenceFrame^,Matrix3x3^>();
	frame->ParentRefFrame=nullptr;
	rotations->Add(frame,modelMatrix(Matrix3d::Identity()));
	ReferenceFrame^ child=frame;
	while(parent!=nullptr)
	{
		ReferenceFrame^ next=parent->ParentRefFrame;
		Matrix3d inverse=readMatrix(parent->matrix3x3).transpose();
		parent->ParentRefFrame=child;
		rotations->Add(parent,modelMatrix(rotation));
		child=parent;
		parent=next;
		rotation=inverse;
	}

	Dictionary<int,ReferenceFrame^>^ renumbered=gcnew Dictionary<int,ReferenceFrame^>();
	int id=1;
	frame->ID=id;
	renumbered->Add(id,frame);
	ReIDRecurcively(renumbered,frame,id);
	aDictionary=renumbered;
	// Setting a rotation notifies observers synchronously; publish only after
	// the new root, all parent links and dictionary IDs are consistent.
	for each(KeyValuePair<ReferenceFrame^,Matrix3x3^> entry in rotations)
		entry.Key->matrix3x3=entry.Value;
}

System::Void ReferenceFramesDictionary::ReIDRecurcively(Dictionary <int, ReferenceFrame^> ^% tempDictionary, ReferenceFrame^ Parent, int & currentID)
{
	for each(int j in aDictionary->Keys)
	{
		if(Parent==((ReferenceFrame^)aDictionary[j])->ParentRefFrame)
		{
			((ReferenceFrame^)aDictionary[j])->ID=++currentID;
			tempDictionary->Add(currentID, ((ReferenceFrame^)aDictionary[j]));
			ReIDRecurcively(tempDictionary, ((ReferenceFrame^)aDictionary[j]), currentID);
		}
	}
}


bool ReferenceFramesDictionary::CheckDescendantRecurcively(ReferenceFrame^ Parent, ReferenceFrame^ Descendant)
{
	if(Descendant==nullptr) return false;

	if(Parent==Descendant) return true;
	else return CheckDescendantRecurcively(Parent, Descendant->ParentRefFrame);

}

System::Void ReferenceFramesDictionary::RemoveAt(int aID) 
{


	if(aID==1 || !aDictionary->ContainsKey(aID)) //Can not delete Root
		throw gcnew System::String("Cannot delete this frame.");
	else
	{
		OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Remove, 
			gcnew KeyValuePair<int, ReferenceFrame^>(aID, ((ReferenceFrame ^)aDictionary[aID]))));

		//Save the parent of the deleted child
		int jj=aDictionary[aID]->ParentRefFrame->getID();


		//Update frame orphaned interactions
		for each (int i in Interactions->Keys)
			if(i>0) 
			{
				for(ReferenceFrame^ next=((Tensor ^)Interactions[i])->Frame;
					next!=nullptr && next!=aDictionary[aID]->ParentRefFrame;
					next=next->ParentRefFrame)
					if(next==aDictionary[aID]) 
					{
						TransformToFrame((ReferenceFrame ^ )this[aID],((Tensor ^)Interactions[i]));
						((Tensor ^)Interactions[i])->Frame=this[jj];
					};
			};


			// Drop the entire removed subtree; a grandchild cannot remain
			// attached to a reference frame no longer present in the model.
			Dictionary<int,ReferenceFrame^>^ surviving=gcnew Dictionary<int,ReferenceFrame^>();
			ReferenceFrame^ removed=aDictionary[aID];
			for each(ReferenceFrame^ candidate in aDictionary->Values)
			{
				if(CheckDescendantRecurcively(removed,candidate))
					candidate->PropertyChanged-=gcnew PropertyChangedEventHandler(this,&ReferenceFramesDictionary::ReferenceFramesDictionaryChange);
				else
				{
					candidate->ID=surviving->Count+1;
					surviving->Add(candidate->ID,candidate);
				}
			}
			aDictionary=surviving;
			newID=aDictionary->Count+1;

	}

};
System::Void ReferenceFramesDictionary::Clear()
{
	aDictionary->Clear();
	newID=2;

	//Create lab-root reference frame
	::Matrix3x3 ^Matrix=gcnew ::Matrix3x3();
	Matrix[0]=1.0; Matrix[4]=1.0; Matrix[8]=1.0;
	Matrix[1]=0.0; Matrix[2]=0.0; Matrix[3]=0.0; 
	Matrix[5]=0.0; Matrix[6]=0.0; Matrix[7]=0.0;
	ReferenceFrame^ frame=gcnew ReferenceFrame(Matrix);
	frame->ID=1;
	frame->ParentRefFrame=nullptr;
	frame->PropertyChanged+=gcnew System::ComponentModel::PropertyChangedEventHandler(
		this,&ReferenceFramesDictionary::ReferenceFramesDictionaryChange);
	aDictionary->Add(1,frame); 

};

System::Void ReferenceFramesDictionary::Add(ReferenceFrame^ frame)
{
	frame->ID=newID;
	aDictionary->Add(newID,frame); 
	newID++;

	//For the new ReferenceFrame add event handler
	frame->PropertyChanged+=gcnew System::ComponentModel::PropertyChangedEventHandler(
		this,&ReferenceFramesDictionary::ReferenceFramesDictionaryChange);


	//Event Handler
	OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Add, 
		gcnew KeyValuePair<int, ReferenceFrame^>(frame->getID(), frame)));


};

System::Void ReferenceFramesDictionary::ReferenceFramesDictionaryChange(System::Object ^ sender,
																		System::ComponentModel::PropertyChangedEventArgs^ e) 
{
	//recalculate for each interactrion the eigen vector and eigenvalues for the new ref frame

	if(Interactions->Count>0) 
		for each(int i in Interactions->Keys)
		{
			if(i>0 && CheckDescendantRecurcively((ReferenceFrame^)sender,((Tensor^)Interactions[i])->Frame)) 
				((Tensor ^)Interactions[i])->Calculate_EigenStates();
		};

		OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Remove, 
			gcnew KeyValuePair<int, ReferenceFrame^>(
			((ReferenceFrame^) sender)->getID(), 
			((ReferenceFrame^) sender))));
};

#pragma endregion

#pragma region Model

/**
* @brief Function for making a specific interaction orientation, the lab frame.
* @param ID The interaction ID.
*/
System::Void Model::MakeInteractionLabFrame(int ID)
{
	Matrix3d A;
	A<<((Tensor ^)InteractionCollection[ID])->matrix3x3[0],
		((Tensor ^)InteractionCollection[ID])->matrix3x3[1],
		((Tensor ^)InteractionCollection[ID])->matrix3x3[2],
		((Tensor ^)InteractionCollection[ID])->matrix3x3[3],
		((Tensor ^)InteractionCollection[ID])->matrix3x3[4],
		((Tensor ^)InteractionCollection[ID])->matrix3x3[5],
		((Tensor ^)InteractionCollection[ID])->matrix3x3[6],
		((Tensor ^)InteractionCollection[ID])->matrix3x3[7],
		((Tensor ^)InteractionCollection[ID])->matrix3x3[8];

	SelfAdjointEigenSolver<Matrix3d> eigensolver(A);
	if (eigensolver.info() != Success) abort();
	Matrix3d Q=eigensolver.eigenvectors();
	Matrix3d QInv=eigensolver.eigenvectors().inverse();


	for each(int i in AtomCollection->Keys)
	{
		Matrix3d XYZ;
		XYZ<<AtomCollection[i]->X, 0.0, 0.0,
			AtomCollection[i]->Y, 0.0, 0.0,
			AtomCollection[i]->Z, 0.0, 0.0;
		XYZ=QInv*XYZ;
		AtomCollection[i]->X=XYZ(0); 
		AtomCollection[i]->Y=XYZ(1);
		AtomCollection[i]->Z=XYZ(2);
	};

	for each(int j in InteractionCollection->Keys)
	{
		if(InteractionCollection[j]->IntKind!=InteractionKind::CBond)
		{
			Matrix3d matrix;
			matrix<<    ((Tensor ^)InteractionCollection[j])->matrix3x3[0],
				((Tensor ^)InteractionCollection[j])->matrix3x3[1],
				((Tensor ^)InteractionCollection[j])->matrix3x3[2],
				((Tensor ^)InteractionCollection[j])->matrix3x3[3],
				((Tensor ^)InteractionCollection[j])->matrix3x3[4],
				((Tensor ^)InteractionCollection[j])->matrix3x3[5],
				((Tensor ^)InteractionCollection[j])->matrix3x3[6],
				((Tensor ^)InteractionCollection[j])->matrix3x3[7],
				((Tensor ^)InteractionCollection[j])->matrix3x3[8];
			matrix=QInv*(matrix*Q);

			for(int i=0;i<9;i++) ((Tensor ^)InteractionCollection[j])->matrix3x3[i]=matrix(i);

		}
	};

	OnPropertyChanged(gcnew PropertyChangedEventArgs("Lab reference frame change"));

}

/**
* @brief Constructor of the model.
*/
Model::Model(System::Void)
{
	//Construct the atoms, interactions and reference frames collections
	InteractionCollection=gcnew InteractionsDictionary();
	AtomCollection=gcnew  AtomsDictionary(InteractionCollection);
	RefFrameCollection=gcnew ReferenceFramesDictionary(InteractionCollection);
	Susceptibility=nullptr;

	//Connect the events of the collections with the event handlers 
	AtomCollection->CollectionChanged+=gcnew System::Collections::Specialized::NotifyCollectionChangedEventHandler(this,&Model::ModelChange);
	InteractionCollection->CollectionChanged+=gcnew System::Collections::Specialized::NotifyCollectionChangedEventHandler(this,&Model::ModelChange);
	RefFrameCollection->CollectionChanged+=gcnew System::Collections::Specialized::NotifyCollectionChangedEventHandler(this,&Model::ModelChange);
}

/**
* @brief Copy Constructor of the model.
* @param aModel The model to be copied.
*/
Model::Model( const Model % aModel )
{
	//Copy each atom, interaction and reference frame

	for each(int i in aModel.RefFrameCollection->Keys)
	{
		RefFrameCollection->Add(aModel.RefFrameCollection[i]);
	};

	for each(int i in aModel.AtomCollection->Keys)
	{
		AtomCollection->Add(aModel.AtomCollection[i]);
	};
	for each(int i in aModel.InteractionCollection->Keys)
	{
		InteractionCollection->Add(aModel.InteractionCollection[i]);
	};
	Susceptibility=aModel.Susceptibility;
}

/**
* @brief Copy assignment Constructor of the model.
* @param aModel The model to be copied.
*/
Model% Model::operator=(Model aModel )
{
	//Swap each collecton
	swap( RefFrameCollection, aModel.RefFrameCollection );
	swap( AtomCollection, aModel.AtomCollection );
	swap( Susceptibility, aModel.Susceptibility );
	swap( InteractionCollection, aModel.InteractionCollection );

	//return the current object
	return *this;
}

/**
* @brief Event Handler when the model is changed to 
* call the handler of model which is the OnPropertyChanged
* @param sender The object that fire this event.
* @param e Arguments for when an property changed. Usually a string.
*/
System::Void Model::ModelChange(System::Object ^ sender,System::Collections::Specialized::NotifyCollectionChangedEventArgs^ e)
{

	OnPropertyChanged(gcnew PropertyChangedEventArgs(sender->ToString())); //Model
}

#pragma endregion



