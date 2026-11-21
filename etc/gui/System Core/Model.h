#pragma once
#include "isotopes.h"
#include "unit.h"
#include <limits>
#include <Dense>

using namespace Eigen;
using namespace System::Collections::Generic;
using namespace System::ComponentModel;
using namespace System::Collections::Specialized;
using namespace System::Collections::ObjectModel;

//Forward declaration
ref class Atom;
ref class ReferenceFrame;

//Definitions
#define epsilo std::numeric_limits<double>::epsilon()
#define dMax   std::numeric_limits<double>::max()

/**
* @brief C++/CLI Swap for exception safe copy.
* @param one First Object.
* @param two Second Object.
*/
template<class T >
Void inline swap( T% one, T% two )
{
	T tmp( one );
	one = two;
	two = tmp;
};

#pragma region "Interaction related classes"

/** \class  InteractionKind
* Collection of Interaction types
*/
public enum class InteractionKind
{
	CBond,        //Chemical Bonding
	HFC,          //Hyperfine Coupling
	CShielding,   //Chemical Shielding
	Shift,        //Chemical Shift
	GTensor,      //G-Tensor
	ZFS,          //Zero Field Splitting
	Exchange,     //Exchange Coupling
	Quadrupolar,  //Quadropolar interaction
	Dipolar,      //Dipolar interaction
	Jcoupling,    //J-Coupling
	spinrotation,  //Spin rotation
	Susceptibility //Molecular magnetic susceptibility: held by Model::Susceptibility, never in InteractionCollection
};

/** \class  Matrix3x3
*  Class for storing a 3x3 observable matrix.
*/
public ref class  Matrix3x3:INotifyPropertyChanged
{
public:
	/**
	* @brief Constructor of Matrix3x3 initalizing ObservableCollection with 9 elements.
	*/
	Matrix3x3(){ matrix3x3=gcnew ObservableCollection<double>(); for (int i=0;i<9;i++) matrix3x3->Add(0.0f);};

	/**
	* @brief Event for when something changes inside the matrix object.
	*/
	virtual event PropertyChangedEventHandler^ PropertyChanged;

	/**
	* @brief Property [] for giving or setting the element of the matrix.
	*/
	property double default[int] {
		double get (int index) {return matrix3x3[index];};
		void set (int index, double value) 
		{
			matrix3x3[index]=value;

			//Call the event handler as element changed
			OnPropertyChanged(gcnew PropertyChangedEventArgs("Matrix3x3"));
		};
};
private:
	ObservableCollection<double>^ matrix3x3;     ///< The observable collection

	/**
	* @brief Event Handler when the collection is changed
	* @param e Arguments for when an property changed. Usually a string.
	*/
	virtual System::Void OnPropertyChanged(PropertyChangedEventArgs^ e) sealed{PropertyChanged(this, e);};
};


/** \class  Interaction
*  Class for storing an observable Interaction.
*/
public ref class Interaction:INotifyPropertyChanged
{
 public:
	 Interaction(Atom^% , Atom^% , InteractionKind );

	 /**
	 * @brief Event for when something changes inside the interaction object
	 */
	 virtual event PropertyChangedEventHandler^ PropertyChanged;

	 int ID;				///< The 1-based ID of interaction
	 String^ reference;		///< The reference substance(only for chemical shifts) of interaction

	 /**
	 * @brief Getter for the interaction ID.
	 * @result Return interaction ID
	 */
	 int getID(){return ID;};

	 //Properties

	 /**
	 * @brief Property dcm for giving or setting the dcm matrix of the interaction.
	 */
	 property ObservableCollection<double>^ dcm{
		 ObservableCollection<double>^ get() {return adcm;}
		 void  set(ObservableCollection<double>^ anotherdcm)
		 {
			 //Change the dcm and fire the change event
			 adcm=anotherdcm;
			 OnPropertyChanged(gcnew PropertyChangedEventArgs("dcm"));
		 }
	 };

	 /**
	 * @brief Property eigen for giving or setting the eigen values of the interaction.
	 */
	 property ObservableCollection<double>^ eigen{
		 ObservableCollection<double>^ get() {return aeigen;}
		 void  set(ObservableCollection<double>^ anothereigen)
		 {
			 //Change the eigen values and fire the change event
			 aeigen=anothereigen;
			 OnPropertyChanged(gcnew PropertyChangedEventArgs("eigen"));
		 }
	 };

	 /**
	 * @brief Property unit for giving or setting the unit of the interaction.
	 */
	 property Unit^ unit{
		 Unit^ get() {return aUnit;}
		 void  set(Unit^ anotherUnit)
		 {
			 //Change the unit and fire the change event
			 aUnit=anotherUnit;
			 OnPropertyChanged(gcnew PropertyChangedEventArgs("Unit"));
		 }
	 };

	 /**
	 * @brief Property IntKind for giving or setting the interaction kind of the interaction.
	 */
	 property InteractionKind IntKind{
		 InteractionKind get() {return intKind;}
		 void  set(InteractionKind aIntKind)
		 {
			 //Try to change kind of interaction, else exception will be thrown
			 try{
				 //Checks the vadility of the interaction using 
				 //new interaction kind, else throws exception.
				 ValidCreationChange(a, b, aIntKind, "change to ");

				 //Save old interaction kind, change interaction kind and fire the change event
				 InteractionKind oldintKind=intKind;
				 intKind=aIntKind;
				 OnPropertyChanged(gcnew PropertyChangedEventArgs("InteractionKind " +oldintKind.ToString()));
			 }
			 catch(String^ e){throw e;};
		 };
	 };

	 /**
	 * @brief Property Label for giving or setting the label of the interaction.
	 */
	 property System::String^ Label{
		 System::String^ get() {return label;}
		 void  set(System::String^ aLabel)
		 {
			 //Change the label and fire the change event
			 label=aLabel;
			 OnPropertyChanged(gcnew PropertyChangedEventArgs("Label"));
		 }
	 };

	 /**
	 * @brief Property A for giving or setting the atom B of the interaction.
	 */
	 property Atom^ A{
		 Atom^ get() {return a;}
		 void  set(Atom^ tA)
		 {
			 //Try to change atom, else exception will be thrown
			 try{
				 //Checks the vadility of the interaction using 
				 //the new A atom and the atom B, else throws exception.
				 ValidCreationChange(tA, b, intKind, "change to ");

				 //Change the atom A and fire the change event
				 a=tA;
				 OnPropertyChanged(gcnew PropertyChangedEventArgs("A"));
			 }
			 catch(String^ e){throw e;};
		 };
	 };

	 /**
	 * @brief Property B for giving or setting the atom B of the interaction.
	 */
	 property Atom^ B{
		 Atom^ get() {return b;}
		 void  set(Atom^ tB)
		 {  
			 //Try to change atom, else exception will be thrown
			 try{
				 //Checks the vadility of the interaction using 
				 //the A atom and the new atom B, else throws exception.
				 ValidCreationChange(a, tB, intKind, "change to ");

				 //Change the atom B and fire the change event
				 b=tB;
				 OnPropertyChanged(gcnew PropertyChangedEventArgs("B"));
			 }
			 catch(String^ e){throw e;};
		 };
	 };

private:
	ObservableCollection<double>^ aeigen;  ///< The matrix with the eigen values.
	ObservableCollection<double>^ adcm;	   ///< The matrix with the eigen vectors (DCM).
	InteractionKind intKind;			   ///< The kind of interaction.
	Unit^ aUnit;						   ///< The units of the interaction.
	String^ label;						   ///< The label of the interaction.
	Atom^ a; Atom ^b;					   ///< The related atoms of the interaction.
	System::Void ValidCreationChange(Atom^ tA, Atom ^tB, InteractionKind aIntKind, String ^ function);

	/**
	* @brief Event Handler when the collections is changed to 
	* call the handler of interaction which is the OnPropertyChanged with eigen string.
	* @param sender The object that fire this event.
	* @param e Arguments for when an property changed. Usually a string.
	*/
	System::Void InteractionChange(System::Object ^ sender, System::Collections::Specialized::NotifyCollectionChangedEventArgs^ e)
	{OnPropertyChanged(gcnew PropertyChangedEventArgs("eigen"));}

protected:
	/**
	* @brief Event Handler when the interaction is changed
	* @param e Arguments for when an property changed. Usually a string.
	*/
	virtual System::Void OnPropertyChanged(PropertyChangedEventArgs^ e) sealed {PropertyChanged(this, e);};
};

/** \class  Tensor
*  Class for storing an observable Tensor which has reference frame.
*/
public ref class Tensor: public Interaction
{
public:
	Tensor(Atom^ A, Atom^ B, InteractionKind Kind, ::Matrix3x3^ anotherMatrix3x3);

	/**
	* @brief Property matrix3x3 for giving or setting the matrix of the tensor of the interaction.
	*/
	property ::Matrix3x3^ matrix3x3{
		::Matrix3x3^ get() {return aMatrix3x3;}
		void  set(::Matrix3x3^ anotherMatrix3x3)
		{
			//Change the matrix of the tensor and fire the change event
			aMatrix3x3->PropertyChanged-=gcnew PropertyChangedEventHandler(this,&Tensor::TensorChange);
			aMatrix3x3=anotherMatrix3x3;
			aMatrix3x3->PropertyChanged+=gcnew PropertyChangedEventHandler(this,&Tensor::TensorChange);
			TensorChange(gcnew System::Object(), gcnew PropertyChangedEventArgs("Matrix3x3"));
		}
	};

	/**
	* @brief Property Frame for giving or setting the reference frame of the tensor of the interaction.
	*/
	property ReferenceFrame^ Frame{
		ReferenceFrame^ get() {return frame;}
		void  set(ReferenceFrame^ aFrame)
		{  
			//Change the reference frame of the tensor and fire the change event
			frame=aFrame;
			TensorChange(gcnew System::Object(), gcnew PropertyChangedEventArgs("Frame"));
		};
	};


internal: 
	System::Void Calculate_EigenStates();

private:
	::Matrix3x3^ aMatrix3x3;	///< The observable matrix of the tensor of the interaction.
	ReferenceFrame^ frame;		///< The reference frame of the tensor of the interaction.

	/**
	* @brief Event Handler when tensor is changed to 
	* call the handler of interaction which is the OnPropertyChanged after refreshing eignstates.
	* @param sender The object that fire this event.
	* @param e Arguments for when an property changed. Usually a string.
	*/
	System::Void TensorChange(System::Object ^ sender,System::ComponentModel::PropertyChangedEventArgs^ e) 
	{Calculate_EigenStates(); OnPropertyChanged(gcnew PropertyChangedEventArgs(e->PropertyName));};

};

/** \class  InteractionsDictionary
*  Class for storing an observable Interactions collection.
*/
public ref class InteractionsDictionary: INotifyCollectionChanged, 
								  INotifyPropertyChanged
{

public:
	InteractionsDictionary();

	/**
	* @brief Event for when something changes inside the collection.
	*/
	virtual event NotifyCollectionChangedEventHandler^ CollectionChanged;

	/**
	* @brief Event for when a property changes. Not used
	*/
	virtual event PropertyChangedEventHandler^ PropertyChanged;

	/**
	* @brief Getter for the number of interactions fo specific kind.
	* @param IntKind The kinf of interaction that the number of which we want.
	* @result Return number of interactions.
	*/
	int GetInterKindCount(InteractionKind IntKind){return InterKindCount[IntKind];};

	//Other functions
	System::Void Add(Interaction^ inter);
	System::Void Clear();
	System::Void RemoveAt(int aID);

	//Properties

	/**
	* @brief Property Count for giving the number of interactions in the collection.
	*/
	property int Count {
		int get ()  {return newID-1;};
	};

	/**
	* @brief Property Keys for giving the keys of the collection.
	*/
	property Dictionary<int, Interaction^>::KeyCollection^ Keys
	{
		Dictionary<int, Interaction^>::KeyCollection^ get ()  {return aDictionary->Keys;};
	};

	/**
	* @brief Property [] for giving or setting an interaction object inside the collection.
	*/
	property Interaction^ default[int] {
		Interaction^ get (int index)  {
			//Check if index is valid else throw exception
			try	{return aDictionary[index];}
			catch(KeyNotFoundException^){throw gcnew System::String("Interaction with ID "+index+" does not exist.");};
		};
		void  set (int index, Interaction^ inter) {
			//Try to change else throw exception
			try{
				Interaction^ old=aDictionary[index];
				inter->ID=old->ID;
				old->PropertyChanged-=gcnew PropertyChangedEventHandler(this,&InteractionsDictionary::InteractionsDictionaryChange);
				InterKindCount[old->IntKind]--;
				InterKindCount[inter->IntKind]++;
				inter->PropertyChanged+=gcnew PropertyChangedEventHandler(this,&InteractionsDictionary::InteractionsDictionaryChange);
				aDictionary[index]=inter;
				UpdateEigenMinMax(nullptr,old->IntKind);
				UpdateEigenMinMax(inter,inter->IntKind);

				//Fire event
				OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Remove, 
					gcnew KeyValuePair<int, Interaction^>(inter->ID, inter)));
			}
			catch(KeyNotFoundException^){throw gcnew System::String("Interaction with ID "+index+" does not exist.");};
		};
	};


	Dictionary<InteractionKind, int>^ InterKindCount;   ///< Collection of the number of each kind of interactions.
	Dictionary<InteractionKind, double>^ InterKindMax;	///< Collection of the maxima eigen value of each kind of interactions.
	Dictionary<InteractionKind, double>^ InterKindMin;	///< Collection of the minima eigen value of each kind of interactions.
	Dictionary<InteractionKind, int>^ iInterKindMax;	///< Collection of the index of interaction with maxima eigen value of each kind of interactions.
	Dictionary<InteractionKind, int>^ iInterKindMin;	///< Collection of the index of interaction with minima eigen value of each kind of interactions.

internal:  //like friend
	int newID;    ///< One-based ID for the new interaction.
	int bondID;	  ///< One-based ID for the new bond interaction.

private:
	Dictionary <int, Interaction^> ^ aDictionary; ///< The observable collection of the interactions.
	System::Void UpdateEigenMinMax(Interaction^ aInteraction, InteractionKind intKind);

	System::Void InteractionsDictionaryChange(System::Object ^ sender, System::ComponentModel::PropertyChangedEventArgs^ e);

	/**
	* @brief Event Handler when the interaction collection is changed
	* @param e Arguments for when an property changed. Usually a string.
	*/
	virtual System::Void OnCollectionChanged(NotifyCollectionChangedEventArgs^ e) sealed {CollectionChanged(this, e);};

};

#pragma endregion

#pragma region "Atom and Reference frame related classes"

/** \class  ReferenceFrame
*  Class for storing an observable reference frame with parent.
*/
public ref class ReferenceFrame:INotifyPropertyChanged
{
public:
	//Constructors
	ReferenceFrame(Matrix3x3^ anotherMatrix3x3)	;
	ReferenceFrame();	

	 /**
	 * @brief Getter for the reference frame ID.
	 * @result Return reference frame ID
	 */
	int getID(){return ID;};

	/**
	* @brief Event for when something changes inside the matrix object.
	*/
	virtual event PropertyChangedEventHandler^ PropertyChanged;

	//Properties

	/**
	* @brief Property matrix3x3 for giving or setting the rotation matrix of the reference frame.
	*/
	property Matrix3x3^ matrix3x3{
		Matrix3x3^ get() {return rotation;}
		void  set(Matrix3x3^ anotherMatrix3x3)
		{
			//Change the rotation matrix and fire the change event
			rotation->PropertyChanged-=gcnew PropertyChangedEventHandler(this,&ReferenceFrame::RotationChange);
			rotation=anotherMatrix3x3;
			rotation->PropertyChanged+=gcnew PropertyChangedEventHandler(this,&ReferenceFrame::RotationChange);
			OnPropertyChanged(gcnew PropertyChangedEventArgs("Rotation"));
		}
	};

	/**
	* @brief Property Label for giving or setting the label of the reference frame.
	*/
	property System::String^ Label{
		System::String^ get() {return label;}
		void  set(System::String^ aLabel)
		{
			//Change the label and fire the change event
			label=aLabel;
			OnPropertyChanged(gcnew PropertyChangedEventArgs("Label"));
		}
	};
	

internal:
	int ID;							 ///< The one-based ID of the reference frame.
	Matrix3x3^ rotation;			 ///< The observable rotation matrix of the reference frame.
	ReferenceFrame ^ ParentRefFrame; ///< The parent reference frame of this reference frame.
	System::String^ label;			 ///< The label of the reference frame.

	/**
	* @brief Event Handler when the reference frame is changed
	* @param e Arguments for when an property changed. Usually a string.
	*/
	virtual System::Void OnPropertyChanged(PropertyChangedEventArgs^ e) sealed{	PropertyChanged(this, e);};

	/**
	* @brief Event Handler when the rotation matrix is changed to 
	* call the handler of reference frame which is the OnPropertyChanged.
	* @param sender The object that fire this event.
	* @param e Arguments for when an property changed. Usually a string.
	*/
	System::Void RotationChange(System::Object ^ sender,System::ComponentModel::PropertyChangedEventArgs^ e) 
	{OnPropertyChanged(gcnew PropertyChangedEventArgs("Rotation"));};
};

/** \class  ReferenceFramesDictionary
*  Class for storing an observable reference frame collection.
*/
public ref class ReferenceFramesDictionary: INotifyCollectionChanged, 
								            INotifyPropertyChanged
{
public:
	ReferenceFramesDictionary(InteractionsDictionary ^%aInteractions);

	/**
	* @brief Event for when something changes inside the collection.
	*/
	virtual event NotifyCollectionChangedEventHandler^ CollectionChanged;

	/**
	* @brief Event for when a property changes. Not used
	*/
	virtual event PropertyChangedEventHandler^ PropertyChanged;

	//Other functions
	System::Void Add(ReferenceFrame^ frame);
	System::Void Clear();
	System::Void RemoveAt(int aID); 
	System::Void MakeFrameLabFrame(ReferenceFrame^ frame);
	
	//Properties

	/**
	* @brief Property Count for giving the number of reference frames in the collection.
	*/
	property int Count {
		int get ()  {return newID-1;};
	};
	
	/**
	* @brief Property Keys for giving the keys of the collection.
	*/
	property Dictionary<int, ReferenceFrame^>::KeyCollection^ Keys{
		Dictionary<int, ReferenceFrame^>::KeyCollection^ get ()  {return aDictionary->Keys;};
	};

	/**
	* @brief Property [] for giving or setting an reference frame object inside the collection.
	*/
	property ReferenceFrame^ default[int] {
		ReferenceFrame^ get (int index)  {
			//Check if index is valid else thro exception
			try	{return aDictionary[index];}
			catch(KeyNotFoundException^){throw gcnew System::String("Reference Frame with ID "+index+" does not exist.");};
		};
		void  set (int index, ReferenceFrame^ frame) {
			//Try to change else throw exception
			try{
				//Set the ID form the old one and set the event handler of the reference frame
				frame->ID=aDictionary[index]->ID;
				frame->PropertyChanged+=gcnew System::ComponentModel::PropertyChangedEventHandler(
					this,&ReferenceFramesDictionary::ReferenceFramesDictionaryChange);

				//Make the change
				aDictionary[index]=frame;

				//Refresh interactions using the new reference frame
				ReferenceFramesDictionaryChange( aDictionary[index],
					gcnew System::ComponentModel::PropertyChangedEventArgs("Reference Frame")); 

				//Fire event
				OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Remove, 
					gcnew KeyValuePair<int, ReferenceFrame^>(frame->ID, frame)));
			}
			catch(KeyNotFoundException^){throw gcnew System::String("Reference Frame with ID "+index+" does not exist.");};
		};
	};



internal:  //like friend
	int newID;    ///< One-based ID for the new reference frame.

private:
	Dictionary <int, ReferenceFrame^> ^ aDictionary; ///< The observable collection of the reference frames.
	InteractionsDictionary ^Interactions;			 ///< The observable collection of the interactions.
	System::Void TransformToFrame(ReferenceFrame^ frame,Tensor ^% inter);
	bool CheckDescendantRecurcively(ReferenceFrame^ Parent, ReferenceFrame^ Descendant);
	System::Void ReIDRecurcively(Dictionary <int, ReferenceFrame^> ^% , ReferenceFrame^ , int & );
	System::Void ReferenceFramesDictionaryChange(System::Object ^ sender, System::ComponentModel::PropertyChangedEventArgs^ e);

	/**
	* @brief Event Handler when the reference frame collection is changed
	* @param e Arguments for when an property changed. Usually a string.
	*/
	virtual System::Void OnCollectionChanged(NotifyCollectionChangedEventArgs^ e) sealed {CollectionChanged(this, e);};
};

/** \class  Atom
*  Class for storing an observable Atom or spin with its properties.
*/
public ref class Atom:INotifyPropertyChanged
{
public:
	/**
	* @brief Constructor of Atom to 0,0,0 postion and H isotope.
	*/
	Atom(System::Void):x(0.0f),y(0.0f),z(0.0f),aIsotope(Isotopes::FindIsotope("H",0, 0)){};

	int ID;  ///< The 1-based ID of atom.

	/**
	* @brief Getter for the atom ID.
	* @result Return atom ID
	*/
	int getID(){return ID;};

	/**
	* @brief Event for when something changes inside the matrix object.
	*/
	virtual event PropertyChangedEventHandler^ PropertyChanged;

	//Properties

	/**
	* @brief Property x coordinate position  of the atom.
	*/
	property double X{
		double get() {return x;}
		void  set(double XX)
		{
			//Change the x coordinate and fire the change event
			x=XX;
			OnPropertyChanged(gcnew PropertyChangedEventArgs("X"));
		}
	};

	/**
	* @brief Property y coordinate position  of the atom.
	*/
	property double Y{
		double get() {return y;}
		void  set(double YY)
		{
			//Change the y coordinate and fire the change event
			y=YY;
			OnPropertyChanged(gcnew PropertyChangedEventArgs("Y"));
		}
	};

	/**
	* @brief Property z coordinate position  of the atom.
	*/
	property double Z{
		double get() {return z;}
		void  set(double ZZ)
		{
			//Change the z coordinate and fire the change event
			z=ZZ;
			OnPropertyChanged(gcnew PropertyChangedEventArgs("Z"));
		}
	};

	/**
	* @brief Property label  of the atom.
	*/
	property System::String^ Label{
		System::String^ get() {return label;}
		void  set(System::String^ aLabel)
		{
			//Change the label and fire the change event
			label=aLabel;
			OnPropertyChanged(gcnew PropertyChangedEventArgs("Label"));
		}
	};

	/**
	* @brief Property isotope  of the atom.
	*/
	property Isotope^ isotope{
		Isotope^ get() {return aIsotope;}
		void  set(Isotope^ anotherIsotope)
		{
			//Remember old element
			String^ iso=aIsotope->Element;

			//Flag for changing from spin different of 0.5 to euqal to zero
			bool chspin=aIsotope->Spin!=0.5 && anotherIsotope->Spin==0.5;

			//Change the isotope and fire the change event
			aIsotope=anotherIsotope;
			OnPropertyChanged(gcnew PropertyChangedEventArgs((chspin? "SpinChange ":"" )+ "Isotope $"+iso));
		}
	};


internal:
	double x,y,z;		   ///< The x,y,z coordinates of the atom.
	Isotope ^ aIsotope;	   ///< The isotope of the atom.
	System::String^ label; ///< The label of the atom.

	/**
	* @brief Event Handler when the atom is changed
	* @param e Arguments for when an property changed. Usually a string.
	*/
	virtual System::Void OnPropertyChanged(PropertyChangedEventArgs^ e) sealed {PropertyChanged(this, e);};
};

/** \class  AtomsDictionary
*  Class for storing an observable Atom or spin collection.
*/
public ref class AtomsDictionary: INotifyCollectionChanged, 
								  INotifyPropertyChanged
{
public:
	AtomsDictionary(InteractionsDictionary ^%aInteractions);

	/**
	* @brief Event for when something changes inside the collection.
	*/
	virtual event NotifyCollectionChangedEventHandler^ CollectionChanged;

	/**
	* @brief Event for when a property changes. Not used
	*/
	virtual event PropertyChangedEventHandler^ PropertyChanged; 

	/**
	* @brief Getter for the number of electrons in the collection.
	* @result Return number of electrons
	*/
	static void CheckAtomCount(int count)
	{
		if(count>4096)
			throw gcnew System::InvalidOperationException("Spinach GUI supports at most 4096 atoms.");
	}

	int get_electron_number(){return electron_number;};

	//Other functions
	System::Void Clear();
	System::Void Add(Atom^ atom);
	System::Void RemoveAt(int aID);

	//Properties

	/**
	* @brief Property refresh for giving or setting the refreshing flag of the atom collection
	* in order to control if the bonds and coordinate maxima are recalculated.
	*/
	property bool Refreshing{
		bool get() {return refresh;}
		void  set(bool Refresh) 
		{ 
			//When the previous state is not refreshing and is change to true
			if(!refresh && Refresh)
			{
				refresh=Refresh;
				//Recalculate bond distances

				//Loop over all bond interactions
				for(int i=-1;i>Interactions->bondID;i--)
				{
					Atom^ A=((Interaction^)Interactions[i])->A;
					Atom^ B=((Interaction^)Interactions[i])->B;

					//Calculating the distance
					double distance=sqrt((A->X-B->X)*(A->X-B->X)+(A->Y-B->Y)*(A->Y-B->Y)+(A->Z-B->Z)*(A->Z-B->Z));

					//Setting to the tensor matrix
					::Matrix3x3 ^Matrix=gcnew ::Matrix3x3();
					Matrix[0]=distance; Matrix[4]=distance; Matrix[8]=distance;
					Matrix[1]=0.0; Matrix[2]=0.0; Matrix[3]=0.0; 
					Matrix[5]=0.0; Matrix[6]=0.0; Matrix[7]=0.0;
					((Tensor^)Interactions[i])->matrix3x3=Matrix;

				};

				MinX=MinY=MinZ=dMax;
				MaxX=MaxY=MaxZ=-dMax;

				//Loop over all atoms to find minim and maximum coordinates
				for each(int i in aDictionary->Keys) 
				{	
					if(aDictionary[i]->isotope->Element!="e")
					{
						//Check if Atom is maximum or minimum and save ID
						if(aDictionary[i]->X < MinX) {MinX=aDictionary[i]->X; iMinX=aDictionary[i]->ID;}
						if(aDictionary[i]->X > MaxX) {MaxX=aDictionary[i]->X; iMaxX=aDictionary[i]->ID;}
						if(aDictionary[i]->Y < MinY) {MinY=aDictionary[i]->Y; iMinY=aDictionary[i]->ID;}
						if(aDictionary[i]->Y > MaxY) {MaxY=aDictionary[i]->Y; iMaxY=aDictionary[i]->ID;}
						if(aDictionary[i]->Z < MinZ) {MinZ=aDictionary[i]->Z; iMinZ=aDictionary[i]->ID;}
						if(aDictionary[i]->Z > MaxZ) {MaxZ=aDictionary[i]->Z; iMaxZ=aDictionary[i]->ID;}
					}
				}
			}
			refresh=Refresh;
		}
	};

	/**
	* @brief Property Count for giving the number of atom in the collection.
	*/
	property int Count {
		int get ()  {return aDictionary->Count;};
	};

	/**
	* @brief Property Keys for giving the keys of the collection.
	*/
	property Dictionary<int, Atom^>::KeyCollection^ Keys{
		Dictionary<int, Atom^>::KeyCollection^ get ()  {return aDictionary->Keys;};
	};

	/**
	* @brief Property [] for giving or setting an Atom object inside the collection.
	*/
	property Atom^ default[int]{
		Atom^ get (int index)  {
			//Check if index is valid else throw exception
			try	{return aDictionary[index];}
			catch(KeyNotFoundException^){throw gcnew System::String("Atom with ID "+index+" does not exist.");};
		};
		void  set (int index, Atom^ atom) {
			//Try to change else throw exception
			try{
				//Remember old element
				String^ iso=aDictionary[index]->aIsotope->Element;

				//Flag for checking spin
				bool chspin=aDictionary[index]->aIsotope->Spin!=0.5 && atom->aIsotope->Spin==0.5;

				//Deep copy
				aDictionary[index]->aIsotope=atom->aIsotope;
				aDictionary[index]->x=atom->x;
				aDictionary[index]->y=atom->y;
				aDictionary[index]->z=atom->z;
				aDictionary[index]->label=atom->label;

				//Check if change is valid
				AtomsDictionaryChange( aDictionary[index],
					gcnew System::ComponentModel::PropertyChangedEventArgs(
					(chspin? "SpinChange ":"" )+ "X Y Z Isotope $"+iso)); 

				//Fire event
				OnCollectionChanged(gcnew NotifyCollectionChangedEventArgs(NotifyCollectionChangedAction::Remove, 
					gcnew KeyValuePair<int, Atom^>(aDictionary[index]->getID(),
					aDictionary[index])));
			}
			catch(KeyNotFoundException^){throw gcnew System::String("Atom with ID "+index+" does not exist.");};
		};
	}
	

	double MinX,MinY,MinZ; ///< Minimum of coordinates of the collection.
	double MaxX,MaxY,MaxZ; ///< Maximum of coordinates of the collection.

private:
	Dictionary <int, Atom^> ^ aDictionary;	///< The observable collection of the atoms.
	int newID;								///< One-based ID for the new element.
	int electron_number;					///< Electron number in the collection.
	bool refresh;							///< Flag for refreshing of the bonds.
	int iMinX,iMinY,iMinZ;					///< Save which element is Minimum of coordinates.
	int iMaxX,iMaxY,iMaxZ;					///< Save which element is Maximum of coordinates.
	InteractionsDictionary ^Interactions;	///< The observable collection of the interactions.

	/**
	* @brief Event Handler when the atom collection is changed
	* @param e Arguments for when an property changed. Usually a string.
	*/
	virtual System::Void OnCollectionChanged(NotifyCollectionChangedEventArgs^ e) sealed{CollectionChanged(this, e);};

	System::Void AtomsDictionaryChange(System::Object ^ sender, System::ComponentModel::PropertyChangedEventArgs^ e);
};

#pragma endregion


/** \class  Model
*  Class for storing an observable Model with all of its properties.
*/
public ref class Model:INotifyPropertyChanged
{
public:
	//Constructors, copy and assigment
	Model(System::Void);
	Model( const Model % aModel);
	Model% operator=( Model aModel);


	AtomsDictionary^ AtomCollection;				///< The collection of atoms in the model.
	InteractionsDictionary^ InteractionCollection;  ///< The collection of interactions in the model.
	ReferenceFramesDictionary^ RefFrameCollection;  ///< The collection of reference frames in the model.
	Tensor^ Susceptibility;  ///< Molecular magnetic susceptibility tensor (cgs-ppm, molecular frame, centre at the origin); nullptr if none. Not a spin interaction.

	/**
	* @brief Event for when something changes inside the model object
	*/
	virtual event PropertyChangedEventHandler^ PropertyChanged;

	System::Void MakeInteractionLabFrame(int ID);

private:
	System::Void ModelChange(System::Object ^ sender, System::Collections::Specialized::NotifyCollectionChangedEventArgs^ e);

	/**
	* @brief Event Handler when the model is changed
	* @param e Arguments for when an property changed. Usually a string.
	*/
	virtual System::Void OnPropertyChanged(PropertyChangedEventArgs^ e) sealed {PropertyChanged(this, e);};
};



