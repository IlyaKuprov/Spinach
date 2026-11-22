#pragma once

#include "FilterFile.h"
#include "ShieldingtoShift.h"
#include "IO Formats/EasySpin.h"
#include "3D Engine/OpenGL.h"
#include "System Core/isotopes.h"
#include "System Core/Model.h"
#include "IO Formats/Gauss.h"
#include "IO Formats/Castep.h"
#include "IO Formats/SpinXML.h"
#include "IO Formats/Spinach.h"
#include "IO Formats/Simpson.h"
#include "IO Formats/XYZ.h"
#include "IO Formats/MOL.h"
#include "IO Formats/Orca.h"
#include "IO Formats/SpinEvolution.h"
#include "OrientationDialog.h"
#include "RefFrameOrientationDialog.h"
#include "IsotopesDialog.h"
#include "About.h"
#include "BalloonToolTip.h"

#include "InteractionTrackBar.h"

namespace SpinachGUI {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace OpenGLForm;


	/// <summary>
	/// Summary for MainForm
	/// </summary>
	public ref class MainForm : public System::Windows::Forms::Form
	{
	public:
	//	MainForm (const MainForm ^ aMainForm){return ((MainForm ^)aMainForm->MemberwiseClone());};
		MainForm(String ^ OpenedFileName)
		{
			
			InitializeComponent();
			
			//
			//TODO: Add the constructor code here
			//
			/*
			this->SetStyle(ControlStyles::UserPaint | 
				           ControlStyles::OptimizedDoubleBuffer |
						   ControlStyles::AllPaintingInWmPaint | 
						   ControlStyles::SupportsTransparentBackColor |
						   ControlStyles::ResizeRedraw, true);
			this->UpdateStyles();
			*/
			    //MaintabControl->ResizeRedraw=false;
			doubleResize=false;
			orientationdialog= gcnew OrientationDialog();
						this->Text+=SpinachGUIver;
			MainForm_Initialization();

				dataGridView1->GetType()->InvokeMember(
				"DoubleBuffered", 
				Reflection::BindingFlags::NonPublic | 
				Reflection::BindingFlags::Instance | 
				Reflection::BindingFlags::SetProperty,
				nullptr, 
				dataGridView1,
				gcnew cli::array< System::Object^  >(1) {true}
				);

				dataGridView2->GetType()->InvokeMember(
				"DoubleBuffered", 
				Reflection::BindingFlags::NonPublic | 
				Reflection::BindingFlags::Instance | 
				Reflection::BindingFlags::SetProperty,
				nullptr, 
				dataGridView2,
				gcnew cli::array< System::Object^  >(1) {true}
				);
				
				
				MaintabControl->GetType()->InvokeMember(
				"DoubleBuffered", 
				Reflection::BindingFlags::NonPublic | 
				Reflection::BindingFlags::Instance | 
				Reflection::BindingFlags::SetProperty,
				nullptr, 
				MaintabControl,
				gcnew cli::array< System::Object^  >(1) {true}
				);

				
				Home->GetType()->InvokeMember(
				"DoubleBuffered", 
				Reflection::BindingFlags::NonPublic | 
				Reflection::BindingFlags::Instance | 
				Reflection::BindingFlags::SetProperty,
				nullptr, 
				Home,
				gcnew cli::array< System::Object^  >(1) {true}
				);

				Visualization->GetType()->InvokeMember(
				"DoubleBuffered", 
				Reflection::BindingFlags::NonPublic | 
				Reflection::BindingFlags::Instance | 
				Reflection::BindingFlags::SetProperty,
				nullptr, 
				Visualization,
				gcnew cli::array< System::Object^  >(1) {true}
				);
				



			if (OpenedFileName!="") OpenSpinXML(OpenedFileName);
			 //this->SetStyle( ControlStyles::OptimizedDoubleBuffer, true );
/*
			this->FinMaintabControl = gcnew myTabControl();
			this->FinMaintabControl->SuspendLayout();
			this->FinMaintabControl->Anchor = this->MaintabControl->Anchor;
			this->MaintabControl->Controls->Remove(this->Home);
			this->MaintabControl->Controls->Remove(this->Visualization);
			this->FinMaintabControl->Controls->Add(this->Home);
			this->FinMaintabControl->Controls->Add(this->Visualization);
			this->FinMaintabControl->Location =this->MaintabControl->Location;
			this->FinMaintabControl->Name =this->MaintabControl->Name;
			this->FinMaintabControl->SelectedIndex = this->MaintabControl->SelectedIndex;
			this->FinMaintabControl->Size =this->MaintabControl->Size;
			this->FinMaintabControl->TabIndex = this->MaintabControl->TabIndex;
			this->Controls->Remove(this->MaintabControl);
			this->Controls->Add(this->FinMaintabControl);
            this->FinMaintabControl->ResumeLayout(false);
	*/		
						
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MainForm()
		{
			if (components)
			{
				delete components;
			}
		}

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
	
		System::ComponentModel::IContainer^ components;
		System::ComponentModel::ComponentResourceManager^  MyResources;
	

private: System::Windows::Forms::ToolTip^  ToolTips;
public: System::Windows::Forms::SaveFileDialog^  SaveFilesDialog;
private: System::Windows::Forms::DataGridViewRow^  Row1;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel4;
private: System::Windows::Forms::Label^  label3;
private: System::Windows::Forms::ColorDialog^  GenColorDialog;
protected: System::Windows::Forms::OpenFileDialog^  OpenFilesDialog;


	public:  OpenGLForm::COpenGL ^ OpenGL;
	private:  int datagridview1_column,datagridview1_row;
	private:  int datagridview2_column,datagridview2_row;
	private: BalloonToolTip^  BalloonToolTip1;
	private: BalloonToolTip^  BalloonToolTip2;
	//public:	  ArrayList^ Isotopes_Collection;
	public:   Model^ SystemModel;
	private:  bool insertAtom;
	private:  bool insertInter;
	private:  int SelectedAtomID,SelectedInterID;
	protected:  String^ CurrentSavedFile; 
	public:   String ^CurrentTitleFile;
	private: System::ComponentModel::ComponentResourceManager^  MyResources2;
	private:  System::EventHandler^ comboevent;
	private: System::Windows::Forms::PrintDialog^  printDialog1;
    private: System::Drawing::Printing::PrintDocument^  DocumenttoPrint;
	

private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel13;
public: System::Windows::Forms::DataGridView^  dataGridView1;
private: System::Windows::Forms::Label^  label10;
private: System::Windows::Forms::SplitContainer^  splitContainer1;
public: System::Windows::Forms::SplitContainer^  splitContainer2;









private: System::Windows::Forms::DataGridViewTextBoxColumn^  Column1;
private: System::Windows::Forms::DataGridViewTextBoxColumn^  Label;
private: System::Windows::Forms::DataGridViewComboBoxColumn^  Column4;
private: System::Windows::Forms::DataGridViewComboBoxColumn^  Column5;
private: System::Windows::Forms::DataGridViewTextBoxColumn^  Column11;
private: System::Windows::Forms::DataGridViewTextBoxColumn^  Column6;
private: System::Windows::Forms::DataGridViewTextBoxColumn^  Column7;
private: System::Windows::Forms::DataGridViewTextBoxColumn^  Column8;
private: System::Windows::Forms::DataGridViewImageColumn^  Column9;
private: System::Windows::Forms::TabControl^  MaintabControl;
private: System::Windows::Forms::TabPage^  Home;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel1;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel3;
private: System::Windows::Forms::Button^  SaveAsbutton;
private: System::Windows::Forms::Button^  OpenFileButton;
private: System::Windows::Forms::Button^  Newbutton;
private: System::Windows::Forms::Button^  AboutButton;
private: System::Windows::Forms::Button^  ExportBitmap;
private: System::Windows::Forms::Button^  Importbutton;
private: System::Windows::Forms::Button^  Printbutton;
private: System::Windows::Forms::Button^  ExitButton;



protected: System::Windows::Forms::Button^  Savebutton;
private: System::Windows::Forms::Label^  label1;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel11;
private: System::Windows::Forms::Label^  label8;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel12;
private: System::Windows::Forms::Button^  ExportSimpson;
private: System::Windows::Forms::Button^  ExportSpinach;
private: System::Windows::Forms::Button^  button2;
private: System::Windows::Forms::Button^  ExportSpinEvolution;

private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel2;
private: System::Windows::Forms::Label^  label2;
private: System::Windows::Forms::Panel^  panel1;
private: System::Windows::Forms::CheckBox^  checkBoxAtomID;
private: System::Windows::Forms::CheckBox^  SelectedInteractionscheckBox;
private: System::Windows::Forms::Label^  BackGroundColorLabel;
private: System::Windows::Forms::Button^  IsotopesButton;
private: System::Windows::Forms::CheckBox^  checkBoxAxes;
public: System::Windows::Forms::PictureBox^  mouseSelectScroll;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel5;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel14;
private: System::Windows::Forms::Button^  button7;
private: System::Windows::Forms::Button^  button5;
private: System::Windows::Forms::Button^  YZbutton;
private: System::Windows::Forms::Button^  XZbutton;
private: System::Windows::Forms::Button^  XYbutton;
private: System::Windows::Forms::Button^  button8;
private: System::Windows::Forms::Label^  label6;
public: System::Windows::Forms::PictureBox^  mouseSelectRight;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel9;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel10;
protected: System::Windows::Forms::Label^  ModelFileName;
public: System::Windows::Forms::PictureBox^  pictureBoxMouse;
public: System::Windows::Forms::PictureBox^  mouseSelectLeft;


private: System::Windows::Forms::Label^  label9;
private: System::Windows::Forms::TabPage^  Visualization;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel6;
private: System::Windows::Forms::Label^  label4;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel7;
private: System::Windows::Forms::Panel^  Bondpanel;
private: System::Windows::Forms::TrackBar^  BondtrackBar;
private: System::Windows::Forms::CheckBox^  BondcheckBox;
private: System::Windows::Forms::Label^  BondColorlabel;
private: System::Windows::Forms::Panel^  CSpanel;
private: System::Windows::Forms::TrackBar^  CStrackBar;
private: System::Windows::Forms::CheckBox^  CScheckBox;
private: System::Windows::Forms::Label^  CSColorlabel;
private: System::Windows::Forms::Panel^  HFCpanel;
private: System::Windows::Forms::TrackBar^  HFCtrackBar;
private: System::Windows::Forms::CheckBox^  HFCcheckBox;
private: System::Windows::Forms::Label^  HFCColorlabel;
private: System::Windows::Forms::Panel^  GTpanel;
private: System::Windows::Forms::TrackBar^  GTtrackBar;
private: System::Windows::Forms::CheckBox^  GTcheckBox;
private: System::Windows::Forms::Label^  GTColorlabel;
private: System::Windows::Forms::Panel^  SUSpanel;
private: System::Windows::Forms::TrackBar^  SUStrackBar;
private: System::Windows::Forms::CheckBox^  SUScheckBox;
private: System::Windows::Forms::Label^  SUSColorlabel;
private: System::Windows::Forms::Panel^  JCpanel;
private: System::Windows::Forms::TrackBar^  JCtrackBar;
private: System::Windows::Forms::CheckBox^  JCcheckBox;
private: System::Windows::Forms::Label^  JCColorlabel;
private: System::Windows::Forms::Panel^  ECpanel;
private: System::Windows::Forms::TrackBar^  ECtrackBar;
private: System::Windows::Forms::CheckBox^  ECcheckBox;
private: System::Windows::Forms::Label^  ECColorlabel;
private: System::Windows::Forms::Panel^  QCpanel;
private: System::Windows::Forms::TrackBar^  QCtrackBar;
private: System::Windows::Forms::CheckBox^  QCcheckBox;
private: System::Windows::Forms::Label^  QCColorlabel;
private: System::Windows::Forms::Panel^  ZFpanel;
private: System::Windows::Forms::TrackBar^  ZFtrackBar;
private: System::Windows::Forms::CheckBox^  ZFcheckBox;
private: System::Windows::Forms::Label^  ZFColorlabel;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel8;
private: System::Windows::Forms::Panel^  panel4;
private: System::Windows::Forms::ToolStrip^  toolStrip2;
private: System::Windows::Forms::ToolStripButton^  toolStripButton1;
private: System::Windows::Forms::ToolStripButton^  toolStripButton2;
private: System::Windows::Forms::ToolStripButton^  toolStripButton3;
private: System::Windows::Forms::CheckBox^  EigenEllipsoidscheckBox;
private: System::Windows::Forms::CheckBox^  EigenAxescheckBox;
private: System::Windows::Forms::Label^  label7;
public: System::Windows::Forms::Timer^  OpenGLtimer;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel15;
public: 
private: System::Windows::Forms::Label^  label5;


private: System::Windows::Forms::DataGridViewTextBoxColumn^  dataGridViewTextBoxColumn1;
private: System::Windows::Forms::DataGridViewTextBoxColumn^  Column10;
private: System::Windows::Forms::DataGridViewComboBoxColumn^  AtomID;
private: System::Windows::Forms::DataGridViewComboBoxColumn^  dataGridViewComboBoxColumn1;
private: System::Windows::Forms::DataGridViewComboBoxColumn^  Column2;
private: System::Windows::Forms::DataGridViewComboBoxColumn^  Column13;
private: System::Windows::Forms::DataGridViewButtonColumn^  Column3;
private: System::Windows::Forms::DataGridViewImageColumn^  Column12;
private: System::Windows::Forms::ImageList^  SpinachGUIimageList;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel16;
private: System::Windows::Forms::TreeView^  RefFramestreeView;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel17;
private: System::Windows::Forms::Button^  AddNodeButton;
private: System::Windows::Forms::Button^  DeleteNodeButton;
private: System::Windows::Forms::Button^  MakeLabButton;



public: System::Windows::Forms::CheckBox^  ReferencecheckBox;
private: System::Windows::Forms::ProgressBar^  GenprogressBar;
private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel18;
private: System::Windows::Forms::StatusStrip^  StatusBar;
public: System::Windows::Forms::ToolStripStatusLabel^  StatusBarLabel;
public: System::Windows::Forms::ContextMenuStrip^  OpenGLcontextMenu;
private: System::Windows::Forms::ToolStripMenuItem^  toolStripMenuItem1;



private: System::ComponentModel::BackgroundWorker^  ImportbackgroundWorker;
private: System::Windows::Forms::CheckBox^  Stereobutton;
private: System::Windows::Forms::CheckBox^  FullScreencheckBox;





	public: int ^testpointer;
	public: System::Void OnMouseDownOpenGL(System::Object^ , System::Windows::Forms::MouseEventArgs^ );
	public: System::Void OnMouseUpOpenGL(System::Object^ , System::Windows::Forms::MouseEventArgs^ );
	private: System::Void MainForm_Initialization();
	private: OrientationDialog^ orientationdialog;
	private: EasySpin^ easyspin;
	private: Spinach^ spinach;
	private: ArrayList^ NodesToDeleted;
	private: Simpson^ simpson;
			 SpinEvolution^ spinEvolution;
			 bool doubleResize;
			 int SplitterDistance1, SplitterDistance2;
			 FormWindowState  FormState;


public: System::Windows::Forms::DataGridView^  dataGridView2;


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle9 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle11 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle12 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle10 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(MainForm::typeid));
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle13 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle15 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle16 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle14 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::TreeNode^  treeNode2 = (gcnew System::Windows::Forms::TreeNode(L"Lab Frame [1]"));








































			this->tableLayoutPanel4 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->dataGridView2 = (gcnew System::Windows::Forms::DataGridView());
			this->dataGridViewTextBoxColumn1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column10 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->AtomID = (gcnew System::Windows::Forms::DataGridViewComboBoxColumn());
			this->dataGridViewComboBoxColumn1 = (gcnew System::Windows::Forms::DataGridViewComboBoxColumn());
			this->Column2 = (gcnew System::Windows::Forms::DataGridViewComboBoxColumn());
			this->Column13 = (gcnew System::Windows::Forms::DataGridViewComboBoxColumn());
			this->Column3 = (gcnew System::Windows::Forms::DataGridViewButtonColumn());
			this->Column12 = (gcnew System::Windows::Forms::DataGridViewImageColumn());
			this->GenColorDialog = (gcnew System::Windows::Forms::ColorDialog());
			this->OpenFilesDialog = (gcnew System::Windows::Forms::OpenFileDialog());
			this->ToolTips = (gcnew System::Windows::Forms::ToolTip(this->components));
			this->SaveAsbutton = (gcnew System::Windows::Forms::Button());
			this->OpenFileButton = (gcnew System::Windows::Forms::Button());
			this->Newbutton = (gcnew System::Windows::Forms::Button());
			this->AboutButton = (gcnew System::Windows::Forms::Button());
			this->ExportBitmap = (gcnew System::Windows::Forms::Button());
			this->Importbutton = (gcnew System::Windows::Forms::Button());
			this->Printbutton = (gcnew System::Windows::Forms::Button());
			this->ExitButton = (gcnew System::Windows::Forms::Button());
			this->Savebutton = (gcnew System::Windows::Forms::Button());
			this->ExportSimpson = (gcnew System::Windows::Forms::Button());
			this->ExportSpinach = (gcnew System::Windows::Forms::Button());
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->ExportSpinEvolution = (gcnew System::Windows::Forms::Button());
			this->BackGroundColorLabel = (gcnew System::Windows::Forms::Label());
			this->IsotopesButton = (gcnew System::Windows::Forms::Button());
			this->BondColorlabel = (gcnew System::Windows::Forms::Label());
			this->CSColorlabel = (gcnew System::Windows::Forms::Label());
			this->HFCColorlabel = (gcnew System::Windows::Forms::Label());
			this->GTColorlabel = (gcnew System::Windows::Forms::Label());
			this->SUSColorlabel = (gcnew System::Windows::Forms::Label());
			this->JCColorlabel = (gcnew System::Windows::Forms::Label());
			this->QCColorlabel = (gcnew System::Windows::Forms::Label());
			this->ZFColorlabel = (gcnew System::Windows::Forms::Label());
			this->ECColorlabel = (gcnew System::Windows::Forms::Label());
			this->AddNodeButton = (gcnew System::Windows::Forms::Button());
			this->DeleteNodeButton = (gcnew System::Windows::Forms::Button());
			this->MakeLabButton = (gcnew System::Windows::Forms::Button());
			this->SaveFilesDialog = (gcnew System::Windows::Forms::SaveFileDialog());
			this->printDialog1 = (gcnew System::Windows::Forms::PrintDialog());
			this->DocumenttoPrint = (gcnew System::Drawing::Printing::PrintDocument());
			this->tableLayoutPanel13 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->dataGridView1 = (gcnew System::Windows::Forms::DataGridView());
			this->Column1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Label = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column4 = (gcnew System::Windows::Forms::DataGridViewComboBoxColumn());
			this->Column5 = (gcnew System::Windows::Forms::DataGridViewComboBoxColumn());
			this->Column11 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column6 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column7 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column8 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column9 = (gcnew System::Windows::Forms::DataGridViewImageColumn());
			this->splitContainer1 = (gcnew System::Windows::Forms::SplitContainer());
			this->splitContainer2 = (gcnew System::Windows::Forms::SplitContainer());
			this->OpenGLcontextMenu = (gcnew System::Windows::Forms::ContextMenuStrip(this->components));
			this->toolStripMenuItem1 = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->MaintabControl = (gcnew System::Windows::Forms::TabControl());
			this->Home = (gcnew System::Windows::Forms::TabPage());
			this->tableLayoutPanel15 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->tableLayoutPanel16 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->RefFramestreeView = (gcnew System::Windows::Forms::TreeView());
			this->SpinachGUIimageList = (gcnew System::Windows::Forms::ImageList(this->components));
			this->tableLayoutPanel17 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->ReferencecheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->tableLayoutPanel1 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->tableLayoutPanel3 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->Stereobutton = (gcnew System::Windows::Forms::CheckBox());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->tableLayoutPanel11 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->tableLayoutPanel12 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->tableLayoutPanel2 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->panel1 = (gcnew System::Windows::Forms::Panel());
			this->FullScreencheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->checkBoxAtomID = (gcnew System::Windows::Forms::CheckBox());
			this->SelectedInteractionscheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->checkBoxAxes = (gcnew System::Windows::Forms::CheckBox());
			this->mouseSelectScroll = (gcnew System::Windows::Forms::PictureBox());
			this->tableLayoutPanel5 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->tableLayoutPanel14 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->button7 = (gcnew System::Windows::Forms::Button());
			this->button5 = (gcnew System::Windows::Forms::Button());
			this->YZbutton = (gcnew System::Windows::Forms::Button());
			this->XZbutton = (gcnew System::Windows::Forms::Button());
			this->XYbutton = (gcnew System::Windows::Forms::Button());
			this->button8 = (gcnew System::Windows::Forms::Button());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->mouseSelectRight = (gcnew System::Windows::Forms::PictureBox());
			this->tableLayoutPanel9 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->tableLayoutPanel10 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->ModelFileName = (gcnew System::Windows::Forms::Label());
			this->pictureBoxMouse = (gcnew System::Windows::Forms::PictureBox());
			this->mouseSelectLeft = (gcnew System::Windows::Forms::PictureBox());
			this->tableLayoutPanel18 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->StatusBar = (gcnew System::Windows::Forms::StatusStrip());
			this->StatusBarLabel = (gcnew System::Windows::Forms::ToolStripStatusLabel());
			this->GenprogressBar = (gcnew System::Windows::Forms::ProgressBar());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->Visualization = (gcnew System::Windows::Forms::TabPage());
			this->tableLayoutPanel6 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->tableLayoutPanel7 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->Bondpanel = (gcnew System::Windows::Forms::Panel());
			this->BondtrackBar = (gcnew InteractionTrackBar());
			this->BondcheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->CSpanel = (gcnew System::Windows::Forms::Panel());
			this->CStrackBar = (gcnew InteractionTrackBar());
			this->CScheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->HFCpanel = (gcnew System::Windows::Forms::Panel());
			this->HFCtrackBar = (gcnew InteractionTrackBar());
			this->HFCcheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->GTpanel = (gcnew System::Windows::Forms::Panel());
			this->GTtrackBar = (gcnew InteractionTrackBar());
			this->GTcheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->SUSpanel = (gcnew System::Windows::Forms::Panel());
			this->SUStrackBar = (gcnew InteractionTrackBar());
			this->SUScheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->JCpanel = (gcnew System::Windows::Forms::Panel());
			this->JCtrackBar = (gcnew InteractionTrackBar());
			this->JCcheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->ECpanel = (gcnew System::Windows::Forms::Panel());
			this->ECtrackBar = (gcnew InteractionTrackBar());
			this->ECcheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->QCpanel = (gcnew System::Windows::Forms::Panel());
			this->QCtrackBar = (gcnew InteractionTrackBar());
			this->QCcheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->ZFpanel = (gcnew System::Windows::Forms::Panel());
			this->ZFtrackBar = (gcnew InteractionTrackBar());
			this->ZFcheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->tableLayoutPanel8 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->panel4 = (gcnew System::Windows::Forms::Panel());
			this->toolStrip2 = (gcnew System::Windows::Forms::ToolStrip());
			this->toolStripButton1 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton2 = (gcnew System::Windows::Forms::ToolStripButton());
			this->toolStripButton3 = (gcnew System::Windows::Forms::ToolStripButton());
			this->EigenEllipsoidscheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->EigenAxescheckBox = (gcnew System::Windows::Forms::CheckBox());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->OpenGLtimer = (gcnew System::Windows::Forms::Timer(this->components));
			this->ImportbackgroundWorker = (gcnew System::ComponentModel::BackgroundWorker());
			this->tableLayoutPanel4->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->dataGridView2))->BeginInit();
			this->tableLayoutPanel13->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->dataGridView1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->splitContainer1))->BeginInit();
			this->splitContainer1->Panel1->SuspendLayout();
			this->splitContainer1->Panel2->SuspendLayout();
			this->splitContainer1->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->splitContainer2))->BeginInit();
			this->splitContainer2->Panel2->SuspendLayout();
			this->splitContainer2->SuspendLayout();
			this->OpenGLcontextMenu->SuspendLayout();
			this->MaintabControl->SuspendLayout();
			this->Home->SuspendLayout();
			this->tableLayoutPanel15->SuspendLayout();
			this->tableLayoutPanel16->SuspendLayout();
			this->tableLayoutPanel17->SuspendLayout();
			this->tableLayoutPanel1->SuspendLayout();
			this->tableLayoutPanel3->SuspendLayout();
			this->tableLayoutPanel11->SuspendLayout();
			this->tableLayoutPanel12->SuspendLayout();
			this->tableLayoutPanel2->SuspendLayout();
			this->panel1->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->mouseSelectScroll))->BeginInit();
			this->tableLayoutPanel5->SuspendLayout();
			this->tableLayoutPanel14->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->mouseSelectRight))->BeginInit();
			this->tableLayoutPanel9->SuspendLayout();
			this->tableLayoutPanel10->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pictureBoxMouse))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->mouseSelectLeft))->BeginInit();
			this->tableLayoutPanel18->SuspendLayout();
			this->StatusBar->SuspendLayout();
			this->Visualization->SuspendLayout();
			this->tableLayoutPanel6->SuspendLayout();
			this->tableLayoutPanel7->SuspendLayout();
			this->Bondpanel->SuspendLayout();
			this->CSpanel->SuspendLayout();
			this->HFCpanel->SuspendLayout();
			this->GTpanel->SuspendLayout();
			this->SUSpanel->SuspendLayout();
			this->JCpanel->SuspendLayout();
			this->ECpanel->SuspendLayout();
			this->QCpanel->SuspendLayout();
			this->ZFpanel->SuspendLayout();
			this->tableLayoutPanel8->SuspendLayout();
			this->panel4->SuspendLayout();
			this->toolStrip2->SuspendLayout();
			this->SuspendLayout();
			// 
			// tableLayoutPanel4
			// 
			this->tableLayoutPanel4->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->tableLayoutPanel4->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel4->ColumnCount = 1;
			this->tableLayoutPanel4->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel4->Controls->Add(this->label3, 0, 0);
			this->tableLayoutPanel4->Controls->Add(this->dataGridView2, 0, 1);
			this->tableLayoutPanel4->Location = System::Drawing::Point(0, 0);
			this->tableLayoutPanel4->Name = L"tableLayoutPanel4";
			this->tableLayoutPanel4->RowCount = 2;
			this->tableLayoutPanel4->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 20)));
			this->tableLayoutPanel4->RowStyles->Add((gcnew System::Windows::Forms::RowStyle()));
			this->tableLayoutPanel4->Size = System::Drawing::Size(408, 627);
			this->tableLayoutPanel4->TabIndex = 2;
			// 
			// label3
			// 
			this->label3->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label3->AutoSize = true;
			this->label3->BackColor = System::Drawing::Color::Gray;
			this->label3->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label3->Location = System::Drawing::Point(1, 1);
			this->label3->Margin = System::Windows::Forms::Padding(0);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(406, 20);
			this->label3->TabIndex = 0;
			this->label3->Text = L"Interactions";
			this->label3->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// dataGridView2
			// 
			this->dataGridView2->AllowUserToResizeColumns = false;
			this->dataGridView2->AllowUserToResizeRows = false;
			this->dataGridView2->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->dataGridView2->AutoSizeColumnsMode = System::Windows::Forms::DataGridViewAutoSizeColumnsMode::Fill;
			dataGridViewCellStyle9->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleCenter;
			dataGridViewCellStyle9->BackColor = System::Drawing::SystemColors::Control;
			dataGridViewCellStyle9->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			dataGridViewCellStyle9->ForeColor = System::Drawing::SystemColors::WindowText;
			dataGridViewCellStyle9->SelectionBackColor = System::Drawing::SystemColors::Highlight;
			dataGridViewCellStyle9->SelectionForeColor = System::Drawing::SystemColors::HighlightText;
			dataGridViewCellStyle9->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->dataGridView2->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle9;
			this->dataGridView2->ColumnHeadersHeight = 45;
			this->dataGridView2->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::DisableResizing;
			this->dataGridView2->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(8) {this->dataGridViewTextBoxColumn1, 
				this->Column10, this->AtomID, this->dataGridViewComboBoxColumn1, this->Column2, this->Column13, this->Column3, this->Column12});
			dataGridViewCellStyle11->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleCenter;
			dataGridViewCellStyle11->BackColor = System::Drawing::SystemColors::Window;
			dataGridViewCellStyle11->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			dataGridViewCellStyle11->ForeColor = System::Drawing::SystemColors::ControlText;
			dataGridViewCellStyle11->SelectionBackColor = System::Drawing::Color::Wheat;
			dataGridViewCellStyle11->SelectionForeColor = System::Drawing::SystemColors::ControlText;
			dataGridViewCellStyle11->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->dataGridView2->DefaultCellStyle = dataGridViewCellStyle11;
			this->dataGridView2->EditMode = System::Windows::Forms::DataGridViewEditMode::EditOnKeystroke;
			this->dataGridView2->Location = System::Drawing::Point(4, 25);
			this->dataGridView2->Name = L"dataGridView2";
			dataGridViewCellStyle12->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle12->BackColor = System::Drawing::SystemColors::Control;
			dataGridViewCellStyle12->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 7.471698F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			dataGridViewCellStyle12->ForeColor = System::Drawing::SystemColors::WindowText;
			dataGridViewCellStyle12->SelectionBackColor = System::Drawing::SystemColors::Highlight;
			dataGridViewCellStyle12->SelectionForeColor = System::Drawing::SystemColors::HighlightText;
			dataGridViewCellStyle12->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->dataGridView2->RowHeadersDefaultCellStyle = dataGridViewCellStyle12;
			this->dataGridView2->RowHeadersWidth = 35;
			this->dataGridView2->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->dataGridView2->ShowCellToolTips = false;
			this->dataGridView2->Size = System::Drawing::Size(400, 600);
			this->dataGridView2->TabIndex = 1;
			this->dataGridView2->CellBeginEdit += gcnew System::Windows::Forms::DataGridViewCellCancelEventHandler(this, &MainForm::dataGridView2_CellBeginEdit);
			this->dataGridView2->CellClick += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MainForm::dataGridView2_CellClick);
			this->dataGridView2->CellEndEdit += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MainForm::dataGridView2_CellEndEdit);
			this->dataGridView2->CellMouseEnter += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MainForm::dataGridView2_CellEnter);
			this->dataGridView2->CurrentCellDirtyStateChanged += gcnew System::EventHandler(this, &MainForm::dataGridView2_CurrentCellDirtyStateChanged);
			this->dataGridView2->RowsAdded += gcnew System::Windows::Forms::DataGridViewRowsAddedEventHandler(this, &MainForm::dataGridView2_RowsAdded);
			this->dataGridView2->RowsRemoved += gcnew System::Windows::Forms::DataGridViewRowsRemovedEventHandler(this, &MainForm::dataGridView2_RowsRemoved);
			this->dataGridView2->SelectionChanged += gcnew System::EventHandler(this, &MainForm::dataGridView2_SelectionChanged);
			this->dataGridView2->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &MainForm::dataGridView2_Paint);
			this->dataGridView2->MouseEnter += gcnew System::EventHandler(this, &MainForm::dataGridView2_MouseEnter);
			// 
			// dataGridViewTextBoxColumn1
			// 
			this->dataGridViewTextBoxColumn1->FillWeight = 70;
			this->dataGridViewTextBoxColumn1->HeaderText = L"ID";
			this->dataGridViewTextBoxColumn1->MinimumWidth = 36;
			this->dataGridViewTextBoxColumn1->Name = L"dataGridViewTextBoxColumn1";
			this->dataGridViewTextBoxColumn1->ReadOnly = true;
			// 
			// Column10
			// 
			this->Column10->FillWeight = 60;
			this->Column10->HeaderText = L"Label";
			this->Column10->MinimumWidth = 33;
			this->Column10->Name = L"Column10";
			// 
			// AtomID
			// 
			this->AtomID->FillWeight = 80;
			this->AtomID->HeaderText = L"Atom ID";
			this->AtomID->MinimumWidth = 39;
			this->AtomID->Name = L"AtomID";
			// 
			// dataGridViewComboBoxColumn1
			// 
			this->dataGridViewComboBoxColumn1->HeaderText = L"Inter. Type";
			this->dataGridViewComboBoxColumn1->Items->AddRange(gcnew cli::array< System::Object^  >(7) {L"Exchange", L"Shift", L"HFC", 
				L"G-tensor", L"Quadrupolar", L"J-coupling", L"Zero-Field"});
			this->dataGridViewComboBoxColumn1->MinimumWidth = 68;
			this->dataGridViewComboBoxColumn1->Name = L"dataGridViewComboBoxColumn1";
			this->dataGridViewComboBoxColumn1->ReadOnly = true;
			// 
			// Column2
			// 
			this->Column2->HeaderText = L"Coupled Atom  ID";
			this->Column2->MinimumWidth = 51;
			this->Column2->Name = L"Column2";
			this->Column2->ReadOnly = true;
			// 
			// Column13
			// 
			this->Column13->HeaderText = L"Ref. Frame ID";
			this->Column13->MinimumWidth = 38;
			this->Column13->Name = L"Column13";
			this->Column13->ReadOnly = true;
			// 
			// Column3
			// 
			this->Column3->HeaderText = L"Magnitude  Orientation";
			this->Column3->MinimumWidth = 65;
			this->Column3->Name = L"Column3";
			this->Column3->ReadOnly = true;
			this->Column3->Text = L"Edit";
			// 
			// Column12
			// 
			this->Column12->AutoSizeMode = System::Windows::Forms::DataGridViewAutoSizeColumnMode::None;
			dataGridViewCellStyle10->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleCenter;
			dataGridViewCellStyle10->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(64)), 
				static_cast<System::Int32>(static_cast<System::Byte>(64)), static_cast<System::Int32>(static_cast<System::Byte>(64)));
			dataGridViewCellStyle10->NullValue = nullptr;
			dataGridViewCellStyle10->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(64)), 
				static_cast<System::Int32>(static_cast<System::Byte>(64)), static_cast<System::Int32>(static_cast<System::Byte>(64)));
			this->Column12->DefaultCellStyle = dataGridViewCellStyle10;
			this->Column12->HeaderText = L"";
			this->Column12->MinimumWidth = 16;
			this->Column12->Name = L"Column12";
			this->Column12->Width = 16;
			// 
			// OpenFilesDialog
			// 
			this->OpenFilesDialog->Filter = L"SpinXML File|*.xml|Gaussian3 File|*.log;*.out|xyz File|*.xyz|All Files|*.*";
			this->OpenFilesDialog->RestoreDirectory = true;
			// 
			// ToolTips
			// 
			this->ToolTips->AutomaticDelay = 50;
			this->ToolTips->AutoPopDelay = 5000;
			this->ToolTips->InitialDelay = 50;
			this->ToolTips->ReshowDelay = 10;
			this->ToolTips->Popup += gcnew System::Windows::Forms::PopupEventHandler(this, &MainForm::ToolTips_Popup);
			// 
			// SaveAsbutton
			// 
			this->SaveAsbutton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"SaveAsbutton.Image")));
			this->SaveAsbutton->Location = System::Drawing::Point(3, 46);
			this->SaveAsbutton->Name = L"SaveAsbutton";
			this->SaveAsbutton->Size = System::Drawing::Size(28, 35);
			this->SaveAsbutton->TabIndex = 7;
			this->ToolTips->SetToolTip(this->SaveAsbutton, L"Save As");
			this->SaveAsbutton->UseVisualStyleBackColor = true;
			this->SaveAsbutton->Click += gcnew System::EventHandler(this, &MainForm::SaveAsButton_Click);
			// 
			// OpenFileButton
			// 
			this->OpenFileButton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"OpenFileButton.Image")));
			this->OpenFileButton->Location = System::Drawing::Point(38, 3);
			this->OpenFileButton->Name = L"OpenFileButton";
			this->OpenFileButton->Size = System::Drawing::Size(28, 36);
			this->OpenFileButton->TabIndex = 1;
			this->ToolTips->SetToolTip(this->OpenFileButton, L"Open File");
			this->OpenFileButton->UseVisualStyleBackColor = true;
			this->OpenFileButton->Click += gcnew System::EventHandler(this, &MainForm::OpenFileButton_Click);
			// 
			// Newbutton
			// 
			this->Newbutton->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 4.754717F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->Newbutton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"Newbutton.Image")));
			this->Newbutton->Location = System::Drawing::Point(3, 3);
			this->Newbutton->Name = L"Newbutton";
			this->Newbutton->Size = System::Drawing::Size(28, 36);
			this->Newbutton->TabIndex = 0;
			this->Newbutton->TextAlign = System::Drawing::ContentAlignment::BottomCenter;
			this->ToolTips->SetToolTip(this->Newbutton, L"New Empty SpinXML");
			this->Newbutton->UseVisualStyleBackColor = true;
			this->Newbutton->Click += gcnew System::EventHandler(this, &MainForm::Newbutton_Click);
			// 
			// AboutButton
			// 
			this->AboutButton->DialogResult = System::Windows::Forms::DialogResult::Cancel;
			this->AboutButton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"AboutButton.Image")));
			this->AboutButton->Location = System::Drawing::Point(143, 3);
			this->AboutButton->Name = L"AboutButton";
			this->AboutButton->Size = System::Drawing::Size(28, 36);
			this->AboutButton->TabIndex = 4;
			this->ToolTips->SetToolTip(this->AboutButton, L"About");
			this->AboutButton->UseVisualStyleBackColor = true;
			this->AboutButton->Click += gcnew System::EventHandler(this, &MainForm::AboutButton_Click);
			// 
			// ExportBitmap
			// 
			this->ExportBitmap->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"ExportBitmap.Image")));
			this->ExportBitmap->Location = System::Drawing::Point(108, 3);
			this->ExportBitmap->Name = L"ExportBitmap";
			this->ExportBitmap->Size = System::Drawing::Size(28, 36);
			this->ExportBitmap->TabIndex = 2;
			this->ToolTips->SetToolTip(this->ExportBitmap, L"Export Bitmap");
			this->ExportBitmap->UseVisualStyleBackColor = true;
			this->ExportBitmap->Click += gcnew System::EventHandler(this, &MainForm::ExportBitmap_Click);
			// 
			// Importbutton
			// 
			this->Importbutton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"Importbutton.Image")));
			this->Importbutton->Location = System::Drawing::Point(73, 3);
			this->Importbutton->Name = L"Importbutton";
			this->Importbutton->Size = System::Drawing::Size(28, 36);
			this->Importbutton->TabIndex = 9;
			this->ToolTips->SetToolTip(this->Importbutton, L"Import");
			this->Importbutton->UseVisualStyleBackColor = true;
			this->Importbutton->Click += gcnew System::EventHandler(this, &MainForm::Importbutton_Click);
			// 
			// Printbutton
			// 
			this->Printbutton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"Printbutton.Image")));
			this->Printbutton->Location = System::Drawing::Point(73, 46);
			this->Printbutton->Name = L"Printbutton";
			this->Printbutton->Size = System::Drawing::Size(28, 35);
			this->Printbutton->TabIndex = 8;
			this->ToolTips->SetToolTip(this->Printbutton, L"Print");
			this->Printbutton->UseVisualStyleBackColor = true;
			this->Printbutton->Click += gcnew System::EventHandler(this, &MainForm::Printbutton_Click);
			// 
			// ExitButton
			// 
			this->ExitButton->BackColor = System::Drawing::Color::White;
			this->ExitButton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"ExitButton.Image")));
			this->ExitButton->Location = System::Drawing::Point(143, 46);
			this->ExitButton->Name = L"ExitButton";
			this->ExitButton->Size = System::Drawing::Size(28, 35);
			this->ExitButton->TabIndex = 5;
			this->ToolTips->SetToolTip(this->ExitButton, L"Exit Spinach Gui");
			this->ExitButton->UseVisualStyleBackColor = true;
			this->ExitButton->Click += gcnew System::EventHandler(this, &MainForm::ExitButton_Click);
			// 
			// Savebutton
			// 
			this->Savebutton->Enabled = false;
			this->Savebutton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"Savebutton.Image")));
			this->Savebutton->Location = System::Drawing::Point(38, 46);
			this->Savebutton->Name = L"Savebutton";
			this->Savebutton->Size = System::Drawing::Size(28, 35);
			this->Savebutton->TabIndex = 3;
			this->ToolTips->SetToolTip(this->Savebutton, L"Save");
			this->Savebutton->UseVisualStyleBackColor = true;
			this->Savebutton->Click += gcnew System::EventHandler(this, &MainForm::Savebutton_Click);
			// 
			// ExportSimpson
			// 
			this->ExportSimpson->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 4.754717F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->ExportSimpson->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"ExportSimpson.Image")));
			this->ExportSimpson->Location = System::Drawing::Point(3, 46);
			this->ExportSimpson->Name = L"ExportSimpson";
			this->ExportSimpson->Size = System::Drawing::Size(28, 35);
			this->ExportSimpson->TabIndex = 5;
			this->ExportSimpson->TextAlign = System::Drawing::ContentAlignment::BottomCenter;
			this->ToolTips->SetToolTip(this->ExportSimpson, L"Export to Simpson");
			this->ExportSimpson->UseVisualStyleBackColor = true;
			this->ExportSimpson->Click += gcnew System::EventHandler(this, &MainForm::ExportSimpson_Click);
			// 
			// ExportSpinach
			// 
			this->ExportSpinach->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 4.754717F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->ExportSpinach->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"ExportSpinach.Image")));
			this->ExportSpinach->Location = System::Drawing::Point(3, 3);
			this->ExportSpinach->Name = L"ExportSpinach";
			this->ExportSpinach->Size = System::Drawing::Size(28, 36);
			this->ExportSpinach->TabIndex = 1;
			this->ExportSpinach->TextAlign = System::Drawing::ContentAlignment::BottomCenter;
			this->ToolTips->SetToolTip(this->ExportSpinach, L"Export to Spinach");
			this->ExportSpinach->UseVisualStyleBackColor = true;
			this->ExportSpinach->Click += gcnew System::EventHandler(this, &MainForm::ExportSpinach_Click);
			// 
			// button2
			// 
			this->button2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 4.754717F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->button2->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"button2.Image")));
			this->button2->Location = System::Drawing::Point(38, 3);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(28, 36);
			this->button2->TabIndex = 3;
			this->button2->TextAlign = System::Drawing::ContentAlignment::BottomCenter;
			this->ToolTips->SetToolTip(this->button2, L"Export to Easyspin");
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &MainForm::ExportEasyspin_Click);
			// 
			// ExportSpinEvolution
			// 
			this->ExportSpinEvolution->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 4.754717F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->ExportSpinEvolution->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"ExportSpinEvolution.Image")));
			this->ExportSpinEvolution->Location = System::Drawing::Point(38, 46);
			this->ExportSpinEvolution->Name = L"ExportSpinEvolution";
			this->ExportSpinEvolution->Size = System::Drawing::Size(28, 35);
			this->ExportSpinEvolution->TabIndex = 6;
			this->ExportSpinEvolution->TextAlign = System::Drawing::ContentAlignment::BottomCenter;
			this->ToolTips->SetToolTip(this->ExportSpinEvolution, L"Export to SpinEvolution");
			this->ExportSpinEvolution->UseVisualStyleBackColor = true;
			this->ExportSpinEvolution->Click += gcnew System::EventHandler(this, &MainForm::ExportSpinEvolution_Click);
			// 
			// BackGroundColorLabel
			// 
			this->BackGroundColorLabel->AutoSize = true;
			this->BackGroundColorLabel->BorderStyle = System::Windows::Forms::BorderStyle::FixedSingle;
			this->BackGroundColorLabel->Location = System::Drawing::Point(6, 60);
			this->BackGroundColorLabel->Name = L"BackGroundColorLabel";
			this->BackGroundColorLabel->Size = System::Drawing::Size(67, 15);
			this->BackGroundColorLabel->TabIndex = 4;
			this->BackGroundColorLabel->Text = L"Background";
			this->ToolTips->SetToolTip(this->BackGroundColorLabel, L"Click to change Background Color");
			this->BackGroundColorLabel->Click += gcnew System::EventHandler(this, &MainForm::BackGroundColorLabel_Click);
			// 
			// IsotopesButton
			// 
			this->IsotopesButton->Location = System::Drawing::Point(81, 57);
			this->IsotopesButton->Name = L"IsotopesButton";
			this->IsotopesButton->Size = System::Drawing::Size(88, 22);
			this->IsotopesButton->TabIndex = 3;
			this->IsotopesButton->Text = L"Isotopes";
			this->ToolTips->SetToolTip(this->IsotopesButton, L"Isotopes Table");
			this->IsotopesButton->UseVisualStyleBackColor = true;
			this->IsotopesButton->Click += gcnew System::EventHandler(this, &MainForm::IsotopesButton_Click);
			// 
			// BondColorlabel
			// 
			this->BondColorlabel->BackColor = System::Drawing::SystemColors::Control;
			this->BondColorlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.830189F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->BondColorlabel->ForeColor = System::Drawing::SystemColors::ControlText;
			this->BondColorlabel->Location = System::Drawing::Point(1, 1);
			this->BondColorlabel->Name = L"BondColorlabel";
			this->BondColorlabel->Size = System::Drawing::Size(191, 13);
			this->BondColorlabel->TabIndex = 3;
			this->BondColorlabel->Text = L"BOND";
			this->BondColorlabel->TextAlign = System::Drawing::ContentAlignment::TopCenter;
			this->ToolTips->SetToolTip(this->BondColorlabel, L"Change Color of Bonds");
			this->BondColorlabel->DoubleClick += gcnew System::EventHandler(this, &MainForm::Colour_DoubleClick);
			// 
			// CSColorlabel
			// 
			this->CSColorlabel->BackColor = System::Drawing::SystemColors::Control;
			this->CSColorlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.830189F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->CSColorlabel->ForeColor = System::Drawing::SystemColors::ControlText;
			this->CSColorlabel->Location = System::Drawing::Point(1, 1);
			this->CSColorlabel->Name = L"CSColorlabel";
			this->CSColorlabel->Size = System::Drawing::Size(191, 13);
			this->CSColorlabel->TabIndex = 3;
			this->CSColorlabel->Text = L"CHEMICAL SHIFT";
			this->CSColorlabel->TextAlign = System::Drawing::ContentAlignment::TopCenter;
			this->ToolTips->SetToolTip(this->CSColorlabel, L"Change Color of Chemical Shift");
			this->CSColorlabel->DoubleClick += gcnew System::EventHandler(this, &MainForm::Colour_DoubleClick);
			// 
			// HFCColorlabel
			// 
			this->HFCColorlabel->BackColor = System::Drawing::SystemColors::Control;
			this->HFCColorlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.830189F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->HFCColorlabel->ForeColor = System::Drawing::SystemColors::ControlText;
			this->HFCColorlabel->Location = System::Drawing::Point(1, 1);
			this->HFCColorlabel->Name = L"HFCColorlabel";
			this->HFCColorlabel->Size = System::Drawing::Size(191, 13);
			this->HFCColorlabel->TabIndex = 3;
			this->HFCColorlabel->Text = L"HYPERFINE COUPLING";
			this->HFCColorlabel->TextAlign = System::Drawing::ContentAlignment::TopCenter;
			this->ToolTips->SetToolTip(this->HFCColorlabel, L"Change Color of Hyperfine Couplings");
			this->HFCColorlabel->DoubleClick += gcnew System::EventHandler(this, &MainForm::Colour_DoubleClick);
			// 
			// GTColorlabel
			// 
			this->GTColorlabel->BackColor = System::Drawing::SystemColors::Control;
			this->GTColorlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.830189F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->GTColorlabel->ForeColor = System::Drawing::SystemColors::ControlText;
			this->GTColorlabel->Location = System::Drawing::Point(1, 1);
			this->GTColorlabel->Name = L"GTColorlabel";
			this->GTColorlabel->Size = System::Drawing::Size(191, 13);
			this->GTColorlabel->TabIndex = 3;
			this->GTColorlabel->Text = L"G-TENSOR";
			this->GTColorlabel->TextAlign = System::Drawing::ContentAlignment::TopCenter;
			this->ToolTips->SetToolTip(this->GTColorlabel, L"Change Color of G-Tensor");
			this->GTColorlabel->DoubleClick += gcnew System::EventHandler(this, &MainForm::Colour_DoubleClick);
			// 
			// SUSColorlabel
			// 
			this->SUSColorlabel->BackColor = System::Drawing::SystemColors::Control;
			this->SUSColorlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.830189F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->SUSColorlabel->ForeColor = System::Drawing::SystemColors::ControlText;
			this->SUSColorlabel->Location = System::Drawing::Point(1, 1);
			this->SUSColorlabel->Name = L"SUSColorlabel";
			this->SUSColorlabel->Size = System::Drawing::Size(191, 13);
			this->SUSColorlabel->TabIndex = 3;
			this->SUSColorlabel->Text = L"SUSCEPTIBILITY";
			this->SUSColorlabel->TextAlign = System::Drawing::ContentAlignment::TopCenter;
			this->ToolTips->SetToolTip(this->SUSColorlabel, L"Change Color of G-Tensor");
			this->SUSColorlabel->DoubleClick += gcnew System::EventHandler(this, &MainForm::Colour_DoubleClick);
			// 
			// JCColorlabel
			// 
			this->JCColorlabel->BackColor = System::Drawing::SystemColors::Control;
			this->JCColorlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.830189F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->JCColorlabel->ForeColor = System::Drawing::SystemColors::ControlText;
			this->JCColorlabel->Location = System::Drawing::Point(1, 1);
			this->JCColorlabel->Name = L"JCColorlabel";
			this->JCColorlabel->Size = System::Drawing::Size(191, 13);
			this->JCColorlabel->TabIndex = 3;
			this->JCColorlabel->Text = L"J-COUPLING";
			this->JCColorlabel->TextAlign = System::Drawing::ContentAlignment::TopCenter;
			this->ToolTips->SetToolTip(this->JCColorlabel, L"Change Color of J-Couplings");
			this->JCColorlabel->MouseDoubleClick += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::Colorlabel_MouseDoubleClick);
			// 
			// QCColorlabel
			// 
			this->QCColorlabel->BackColor = System::Drawing::SystemColors::Control;
			this->QCColorlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.830189F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->QCColorlabel->ForeColor = System::Drawing::SystemColors::ControlText;
			this->QCColorlabel->Location = System::Drawing::Point(1, 1);
			this->QCColorlabel->Name = L"QCColorlabel";
			this->QCColorlabel->Size = System::Drawing::Size(191, 13);
			this->QCColorlabel->TabIndex = 3;
			this->QCColorlabel->Text = L"QUADRUPOLAR COUPLING";
			this->QCColorlabel->TextAlign = System::Drawing::ContentAlignment::TopCenter;
			this->ToolTips->SetToolTip(this->QCColorlabel, L"Change Color of Hyperfine Couplings");
			this->QCColorlabel->DoubleClick += gcnew System::EventHandler(this, &MainForm::Colour_DoubleClick);
			// 
			// ZFColorlabel
			// 
			this->ZFColorlabel->BackColor = System::Drawing::SystemColors::Control;
			this->ZFColorlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.830189F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->ZFColorlabel->ForeColor = System::Drawing::SystemColors::ControlText;
			this->ZFColorlabel->Location = System::Drawing::Point(1, 1);
			this->ZFColorlabel->Name = L"ZFColorlabel";
			this->ZFColorlabel->Size = System::Drawing::Size(191, 13);
			this->ZFColorlabel->TabIndex = 3;
			this->ZFColorlabel->Text = L"ZERO-FIELD SPLITTING";
			this->ZFColorlabel->TextAlign = System::Drawing::ContentAlignment::TopCenter;
			this->ToolTips->SetToolTip(this->ZFColorlabel, L"Change Color of Zero-Field Splitting");
			this->ZFColorlabel->DoubleClick += gcnew System::EventHandler(this, &MainForm::Colour_DoubleClick);
			// 
			// ECColorlabel
			// 
			this->ECColorlabel->BackColor = System::Drawing::SystemColors::Control;
			this->ECColorlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.830189F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->ECColorlabel->ForeColor = System::Drawing::SystemColors::ControlText;
			this->ECColorlabel->Location = System::Drawing::Point(1, 1);
			this->ECColorlabel->Name = L"ECColorlabel";
			this->ECColorlabel->Size = System::Drawing::Size(191, 13);
			this->ECColorlabel->TabIndex = 3;
			this->ECColorlabel->Text = L"EXCHANGE COUPLING";
			this->ToolTips->SetToolTip(this->ECColorlabel, L"Change Color of  Exchange Couplings");
			this->ECColorlabel->MouseDoubleClick += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::Colorlabel_MouseDoubleClick);
			// 
			// AddNodeButton
			// 
			this->AddNodeButton->Font = (gcnew System::Drawing::Font(L"Arial Black", 5.25F, System::Drawing::FontStyle::Bold));
			this->AddNodeButton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"AddNodeButton.Image")));
			this->AddNodeButton->ImageAlign = System::Drawing::ContentAlignment::TopRight;
			this->AddNodeButton->Location = System::Drawing::Point(1, 1);
			this->AddNodeButton->Margin = System::Windows::Forms::Padding(1);
			this->AddNodeButton->Name = L"AddNodeButton";
			this->AddNodeButton->Size = System::Drawing::Size(23, 23);
			this->AddNodeButton->TabIndex = 0;
			this->ToolTips->SetToolTip(this->AddNodeButton, L"Add Reference Frame Child \r\nto Selected Reference Frame");
			this->AddNodeButton->UseVisualStyleBackColor = true;
			this->AddNodeButton->Click += gcnew System::EventHandler(this, &MainForm::AddNodeButton_Click);
			// 
			// DeleteNodeButton
			// 
			this->DeleteNodeButton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"DeleteNodeButton.Image")));
			this->DeleteNodeButton->Location = System::Drawing::Point(1, 27);
			this->DeleteNodeButton->Margin = System::Windows::Forms::Padding(1);
			this->DeleteNodeButton->Name = L"DeleteNodeButton";
			this->DeleteNodeButton->Size = System::Drawing::Size(23, 23);
			this->DeleteNodeButton->TabIndex = 1;
			this->ToolTips->SetToolTip(this->DeleteNodeButton, L"Delete Selected Reference Frame");
			this->DeleteNodeButton->UseVisualStyleBackColor = true;
			this->DeleteNodeButton->Click += gcnew System::EventHandler(this, &MainForm::DeleteNodeButton_Click);
			// 
			// MakeLabButton
			// 
			this->MakeLabButton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"MakeLabButton.Image")));
			this->MakeLabButton->Location = System::Drawing::Point(1, 53);
			this->MakeLabButton->Margin = System::Windows::Forms::Padding(1);
			this->MakeLabButton->Name = L"MakeLabButton";
			this->MakeLabButton->Size = System::Drawing::Size(23, 23);
			this->MakeLabButton->TabIndex = 2;
			this->ToolTips->SetToolTip(this->MakeLabButton, L"Make selected frame lab frame.");
			this->MakeLabButton->UseVisualStyleBackColor = true;
			this->MakeLabButton->Click += gcnew System::EventHandler(this, &MainForm::MakeLabButton_Click);
			// 
			// SaveFilesDialog
			// 
			this->SaveFilesDialog->RestoreDirectory = true;
			// 
			// printDialog1
			// 
			this->printDialog1->Document = this->DocumenttoPrint;
			this->printDialog1->UseEXDialog = true;
			// 
			// DocumenttoPrint
			// 
			this->DocumenttoPrint->PrintPage += gcnew System::Drawing::Printing::PrintPageEventHandler(this, &MainForm::DocumenttoPrint_PrintPage);
			// 
			// tableLayoutPanel13
			// 
			this->tableLayoutPanel13->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->tableLayoutPanel13->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel13->ColumnCount = 1;
			this->tableLayoutPanel13->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel13->Controls->Add(this->label10, 0, 0);
			this->tableLayoutPanel13->Controls->Add(this->dataGridView1, 0, 1);
			this->tableLayoutPanel13->Location = System::Drawing::Point(3, 0);
			this->tableLayoutPanel13->Name = L"tableLayoutPanel13";
			this->tableLayoutPanel13->RowCount = 2;
			this->tableLayoutPanel13->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 
				20)));
			this->tableLayoutPanel13->RowStyles->Add((gcnew System::Windows::Forms::RowStyle()));
			this->tableLayoutPanel13->Size = System::Drawing::Size(409, 627);
			this->tableLayoutPanel13->TabIndex = 5;
			this->tableLayoutPanel13->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &MainForm::dataGridView1_Paint);
			// 
			// label10
			// 
			this->label10->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label10->AutoSize = true;
			this->label10->BackColor = System::Drawing::Color::Gray;
			this->label10->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label10->Location = System::Drawing::Point(1, 1);
			this->label10->Margin = System::Windows::Forms::Padding(0);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(407, 20);
			this->label10->TabIndex = 0;
			this->label10->Text = L"Atoms";
			this->label10->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// dataGridView1
			// 
			this->dataGridView1->AllowUserToResizeColumns = false;
			this->dataGridView1->AllowUserToResizeRows = false;
			this->dataGridView1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->dataGridView1->AutoSizeColumnsMode = System::Windows::Forms::DataGridViewAutoSizeColumnsMode::Fill;
			dataGridViewCellStyle13->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleCenter;
			dataGridViewCellStyle13->BackColor = System::Drawing::SystemColors::Control;
			dataGridViewCellStyle13->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			dataGridViewCellStyle13->ForeColor = System::Drawing::SystemColors::WindowText;
			dataGridViewCellStyle13->SelectionBackColor = System::Drawing::SystemColors::Highlight;
			dataGridViewCellStyle13->SelectionForeColor = System::Drawing::SystemColors::HighlightText;
			dataGridViewCellStyle13->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->dataGridView1->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle13;
			this->dataGridView1->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView1->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(9) {this->Column1, 
				this->Label, this->Column4, this->Column5, this->Column11, this->Column6, this->Column7, this->Column8, this->Column9});
			dataGridViewCellStyle15->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleCenter;
			dataGridViewCellStyle15->BackColor = System::Drawing::SystemColors::Window;
			dataGridViewCellStyle15->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			dataGridViewCellStyle15->ForeColor = System::Drawing::SystemColors::ControlText;
			dataGridViewCellStyle15->SelectionBackColor = System::Drawing::Color::Pink;
			dataGridViewCellStyle15->SelectionForeColor = System::Drawing::SystemColors::ControlText;
			dataGridViewCellStyle15->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->dataGridView1->DefaultCellStyle = dataGridViewCellStyle15;
			this->dataGridView1->EditMode = System::Windows::Forms::DataGridViewEditMode::EditOnKeystroke;
			this->dataGridView1->Location = System::Drawing::Point(4, 25);
			this->dataGridView1->Name = L"dataGridView1";
			dataGridViewCellStyle16->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle16->BackColor = System::Drawing::SystemColors::Control;
			dataGridViewCellStyle16->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 4.471698F));
			dataGridViewCellStyle16->ForeColor = System::Drawing::SystemColors::WindowText;
			dataGridViewCellStyle16->SelectionBackColor = System::Drawing::SystemColors::Highlight;
			dataGridViewCellStyle16->SelectionForeColor = System::Drawing::SystemColors::HighlightText;
			dataGridViewCellStyle16->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->dataGridView1->RowHeadersDefaultCellStyle = dataGridViewCellStyle16;
			this->dataGridView1->RowHeadersWidth = 33;
			this->dataGridView1->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->dataGridView1->ShowCellToolTips = false;
			this->dataGridView1->Size = System::Drawing::Size(401, 600);
			this->dataGridView1->TabIndex = 0;
			this->dataGridView1->CellBeginEdit += gcnew System::Windows::Forms::DataGridViewCellCancelEventHandler(this, &MainForm::dataGridView1_CellBeginEdit);
			this->dataGridView1->CellClick += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MainForm::dataGridView1_CellClick);
			this->dataGridView1->CellEndEdit += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MainForm::dataGridView1_CellEndEdit);
			this->dataGridView1->CellLeave += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MainForm::dataGridView1_CellLeave);
			this->dataGridView1->CellMouseEnter += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MainForm::dataGridView1_CellEnter);
			this->dataGridView1->CurrentCellDirtyStateChanged += gcnew System::EventHandler(this, &MainForm::dataGridView1_CurrentCellDirtyStateChanged);
			this->dataGridView1->EditingControlShowing += gcnew System::Windows::Forms::DataGridViewEditingControlShowingEventHandler(this, &MainForm::dataGridView1_EditingControlShowing);
			this->dataGridView1->RowsAdded += gcnew System::Windows::Forms::DataGridViewRowsAddedEventHandler(this, &MainForm::dataGridView1_RowsAdded);
			this->dataGridView1->RowsRemoved += gcnew System::Windows::Forms::DataGridViewRowsRemovedEventHandler(this, &MainForm::dataGridView1_RowsRemoved);
			this->dataGridView1->SelectionChanged += gcnew System::EventHandler(this, &MainForm::dataGridView1_SelectionChanged);
			this->dataGridView1->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &MainForm::dataGridView1_Paint);
			this->dataGridView1->MouseEnter += gcnew System::EventHandler(this, &MainForm::dataGridView1_MouseEnter);
			// 
			// Column1
			// 
			this->Column1->FillWeight = 87.25761F;
			this->Column1->HeaderText = L"ID";
			this->Column1->MinimumWidth = 33;
			this->Column1->Name = L"Column1";
			this->Column1->ReadOnly = true;
			// 
			// Label
			// 
			this->Label->FillWeight = 119.2028F;
			this->Label->HeaderText = L"Label";
			this->Label->MinimumWidth = 42;
			this->Label->Name = L"Label";
			this->Label->ReadOnly = true;
			// 
			// Column4
			// 
			this->Column4->FillWeight = 108.429F;
			this->Column4->HeaderText = L"Elem.";
			this->Column4->MinimumWidth = 42;
			this->Column4->Name = L"Column4";
			this->Column4->ToolTipText = L"add atom";
			// 
			// Column5
			// 
			this->Column5->FillWeight = 111.9439F;
			this->Column5->HeaderText = L"Mass Number";
			this->Column5->MinimumWidth = 48;
			this->Column5->Name = L"Column5";
			this->Column5->ReadOnly = true;
			this->Column5->SortMode = System::Windows::Forms::DataGridViewColumnSortMode::Automatic;
			// 
			// Column11
			// 
			this->Column11->FillWeight = 104.8731F;
			this->Column11->HeaderText = L"Color";
			this->Column11->MinimumWidth = 40;
			this->Column11->Name = L"Column11";
			this->Column11->ReadOnly = true;
			this->Column11->SortMode = System::Windows::Forms::DataGridViewColumnSortMode::NotSortable;
			// 
			// Column6
			// 
			this->Column6->FillWeight = 103.7134F;
			this->Column6->HeaderText = L"x(Å)";
			this->Column6->MinimumWidth = 42;
			this->Column6->Name = L"Column6";
			this->Column6->ReadOnly = true;
			// 
			// Column7
			// 
			this->Column7->FillWeight = 102.6886F;
			this->Column7->HeaderText = L"y(Å)";
			this->Column7->MinimumWidth = 42;
			this->Column7->Name = L"Column7";
			this->Column7->ReadOnly = true;
			// 
			// Column8
			// 
			this->Column8->FillWeight = 101.7829F;
			this->Column8->HeaderText = L"z(Å)";
			this->Column8->MinimumWidth = 42;
			this->Column8->Name = L"Column8";
			this->Column8->ReadOnly = true;
			// 
			// Column9
			// 
			this->Column9->AutoSizeMode = System::Windows::Forms::DataGridViewAutoSizeColumnMode::None;
			dataGridViewCellStyle14->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleCenter;
			dataGridViewCellStyle14->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(64)), 
				static_cast<System::Int32>(static_cast<System::Byte>(64)), static_cast<System::Int32>(static_cast<System::Byte>(64)));
			dataGridViewCellStyle14->NullValue = nullptr;
			dataGridViewCellStyle14->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(64)), 
				static_cast<System::Int32>(static_cast<System::Byte>(64)), static_cast<System::Int32>(static_cast<System::Byte>(64)));
			this->Column9->DefaultCellStyle = dataGridViewCellStyle14;
			this->Column9->HeaderText = L"";
			this->Column9->ImageLayout = System::Windows::Forms::DataGridViewImageCellLayout::Zoom;
			this->Column9->MinimumWidth = 16;
			this->Column9->Name = L"Column9";
			this->Column9->Resizable = System::Windows::Forms::DataGridViewTriState::True;
			this->Column9->Width = 16;
			// 
			// splitContainer1
			// 
			this->splitContainer1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->splitContainer1->FixedPanel = System::Windows::Forms::FixedPanel::Panel1;
			this->splitContainer1->Location = System::Drawing::Point(0, 141);
			this->splitContainer1->Name = L"splitContainer1";
			// 
			// splitContainer1.Panel1
			// 
			this->splitContainer1->Panel1->Controls->Add(this->tableLayoutPanel13);
			this->splitContainer1->Panel1MinSize = 0;
			// 
			// splitContainer1.Panel2
			// 
			this->splitContainer1->Panel2->Controls->Add(this->splitContainer2);
			this->splitContainer1->Panel2MinSize = 200;
			this->splitContainer1->Size = System::Drawing::Size(1272, 627);
			this->splitContainer1->SplitterDistance = 408;
			this->splitContainer1->TabIndex = 6;
			// 
			// splitContainer2
			// 
			this->splitContainer2->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->splitContainer2->FixedPanel = System::Windows::Forms::FixedPanel::Panel2;
			this->splitContainer2->Location = System::Drawing::Point(0, 0);
			this->splitContainer2->Name = L"splitContainer2";
			// 
			// splitContainer2.Panel1
			// 
			this->splitContainer2->Panel1->Resize += gcnew System::EventHandler(this, &MainForm::splitContainer2_Panel1_Resize);
			this->splitContainer2->Panel1MinSize = 200;
			// 
			// splitContainer2.Panel2
			// 
			this->splitContainer2->Panel2->Controls->Add(this->tableLayoutPanel4);
			this->splitContainer2->Panel2MinSize = 0;
			this->splitContainer2->Size = System::Drawing::Size(860, 627);
			this->splitContainer2->SplitterDistance = 449;
			this->splitContainer2->TabIndex = 0;
			// 
			// OpenGLcontextMenu
			// 
			this->OpenGLcontextMenu->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) {this->toolStripMenuItem1});
			this->OpenGLcontextMenu->Name = L"OpenGLcontextMenu";
			this->OpenGLcontextMenu->Size = System::Drawing::Size(132, 26);
			this->OpenGLcontextMenu->Closed += gcnew System::Windows::Forms::ToolStripDropDownClosedEventHandler(this, &MainForm::OpenGLcontextMenu_Closed);
			this->OpenGLcontextMenu->ItemClicked += gcnew System::Windows::Forms::ToolStripItemClickedEventHandler(this, &MainForm::OpenGLcontextMenu_ItemClicked);
			// 
			// toolStripMenuItem1
			// 
			this->toolStripMenuItem1->Name = L"toolStripMenuItem1";
			this->toolStripMenuItem1->Size = System::Drawing::Size(131, 22);
			this->toolStripMenuItem1->Text = L"something";
			// 
			// MaintabControl
			// 
			this->MaintabControl->Anchor = static_cast<System::Windows::Forms::AnchorStyles>(((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->MaintabControl->Controls->Add(this->Home);
			this->MaintabControl->Controls->Add(this->Visualization);
			this->MaintabControl->Location = System::Drawing::Point(1, 0);
			this->MaintabControl->Name = L"MaintabControl";
			this->MaintabControl->SelectedIndex = 0;
			this->MaintabControl->Size = System::Drawing::Size(1271, 141);
			this->MaintabControl->TabIndex = 3;
			// 
			// Home
			// 
			this->Home->BackColor = System::Drawing::Color::White;
			this->Home->Controls->Add(this->tableLayoutPanel15);
			this->Home->Controls->Add(this->tableLayoutPanel1);
			this->Home->Controls->Add(this->tableLayoutPanel11);
			this->Home->Controls->Add(this->tableLayoutPanel2);
			this->Home->Controls->Add(this->mouseSelectScroll);
			this->Home->Controls->Add(this->tableLayoutPanel5);
			this->Home->Controls->Add(this->mouseSelectRight);
			this->Home->Controls->Add(this->tableLayoutPanel9);
			this->Home->Location = System::Drawing::Point(4, 22);
			this->Home->Name = L"Home";
			this->Home->Padding = System::Windows::Forms::Padding(3);
			this->Home->Size = System::Drawing::Size(1263, 115);
			this->Home->TabIndex = 0;
			this->Home->Text = L"Home";
			// 
			// tableLayoutPanel15
			// 
			this->tableLayoutPanel15->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel15->ColumnCount = 1;
			this->tableLayoutPanel15->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel15->Controls->Add(this->label5, 0, 1);
			this->tableLayoutPanel15->Controls->Add(this->tableLayoutPanel16, 0, 0);
			this->tableLayoutPanel15->Location = System::Drawing::Point(773, 0);
			this->tableLayoutPanel15->Name = L"tableLayoutPanel15";
			this->tableLayoutPanel15->RowCount = 2;
			this->tableLayoutPanel15->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 80.67633F)));
			this->tableLayoutPanel15->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 19.32367F)));
			this->tableLayoutPanel15->Size = System::Drawing::Size(347, 116);
			this->tableLayoutPanel15->TabIndex = 8;
			// 
			// label5
			// 
			this->label5->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label5->AutoSize = true;
			this->label5->BackColor = System::Drawing::Color::Gray;
			this->label5->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label5->Location = System::Drawing::Point(1, 93);
			this->label5->Margin = System::Windows::Forms::Padding(0);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(345, 22);
			this->label5->TabIndex = 0;
			this->label5->Text = L"Reference Frames";
			this->label5->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// tableLayoutPanel16
			// 
			this->tableLayoutPanel16->ColumnCount = 3;
			this->tableLayoutPanel16->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				88.19188F)));
			this->tableLayoutPanel16->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				11.80812F)));
			this->tableLayoutPanel16->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Absolute, 
				73)));
			this->tableLayoutPanel16->Controls->Add(this->RefFramestreeView, 0, 0);
			this->tableLayoutPanel16->Controls->Add(this->tableLayoutPanel17, 1, 0);
			this->tableLayoutPanel16->Controls->Add(this->ReferencecheckBox, 2, 0);
			this->tableLayoutPanel16->Location = System::Drawing::Point(4, 4);
			this->tableLayoutPanel16->Name = L"tableLayoutPanel16";
			this->tableLayoutPanel16->RowCount = 1;
			this->tableLayoutPanel16->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 50)));
			this->tableLayoutPanel16->Size = System::Drawing::Size(339, 85);
			this->tableLayoutPanel16->TabIndex = 1;
			// 
			// RefFramestreeView
			// 
			this->RefFramestreeView->HideSelection = false;
			this->RefFramestreeView->ImageIndex = 1;
			this->RefFramestreeView->ImageList = this->SpinachGUIimageList;
			this->RefFramestreeView->LabelEdit = true;
			this->RefFramestreeView->Location = System::Drawing::Point(3, 3);
			this->RefFramestreeView->Name = L"RefFramestreeView";
			treeNode2->Name = L"1";
			treeNode2->Text = L"Lab Frame [1]";
			this->RefFramestreeView->Nodes->AddRange(gcnew cli::array< System::Windows::Forms::TreeNode^  >(1) {treeNode2});
			this->RefFramestreeView->SelectedImageIndex = 0;
			this->RefFramestreeView->Size = System::Drawing::Size(228, 79);
			this->RefFramestreeView->TabIndex = 2;
			this->RefFramestreeView->BeforeLabelEdit += gcnew System::Windows::Forms::NodeLabelEditEventHandler(this, &MainForm::RefFramestreeView_BeforeLabelEdit);
			this->RefFramestreeView->AfterLabelEdit += gcnew System::Windows::Forms::NodeLabelEditEventHandler(this, &MainForm::RefFramestreeView_AfterLabelEdit);
			this->RefFramestreeView->NodeMouseDoubleClick += gcnew System::Windows::Forms::TreeNodeMouseClickEventHandler(this, &MainForm::RefFramestreeView_NodeMouseDoubleClick);
			// 
			// SpinachGUIimageList
			// 
			this->SpinachGUIimageList->ImageStream = (cli::safe_cast<System::Windows::Forms::ImageListStreamer^  >(resources->GetObject(L"SpinachGUIimageList.ImageStream")));
			this->SpinachGUIimageList->TransparentColor = System::Drawing::Color::White;
			this->SpinachGUIimageList->Images->SetKeyName(0, L"RefFrame.bmp");
			this->SpinachGUIimageList->Images->SetKeyName(1, L"RefFrameico.bmp");
			// 
			// tableLayoutPanel17
			// 
			this->tableLayoutPanel17->ColumnCount = 1;
			this->tableLayoutPanel17->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel17->Controls->Add(this->AddNodeButton, 0, 0);
			this->tableLayoutPanel17->Controls->Add(this->DeleteNodeButton, 0, 1);
			this->tableLayoutPanel17->Controls->Add(this->MakeLabButton, 0, 2);
			this->tableLayoutPanel17->Location = System::Drawing::Point(237, 3);
			this->tableLayoutPanel17->Name = L"tableLayoutPanel17";
			this->tableLayoutPanel17->RowCount = 3;
			this->tableLayoutPanel17->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 33.33333F)));
			this->tableLayoutPanel17->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 33.33333F)));
			this->tableLayoutPanel17->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 33.33334F)));
			this->tableLayoutPanel17->Size = System::Drawing::Size(25, 79);
			this->tableLayoutPanel17->TabIndex = 3;
			// 
			// ReferencecheckBox
			// 
			this->ReferencecheckBox->AutoSize = true;
			this->ReferencecheckBox->Location = System::Drawing::Point(266, 1);
			this->ReferencecheckBox->Margin = System::Windows::Forms::Padding(1);
			this->ReferencecheckBox->Name = L"ReferencecheckBox";
			this->ReferencecheckBox->Size = System::Drawing::Size(72, 82);
			this->ReferencecheckBox->TabIndex = 4;
			this->ReferencecheckBox->Text = L"\r\nShow \r\naxes of \r\nselected \r\nreference \r\nframe";
			this->ReferencecheckBox->UseVisualStyleBackColor = true;
			this->ReferencecheckBox->Click += gcnew System::EventHandler(this, &MainForm::ReferencecheckBox_Click);
			// 
			// tableLayoutPanel1
			// 
			this->tableLayoutPanel1->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel1->ColumnCount = 1;
			this->tableLayoutPanel1->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel1->Controls->Add(this->tableLayoutPanel3, 0, 0);
			this->tableLayoutPanel1->Controls->Add(this->label1, 0, 1);
			this->tableLayoutPanel1->Location = System::Drawing::Point(1, 0);
			this->tableLayoutPanel1->Name = L"tableLayoutPanel1";
			this->tableLayoutPanel1->RowCount = 2;
			this->tableLayoutPanel1->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 80.06403F)));
			this->tableLayoutPanel1->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 19.93596F)));
			this->tableLayoutPanel1->Size = System::Drawing::Size(185, 116);
			this->tableLayoutPanel1->TabIndex = 0;
			// 
			// tableLayoutPanel3
			// 
			this->tableLayoutPanel3->ColumnCount = 5;
			this->tableLayoutPanel3->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				20)));
			this->tableLayoutPanel3->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				20)));
			this->tableLayoutPanel3->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				20)));
			this->tableLayoutPanel3->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				20)));
			this->tableLayoutPanel3->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				20)));
			this->tableLayoutPanel3->Controls->Add(this->SaveAsbutton, 0, 1);
			this->tableLayoutPanel3->Controls->Add(this->OpenFileButton, 1, 0);
			this->tableLayoutPanel3->Controls->Add(this->Newbutton, 0, 0);
			this->tableLayoutPanel3->Controls->Add(this->AboutButton, 4, 0);
			this->tableLayoutPanel3->Controls->Add(this->ExportBitmap, 3, 0);
			this->tableLayoutPanel3->Controls->Add(this->Importbutton, 2, 0);
			this->tableLayoutPanel3->Controls->Add(this->Printbutton, 2, 1);
			this->tableLayoutPanel3->Controls->Add(this->Savebutton, 1, 1);
			this->tableLayoutPanel3->Controls->Add(this->ExitButton, 4, 1);
			this->tableLayoutPanel3->Controls->Add(this->Stereobutton, 3, 1);
			this->tableLayoutPanel3->Location = System::Drawing::Point(4, 4);
			this->tableLayoutPanel3->Name = L"tableLayoutPanel3";
			this->tableLayoutPanel3->RowCount = 2;
			this->tableLayoutPanel3->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 52.38095F)));
			this->tableLayoutPanel3->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 47.61905F)));
			this->tableLayoutPanel3->Size = System::Drawing::Size(177, 84);
			this->tableLayoutPanel3->TabIndex = 1;
			// 
			// Stereobutton
			// 
			this->Stereobutton->Appearance = System::Windows::Forms::Appearance::Button;
			this->Stereobutton->BackgroundImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"Stereobutton.BackgroundImage")));
			this->Stereobutton->Location = System::Drawing::Point(108, 46);
			this->Stereobutton->Name = L"Stereobutton";
			this->Stereobutton->Size = System::Drawing::Size(29, 35);
			this->Stereobutton->TabIndex = 10;
			this->Stereobutton->UseVisualStyleBackColor = true;
			this->Stereobutton->Click += gcnew System::EventHandler(this, &MainForm::Stereobutton_Click);
			// 
			// label1
			// 
			this->label1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label1->AutoSize = true;
			this->label1->BackColor = System::Drawing::Color::Gray;
			this->label1->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label1->Location = System::Drawing::Point(1, 92);
			this->label1->Margin = System::Windows::Forms::Padding(0);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(183, 23);
			this->label1->TabIndex = 0;
			this->label1->Text = L"File";
			this->label1->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// tableLayoutPanel11
			// 
			this->tableLayoutPanel11->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel11->ColumnCount = 1;
			this->tableLayoutPanel11->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel11->Controls->Add(this->label8, 0, 1);
			this->tableLayoutPanel11->Controls->Add(this->tableLayoutPanel12, 0, 0);
			this->tableLayoutPanel11->Location = System::Drawing::Point(189, 0);
			this->tableLayoutPanel11->Name = L"tableLayoutPanel11";
			this->tableLayoutPanel11->RowCount = 2;
			this->tableLayoutPanel11->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 80.06403F)));
			this->tableLayoutPanel11->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 19.93596F)));
			this->tableLayoutPanel11->Size = System::Drawing::Size(80, 116);
			this->tableLayoutPanel11->TabIndex = 7;
			// 
			// label8
			// 
			this->label8->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label8->AutoSize = true;
			this->label8->BackColor = System::Drawing::Color::Gray;
			this->label8->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label8->Location = System::Drawing::Point(1, 92);
			this->label8->Margin = System::Windows::Forms::Padding(0);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(78, 23);
			this->label8->TabIndex = 0;
			this->label8->Text = L"Export";
			this->label8->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// tableLayoutPanel12
			// 
			this->tableLayoutPanel12->ColumnCount = 2;
			this->tableLayoutPanel12->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				20)));
			this->tableLayoutPanel12->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				20)));
			this->tableLayoutPanel12->Controls->Add(this->ExportSimpson, 0, 1);
			this->tableLayoutPanel12->Controls->Add(this->ExportSpinach, 0, 0);
			this->tableLayoutPanel12->Controls->Add(this->button2, 1, 0);
			this->tableLayoutPanel12->Controls->Add(this->ExportSpinEvolution, 1, 1);
			this->tableLayoutPanel12->Location = System::Drawing::Point(4, 4);
			this->tableLayoutPanel12->Name = L"tableLayoutPanel12";
			this->tableLayoutPanel12->RowCount = 2;
			this->tableLayoutPanel12->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 52.38095F)));
			this->tableLayoutPanel12->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 47.61905F)));
			this->tableLayoutPanel12->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 
				20)));
			this->tableLayoutPanel12->Size = System::Drawing::Size(71, 84);
			this->tableLayoutPanel12->TabIndex = 1;
			// 
			// tableLayoutPanel2
			// 
			this->tableLayoutPanel2->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel2->ColumnCount = 1;
			this->tableLayoutPanel2->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel2->Controls->Add(this->label2, 0, 1);
			this->tableLayoutPanel2->Controls->Add(this->panel1, 0, 0);
			this->tableLayoutPanel2->Location = System::Drawing::Point(273, 0);
			this->tableLayoutPanel2->Name = L"tableLayoutPanel2";
			this->tableLayoutPanel2->RowCount = 2;
			this->tableLayoutPanel2->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 80.17153F)));
			this->tableLayoutPanel2->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 19.82847F)));
			this->tableLayoutPanel2->Size = System::Drawing::Size(185, 116);
			this->tableLayoutPanel2->TabIndex = 1;
			// 
			// label2
			// 
			this->label2->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label2->AutoSize = true;
			this->label2->BackColor = System::Drawing::Color::Gray;
			this->label2->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label2->Location = System::Drawing::Point(1, 92);
			this->label2->Margin = System::Windows::Forms::Padding(0);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(183, 23);
			this->label2->TabIndex = 0;
			this->label2->Text = L"3D View Properties";
			this->label2->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// panel1
			// 
			this->panel1->Controls->Add(this->FullScreencheckBox);
			this->panel1->Controls->Add(this->checkBoxAtomID);
			this->panel1->Controls->Add(this->SelectedInteractionscheckBox);
			this->panel1->Controls->Add(this->BackGroundColorLabel);
			this->panel1->Controls->Add(this->IsotopesButton);
			this->panel1->Controls->Add(this->checkBoxAxes);
			this->panel1->Location = System::Drawing::Point(4, 4);
			this->panel1->Name = L"panel1";
			this->panel1->Size = System::Drawing::Size(177, 84);
			this->panel1->TabIndex = 2;
			// 
			// FullScreencheckBox
			// 
			this->FullScreencheckBox->AutoSize = true;
			this->FullScreencheckBox->Location = System::Drawing::Point(3, 40);
			this->FullScreencheckBox->Name = L"FullScreencheckBox";
			this->FullScreencheckBox->Size = System::Drawing::Size(74, 17);
			this->FullScreencheckBox->TabIndex = 6;
			this->FullScreencheckBox->Text = L"Fullscreen";
			this->FullScreencheckBox->UseVisualStyleBackColor = true;
			this->FullScreencheckBox->Click += gcnew System::EventHandler(this, &MainForm::FullScreencheckBox_Click);
			// 
			// checkBoxAtomID
			// 
			this->checkBoxAtomID->AutoSize = true;
			this->checkBoxAtomID->Location = System::Drawing::Point(3, 21);
			this->checkBoxAtomID->Name = L"checkBoxAtomID";
			this->checkBoxAtomID->Size = System::Drawing::Size(69, 17);
			this->checkBoxAtomID->TabIndex = 5;
			this->checkBoxAtomID->Text = L"Atom IDs";
			this->checkBoxAtomID->UseVisualStyleBackColor = true;
			this->checkBoxAtomID->Click += gcnew System::EventHandler(this, &MainForm::checkBoxAtomID_Click);
			// 
			// SelectedInteractionscheckBox
			// 
			this->SelectedInteractionscheckBox->Checked = true;
			this->SelectedInteractionscheckBox->CheckState = System::Windows::Forms::CheckState::Checked;
			this->SelectedInteractionscheckBox->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 7.4F, System::Drawing::FontStyle::Regular, 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->SelectedInteractionscheckBox->Location = System::Drawing::Point(73, 3);
			this->SelectedInteractionscheckBox->Margin = System::Windows::Forms::Padding(3, 3, 3, 2);
			this->SelectedInteractionscheckBox->Name = L"SelectedInteractionscheckBox";
			this->SelectedInteractionscheckBox->Size = System::Drawing::Size(99, 45);
			this->SelectedInteractionscheckBox->TabIndex = 2;
			this->SelectedInteractionscheckBox->Text = L"Show Interactions in selected atoms";
			this->SelectedInteractionscheckBox->UseVisualStyleBackColor = true;
			this->SelectedInteractionscheckBox->Click += gcnew System::EventHandler(this, &MainForm::SelectedInteractionscheckBox_Click);
			// 
			// checkBoxAxes
			// 
			this->checkBoxAxes->AutoSize = true;
			this->checkBoxAxes->Checked = true;
			this->checkBoxAxes->CheckState = System::Windows::Forms::CheckState::Checked;
			this->checkBoxAxes->Location = System::Drawing::Point(3, 3);
			this->checkBoxAxes->Name = L"checkBoxAxes";
			this->checkBoxAxes->Size = System::Drawing::Size(49, 17);
			this->checkBoxAxes->TabIndex = 0;
			this->checkBoxAxes->Text = L"Axes";
			this->checkBoxAxes->UseVisualStyleBackColor = true;
			this->checkBoxAxes->Click += gcnew System::EventHandler(this, &MainForm::checkBoxAxes_CheckStateChanged);
			// 
			// mouseSelectScroll
			// 
			this->mouseSelectScroll->Enabled = false;
			this->mouseSelectScroll->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"mouseSelectScroll.Image")));
			this->mouseSelectScroll->InitialImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"mouseSelectScroll.InitialImage")));
			this->mouseSelectScroll->Location = System::Drawing::Point(728, 21);
			this->mouseSelectScroll->Name = L"mouseSelectScroll";
			this->mouseSelectScroll->Size = System::Drawing::Size(32, 10);
			this->mouseSelectScroll->SizeMode = System::Windows::Forms::PictureBoxSizeMode::StretchImage;
			this->mouseSelectScroll->TabIndex = 6;
			this->mouseSelectScroll->TabStop = false;
			this->mouseSelectScroll->Visible = false;
			// 
			// tableLayoutPanel5
			// 
			this->tableLayoutPanel5->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel5->ColumnCount = 1;
			this->tableLayoutPanel5->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel5->Controls->Add(this->tableLayoutPanel14, 0, 0);
			this->tableLayoutPanel5->Controls->Add(this->label6, 0, 2);
			this->tableLayoutPanel5->Location = System::Drawing::Point(461, 0);
			this->tableLayoutPanel5->Name = L"tableLayoutPanel5";
			this->tableLayoutPanel5->RowCount = 3;
			this->tableLayoutPanel5->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 65.21739F)));
			this->tableLayoutPanel5->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 13.91304F)));
			this->tableLayoutPanel5->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 19.19994F)));
			this->tableLayoutPanel5->Size = System::Drawing::Size(117, 116);
			this->tableLayoutPanel5->TabIndex = 2;
			// 
			// tableLayoutPanel14
			// 
			this->tableLayoutPanel14->ColumnCount = 3;
			this->tableLayoutPanel14->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				33.33333F)));
			this->tableLayoutPanel14->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				33.33334F)));
			this->tableLayoutPanel14->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				33.33334F)));
			this->tableLayoutPanel14->Controls->Add(this->button7, 1, 1);
			this->tableLayoutPanel14->Controls->Add(this->button5, 0, 1);
			this->tableLayoutPanel14->Controls->Add(this->YZbutton, 2, 0);
			this->tableLayoutPanel14->Controls->Add(this->XZbutton, 1, 0);
			this->tableLayoutPanel14->Controls->Add(this->XYbutton, 0, 0);
			this->tableLayoutPanel14->Controls->Add(this->button8, 2, 1);
			this->tableLayoutPanel14->Location = System::Drawing::Point(2, 1);
			this->tableLayoutPanel14->Margin = System::Windows::Forms::Padding(1, 0, 0, 0);
			this->tableLayoutPanel14->Name = L"tableLayoutPanel14";
			this->tableLayoutPanel14->RowCount = 2;
			this->tableLayoutPanel14->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 50)));
			this->tableLayoutPanel14->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 50)));
			this->tableLayoutPanel14->Size = System::Drawing::Size(114, 74);
			this->tableLayoutPanel14->TabIndex = 8;
			// 
			// button7
			// 
			this->button7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10.85F));
			this->button7->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"button7.Image")));
			this->button7->Location = System::Drawing::Point(38, 38);
			this->button7->Margin = System::Windows::Forms::Padding(1);
			this->button7->Name = L"button7";
			this->button7->Size = System::Drawing::Size(35, 35);
			this->button7->TabIndex = 14;
			this->button7->Tag = L"5";
			this->button7->Text = L"ZX";
			this->button7->UseVisualStyleBackColor = true;
			this->button7->Click += gcnew System::EventHandler(this, &MainForm::RotateToPlaneClick);
			// 
			// button5
			// 
			this->button5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10.85F));
			this->button5->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"button5.Image")));
			this->button5->Location = System::Drawing::Point(1, 38);
			this->button5->Margin = System::Windows::Forms::Padding(1);
			this->button5->Name = L"button5";
			this->button5->Size = System::Drawing::Size(35, 35);
			this->button5->TabIndex = 13;
			this->button5->Tag = L"4";
			this->button5->Text = L"YX";
			this->button5->UseVisualStyleBackColor = true;
			this->button5->Click += gcnew System::EventHandler(this, &MainForm::RotateToPlaneClick);
			// 
			// YZbutton
			// 
			this->YZbutton->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10.85F));
			this->YZbutton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"YZbutton.Image")));
			this->YZbutton->Location = System::Drawing::Point(76, 1);
			this->YZbutton->Margin = System::Windows::Forms::Padding(1);
			this->YZbutton->Name = L"YZbutton";
			this->YZbutton->Size = System::Drawing::Size(35, 35);
			this->YZbutton->TabIndex = 12;
			this->YZbutton->Tag = L"3";
			this->YZbutton->Text = L"ZY";
			this->YZbutton->UseVisualStyleBackColor = true;
			this->YZbutton->Click += gcnew System::EventHandler(this, &MainForm::RotateToPlaneClick);
			// 
			// XZbutton
			// 
			this->XZbutton->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10.85F));
			this->XZbutton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"XZbutton.Image")));
			this->XZbutton->Location = System::Drawing::Point(38, 1);
			this->XZbutton->Margin = System::Windows::Forms::Padding(1);
			this->XZbutton->Name = L"XZbutton";
			this->XZbutton->Size = System::Drawing::Size(35, 35);
			this->XZbutton->TabIndex = 11;
			this->XZbutton->Tag = L"2";
			this->XZbutton->Text = L"XZ";
			this->XZbutton->UseVisualStyleBackColor = true;
			this->XZbutton->Click += gcnew System::EventHandler(this, &MainForm::RotateToPlaneClick);
			// 
			// XYbutton
			// 
			this->XYbutton->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10.85F));
			this->XYbutton->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"XYbutton.Image")));
			this->XYbutton->Location = System::Drawing::Point(1, 1);
			this->XYbutton->Margin = System::Windows::Forms::Padding(1);
			this->XYbutton->Name = L"XYbutton";
			this->XYbutton->Size = System::Drawing::Size(35, 35);
			this->XYbutton->TabIndex = 10;
			this->XYbutton->Tag = L"1";
			this->XYbutton->Text = L"XY";
			this->XYbutton->UseVisualStyleBackColor = true;
			this->XYbutton->Click += gcnew System::EventHandler(this, &MainForm::HomeButton_Click);
			// 
			// button8
			// 
			this->button8->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10.85F));
			this->button8->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"button8.Image")));
			this->button8->Location = System::Drawing::Point(76, 38);
			this->button8->Margin = System::Windows::Forms::Padding(1);
			this->button8->Name = L"button8";
			this->button8->Size = System::Drawing::Size(35, 35);
			this->button8->TabIndex = 9;
			this->button8->Tag = L"6";
			this->button8->Text = L"YZ";
			this->button8->UseVisualStyleBackColor = true;
			this->button8->Click += gcnew System::EventHandler(this, &MainForm::RotateToPlaneClick);
			// 
			// label6
			// 
			this->label6->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label6->AutoSize = true;
			this->label6->BackColor = System::Drawing::Color::Gray;
			this->label6->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label6->Location = System::Drawing::Point(1, 92);
			this->label6->Margin = System::Windows::Forms::Padding(0);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(115, 23);
			this->label6->TabIndex = 0;
			this->label6->Text = L"View Planes";
			this->label6->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// mouseSelectRight
			// 
			this->mouseSelectRight->Enabled = false;
			this->mouseSelectRight->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"mouseSelectRight.Image")));
			this->mouseSelectRight->InitialImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"mouseSelectRight.InitialImage")));
			this->mouseSelectRight->Location = System::Drawing::Point(730, 37);
			this->mouseSelectRight->Name = L"mouseSelectRight";
			this->mouseSelectRight->Size = System::Drawing::Size(32, 10);
			this->mouseSelectRight->SizeMode = System::Windows::Forms::PictureBoxSizeMode::StretchImage;
			this->mouseSelectRight->TabIndex = 5;
			this->mouseSelectRight->TabStop = false;
			this->mouseSelectRight->Visible = false;
			// 
			// tableLayoutPanel9
			// 
			this->tableLayoutPanel9->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel9->ColumnCount = 1;
			this->tableLayoutPanel9->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel9->Controls->Add(this->tableLayoutPanel10, 0, 0);
			this->tableLayoutPanel9->Controls->Add(this->label9, 0, 1);
			this->tableLayoutPanel9->Location = System::Drawing::Point(582, 0);
			this->tableLayoutPanel9->Name = L"tableLayoutPanel9";
			this->tableLayoutPanel9->RowCount = 2;
			this->tableLayoutPanel9->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 80.06403F)));
			this->tableLayoutPanel9->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 19.93596F)));
			this->tableLayoutPanel9->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 20)));
			this->tableLayoutPanel9->Size = System::Drawing::Size(187, 116);
			this->tableLayoutPanel9->TabIndex = 3;
			// 
			// tableLayoutPanel10
			// 
			this->tableLayoutPanel10->ColumnCount = 2;
			this->tableLayoutPanel10->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				77.65363F)));
			this->tableLayoutPanel10->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				22.34637F)));
			this->tableLayoutPanel10->Controls->Add(this->ModelFileName, 0, 0);
			this->tableLayoutPanel10->Controls->Add(this->pictureBoxMouse, 1, 1);
			this->tableLayoutPanel10->Controls->Add(this->mouseSelectLeft, 1, 0);
			this->tableLayoutPanel10->Controls->Add(this->tableLayoutPanel18, 0, 1);
			this->tableLayoutPanel10->Location = System::Drawing::Point(4, 4);
			this->tableLayoutPanel10->Name = L"tableLayoutPanel10";
			this->tableLayoutPanel10->RowCount = 2;
			this->tableLayoutPanel10->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 53.42832F)));
			this->tableLayoutPanel10->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 46.57168F)));
			this->tableLayoutPanel10->Size = System::Drawing::Size(179, 84);
			this->tableLayoutPanel10->TabIndex = 2;
			// 
			// ModelFileName
			// 
			this->ModelFileName->AutoSize = true;
			this->ModelFileName->Location = System::Drawing::Point(3, 0);
			this->ModelFileName->Name = L"ModelFileName";
			this->ModelFileName->Size = System::Drawing::Size(105, 13);
			this->ModelFileName->TabIndex = 1;
			this->ModelFileName->Text = L"Imported File: (None)";
			// 
			// pictureBoxMouse
			// 
			this->pictureBoxMouse->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"pictureBoxMouse.Image")));
			this->pictureBoxMouse->InitialImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"pictureBoxMouse.InitialImage")));
			this->pictureBoxMouse->Location = System::Drawing::Point(142, 47);
			this->pictureBoxMouse->Name = L"pictureBoxMouse";
			this->pictureBoxMouse->Size = System::Drawing::Size(32, 32);
			this->pictureBoxMouse->SizeMode = System::Windows::Forms::PictureBoxSizeMode::AutoSize;
			this->pictureBoxMouse->TabIndex = 2;
			this->pictureBoxMouse->TabStop = false;
			// 
			// mouseSelectLeft
			// 
			this->mouseSelectLeft->Enabled = false;
			this->mouseSelectLeft->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"mouseSelectLeft.Image")));
			this->mouseSelectLeft->InitialImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"mouseSelectLeft.InitialImage")));
			this->mouseSelectLeft->Location = System::Drawing::Point(142, 3);
			this->mouseSelectLeft->Name = L"mouseSelectLeft";
			this->mouseSelectLeft->Size = System::Drawing::Size(32, 10);
			this->mouseSelectLeft->SizeMode = System::Windows::Forms::PictureBoxSizeMode::StretchImage;
			this->mouseSelectLeft->TabIndex = 4;
			this->mouseSelectLeft->TabStop = false;
			this->mouseSelectLeft->Visible = false;
			// 
			// tableLayoutPanel18
			// 
			this->tableLayoutPanel18->ColumnCount = 1;
			this->tableLayoutPanel18->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel18->Controls->Add(this->StatusBar, 0, 1);
			this->tableLayoutPanel18->Controls->Add(this->GenprogressBar, 0, 0);
			this->tableLayoutPanel18->Location = System::Drawing::Point(3, 44);
			this->tableLayoutPanel18->Margin = System::Windows::Forms::Padding(3, 0, 3, 0);
			this->tableLayoutPanel18->Name = L"tableLayoutPanel18";
			this->tableLayoutPanel18->RowCount = 2;
			this->tableLayoutPanel18->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 
				10)));
			this->tableLayoutPanel18->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 
				25)));
			this->tableLayoutPanel18->Size = System::Drawing::Size(133, 34);
			this->tableLayoutPanel18->TabIndex = 6;
			// 
			// StatusBar
			// 
			this->StatusBar->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) {this->StatusBarLabel});
			this->StatusBar->Location = System::Drawing::Point(0, 13);
			this->StatusBar->Name = L"StatusBar";
			this->StatusBar->Size = System::Drawing::Size(133, 22);
			this->StatusBar->TabIndex = 10;
			this->StatusBar->Text = L"status";
			// 
			// StatusBarLabel
			// 
			this->StatusBarLabel->Name = L"StatusBarLabel";
			this->StatusBarLabel->Size = System::Drawing::Size(67, 17);
			this->StatusBarLabel->Text = L"SpinachGUI";
			// 
			// GenprogressBar
			// 
			this->GenprogressBar->Location = System::Drawing::Point(3, 1);
			this->GenprogressBar->Margin = System::Windows::Forms::Padding(3, 1, 3, 0);
			this->GenprogressBar->MarqueeAnimationSpeed = 10;
			this->GenprogressBar->Name = L"GenprogressBar";
			this->GenprogressBar->Size = System::Drawing::Size(127, 9);
			this->GenprogressBar->Style = System::Windows::Forms::ProgressBarStyle::Marquee;
			this->GenprogressBar->TabIndex = 9;
			this->GenprogressBar->Visible = false;
			// 
			// label9
			// 
			this->label9->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label9->AutoSize = true;
			this->label9->BackColor = System::Drawing::Color::Gray;
			this->label9->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label9->Location = System::Drawing::Point(1, 92);
			this->label9->Margin = System::Windows::Forms::Padding(0);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(185, 23);
			this->label9->TabIndex = 0;
			this->label9->Text = L"Status";
			this->label9->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// Visualization
			// 
			this->Visualization->BackColor = System::Drawing::Color::White;
			this->Visualization->Controls->Add(this->tableLayoutPanel6);
			this->Visualization->Controls->Add(this->tableLayoutPanel8);
			this->Visualization->Location = System::Drawing::Point(4, 22);
			this->Visualization->Name = L"Visualization";
			this->Visualization->Padding = System::Windows::Forms::Padding(3);
			this->Visualization->Size = System::Drawing::Size(1263, 115);
			this->Visualization->TabIndex = 1;
			this->Visualization->Text = L"Visualization";
			// 
			// tableLayoutPanel6
			// 
			this->tableLayoutPanel6->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel6->ColumnCount = 1;
			this->tableLayoutPanel6->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel6->Controls->Add(this->label4, 0, 1);
			this->tableLayoutPanel6->Controls->Add(this->tableLayoutPanel7, 0, 0);
			this->tableLayoutPanel6->Location = System::Drawing::Point(1, 0);
			this->tableLayoutPanel6->Name = L"tableLayoutPanel6";
			this->tableLayoutPanel6->RowCount = 2;
			this->tableLayoutPanel6->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 80.70175F)));
			this->tableLayoutPanel6->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 19.29825F)));
			this->tableLayoutPanel6->Size = System::Drawing::Size(959, 115);
			this->tableLayoutPanel6->TabIndex = 1;
			// 
			// label4
			// 
			this->label4->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label4->AutoSize = true;
			this->label4->BackColor = System::Drawing::Color::Gray;
			this->label4->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label4->Location = System::Drawing::Point(1, 92);
			this->label4->Margin = System::Windows::Forms::Padding(0);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(957, 22);
			this->label4->TabIndex = 0;
			this->label4->Text = L"Tensor";
			this->label4->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// tableLayoutPanel7
			// 
			this->tableLayoutPanel7->BackColor = System::Drawing::Color::White;
			this->tableLayoutPanel7->ColumnCount = 5;
			this->tableLayoutPanel7->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 20)));
			this->tableLayoutPanel7->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 20)));
			this->tableLayoutPanel7->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 20)));
			this->tableLayoutPanel7->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 20)));
			this->tableLayoutPanel7->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 20)));
			this->tableLayoutPanel7->Controls->Add(this->SUSpanel, 4, 0);
			this->tableLayoutPanel7->Controls->Add(this->Bondpanel, 0, 0);
			this->tableLayoutPanel7->Controls->Add(this->CSpanel, 0, 1);
			this->tableLayoutPanel7->Controls->Add(this->HFCpanel, 1, 0);
			this->tableLayoutPanel7->Controls->Add(this->GTpanel, 1, 1);
			this->tableLayoutPanel7->Controls->Add(this->JCpanel, 2, 0);
			this->tableLayoutPanel7->Controls->Add(this->ECpanel, 2, 1);
			this->tableLayoutPanel7->Controls->Add(this->QCpanel, 3, 0);
			this->tableLayoutPanel7->Controls->Add(this->ZFpanel, 3, 1);
			this->tableLayoutPanel7->ForeColor = System::Drawing::Color::LightGray;
			this->tableLayoutPanel7->Location = System::Drawing::Point(1, 1);
			this->tableLayoutPanel7->Margin = System::Windows::Forms::Padding(0);
			this->tableLayoutPanel7->Name = L"tableLayoutPanel7";
			this->tableLayoutPanel7->RowCount = 2;
			this->tableLayoutPanel7->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 50)));
			this->tableLayoutPanel7->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 50)));
			this->tableLayoutPanel7->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 20)));
			this->tableLayoutPanel7->Size = System::Drawing::Size(957, 90);
			this->tableLayoutPanel7->TabIndex = 1;
			// 
			// Bondpanel
			// 
			this->Bondpanel->BackColor = System::Drawing::SystemColors::Control;
			this->Bondpanel->Controls->Add(this->BondtrackBar);
			this->Bondpanel->Controls->Add(this->BondcheckBox);
			this->Bondpanel->Controls->Add(this->BondColorlabel);
			this->Bondpanel->Enabled = false;
			this->Bondpanel->Location = System::Drawing::Point(0, 0);
			this->Bondpanel->Margin = System::Windows::Forms::Padding(0);
			this->Bondpanel->Name = L"Bondpanel";
			this->Bondpanel->Size = System::Drawing::Size(190, 44);
			this->Bondpanel->TabIndex = 5;
			// 
			// BondtrackBar
			// 












			this->BondtrackBar->AutoSize = false;
			this->BondtrackBar->Location = System::Drawing::Point(15, 18);
			this->BondtrackBar->Margin = System::Windows::Forms::Padding(2, 3, 2, 3);
			this->BondtrackBar->Maximum = 100;
			this->BondtrackBar->Name = L"BondtrackBar";
			this->BondtrackBar->Size = System::Drawing::Size(177, 23);


			this->BondtrackBar->TabIndex = 4;
			this->BondtrackBar->TickFrequency = 4;
			this->BondtrackBar->Value = 50;
			// 
			// BondcheckBox
			// 
			this->BondcheckBox->AutoSize = true;
			this->BondcheckBox->BackColor = System::Drawing::Color::Transparent;
			this->BondcheckBox->Location = System::Drawing::Point(1, 18);
			this->BondcheckBox->Name = L"BondcheckBox";
			this->BondcheckBox->Size = System::Drawing::Size(15, 14);
			this->BondcheckBox->TabIndex = 2;
			this->BondcheckBox->UseVisualStyleBackColor = false;
			this->BondcheckBox->Click += gcnew System::EventHandler(this, &MainForm::BondcheckBox_Click);
			// 
			// CSpanel
			// 
			this->CSpanel->BackColor = System::Drawing::SystemColors::Control;
			this->CSpanel->Controls->Add(this->CStrackBar);
			this->CSpanel->Controls->Add(this->CScheckBox);
			this->CSpanel->Controls->Add(this->CSColorlabel);
			this->CSpanel->Enabled = false;
			this->CSpanel->Location = System::Drawing::Point(0, 45);
			this->CSpanel->Margin = System::Windows::Forms::Padding(0);
			this->CSpanel->Name = L"CSpanel";
			this->CSpanel->Size = System::Drawing::Size(190, 44);
			this->CSpanel->TabIndex = 2;
			// 
			// CStrackBar
			// 












			this->CStrackBar->AutoSize = false;
			this->CStrackBar->Location = System::Drawing::Point(15, 18);
			this->CStrackBar->Margin = System::Windows::Forms::Padding(2, 3, 2, 3);
			this->CStrackBar->Maximum = 30;
			this->CStrackBar->Minimum = -30;
			this->CStrackBar->Name = L"CStrackBar";
			this->CStrackBar->Size = System::Drawing::Size(177, 23);


			this->CStrackBar->TabIndex = 5;
			this->CStrackBar->TickFrequency = 4;
			this->CStrackBar->Value = 0;
			// 
			// CScheckBox
			// 
			this->CScheckBox->AutoSize = true;
			this->CScheckBox->BackColor = System::Drawing::Color::Transparent;
			this->CScheckBox->Location = System::Drawing::Point(1, 18);
			this->CScheckBox->Name = L"CScheckBox";
			this->CScheckBox->Size = System::Drawing::Size(15, 14);
			this->CScheckBox->TabIndex = 2;
			this->CScheckBox->UseVisualStyleBackColor = false;
			this->CScheckBox->Click += gcnew System::EventHandler(this, &MainForm::CScheckBox_Click);
			// 
			// HFCpanel
			// 
			this->HFCpanel->BackColor = System::Drawing::SystemColors::Control;
			this->HFCpanel->Controls->Add(this->HFCtrackBar);
			this->HFCpanel->Controls->Add(this->HFCcheckBox);
			this->HFCpanel->Controls->Add(this->HFCColorlabel);
			this->HFCpanel->Enabled = false;
			this->HFCpanel->ForeColor = System::Drawing::Color::LightGray;
			this->HFCpanel->Location = System::Drawing::Point(191, 0);
			this->HFCpanel->Margin = System::Windows::Forms::Padding(0);
			this->HFCpanel->Name = L"HFCpanel";
			this->HFCpanel->Size = System::Drawing::Size(190, 44);
			this->HFCpanel->TabIndex = 0;
			// 
			// HFCtrackBar
			// 












			this->HFCtrackBar->AutoSize = false;
			this->HFCtrackBar->Location = System::Drawing::Point(15, 18);
			this->HFCtrackBar->Margin = System::Windows::Forms::Padding(2, 3, 2, 3);
			this->HFCtrackBar->Maximum = 30;
			this->HFCtrackBar->Minimum = -30;
			this->HFCtrackBar->Name = L"HFCtrackBar";
			this->HFCtrackBar->Size = System::Drawing::Size(177, 23);


			this->HFCtrackBar->TabIndex = 7;
			this->HFCtrackBar->TickFrequency = 4;
			this->HFCtrackBar->Value = 0;
			// 
			// HFCcheckBox
			// 
			this->HFCcheckBox->AutoSize = true;
			this->HFCcheckBox->BackColor = System::Drawing::Color::Transparent;
			this->HFCcheckBox->Location = System::Drawing::Point(1, 18);
			this->HFCcheckBox->Name = L"HFCcheckBox";
			this->HFCcheckBox->Size = System::Drawing::Size(15, 14);
			this->HFCcheckBox->TabIndex = 2;
			this->HFCcheckBox->UseVisualStyleBackColor = false;
			this->HFCcheckBox->Click += gcnew System::EventHandler(this, &MainForm::HFCcheckBox_Click);
			// 
			// GTpanel
			// 
			this->GTpanel->BackColor = System::Drawing::SystemColors::Control;
			this->GTpanel->Controls->Add(this->GTtrackBar);
			this->GTpanel->Controls->Add(this->GTcheckBox);
			this->GTpanel->Controls->Add(this->GTColorlabel);
			this->GTpanel->Enabled = false;
			this->GTpanel->Location = System::Drawing::Point(191, 45);
			this->GTpanel->Margin = System::Windows::Forms::Padding(0);
			this->GTpanel->Name = L"GTpanel";
			this->GTpanel->Size = System::Drawing::Size(190, 44);
			this->GTpanel->TabIndex = 10;
			// 
			// SUSpanel
			// 
			this->SUSpanel->BackColor = System::Drawing::SystemColors::Control;
			this->SUSpanel->Controls->Add(this->SUStrackBar);
			this->SUSpanel->Controls->Add(this->SUScheckBox);
			this->SUSpanel->Controls->Add(this->SUSColorlabel);
			this->SUSpanel->Enabled = false;
			this->SUSpanel->Location = System::Drawing::Point(766, 0);
			this->SUSpanel->Margin = System::Windows::Forms::Padding(0);
			this->SUSpanel->Name = L"SUSpanel";
			this->SUSpanel->Size = System::Drawing::Size(190, 44);
			this->SUSpanel->TabIndex = 10;
			// 
			// GTtrackBar
			// 












			this->GTtrackBar->AutoSize = false;
			this->GTtrackBar->Location = System::Drawing::Point(15, 18);
			this->GTtrackBar->Margin = System::Windows::Forms::Padding(2, 3, 2, 3);
			this->GTtrackBar->Maximum = 30;
			this->GTtrackBar->Minimum = -30;
			this->GTtrackBar->Name = L"GTtrackBar";
			this->GTtrackBar->Size = System::Drawing::Size(177, 23);


			this->GTtrackBar->TabIndex = 7;
			this->GTtrackBar->TickFrequency = 4;
			this->GTtrackBar->Value = 0;
			// 
			// SUStrackBar
			// 












			this->SUStrackBar->AutoSize = false;
			this->SUStrackBar->Location = System::Drawing::Point(15, 18);
			this->SUStrackBar->Margin = System::Windows::Forms::Padding(2, 3, 2, 3);
			this->SUStrackBar->Maximum = 30;
			this->SUStrackBar->Minimum = -30;
			this->SUStrackBar->Name = L"SUStrackBar";
			this->SUStrackBar->Size = System::Drawing::Size(177, 23);


			this->SUStrackBar->TabIndex = 7;
			this->SUStrackBar->TickFrequency = 4;
			this->SUStrackBar->Value = 0;
			// 
			// GTcheckBox
			// 
			this->GTcheckBox->AutoSize = true;
			this->GTcheckBox->BackColor = System::Drawing::Color::Transparent;
			this->GTcheckBox->Location = System::Drawing::Point(1, 18);
			this->GTcheckBox->Name = L"GTcheckBox";
			this->GTcheckBox->Size = System::Drawing::Size(15, 14);
			this->GTcheckBox->TabIndex = 2;
			this->GTcheckBox->UseVisualStyleBackColor = false;
			this->GTcheckBox->Click += gcnew System::EventHandler(this, &MainForm::GTcheckBox_Click);
			// 
			// SUScheckBox
			// 
			this->SUScheckBox->AutoSize = true;
			this->SUScheckBox->BackColor = System::Drawing::Color::Transparent;
			this->SUScheckBox->Location = System::Drawing::Point(1, 18);
			this->SUScheckBox->Name = L"SUScheckBox";
			this->SUScheckBox->Size = System::Drawing::Size(15, 14);
			this->SUScheckBox->TabIndex = 2;
			this->SUScheckBox->UseVisualStyleBackColor = false;
			this->SUScheckBox->Click += gcnew System::EventHandler(this, &MainForm::SUScheckBox_Click);
			// 
			// JCpanel
			// 
			this->JCpanel->BackColor = System::Drawing::SystemColors::Control;
			this->JCpanel->Controls->Add(this->JCColorlabel);
			this->JCpanel->Controls->Add(this->JCtrackBar);
			this->JCpanel->Controls->Add(this->JCcheckBox);
			this->JCpanel->Enabled = false;
			this->JCpanel->ForeColor = System::Drawing::Color::LightGray;
			this->JCpanel->Location = System::Drawing::Point(382, 0);
			this->JCpanel->Margin = System::Windows::Forms::Padding(0);
			this->JCpanel->Name = L"JCpanel";
			this->JCpanel->Size = System::Drawing::Size(190, 44);
			this->JCpanel->TabIndex = 3;
			// 
			// JCtrackBar
			// 












			this->JCtrackBar->AutoSize = false;
			this->JCtrackBar->Location = System::Drawing::Point(15, 18);
			this->JCtrackBar->Margin = System::Windows::Forms::Padding(2, 3, 2, 3);
			this->JCtrackBar->Maximum = 100;
			this->JCtrackBar->Name = L"JCtrackBar";
			this->JCtrackBar->Size = System::Drawing::Size(177, 23);


			this->JCtrackBar->TabIndex = 5;
			this->JCtrackBar->TickFrequency = 4;
			this->JCtrackBar->Value = 50;
			// 
			// JCcheckBox
			// 
			this->JCcheckBox->AutoSize = true;
			this->JCcheckBox->BackColor = System::Drawing::Color::Transparent;
			this->JCcheckBox->Location = System::Drawing::Point(1, 18);
			this->JCcheckBox->Name = L"JCcheckBox";
			this->JCcheckBox->Size = System::Drawing::Size(15, 14);
			this->JCcheckBox->TabIndex = 2;
			this->JCcheckBox->UseVisualStyleBackColor = false;
			this->JCcheckBox->Click += gcnew System::EventHandler(this, &MainForm::JCcheckBox_Click);
			// 
			// ECpanel
			// 
			this->ECpanel->BackColor = System::Drawing::SystemColors::Control;
			this->ECpanel->Controls->Add(this->ECtrackBar);
			this->ECpanel->Controls->Add(this->ECcheckBox);
			this->ECpanel->Controls->Add(this->ECColorlabel);
			this->ECpanel->Enabled = false;
			this->ECpanel->ForeColor = System::Drawing::Color::LightGray;
			this->ECpanel->Location = System::Drawing::Point(382, 45);
			this->ECpanel->Margin = System::Windows::Forms::Padding(0);
			this->ECpanel->Name = L"ECpanel";
			this->ECpanel->Size = System::Drawing::Size(190, 44);
			this->ECpanel->TabIndex = 9;
			// 
			// ECtrackBar
			// 












			this->ECtrackBar->AutoSize = false;
			this->ECtrackBar->Location = System::Drawing::Point(15, 18);
			this->ECtrackBar->Margin = System::Windows::Forms::Padding(2, 3, 2, 3);
			this->ECtrackBar->Maximum = 100;
			this->ECtrackBar->Name = L"ECtrackBar";
			this->ECtrackBar->Size = System::Drawing::Size(177, 23);


			this->ECtrackBar->TabIndex = 6;
			this->ECtrackBar->TickFrequency = 4;
			this->ECtrackBar->Value = 50;
			// 
			// ECcheckBox
			// 
			this->ECcheckBox->AutoSize = true;
			this->ECcheckBox->BackColor = System::Drawing::Color::Transparent;
			this->ECcheckBox->Location = System::Drawing::Point(1, 18);
			this->ECcheckBox->Name = L"ECcheckBox";
			this->ECcheckBox->Size = System::Drawing::Size(15, 14);
			this->ECcheckBox->TabIndex = 2;
			this->ECcheckBox->UseVisualStyleBackColor = false;
			this->ECcheckBox->Click += gcnew System::EventHandler(this, &MainForm::ECcheckBox_Click);
			// 
			// QCpanel
			// 
			this->QCpanel->BackColor = System::Drawing::SystemColors::Control;
			this->QCpanel->Controls->Add(this->QCtrackBar);
			this->QCpanel->Controls->Add(this->QCcheckBox);
			this->QCpanel->Controls->Add(this->QCColorlabel);
			this->QCpanel->Enabled = false;
			this->QCpanel->Location = System::Drawing::Point(573, 0);
			this->QCpanel->Margin = System::Windows::Forms::Padding(0);
			this->QCpanel->Name = L"QCpanel";
			this->QCpanel->Size = System::Drawing::Size(190, 44);
			this->QCpanel->TabIndex = 13;
			// 
			// QCtrackBar
			// 












			this->QCtrackBar->AutoSize = false;
			this->QCtrackBar->Location = System::Drawing::Point(15, 18);
			this->QCtrackBar->Margin = System::Windows::Forms::Padding(2, 3, 2, 3);
			this->QCtrackBar->Maximum = 30;
			this->QCtrackBar->Minimum = -30;
			this->QCtrackBar->Name = L"QCtrackBar";
			this->QCtrackBar->Size = System::Drawing::Size(177, 23);


			this->QCtrackBar->TabIndex = 7;
			this->QCtrackBar->TickFrequency = 4;
			this->QCtrackBar->Value = 0;
			// 
			// QCcheckBox
			// 
			this->QCcheckBox->AutoSize = true;
			this->QCcheckBox->BackColor = System::Drawing::Color::Transparent;
			this->QCcheckBox->Location = System::Drawing::Point(1, 18);
			this->QCcheckBox->Name = L"QCcheckBox";
			this->QCcheckBox->Size = System::Drawing::Size(15, 14);
			this->QCcheckBox->TabIndex = 2;
			this->QCcheckBox->UseVisualStyleBackColor = false;
			this->QCcheckBox->Click += gcnew System::EventHandler(this, &MainForm::QCcheckBox_Click);
			// 
			// ZFpanel
			// 
			this->ZFpanel->BackColor = System::Drawing::SystemColors::Control;
			this->ZFpanel->Controls->Add(this->ZFtrackBar);
			this->ZFpanel->Controls->Add(this->ZFcheckBox);
			this->ZFpanel->Controls->Add(this->ZFColorlabel);
			this->ZFpanel->Enabled = false;
			this->ZFpanel->Location = System::Drawing::Point(573, 45);
			this->ZFpanel->Margin = System::Windows::Forms::Padding(0);
			this->ZFpanel->Name = L"ZFpanel";
			this->ZFpanel->Size = System::Drawing::Size(190, 44);
			this->ZFpanel->TabIndex = 6;
			// 
			// ZFtrackBar
			// 












			this->ZFtrackBar->AutoSize = false;
			this->ZFtrackBar->Location = System::Drawing::Point(15, 18);
			this->ZFtrackBar->Margin = System::Windows::Forms::Padding(2, 3, 2, 3);
			this->ZFtrackBar->Maximum = 30;
			this->ZFtrackBar->Minimum = -30;
			this->ZFtrackBar->Name = L"ZFtrackBar";
			this->ZFtrackBar->Size = System::Drawing::Size(177, 23);


			this->ZFtrackBar->TabIndex = 6;
			this->ZFtrackBar->TickFrequency = 4;
			this->ZFtrackBar->Value = 0;
			// 
			// ZFcheckBox
			// 
			this->ZFcheckBox->AutoSize = true;
			this->ZFcheckBox->BackColor = System::Drawing::Color::Transparent;
			this->ZFcheckBox->Location = System::Drawing::Point(1, 18);
			this->ZFcheckBox->Name = L"ZFcheckBox";
			this->ZFcheckBox->Size = System::Drawing::Size(15, 14);
			this->ZFcheckBox->TabIndex = 2;
			this->ZFcheckBox->UseVisualStyleBackColor = false;
			this->ZFcheckBox->Click += gcnew System::EventHandler(this, &MainForm::ZFcheckBox_Click);
			// 
			// tableLayoutPanel8
			// 
			this->tableLayoutPanel8->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel8->ColumnCount = 1;
			this->tableLayoutPanel8->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent, 
				100)));
			this->tableLayoutPanel8->Controls->Add(this->panel4, 0, 0);
			this->tableLayoutPanel8->Controls->Add(this->label7, 0, 1);
			this->tableLayoutPanel8->Location = System::Drawing::Point(963, 0);
			this->tableLayoutPanel8->Name = L"tableLayoutPanel8";
			this->tableLayoutPanel8->RowCount = 2;
			this->tableLayoutPanel8->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 100)));
			this->tableLayoutPanel8->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 22)));
			this->tableLayoutPanel8->Size = System::Drawing::Size(163, 115);
			this->tableLayoutPanel8->TabIndex = 2;
			// 
			// panel4
			// 
			this->panel4->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->panel4->BackColor = System::Drawing::Color::Silver;
			this->panel4->Controls->Add(this->toolStrip2);
			this->panel4->Controls->Add(this->EigenEllipsoidscheckBox);
			this->panel4->Controls->Add(this->EigenAxescheckBox);
			this->panel4->Location = System::Drawing::Point(4, 4);
			this->panel4->Name = L"panel4";
			this->panel4->Size = System::Drawing::Size(155, 84);
			this->panel4->TabIndex = 3;
			// 
			// toolStrip2
			// 
			this->toolStrip2->GripStyle = System::Windows::Forms::ToolStripGripStyle::Hidden;
			this->toolStrip2->ImageScalingSize = System::Drawing::Size(36, 36);
			this->toolStrip2->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(3) {this->toolStripButton1, 
				this->toolStripButton2, this->toolStripButton3});
			this->toolStrip2->Location = System::Drawing::Point(0, 0);
			this->toolStrip2->Name = L"toolStrip2";
			this->toolStrip2->Size = System::Drawing::Size(155, 43);
			this->toolStrip2->Stretch = true;
			this->toolStrip2->TabIndex = 8;
			this->toolStrip2->Text = L"toolStrip2";
			// 
			// toolStripButton1
			// 
			this->toolStripButton1->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->toolStripButton1->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"toolStripButton1.Image")));
			this->toolStripButton1->ImageScaling = System::Windows::Forms::ToolStripItemImageScaling::None;
			this->toolStripButton1->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton1->Name = L"toolStripButton1";
			this->toolStripButton1->Size = System::Drawing::Size(40, 40);
			this->toolStripButton1->Tag = L"1";
			this->toolStripButton1->Text = L"NMR";
			this->toolStripButton1->Click += gcnew System::EventHandler(this, &MainForm::TensorButton_Click);
			// 
			// toolStripButton2
			// 
			this->toolStripButton2->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->toolStripButton2->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"toolStripButton2.Image")));
			this->toolStripButton2->ImageScaling = System::Windows::Forms::ToolStripItemImageScaling::None;
			this->toolStripButton2->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton2->Name = L"toolStripButton2";
			this->toolStripButton2->Size = System::Drawing::Size(40, 40);
			this->toolStripButton2->Tag = L"2";
			this->toolStripButton2->Text = L"toolStripButton2";
			this->toolStripButton2->Click += gcnew System::EventHandler(this, &MainForm::TensorButton_Click);
			// 
			// toolStripButton3
			// 
			this->toolStripButton3->DisplayStyle = System::Windows::Forms::ToolStripItemDisplayStyle::Image;
			this->toolStripButton3->Image = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"toolStripButton3.Image")));
			this->toolStripButton3->ImageScaling = System::Windows::Forms::ToolStripItemImageScaling::None;
			this->toolStripButton3->ImageTransparentColor = System::Drawing::Color::Magenta;
			this->toolStripButton3->Name = L"toolStripButton3";
			this->toolStripButton3->Size = System::Drawing::Size(40, 40);
			this->toolStripButton3->Tag = L"3";
			this->toolStripButton3->Text = L"toolStripButton3";
			this->toolStripButton3->Click += gcnew System::EventHandler(this, &MainForm::TensorButton_Click);
			// 
			// EigenEllipsoidscheckBox
			// 
			this->EigenEllipsoidscheckBox->AutoSize = true;
			this->EigenEllipsoidscheckBox->BackColor = System::Drawing::Color::Silver;
			this->EigenEllipsoidscheckBox->Checked = true;
			this->EigenEllipsoidscheckBox->CheckState = System::Windows::Forms::CheckState::Checked;
			this->EigenEllipsoidscheckBox->Location = System::Drawing::Point(3, 68);
			this->EigenEllipsoidscheckBox->Name = L"EigenEllipsoidscheckBox";
			this->EigenEllipsoidscheckBox->Size = System::Drawing::Size(158, 17);
			this->EigenEllipsoidscheckBox->TabIndex = 5;
			this->EigenEllipsoidscheckBox->Text = L"Interaction Tensor Ellipsoids";
			this->EigenEllipsoidscheckBox->UseVisualStyleBackColor = false;
			this->EigenEllipsoidscheckBox->Click += gcnew System::EventHandler(this, &MainForm::EigenEllipsoidscheckBox_Click);
			// 
			// EigenAxescheckBox
			// 
			this->EigenAxescheckBox->AutoSize = true;
			this->EigenAxescheckBox->BackColor = System::Drawing::Color::Silver;
			this->EigenAxescheckBox->Checked = true;
			this->EigenAxescheckBox->CheckState = System::Windows::Forms::CheckState::Checked;
			this->EigenAxescheckBox->Location = System::Drawing::Point(3, 50);
			this->EigenAxescheckBox->Name = L"EigenAxescheckBox";
			this->EigenAxescheckBox->Size = System::Drawing::Size(141, 17);
			this->EigenAxescheckBox->TabIndex = 0;
			this->EigenAxescheckBox->Text = L"Interaction Tensor Axes ";
			this->EigenAxescheckBox->UseVisualStyleBackColor = false;
			this->EigenAxescheckBox->Click += gcnew System::EventHandler(this, &MainForm::EigenAxescheckBox_Click);
			// 
			// label7
			// 
			this->label7->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom) 
				| System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->label7->AutoSize = true;
			this->label7->BackColor = System::Drawing::Color::Gray;
			this->label7->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label7->Location = System::Drawing::Point(1, 92);
			this->label7->Margin = System::Windows::Forms::Padding(0);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(161, 22);
			this->label7->TabIndex = 0;
			this->label7->Text = L"View";
			this->label7->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			// 
			// OpenGLtimer
			// 
			this->OpenGLtimer->Interval = 15;
			this->OpenGLtimer->Tick += gcnew System::EventHandler(this, &MainForm::OpenGLtimer_Tick);
			// 
			// ImportbackgroundWorker
			// 
			this->ImportbackgroundWorker->DoWork += gcnew System::ComponentModel::DoWorkEventHandler(this, &MainForm::ImportbackgroundWorker_DoWork);
			this->ImportbackgroundWorker->RunWorkerCompleted += gcnew System::ComponentModel::RunWorkerCompletedEventHandler(this, &MainForm::ImportbackgroundWorker_RunWorkerCompleted);
			// 
			// MainForm
			// 
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::None;
			this->BackColor = System::Drawing::SystemColors::Control;
			this->ClientSize = System::Drawing::Size(1271, 768);
			this->Controls->Add(this->MaintabControl);
			this->Controls->Add(this->splitContainer1);
			this->HelpButton = true;
			this->Icon = (cli::safe_cast<System::Drawing::Icon^  >(resources->GetObject(L"$this.Icon")));
			this->Name = L"MainForm";
			this->Text = L"SpinXML1 - Spinach GUI ";
			this->Shown += gcnew System::EventHandler(this, &MainForm::MainForm_Shown);
			this->Move += gcnew System::EventHandler(this, &MainForm::MainForm_Move);
			this->Resize += gcnew System::EventHandler(this, &MainForm::MainForm_Resize);
			this->tableLayoutPanel4->ResumeLayout(false);
			this->tableLayoutPanel4->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->dataGridView2))->EndInit();
			this->tableLayoutPanel13->ResumeLayout(false);
			this->tableLayoutPanel13->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->dataGridView1))->EndInit();
			this->splitContainer1->Panel1->ResumeLayout(false);
			this->splitContainer1->Panel2->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->splitContainer1))->EndInit();
			this->splitContainer1->ResumeLayout(false);
			this->splitContainer2->Panel2->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->splitContainer2))->EndInit();
			this->splitContainer2->ResumeLayout(false);
			this->OpenGLcontextMenu->ResumeLayout(false);
			this->MaintabControl->ResumeLayout(false);
			this->Home->ResumeLayout(false);
			this->tableLayoutPanel15->ResumeLayout(false);
			this->tableLayoutPanel15->PerformLayout();
			this->tableLayoutPanel16->ResumeLayout(false);
			this->tableLayoutPanel16->PerformLayout();
			this->tableLayoutPanel17->ResumeLayout(false);
			this->tableLayoutPanel1->ResumeLayout(false);
			this->tableLayoutPanel1->PerformLayout();
			this->tableLayoutPanel3->ResumeLayout(false);
			this->tableLayoutPanel11->ResumeLayout(false);
			this->tableLayoutPanel11->PerformLayout();
			this->tableLayoutPanel12->ResumeLayout(false);
			this->tableLayoutPanel2->ResumeLayout(false);
			this->tableLayoutPanel2->PerformLayout();
			this->panel1->ResumeLayout(false);
			this->panel1->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->mouseSelectScroll))->EndInit();
			this->tableLayoutPanel5->ResumeLayout(false);
			this->tableLayoutPanel5->PerformLayout();
			this->tableLayoutPanel14->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->mouseSelectRight))->EndInit();
			this->tableLayoutPanel9->ResumeLayout(false);
			this->tableLayoutPanel9->PerformLayout();
			this->tableLayoutPanel10->ResumeLayout(false);
			this->tableLayoutPanel10->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pictureBoxMouse))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->mouseSelectLeft))->EndInit();
			this->tableLayoutPanel18->ResumeLayout(false);
			this->tableLayoutPanel18->PerformLayout();
			this->StatusBar->ResumeLayout(false);
			this->StatusBar->PerformLayout();
			this->Visualization->ResumeLayout(false);
			this->tableLayoutPanel6->ResumeLayout(false);
			this->tableLayoutPanel6->PerformLayout();
			this->tableLayoutPanel7->ResumeLayout(false);
			this->Bondpanel->ResumeLayout(false);
			this->Bondpanel->PerformLayout();
			this->CSpanel->ResumeLayout(false);
			this->CSpanel->PerformLayout();
			this->HFCpanel->ResumeLayout(false);
			this->HFCpanel->PerformLayout();
			this->GTpanel->ResumeLayout(false);
			this->GTpanel->PerformLayout();
			this->SUSpanel->ResumeLayout(false);
			this->SUSpanel->PerformLayout();
			this->JCpanel->ResumeLayout(false);
			this->JCpanel->PerformLayout();
			this->ECpanel->ResumeLayout(false);
			this->ECpanel->PerformLayout();
			this->QCpanel->ResumeLayout(false);
			this->QCpanel->PerformLayout();
			this->ZFpanel->ResumeLayout(false);
			this->ZFpanel->PerformLayout();
			this->tableLayoutPanel8->ResumeLayout(false);
			this->tableLayoutPanel8->PerformLayout();
			this->panel4->ResumeLayout(false);
			this->panel4->PerformLayout();
			this->toolStrip2->ResumeLayout(false);
			this->toolStrip2->PerformLayout();
			this->ResumeLayout(false);

		}
#pragma endregion
private: System::Void MainForm_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) ;
private: System::Void MainForm_Resize(System::Object^  sender, System::EventArgs^  e) ;
private: System::Void OpenFileButton_Click(System::Object^  sender, System::EventArgs^  e) ;
private: System::Void dataGridView1_CellEndEdit(System::Object^  sender, System::Windows::Forms::DataGridViewCellEventArgs^  e) ;
private: System::Void dataGridView1_CellBeginEdit(System::Object^  sender, System::Windows::Forms::DataGridViewCellCancelEventArgs^  e) ;
private: System::Void dataGridView1_EditingControlShowing(System::Object^  sender, System::Windows::Forms::DataGridViewEditingControlShowingEventArgs^  e) ;
private: System::Void txtbox_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) ;
private: System::Void cmbbox_SelectedIndexChanged(System::Object^  sender, EventArgs^ e) ;
private: System::Void dataGridView1_RowsAdded(System::Object^  sender, System::Windows::Forms::DataGridViewRowsAddedEventArgs^  e);
private: System::Void HFCtrackBar_Scroll(System::Object^  sender, System::EventArgs^  e) ;
private: System::Void CStrackBar_Scroll(System::Object^  sender, System::EventArgs^  e) ;
private: System::Void BondtrackBar_Scroll(System::Object^  , System::EventArgs^  );
private: System::Void SaveAsButton_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void Colour_DoubleClick(System::Object^  sender, System::EventArgs^  e) ;
private: System::Void ExitButton_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void dataGridView1_CellClick(System::Object^  sender, System::Windows::Forms::DataGridViewCellEventArgs^  e) ;
private: System::Void checkBoxAxes_CheckStateChanged(System::Object^  sender, System::EventArgs^  e) ;
private: System::Void dataGridView2_CellClick(System::Object^  sender, System::Windows::Forms::DataGridViewCellEventArgs^  e) ;			 				
private: System::Void IsotopesButton_Click(System::Object^  sender, System::EventArgs^  e) ;
private: System::Void dataGridView1_SelectionChanged(System::Object^  sender, System::EventArgs^  e);
private: System::Void HomeButton_Click(System::Object^  , System::EventArgs^ ) ;
private: System::Void ToolTips_Popup(System::Object^ , System::Windows::Forms::PopupEventArgs^ );
private: System::Void BackGroundColorLabel_Click(System::Object^ , System::EventArgs^  );
private: System::Void CScheckBox_Click(System::Object^  , System::EventArgs^  ) ;
private: System::Void HFCcheckBox_Click(System::Object^  , System::EventArgs^  ) ;
private: System::Void dataGridView1_RowsRemoved(System::Object^  , System::Windows::Forms::DataGridViewRowsRemovedEventArgs^  );
private: System::Void BondcheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void JCcheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void JCtrackBar_Scroll(System::Object^  , System::EventArgs^  );
private: System::Void GTcheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void GTtrackBar_Scroll(System::Object^  , System::EventArgs^  );
private: System::Void SUScheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void SUStrackBar_Scroll(System::Object^  , System::EventArgs^  );
private: System::Void QCcheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void QCtrackBar_Scroll(System::Object^  , System::EventArgs^  );
private: System::Void ECcheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void ZFcheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void ECtrackBar_Scroll(System::Object^ , System::EventArgs^  );
private: System::Void ZFtrackBar_Scroll(System::Object^  , System::EventArgs^  );
private: System::Void dataGridView1_CellEnter(System::Object^  , System::Windows::Forms::DataGridViewCellEventArgs^  );
private: System::Void dataGridView1_CurrentCellDirtyStateChanged(System::Object^  , System::EventArgs^  );
private: System::Void DataGridView_Interactions_Update();
private: System::Void treeView_RefFrames_Update();
private: System::Void dataGridView2_RowsRemoved(System::Object^  , System::Windows::Forms::DataGridViewRowsRemovedEventArgs^  );
private: System::Void dataGridView2_CellBeginEdit(System::Object^  , System::Windows::Forms::DataGridViewCellCancelEventArgs^  );
private: System::Void dataGridView2_CellEnter(System::Object^  , System::Windows::Forms::DataGridViewCellEventArgs^ );
private: System::Void dataGridView2_CurrentCellDirtyStateChanged(System::Object^  , System::EventArgs^  );
private: System::Void dataGridView2_RowsAdded(System::Object^  , System::Windows::Forms::DataGridViewRowsAddedEventArgs^  );
private: System::Void dataGridView2_CellEndEdit(System::Object^  , System::Windows::Forms::DataGridViewCellEventArgs^  );
protected: System::Void LoadModel();
private: System::Void ExportEasyspin_Click(System::Object^  , System::EventArgs^ );
private: System::Void ExportBitmap_Click(System::Object^  , System::EventArgs^  ); 
private: System::Void Newbutton_Click(System::Object^  , System::EventArgs^  );
private: System::Void Importbutton_Click(System::Object^  , System::EventArgs^  ); 
private: System::Void AboutButton_Click(System::Object^  , System::EventArgs^  ); 
private: System::Void Printbutton_Click(System::Object^  , System::EventArgs^  ) ;
private: System::Void checkBoxAtomID_Click(System::Object^ , System::EventArgs^  ); 
private: System::Void Savebutton_Click(System::Object^  , System::EventArgs^  );
private: System::Void Update_OpenGL_Color(String^, Color^, bool );
private: System::Void UpdateCoordinates();
private: System::Void RefreshThresholds();
private: DataGridViewComboBoxCell^ UpdateCoupledAtomIDComboBox(int , String^ );
private: DataGridViewComboBoxCell^ UpdateInteractionTypeComboBox(int );
private: System::Void dataGridView2_SelectionChanged(System::Object^ , System::EventArgs^  ); 
private: System::Void MainForm_Shown(System::Object^  , System::EventArgs^  ); 
private: System::Void DocumenttoPrint_PrintPage(System::Object^  , System::Drawing::Printing::PrintPageEventArgs^  ) ;
private: System::Void Reset();
private: System::Void CheckExistentInteractions(bool);
private: System::Void OpenSpinXML(String ^ ); 
private: System::Void DelegateMethod(int , int );
private: System::Void Delegate2Method(int , int );
private: delegate System::Void SetColumnIndex(int i, int j);
private: System::Void OnMouseWheel(System::Object^ , System::Windows::Forms::MouseEventArgs^  );
private: System::Void TensorButton_Click(System::Object^  , System::EventArgs^  );
private: System::Void EigenAxescheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void EigenEllipsoidscheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void SelectedInteractionscheckBox_Click(System::Object^  , System::EventArgs^  );
private: System::Void dataGridView1_CellLeave(System::Object^ , System::Windows::Forms::DataGridViewCellEventArgs^  );
private: System::Void dataGridView2_MouseEnter(System::Object^  , System::EventArgs^  );
private: System::Void dataGridView1_MouseEnter(System::Object^  , System::EventArgs^  );
private: System::Void RotateToPlaneClick(System::Object^ , System::EventArgs^  ) ;
private: System::Void dataGridView1_Paint(System::Object^  , System::Windows::Forms::PaintEventArgs^  );
private: System::Void dataGridView2_Paint(System::Object^  , System::Windows::Forms::PaintEventArgs^  );
private: System::Void splitContainer2_Panel1_Resize(System::Object^  sender, System::EventArgs^  e) {
			 if(OpenGL!=nullptr)
			 {
				 OpenGL->ReSizeWindow(this,System::Drawing::Point(0,0),
					 splitContainer2->Panel1->Width,
					 splitContainer2->Panel1->Height);
			 }
		 }

private: System::Void dataGridView_Helpers(System::Void)
		 {

			 int AX=this->Location.X+
				 dataGridView1->RowHeadersWidth+dataGridView1->Columns[0]->Width+
				 dataGridView1->Columns[1]->Width+dataGridView1->Columns[2]->Width-5;
			 int AY=int(this->Location.Y+
				 this->Height-this->ClientSize.Height+
				 splitContainer1->Location.Y+
				 tableLayoutPanel13->RowStyles[0]->Height+
				 dataGridView1->ColumnHeadersHeight+
				 dataGridView1->RowTemplate->Height/1.25);
			 if(dataGridView1->Rows->Count==1 && dataGridView1->Size.Width>150) 
			 {   BalloonToolTip1->Location=Point(AX,AY);
			 BalloonToolTip1->Show();

			 }
			 else BalloonToolTip1->Hide();

			 int IX=this->Location.X+
				 splitContainer1->Panel2->Location.X+
				 splitContainer2->Panel2->Location.X+
				 dataGridView2->RowHeadersWidth+dataGridView2->Columns[0]->Width+
				 dataGridView2->Columns[1]->Width+dataGridView2->Columns[2]->Width-5;
			 int IY=int(this->Location.Y+
				 this->Height-this->ClientSize.Height+
				 splitContainer1->Location.Y+
				 tableLayoutPanel4->RowStyles[0]->Height+dataGridView2->ColumnHeadersHeight+
				 dataGridView2->Rows[0]->Height/1.25);
			 if(dataGridView2->Rows->Count==1 && dataGridView2->Size.Width>150) 
			 {
				 BalloonToolTip2->Location=Point(IX,IY);
				 BalloonToolTip2->Show();
			 }
			 else BalloonToolTip2->Hide();
		 }

private: System::Void MainForm_Move(System::Object^  sender, System::EventArgs^  e) {

			 dataGridView_Helpers();
		 }

private: System::Void OpenGLtimer_Tick(System::Object^  sender, System::EventArgs^  e) {
			 OpenGL->Travelling();
		 }

private: System::Void ModelChange(System::Object ^ sender,PropertyChangedEventArgs^ e) 
		 {

			 if(!Text->Contains("*")) this->Text+=" *";
			 Savebutton->Enabled=true;
			 if(e->PropertyName->Contains("Frame")) UpdateCoordinates();
			
		 }

private: System::Void ExportSpinEvolution_Click(System::Object^  , System::EventArgs^  );
private: System::Void ExportSpinach_Click(System::Object^  , System::EventArgs^  );
private: System::Void ExportSimpson_Click(System::Object^  , System::EventArgs^  );
private: System::Void RefFramestreeView_NodeMouseDoubleClick(System::Object^ , System::Windows::Forms::TreeNodeMouseClickEventArgs^ ); 
private: TreeNode^ FindNodeRecurcively(TreeNode ^Root, String ^Name);
private: System::Void Delegate3Method(System::Object^ , System::Windows::Forms::NodeLabelEditEventArgs^  );
private: delegate System::Void SetSenders(System::Object^  sender, System::Windows::Forms::NodeLabelEditEventArgs^  e);
private: bool labelflag;
private: bool UseTrackBarEventHandler;
private: System::Void RefFramestreeView_BeforeLabelEdit(System::Object^ , System::Windows::Forms::NodeLabelEditEventArgs^ ); 
private: System::Void RefFramestreeView_AfterLabelEdit(System::Object^  , System::Windows::Forms::NodeLabelEditEventArgs^ ) ;
private: System::Void AddNodeButton_Click(System::Object^ , System::EventArgs^  );
private: System::Void DeleteNodeButton_Click(System::Object^ , System::EventArgs^  );
private: System::Void ReferencecheckBox_Click(System::Object^  , System::EventArgs^ );
private: System::Void MakeLabButton_Click(System::Object^ , System::EventArgs^  );
private: System::Void OpenGLcontextMenu_ItemClicked(System::Object^  , System::Windows::Forms::ToolStripItemClickedEventArgs^  );
private: System::Void OpenGLcontextMenu_Closed(System::Object^  , System::Windows::Forms::ToolStripDropDownClosedEventArgs^  );
private: System::Void ImportbackgroundWorker_DoWork(System::Object^  , System::ComponentModel::DoWorkEventArgs^  );
private: System::Void ImportbackgroundWorker_RunWorkerCompleted(System::Object^  , System::ComponentModel::RunWorkerCompletedEventArgs^  );



private: System::Void Stereobutton_Click(System::Object^  sender, System::EventArgs^  e) {
			 OpenGL->stereoview=Stereobutton->Checked;
			 if(Stereobutton->Checked)
			 {
			 OpenGL->BackgroundColor[0]=OpenGL->BackgroundColor[1]=
			 OpenGL->BackgroundColor[2]=0.80;
			 }
			 else
				  {
			 OpenGL->BackgroundColor[0]=OpenGL->BackgroundColor[1]=
			 OpenGL->BackgroundColor[2]=1.0;
			 }

			 OpenGL->Render();
		 }
private: System::Void FullScreencheckBox_Click(System::Object^  sender, System::EventArgs^  e) {
			 if(FullScreencheckBox->Checked)
			 {
			 //Save  previous state
			 SplitterDistance1=splitContainer1->SplitterDistance;
			 SplitterDistance2=splitContainer2->SplitterDistance;
			 FormState=this->WindowState;


			//Change
			 this->FormBorderStyle=::FormBorderStyle::None;
             this->WindowState = FormWindowState::Maximized;
			 //MaintabControl->Visible=false;
			 //delete  MaintabControl;
			 //tableLayoutPanel13->Visible=false;
			 //tableLayoutPanel4->Visible=false;
			 splitContainer1->SplitterDistance=0;
			 splitContainer2->SplitterDistance=5000;
			 }
			 else
			 {
			 this->WindowState = FormState;
			 this->FormBorderStyle=::FormBorderStyle::Sizable;
			 splitContainer1->SplitterDistance=SplitterDistance1;
			 splitContainer2->SplitterDistance=SplitterDistance2;
			 }
		 }

private: System::Void Colorlabel_MouseDoubleClick(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) ;



};
};

