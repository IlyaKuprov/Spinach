#include "../StdAfx.h"
#include <windows.h>
#include <iostream>
#include "OpenGL.h"
#include "../MainForm.h"

/**
*\file OpenGL.cpp
*/

#pragma region Namespaces
using namespace std;
using namespace System::IO;
using namespace System::Windows::Forms;
#pragma endregion

#pragma region Definitions
//Definitions
#define MainICollection  SystemModel->InteractionCollection
#define MainCollection  SystemModel->AtomCollection
#define epsilo std::numeric_limits<double>::epsilon()
#define fun(alpha,factor,check) check?alpha*factor:(alpha-alpha*factor)/(1.0f-alpha*factor)
#pragma endregion


///// Initialization, reset and Update/////
#pragma region Initialization, reset and Update

/**
* @brief Constructor of COpenGL.
* @param parentForm The parent winodws form that the opengl will bind. 
* @param origin The origin point where the OpenGL panel will lie. 
* @param iWidth Width of OpenGL panel.
* @param iHeight Heigh of OpenGL panel.
*/
OpenGLForm::COpenGL::COpenGL(Form ^ parentForm, System::Drawing::Point origin, GLsizei iWidth, GLsizei iHeight)
{
	//Initialization for function members
    SystemModel=((SpinachGUI::MainForm^)parentForm)->SystemModel;
	mainForm=((SpinachGUI::MainForm^)parentForm);
	this->Location= origin;
	this->Height = iHeight;
	this->Width = iWidth;
	mWidth=iWidth;mHeight=iHeight;

	//ArcBall initialization
	ArcMargin=gcnew Point2fT();
	ArcMargin->X=(GLfloat)(-0.05*mWidth);// For 10% margin-no arcball, 5% each side
	ArcMargin->Y=(GLfloat)(-0.05*mHeight);
	ArcBall=gcnew ArcBallT((GLfloat)iWidth+2*ArcMargin->X,(GLfloat) iHeight+2*ArcMargin->Y);

	//Initialization of cubegrid for Metaballs
	cubeGrid=new CUBE_GRID();

	//Initalization od resolution parameters
	LowRes=false;
	focused=false;
	Spherediv=15;

	//Initialization of FPS counter
	last_update = GetTickCount()/1000.0; // setup time of day with initial value
	loop_count = 0;

	//Colors construction
	CScolor=new double[3];
	HCcolor=new double[3];
	GTcolor=new double[3];
	SUScolor=new double[3];
	CBcolor=new double[3];
	JCcolor=new double[6];
	ECcolor=new double[6];
	QCcolor=new double[3];
	ZFcolor=new double[3];
	BackgroundColor= new double[3];

	//Matrix for Rotations
	Transform=gcnew Matrix4fT();
	LastRot=gcnew Matrix3fT();
    ThisRot=gcnew Matrix3fT();
	MousePt=gcnew Point2fT();
	TransformMat=new GLfloat[16];
	SaveMat=new GLfloat[16];

	//Set defaults		
	this->Reset();
	this->Home();

	//Initialization of cubegrid for Metaballs
	this->Init_Metaballs();

	//Setup the Event handlers
	this->MouseEnter += gcnew System::EventHandler(this,&OpenGLForm::COpenGL::OnMouseEnter);
	this->MouseLeave += gcnew System::EventHandler(this,&OpenGLForm::COpenGL::OnMouseLeave);
	this->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this,&OpenGLForm::COpenGL::OnMouseMove);
	this->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this,&OpenGLForm::COpenGL::OnMouseDown);
	this->MouseUp += gcnew System::Windows::Forms::MouseEventHandler(this,&OpenGLForm::COpenGL::OnMouseUp);
	this->MouseWheel += gcnew System::Windows::Forms::MouseEventHandler(this,&OpenGLForm::COpenGL::OnMouseWheel);
	this->KeyDown+=gcnew System::Windows::Forms::KeyEventHandler(this,&OpenGLForm::COpenGL::Del_KeyDown);

	// Create as a child of the specified parent and make OpenGL compliant (no clipping)
	this->CreateParams->Style = WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_CLIPCHILDREN | WS_DISABLED;

	// Create the actual window
	this->Parent = parentForm;	
	m_hDC = GetDC((HWND)this->Handle.ToPointer());

	if(m_hDC)
		{      
		wglMakeCurrent(m_hDC, NULL);
		MySetPixelFormat(m_hDC);
		MySetFontFormat(m_hDC);
		ReSizeGLScene(iWidth, iHeight);
		InitGL();
 		}

		
}

/**
* @brief Build of Display list for fast Rendering.
*/
System::Void OpenGLForm::COpenGL:: BuildDisplayLists(System::Void)
{
	//----sphere list
  MySphere =glGenLists(15);
  glNewList(MySphere,GL_COMPILE);
  // counters for drawing the vertices in 3D space
  // calculate PI using the asin trig function
  int m = 15, p = 15;
  float pi = asin(1.0)*2;
  float rad = 1.0;
  glBegin(GL_QUADS);
  // for each of the latitude lines
  for(int a = 0; a != m; a++) {
	// for each of the longitude intersections with a latitude line
    for(int c = 0; c != p; c++) {

      /* for every set of 4 vertices, first calculate the texture coords for materials
		 then calculate the normals so that lighting is applied correctly to the shape
		 and lastly calculate where to put each vertex based on the radius of the sphere
		 the vertices are calculated clockwise as each corner of the quad */
	  // 1st vertex of the quad
	  glTexCoord2f(-a*1.0f/m, -c/(p+1.0f));
	  glNormal3f(cos(2*pi/m*a)*cos(pi/p*c-pi/2), sin(2*pi/m*a)*cos(pi/p*c-pi/2),
				 sin(pi/p*c-pi/2));
	  glVertex3f(cos(2*pi/m*a)*cos(pi/p*c-pi/2)*rad,
				 sin(2*pi/m*a)*cos(pi/p*c-pi/2)*rad,sin(pi/p*c-pi/2)*rad);
	  // 2nd vertex of the quad
	  glTexCoord2f(-(a+1.0f)/m, -c/(p+1.0f));
	  glNormal3f(cos(2*pi/m*(a+1))*cos(pi/p*c-pi/2), sin(2*pi/m*(a+1))*cos(pi/p*c-pi/2),
				 sin(pi/p*c-pi/2));
	  glVertex3f(cos(2*pi/m*(a+1))*cos(pi/p*c-pi/2)*rad, sin(2*pi/m*(a+1))*cos(pi/p*c-pi/2)*rad,
				 sin(pi/p*c-pi/2)*rad);
	  // 3rd vertex of the quad
	  glTexCoord2f(-a*1.0f/m, -(c+1.0f)/(p+1.0f));
	  glNormal3f(cos(2*pi/m*(a+1))*cos(pi/p*(c+1)-pi/2), sin(2*pi/m*(a+1))*cos(pi/p*(c+1)-pi/2),
				 sin(pi/p*(c+1)-pi/2));

	  glVertex3f(cos(2*pi/m*(a+1))*cos(pi/p*(c+1)-pi/2)*rad, sin(2*pi/m*(a+1))*cos(pi/p*(c+1)-pi/2)*rad,
				 sin(pi/p*(c+1)-pi/2)*rad);
	  // 4th and last vertex of the quad
	  glTexCoord2f(-(a+1.0f)/m, -(c+1.0f)/(p+1.0f));
	  glNormal3f(cos(2*pi/m*a)*cos(pi/p*(c+1)-pi/2), sin(2*pi/m*a)*cos(pi/p*(c+1)-pi/2),
				 sin(pi/p*(c+1)-pi/2));
	  glVertex3f(cos(2*pi/m*a)*cos(pi/p*(c+1)-pi/2)*rad, sin(2*pi/m*a)*cos(pi/p*(c+1)-pi/2)*rad,
				 sin(pi/p*(c+1)-pi/2)*rad);
	}//c
  }//a
  // end the openGL definition of the shape
  glEnd();
  glEndList();

 

  SphereHighRes=MySphere+2;
  glNewList(SphereHighRes,GL_COMPILE);
  gluQuadricDrawStyle(quadratic, GLU_FILL);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  gluSphere(quadratic,1.0,15,15);
  glEndList();

  SphereLowRes=MySphere+3;
  gluQuadricDrawStyle(quadratic, GLU_FILL);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  glNewList(SphereLowRes,GL_COMPILE);
  gluSphere(quadratic,1.0,12,12);
  glEndList();

  SelectedSphereHighRes=MySphere+4;
  glNewList(SelectedSphereHighRes,GL_COMPILE);
  gluQuadricDrawStyle(quadratic, GLU_FILL);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  gluSphere(quadratic,1.0,25,25);
  glEndList();

  SelectedSphereLowRes=MySphere+5;
  glNewList(SelectedSphereLowRes,GL_COMPILE);
  gluQuadricDrawStyle(quadratic, GLU_FILL);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  gluSphere(quadratic,1.0,22,22);
  glEndList();

  WiredSphereLowRes=MySphere+6;
  glNewList(WiredSphereLowRes,GL_COMPILE);
  gluQuadricDrawStyle(quadratic, GLU_LINE);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  gluSphere(quadratic,1.0,22,22);
  glEndList();

  WiredSphereHighRes=MySphere+7;
  glNewList(WiredSphereHighRes,GL_COMPILE);
  gluQuadricDrawStyle(quadratic, GLU_LINE);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  gluSphere(quadratic,1.0,25,25);
  glEndList();
    
  float cylindersize=0.03f;
  CylinderHighRes=MySphere+8;
  glNewList(CylinderHighRes,GL_COMPILE);
  gluQuadricDrawStyle(quadratic, GLU_FILL);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  gluCylinder(quadratic,cylindersize,cylindersize,1.0,15,15);
  glEnd();
  glEndList();

  CylinderLowRes=MySphere+9;
  glNewList(CylinderLowRes,GL_COMPILE);
  gluQuadricDrawStyle(quadratic, GLU_FILL);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  gluCylinder(quadratic,cylindersize,cylindersize,1.0,12,12);
  glEndList();

  cylindersize+=0.005f;
  FullCylinderHighRes=MySphere+10;
  glNewList(FullCylinderHighRes,GL_COMPILE);
  gluQuadricDrawStyle(quadratic, GLU_FILL);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  gluCylinder(quadratic,cylindersize,cylindersize,0.15,15,15);
  glRotatef(180, 1,0,0);
  gluDisk(quadratic,0.0,cylindersize,15,15);
  glRotatef(180, 1,0,0); 
  glTranslatef(0.0f, 0.0f, 0.15f); 
  gluDisk(quadratic,0.0,cylindersize,15,15);
  glEndList();

  FullCylinderLowRes=MySphere+11;
  glNewList(FullCylinderLowRes,GL_COMPILE);
  gluQuadricDrawStyle(quadratic, GLU_FILL);
  gluQuadricNormals(quadratic, GLU_SMOOTH );
  gluCylinder(quadratic,cylindersize,cylindersize,0.15,12,12);
  glRotatef(180, 1,0,0);
  gluDisk(quadratic,0.0,cylindersize,12,12);
  glRotatef(180, 1,0,0); 
  glTranslatef(0.0f, 0.0f, 0.15f); 
  gluDisk(quadratic,0.0,cylindersize,12,12);
  glEndList();


  Navigation=MySphere+12;
  glNewList(Navigation,GL_COMPILE);
  Render_Navigation();
  glEndList();

}

/**
* @brief  Destructor of COpenGL.
*/
OpenGLForm::COpenGL::~COpenGL(System::Void)
{
	if(wglGetCurrentContext()==m_hglrc) wglMakeCurrent(NULL, NULL);
	if(m_hglrc) wglDeleteContext(m_hglrc);
	if(m_hDC) ReleaseDC((HWND)this->Handle.ToPointer(), m_hDC);
	gluDeleteQuadric(quadratic);
	delete cubeGrid;
	delete [] CScolor; delete [] HCcolor; delete [] GTcolor; delete [] SUScolor;
	delete [] CBcolor; delete [] JCcolor; delete [] ECcolor; delete [] QCcolor; delete [] ZFcolor;
	delete [] TransformMat; delete [] SaveMat; delete [] BackgroundColor;
	this->DestroyHandle();
}

/**
* @brief Resetting viewing parameters, colors, flags.
*/
System::Void OpenGLForm::COpenGL::Reset()
{	
	//Default Background, view and select
	BackgroundColor[0]=1.0f;BackgroundColor[1]=1.0f;BackgroundColor[2]=1.0f;
	AxesFlag=true;
	AtomIDFlag=false;
	SelectFlag=true;
	EigenAxesFlag=true;
	EigenEllipsoidsFlag=true;
	SelectedAtomsInteractionsFlag=true;
	stereoview=false;
	ReferenceFrameAxesFlag=false;

	//Default trackbars flags and values
	CSFlag=true;
	HCFlag=true;
	GTFlag=true;
	SUSFlag=true;
	CBFlag=true;
	JCFlag=true;
	DCFlag=true;
	QCFlag=true;
	ECFlag=true;
	ZFFlag=true;
	//
	CSfactor=1.0;
	HCfactor=1.0;
	GTfactor=1.0;
	SUSfactor=1.0;
	CBfactor=0.5f;
	CBfinfactor=1.6f;
	JCfactor=0.5;
	QCfactor=1.0;
	ECfactor=0.5;
	ZFfactor=1.0;

	//MinX=MinY=MinZ=100000;
	//MaxX=MaxY=MaxZ=-100000;

	//Default Minima and maxima trackbar
//	MinCS=MinAC=MinGT=MinCB=MinJC=MinDC=MinQC=MinEC=MinZF=100000;
//	MaxCS=MaxAC=MaxGT=MaxCB=MaxJC=MaxDC=MaxQC=MaxEC=MaxZF=-100000;

	//Default colors
	CScolor[0]=0.0f;CScolor[1]=0.0f;CScolor[2]=0.9f;  //Blue
	HCcolor[0]=1.0f;HCcolor[1]=0.5f;HCcolor[2]=0.25f; //Oramge
	GTcolor[0]=1.0f;GTcolor[1]=0.0f;GTcolor[2]=0.0f;  //Red
	SUScolor[0]=0.0f;SUScolor[1]=0.7f;SUScolor[2]=0.8f; //Cyan: magnetic susceptibility
	CBcolor[0]=0.0f;CBcolor[1]=0.0f;CBcolor[2]=0.0f;  //Black
	QCcolor[0]=0.25f;QCcolor[1]=0.0f;QCcolor[2]=0.25f;//Purple 
	ZFcolor[0]=0.0f;ZFcolor[1]=0.5f;ZFcolor[2]=0.25f; //Green

	ECcolor[0]=1.0f;ECcolor[1]=0.0f;ECcolor[2]=0.0f; //Red to--
	ECcolor[3]=0.0f;ECcolor[4]=0.5f;ECcolor[5]=0.25f;//Green
	JCcolor[0]=0.0f;JCcolor[1]=0.0f;JCcolor[2]=0.9f; //Blue to--
	JCcolor[3]=1.0f;JCcolor[4]=0.0f;JCcolor[5]=0.0f; //Red


}

/**
* @brief Setting 3D system at original rotation, origin, zoom.
*/
System::Void OpenGLForm::COpenGL::Home()
{
	//Set the transformation matrix to origin
	Transform->XX=1.0f;Transform->XY=0.0f;Transform->XZ=0.0f;Transform->XW=0.0f;
	Transform->YX=0.0f;Transform->YY=1.0f;Transform->YZ=0.0f;Transform->YW=0.0f;
	Transform->ZX=0.0f;Transform->ZY=0.0f;Transform->ZZ=1.0f;Transform->ZW=0.0f;
	Transform->TX=0.0f;Transform->TY=0.0f;Transform->TZ=0.0f;Transform->TW=1.0f;
	ConMat();
	Matrix3fSetIdentity(LastRot);								// Reset Rotation
	Matrix3fSetIdentity(ThisRot);								// Reset Rotation
    ThisRot = LastRot;
	Zoom=-30.0;
	TranslX=TranslY=0.0f;

	//Not selecting any thing
	SelectedAtom=-1;
	ToBeSelectedAtom=-1;

	
}

/**
* @brief Setting the initial zoom of the system in order to fit perfctly in screen.
*/
System::Void OpenGLForm::COpenGL::ZoomtoFillScreen()
{
	double ratio=(double)System::Math::Max(1,mHeight)/System::Math::Max(1,mWidth);

	//Fit both projected dimensions; only the horizontal extent needs the aspect ratio.
	double extent=System::Math::Max(ratio*(MainCollection->MaxX-MainCollection->MinX),
	                              MainCollection->MaxY-MainCollection->MinY);
	Zoom=extent>epsilo ? -extent/(2.0*tan(PI/8.0))-5.0 : -30.0;
	if(MainCollection->MaxZ>MainCollection->MinZ)
		Zoom-=0.5*(MainCollection->MaxZ-MainCollection->MinZ);
}

/**
* @brief Converts Matrix4fT form to C++ native matrix form for OpenGL use 
*/
System::Void OpenGLForm::COpenGL::ConMat(System::Void)
{
	TransformMat[0]=Transform->XX;TransformMat[1]=Transform->XY;
	TransformMat[2]=Transform->XZ;TransformMat[3]=Transform->XW;
	TransformMat[4]=Transform->YX;TransformMat[5]=Transform->YY;
	TransformMat[6]=Transform->YZ;TransformMat[7]=Transform->YW;
	TransformMat[8]=Transform->ZX;TransformMat[9]=Transform->ZY;
	TransformMat[10]=Transform->ZZ;TransformMat[11]=Transform->ZW;
	TransformMat[12]=Transform->TX;TransformMat[13]=Transform->TY;
	TransformMat[14]=Transform->TZ;TransformMat[15]=Transform->TW;
	
}

#pragma endregion

///// Travelling to Plane/////////////////
#pragma region Travelling to Plane
/**
* @brief Rotates the 3D model to a plane.
* @param axis1 Name of the first axis of the plane to rotate.
* @param axis2 Name of the second axis of the plane to rotate.
*/
System::Void OpenGLForm::COpenGL::RotateToPlane(char axis1, char axis2)
{
	//Make identity
	 Transform->XX=1.0; Transform->XY=Transform->XZ=Transform->XW=0.0;
	 Transform->YX=0.0; Transform->YY=1.0; Transform->YZ=Transform->YW=0.0; 
	 Transform->ZX=Transform->ZY=0.0; Transform->ZZ=1.0; Transform->ZW=0.0; 
	 Transform->TX=Transform->TY=Transform->TZ=0.0; Transform->TW=1.0;
	  Matrix3fSetIdentity(ThisRot);

	
	if(axis1==axis2) return;

	//Set up the rotattion matrix according to selected plane
	if(axis1=='X' && axis2=='Z')
	{
	Transform->YY=ThisRot->YY=0.0;Transform->ZZ=ThisRot->ZZ=0.0;
	Transform->YZ=ThisRot->YZ=-1.0;Transform->ZY=ThisRot->ZY=1.0;
	};

	if(axis1=='Z' && axis2=='Y')
	{
	Transform->XX=ThisRot->XX=0.0;Transform->ZZ=ThisRot->ZZ=0.0;
	Transform->XZ=ThisRot->XZ=-1.0;Transform->ZX=ThisRot->ZX=1.0;
	};

	if(axis1=='Y' && axis2=='X')
	{
	Transform->XX=ThisRot->XX=0.0;Transform->YY=ThisRot->YY=0.0;Transform->ZZ=ThisRot->ZZ=0.0;
	Transform->XY=ThisRot->XY=1.0;Transform->YX=ThisRot->YX=1.0;Transform->ZZ=ThisRot->ZZ=-1.0;
	};

	if(axis1=='Z' && axis2=='X')
	{
	Transform->XX=ThisRot->XX=0.0;Transform->YY=ThisRot->YY=0.0;Transform->ZZ=ThisRot->ZZ=0.0;
	Transform->XY=ThisRot->XY=1.0;Transform->YZ=ThisRot->YZ=1.0;Transform->ZX=ThisRot->ZX=1.0;
	};

	if(axis1=='Y' && axis2=='Z')
	{
	Transform->XX=ThisRot->XX=0.0;Transform->YY=ThisRot->YY=0.0;Transform->ZZ=ThisRot->ZZ=0.0;
	Transform->XZ=ThisRot->XZ=1.0;Transform->YX=ThisRot->YX=1.0;Transform->ZY=ThisRot->ZY=1.0;
	};
	
    //Extract the current rotation matrix
	GLdouble modelview[16]; //var to hold the modelview info
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview ); //get the modelview info

	//Migrate the final matrix
	GLdouble temp[16];
	temp[0]=Transform->XX;temp[1]=Transform->XY;temp[2]=Transform->XZ;temp[3]=Transform->XW;
	temp[4]=Transform->YX;temp[5]=Transform->YY;temp[6]=Transform->YZ;temp[7]=Transform->YW;
	temp[8]=Transform->ZX;temp[9]=Transform->ZY;temp[10]=Transform->ZZ;temp[11]=Transform->ZW;
	temp[12]=Transform->TX;temp[13]=Transform->TY;temp[14]=Transform->TZ;temp[15]=Transform->TW;

	//Create the incrementing step matrix
	double steps=15.0f;
	for(int i=0;i<12;i++)  SaveMat[i]=(temp[i]-modelview[i])/steps;

	//Start the Slow travelling
	((SpinachGUI::MainForm^)mainForm)->OpenGLtimer->Start();

}

/**
* @brief Creates the transformation matrix for every frame during travelling
*/
System::Void OpenGLForm::COpenGL::Travelling()
{

	//Extract the current rotation matrix
	GLdouble modelview[16]; //var to hold the modelview info
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview ); //get the modelview info

	//Migrate the final matrix
	GLdouble temp[16];
	temp[0]=Transform->XX;temp[1]=Transform->XY;temp[2]=Transform->XZ;temp[3]=Transform->XW;
	temp[4]=Transform->YX;temp[5]=Transform->YY;temp[6]=Transform->YZ;temp[7]=Transform->YW;
	temp[8]=Transform->ZX;temp[9]=Transform->ZY;temp[10]=Transform->ZZ;temp[11]=Transform->ZW;
	temp[12]=Transform->TX;temp[13]=Transform->TY;temp[14]=Transform->TZ;temp[15]=Transform->TW;

	//Create the transformation matrix for the current frame of travelling
	// and calculate the difference from the final transformation
	double sum=0.0f;
	for(int i=0;i<12;i++)
	{
		TransformMat[i]=modelview[i]+SaveMat[i];
		sum+=abs(temp[i]-TransformMat[i]);
	}

	// Travel termination critirion
	if(sum<0.001) 	((SpinachGUI::MainForm^)mainForm)->OpenGLtimer->Stop();
	Render();
}
#pragma endregion

///// Total Rendering Routines////////////
#pragma region Total Rendering Routines
/**
* @brief Updates the display list of the components of the system and renders at the end
* @param type The components type to update. 
*/
System::Void OpenGLForm::COpenGL::UpdateDisplay(UpdateType type)
{
	//Updates only atoms
	if(type==UpdateType::ATOMS || type==UpdateType::ALL)
	{
		AllAtoms=MySphere+13;
		glNewList(AllAtoms,GL_COMPILE);
		Render_Atoms();
		glEndList();
	}
	
	//Updates only bonds and j-couplings
	if(type==UpdateType::BONDS_LIKE || type==UpdateType::ALL)
	{
		AllBonds_Like=MySphere+14;
		glNewList(AllBonds_Like,GL_COMPILE);
		Render_Bonds();
		glEndList();
	}

	//Updates only ellpsoids
	if(type==UpdateType::ELLIPSOIDS || type==UpdateType::ALL)
	{
		All_Ellipsoids=MySphere+15;
		glNewList(All_Ellipsoids,GL_COMPILE);
		//Using second method for rendering tranlucent ellipsoids
		Render_Transparent_Atom_Properties(2);
		Render_Susceptibility();
		glEndList();
		
	}

	Render();
	
}

/**
* @brief Rendering the all of the components of the 3d model
*/
System::Void OpenGLForm::COpenGL::TotalRender(System::Void)
{
	//Render "Nuclei" name
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, (GLfloat)mWidth, 0.0,(GLfloat)mHeight);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	
	glListBase(nFontList);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	char buf[23];
	glColor3f(0.7f, 0.7f, 0.7f);
	glRasterPos2i((GLint)10, (GLint)mHeight-25);

	//Measuring FPS in debug mode
#if defined(_DEBUG)
	sprintf_s(buf,"NUCLEI FPS: %d    ", (int)FPSrate);
	glCallLists (17, GL_UNSIGNED_BYTE, &buf);
#else
	sprintf_s(buf,"NUCLEI");
	glCallLists (6, GL_UNSIGNED_BYTE, &buf);
#endif

	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);

	glMatrixMode(GL_PROJECTION); //swapped this with...
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW); //...this
	glPopMatrix();

	//Render Electron cloud Plane
	if(MainCollection->get_electron_number()>0) Render_Electrons();
	
	
	//Zoom & Translation, Rotation and Centering
	glTranslatef((GLfloat)(TranslX/abs(Zoom)), (GLfloat)(TranslY/abs(Zoom)),(GLfloat) Zoom);
	glMultMatrixf(TransformMat);
	glTranslatef((GLfloat)(-(MainCollection->MaxX+MainCollection->MinX)/2.0f),
		(GLfloat)(-(MainCollection->MaxY+MainCollection->MinY)/2.0f),
		(GLfloat)(-(MainCollection->MaxZ+MainCollection->MinZ)/2.0f));

  		
	//Render Navigation Axes
	if (AxesFlag)  Render_Navigation();
	//glCallList(Navigation);

	//Render Atoms and selected Atoms
	glCallList(AllAtoms);
	Render_SelectedAtoms();
	
	//Render Bonds
	glCallList(AllBonds_Like);

	//Render translucent Atom properties
	glCallList(All_Ellipsoids);
}

/**
* @brief Final Render routine
*/
System::Void OpenGLForm::COpenGL::Render(System::Void)
{
	//For measuring frame rate
	++loop_count;
    double curtime = GetTickCount()/1000.0;

	//Initializing Rendering
	// Background color
	glClearColor((GLfloat)BackgroundColor[0],  
		         (GLfloat)BackgroundColor[1], 
		         (GLfloat)BackgroundColor[2], 
				 1.0f); 
	// Fog-culling color
	GLfloat fogcol[4];
	fogcol[0] = (GLfloat)BackgroundColor[0];  
	fogcol[1] = (GLfloat)BackgroundColor[0];
	fogcol[2] = (GLfloat)BackgroundColor[0];  
	fogcol[3] = 1.0;
	glFogfv(GL_FOG_COLOR, fogcol);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);	// Clear screen and depth buffer
	
	//Setting the perspective
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f,(GLfloat)mWidth/(GLfloat)mHeight,0.1f,300.0f); //gluPerspective(fov, Aspect, nplane, fplane);

	//Back to modelview for rotation
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//Quadbuffering for stereo view or not
	if(stereoview)
	{     
	//Left Eye
    glDrawBuffer(GL_BACK_LEFT);
	Calc_eye(1, true);
	TotalRender();

	//Right Eye
	glDepthMask(GL_TRUE);
	glClear(GL_DEPTH_BUFFER_BIT);
	glDrawBuffer(GL_BACK_RIGHT);
	Calc_eye(2, true);
	TotalRender();
	}
	else 
	{
	Calc_eye(0, true);
	TotalRender();
	}
	
	//For FPS measure
	if (curtime - last_update > 0.5)
	{
		FPSrate = loop_count / (curtime - last_update);
		last_update = curtime;
		loop_count = 0;
	}

	//OpenGL exception catching
	GLenum x = glGetError() ;
	if ( x != GL_NO_ERROR ) cout << "OpenGL error:" << gluErrorString(x) << endl;
		
	
	//Swap the Buffers and send it to the graphics card(flush)
	SwapBuffers(m_hDC);	
	glFlush();
	glDrawBuffer(GL_BACK);
}

/**
* @brief Calculates the transformation matrix and parameters for double or quad buffer display
* @param whicheye For which eye to create the transformation matrix, 1:Left, 2:Right, 0: Nothing.
* @param calc Flag used for rendering of selection buffer.
*/
System::Void OpenGLForm::COpenGL::Calc_eye(int whicheye, bool calc)
{
  float Aspect=(GLfloat)System::Math::Max(1,mWidth)/System::Math::Max(1,mHeight);
  float d,cpUp, cpDown, cpLeft, cpRight;
  float nearClip=1.0f, farClip=6000.0f, vSize=2.0f, zDist=0.0f;
  //float nearClip=0.5f, farClip=10.0f, vSize=4.0f, zDist=0.0f;
  float halfvsize = 0.5f * vSize; //OK
  float halfhsize = Aspect * halfvsize;//OK // width = aspect * height
  float eyePos[3]={0,0,2}; //OK
  float upDir[3]={0,1,0}; //OK
  float eyeDir[3]={0,0,-2}; //OK
  float eyeSep=0.080f; //OK
  float eyeDist=eyePos[2]; //OK
  float flen=3.0; //3
  float eyeSepDir[3];

  if (!eyeDist) return;
  for (int i=0; i<3; i++) eyeDir[i] *= flen/eyeDist;
  eyeDist = flen;

  if(eyePos[2] - zDist != 0.0f)
    {
	 d = nearClip / (eyePos[2] - zDist);
	 cpRight = d * halfhsize;     // right side is at half width
     cpLeft = -cpRight;           // left side is at negative half width
     cpUp = d * halfvsize;        // top side is at half height
     cpDown = -cpUp;              // bottom is at negative half height
     };

  float *L = eyeDir;
  float *U = upDir;
  float m, A = 0.5f * eyeSep;
  eyeSepDir[0] = L[1] * U[2] - L[2] * U[1];
  eyeSepDir[1] = L[2] * U[0] - L[0] * U[2];
  eyeSepDir[2] = L[0] * U[1] - L[1] * U[0];
  m = sqrt(eyeSepDir[0] * eyeSepDir[0] + eyeSepDir[1] * eyeSepDir[1] +
            eyeSepDir[2] * eyeSepDir[2]);

  if(m > 0.0)  A /= m;
  else  A = 0.0;

  eyeSepDir[0] *= A;
  eyeSepDir[1] *= A;
  eyeSepDir[2] *= A;

  int i;
  GLdouble ep[3];

  

  switch (whicheye) {
  case 1:  //Left eye
	  // draw to BACK LEFT only
	  // glDepthMask(GL_TRUE);
	  glDrawBuffer(GL_BACK_LEFT);
	  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  for (i=0; i<3; i++)  ep[i] = eyePos[i] - eyeSepDir[i];
	  break;

  case 2: //Right eye
	  // draw to BACK RIGHT only
	  if(calc)
	  {
	  glDepthMask(GL_TRUE);
	  //glClear(GL_DEPTH_BUFFER_BIT); //normal

	  glDrawBuffer(GL_BACK_RIGHT);//alternative
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);//alternative

	  glDrawBuffer(GL_BACK_RIGHT);
	  }
	  for (i=0; i<3; i++) ep[i] = eyePos[i] + eyeSepDir[i];
	  break;

  case 0://No stereo
	  for (i=0; i<3; i++) ep[i] = eyePos[i];
	  break;
  };

 
  float eyeshift = float(ep[0] - eyePos[0]);
  if(calc) glMatrixMode(GL_PROJECTION);
  if(calc) glLoadIdentity();
  glFrustum((GLdouble)cpLeft  + (eyeshift * nearClip / eyeDist),
	  (GLdouble)cpRight + (eyeshift * nearClip / eyeDist),
	  (GLdouble)cpDown,
	  (GLdouble)cpUp,
	  (GLdouble)nearClip,
	  (GLdouble)farClip);


  glTranslatef(-eyeshift, 0.0, 0.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(eyePos[0], eyePos[1], eyePos[2],
	  (GLdouble)(eyePos[0] + eyeDir[0]),
	  (GLdouble)(eyePos[1] + eyeDir[1]),
	  (GLdouble)(eyePos[2] + eyeDir[2]),
	  (GLdouble)(upDir[0]),
	  (GLdouble)(upDir[1]),
	  (GLdouble)(upDir[2]));
/*

      glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum((GLdouble)cpLeft,  (GLdouble)cpRight,  (GLdouble)cpDown,
              (GLdouble)cpUp,   (GLdouble)nearClip, (GLdouble)farClip);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(ep[0], ep[1], ep[2],
              (GLdouble)(eyePos[0] + eyeDir[0]),
              (GLdouble)(eyePos[1] + eyeDir[1]),
              (GLdouble)(eyePos[2] + eyeDir[2]),
              (GLdouble)(upDir[0]),
              (GLdouble)(upDir[1]),
              (GLdouble)(upDir[2]));

*/

}
#pragma endregion

///// Rendering Nuclei components/////////
#pragma region Rendering Nuclei components
/**
* @brief Render in total the transparent ellispoids of all atoms nuclei selecting render method.
* @param RenderKind Rendering method for translucent ellipsoids 
* 1:Nothing just blending, 2:No depthing, 3:Using multiple passes.
*/
System::Void OpenGLForm::COpenGL::Render_Transparent_Atom_Properties(int RenderKind)
{
	if(RenderKind==1) //Using nothing
	{
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable (GL_BLEND);
		Render_Transparency(1.0f,true);
	}
	else if(RenderKind==2) //No depthing
	{
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable (GL_BLEND);
		glDepthMask(GL_FALSE);
		Render_Transparency(0.8f,true);
		glDepthMask(GL_TRUE);

	}
	else if(RenderKind==3) //Using multiple passes
	{
		double fac=0.75f;
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable (GL_BLEND);
		///////////////////////////////////////
		glDisable(GL_CULL_FACE);
		glDepthFunc(GL_LESS);
		Render_Atoms();
		Render_Transparency(0.0f,true);
		// render with alpha = 0, to prime the depth buffer

		glEnable(GL_CULL_FACE);
		glCullFace(GL_FRONT);
		glDepthFunc(GL_ALWAYS);
		Render_Transparency(fac,true);
		// render with alpha = f*alpha

		glEnable(GL_CULL_FACE);
		glCullFace(GL_FRONT);
		glDepthFunc(GL_LEQUAL);
		Render_Transparency(fac,false);
		// render with alpha = (alpha-f*alpha)/(1.0-f*alpha)


		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		glDepthFunc(GL_ALWAYS);
		Render_Transparency(fac,true);
		// render with alpha = f*alpha

		// There's a trade off here. With culling enabled then a perfectly
		// opaque object (alpha=1) may be wrong. With it disabled, ordering
		// artifacts may appear
		// glEnable(GL_CULL_FACE);
		// glCullFace(GL_BACK);
		glDisable(GL_CULL_FACE);
		glDepthFunc(GL_LEQUAL);
		Render_Transparency(fac,false);

		//back to original
		glDisable(GL_CULL_FACE);
		glDisable (GL_BLEND);
		glDepthFunc(GL_LEQUAL);	
	};



}

/**
* @brief Render an translucent outline of the atom to be selected and an fainted outline
* of the atom selected. Also show the ID of the atom. If proper flag is true show
* the ID of all atoms.
*/
System::Void OpenGLForm::COpenGL::Render_SelectedAtoms(System::Void)
{
	float red,green,blue,scaling; 
	float num=1.5f/300.f;			//scaling factor
	
	//Enable culling
	glFrontFace(GL_CCW);              // set CCW as fron face direction
	glPolygonMode(GL_FRONT, GL_FILL); // set front face fill mode
	glPolygonMode(GL_BACK,  GL_LINE); // set back face fill mode
	glCullFace(GL_BACK);              // set for culling back faces
	glEnable(GL_CULL_FACE);           // enable face culling
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

//Loop over all atoms for drawing all the IDs when needed
 for each(int i in MainCollection->Keys)
	{
	if (((Atom^)MainCollection[i])->isotope->Element!="e")
	   {
		//Draw the fainted outline for the atom selected
		if(SelectFlag && (SelectedAtom==i-1 ))
		{
			glPushMatrix();
			//The same position as the atom
			glTranslatef((GLfloat)((Atom^)MainCollection[i])->X,
						 (GLfloat)((Atom^)MainCollection[i])->Y,
						 (GLfloat)((Atom^)MainCollection[i])->Z);

			//Color a litle bit more fainted than  the  color of the atom
			red=((Atom^)MainCollection[i])->isotope->Red;
			green=((Atom^)MainCollection[i])->isotope->Green;
			blue=((Atom^)MainCollection[i])->isotope->Blue;
			glColor4f(red+0.1f,green+0.1f,blue+0.1f,0.6f);

			//Large than the atom by a factor num
			scaling=num*((Atom^)MainCollection[i])->isotope->Radius;
			glScalef((GLfloat)scaling,(GLfloat)scaling,(GLfloat)scaling);

			//Draw
			if (LowRes)glCallList(SelectedSphereLowRes); else glCallList(SelectedSphereHighRes);
			glPopMatrix();
		};
		
		//Draw the translucent outline for the atom to be selected
		if(SelectFlag && (ToBeSelectedAtom==i-1))
		{
			glPushMatrix();
			//The same position as the atom
			glTranslatef((GLfloat)((Atom^)MainCollection[i])->X,
						 (GLfloat)((Atom^)MainCollection[i])->Y,
						 (GLfloat)((Atom^)MainCollection[i])->Z);

			//Color a litle bit more translucent than  the  color of the atom
			red=((Atom^)MainCollection[i])->isotope->Red;
			green=((Atom^)MainCollection[i])->isotope->Green;
			blue=((Atom^)MainCollection[i])->isotope->Blue;
			glColor4f(red+0.1f,green+0.1f,blue+0.1f,0.2f);

			//Large than the atom by a factor num
			scaling=num*((Atom^)MainCollection[i])->isotope->Radius;
			glScalef((GLfloat)scaling,(GLfloat)scaling,(GLfloat)scaling);
			//Draw
			if (LowRes)glCallList(SelectedSphereLowRes); else glCallList(SelectedSphereHighRes);
			glPopMatrix();
		};

		//Label atoms with their ID to be selected or when flag of drawing all atoms ID is true
		if(AtomIDFlag || ToBeSelectedAtom==i-1) 
			Render_IDText(((Atom^)MainCollection[i])->X, 
				          ((Atom^)MainCollection[i])->Y, 
				          ((Atom^)MainCollection[i])->Z,
						  (GLfloat)((Atom^)MainCollection[i])->isotope->Radius/300.0f,
				          ((Atom^)MainCollection[i])->getID());
		
	    };
	};

		//Disable culling
 	    //glEnable(GL_CULL_FACE);
	    glPolygonMode(GL_FRONT, GL_FILL); // set front face fill mode
        glPolygonMode(GL_BACK,  GL_FILL); // set back face fill mode
        glCullFace(GL_BACK);              // set for culling back faces
        glDisable(GL_CULL_FACE); 
		glDisable (GL_BLEND);

}

/**
* @brief Renders OpenGL text of the integer number ID in distance radius from a point lied in (x,y,z). 
* @param x Coordinate x.
* @param y Coordinate y.
* @param z Coordinate z.
* @param radius Distance form the point where the text will be printed at the right.
* @param ID Integer number which will be printed as a text.
*/
System::Void OpenGLForm::COpenGL::Render_IDText(double x, double y, double z, double radius, int ID)
{
	//Extract coordinates
	glPushMatrix();
	glTranslatef((GLfloat)x,(GLfloat)y,	(GLfloat)z);
	GLdouble modelview[16];							//var to hold the modelview info
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview );  //get the modelview info
	glPopMatrix();


	glPushMatrix();				 //save rotation
	glMatrixMode(GL_MODELVIEW);  // reset to initial  reference frame
	glLoadIdentity();

	//Disable Lighting and depth test
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glColor3f(0.0f, 0.0f, 0.0f);

	//Draw to world coordinates
	glRasterPos3f((GLfloat)(modelview[12]+radius),
				  (GLfloat)(modelview[13]+radius),
		  		  (GLfloat)(modelview[14]+radius));
	
	//Set the font
	glListBase(nFontList);

	//Write in a buffer and then to OpenGL machinery
	char buf[std::numeric_limits<int>::digits10+3];
	sprintf_s(buf,"%d", ID);
	string num=(string)buf;
	glCallLists((GLsizei)num.length(), GL_UNSIGNED_BYTE, buf);

	//Enable  Lighting and depth test
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	glPopMatrix();// back to original orientation


}

/**
* @brief Renders all atoms and compiles the selections ids for The OpenGL machinery. 
*/
System::Void OpenGLForm::COpenGL::Render_Atoms(System::Void)
{
    // Opaque sphere facets must not blend their shared edges with the background.
    glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
    glFrontFace(GL_CCW);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);
    glDisable(GL_POLYGON_SMOOTH);
    glDisable(GL_BLEND);

	float red,green,blue,scaling; 
	float num=1.0f/300.f;			//scaling factor

	//Draw all atoms with color and prepare for selection
	for each(int i in MainCollection->Keys)
	{
		if (((Atom^)MainCollection[i])->isotope->Element!="e")
		{
			glPushMatrix();

			//Set the ID in OpenGL machinery for selection
			glLoadName (i-1);
			//The position of the atom
			glTranslatef((GLfloat)((Atom^)MainCollection[i])->X,
						 (GLfloat)((Atom^)MainCollection[i])->Y,
						 (GLfloat)((Atom^)MainCollection[i])->Z);
				
			//Scaling with a num factor
			scaling=num*((Atom^)MainCollection[i])->isotope->Radius;
			glScalef((GLfloat)scaling,(GLfloat)scaling,(GLfloat)scaling);

			//Color the atom
			red=((Atom^)MainCollection[i])->isotope->Red;
			green=((Atom^)MainCollection[i])->isotope->Green;
			blue=((Atom^)MainCollection[i])->isotope->Blue;
			glColor3f(red,green,blue);

			//Draw
			if (LowRes) glCallList(SphereLowRes); else glCallList(SphereHighRes);
			glPopMatrix();
		};
	};

    glPopAttrib();
}

/**
* @brief Render the transparent ellispoids of all atoms nuclei with the specified translucency
* @param factor A factor multiplied  with the alpha channel of the ellipsoids.
* @param check A flag for how to use facotr and alpha in fun function.
*/
System::Void OpenGLForm::COpenGL::Render_Transparency(double factor, bool check)
{
	//Setup the minimum chemical bond, if zero put it sero
	double validminCB=MainICollection->InterKindMin[InteractionKind::CBond]> epsilo && 
		              MainICollection->InterKindMin[InteractionKind::CBond]< dMax ? 
					  MainICollection->InterKindMin[InteractionKind::CBond]: 1.0f;

	//Calculate the factors of Chemical shifiting, Hyperfine and Quadrupolar according to
	// the chemical bond minimum and their maxima.
	CSfinfactor= CSfactor*validminCB/(abs(MainICollection->InterKindMax[InteractionKind::Shift])>
		epsilo ?
		(MainICollection->InterKindMax[InteractionKind::Shift]) 
		: 1.0f);
	HCfinfactor=HCfactor*validminCB/(abs(MainICollection->InterKindMax[InteractionKind::HFC])>
		epsilo ?
		(MainICollection->InterKindMax[InteractionKind::HFC]) 
		: 1.0f);
	QCfinfactor= QCfactor*validminCB/(abs(MainICollection->InterKindMax[InteractionKind::Quadrupolar])>
		epsilo ?
		(MainICollection->InterKindMax[InteractionKind::Quadrupolar]) 
		: 1.0f);

	//loop over all interaction to draw the interaction ellipsoids
	for each(int i in MainICollection->Keys)
		{
		//Check to draw the interactions in the selected atom or not
		if(!(((Interaction^)MainICollection[i])->A->getID()==(SelectedAtom+1)) ||
		   SelectedAtomsInteractionsFlag)	
		{
			//Drawing CS (Chemical shift) ellipsoids
			if (((Interaction^)MainICollection[i])->IntKind==InteractionKind::Shift && (CSFlag))
				Render_Atom_Properties(((Interaction^)MainICollection[i])->A->X,
				((Interaction^)MainICollection[i])->A->Y,
				((Interaction^)MainICollection[i])->A->Z,
				((Tensor^)MainICollection[i])->dcm, 
				((Tensor^)MainICollection[i])->eigen, 
				CSfinfactor, 
				CScolor, fun(0.4f,factor,check),
				i);

			//Drawing HFC (Hyperfine) or AC (Anisotropic) ellipsoids
			if (((Interaction^)MainICollection[i])->IntKind==InteractionKind::HFC && (HCFlag))
				Render_Atom_Properties(((Interaction^)MainICollection[i])->A->X,
				((Interaction^)MainICollection[i])->A->Y,
				((Interaction^)MainICollection[i])->A->Z,
				((Tensor^)MainICollection[i])->dcm, 
				((Tensor^)MainICollection[i])->eigen, 
				HCfinfactor,
				HCcolor, fun(0.35f,factor,check),
				i);

			//Drawing Zero-field ellipsoids
			if (((Interaction^)MainICollection[i])->IntKind==InteractionKind::Quadrupolar && (QCFlag))
				Render_Atom_Properties(((Interaction^)MainICollection[i])->A->X,
				((Interaction^)MainICollection[i])->A->Y,
				((Interaction^)MainICollection[i])->A->Z,
				((Tensor^)MainICollection[i])->dcm, 
				((Tensor^)MainICollection[i])->eigen, 
				QCfinfactor,  
				QCcolor, fun(0.35f,factor,check),
				i);
		};
	
		};

}

/**
* @brief Render translucent ellispoid and eigen values axes
* @param X Coordinate x.
* @param Y Coordinate y.
* @param Z Coordinate z.
* @param eigenvec Matrix with the eigenvectors(DCM).
* @param eigen Matrix with the eigenvalues.
* @param factor A factor multiplied  with the alpha channel of the ellipsoids.
* @param color Matrix with RGB color.
* @param alpha Alpha channel.
* @param ID ID of the atom to draw.
*/
System::Void OpenGLForm::COpenGL::Render_Atom_Properties(double X, double Y, double Z, 
														 ObservableCollection<double>^ eigenvec, 
	                                                     ObservableCollection<double>^ eigen,  
														 double factor, 
														 double *color, 
														 double alpha, 
														 int ID)
{
	float TransformMatrix[16];

	//Draw Axes of reference frame	
	glPushMatrix();
	float scalesize=0.75;
	//Draw when the flag is true and when an interation is selected
	if(ReferenceFrameAxesFlag &&
		((SpinachGUI::MainForm^)mainForm)->dataGridView2->SelectedRows->Count>0 &&
		Convert::ToInt32(((SpinachGUI::MainForm^)mainForm)->dataGridView2->SelectedRows[0]->Cells[0]->Value)==ID)
	{
		glDepthMask(GL_TRUE);
		float tX=0, tY=0, tZ=0;
		//Position in the same point as ellispoid
		glTranslatef((GLfloat)X, (GLfloat)Y,(GLfloat)Z);

	   //Rotation according to eigen vector of the reference frame
	   TransformMatrix[0]=(GLfloat)((Tensor^)MainICollection[ID])->Frame->matrix3x3[0];
	   TransformMatrix[1]=(GLfloat)((Tensor^)MainICollection[ID])->Frame->matrix3x3[1];
	   TransformMatrix[2]=(GLfloat)((Tensor^)MainICollection[ID])->Frame->matrix3x3[2];
	   TransformMatrix[3]=(GLfloat)0.0;

	   TransformMatrix[4]=(GLfloat)((Tensor^)MainICollection[ID])->Frame->matrix3x3[3];
	   TransformMatrix[5]=(GLfloat)((Tensor^)MainICollection[ID])->Frame->matrix3x3[4];
	   TransformMatrix[6]=(GLfloat)((Tensor^)MainICollection[ID])->Frame->matrix3x3[5];
	   TransformMatrix[7]=(GLfloat)0.0;

	   TransformMatrix[8]=(GLfloat)((Tensor^)MainICollection[ID])->Frame->matrix3x3[6];
	   TransformMatrix[9]=(GLfloat)((Tensor^)MainICollection[ID])->Frame->matrix3x3[7];
	   TransformMatrix[10]=(GLfloat)((Tensor^)MainICollection[ID])->Frame->matrix3x3[8];
	   TransformMatrix[11]=(GLfloat)0.0;
	
	   TransformMatrix[12]=TransformMatrix[13]=TransformMatrix[14]=0.0;TransformMatrix[15]=(GLfloat)1.0;
	   glMultMatrixf(TransformMatrix);

		// DRAW WORLD ORIGIN AXIS
		glShadeModel(GL_SMOOTH);

		//X-Axis
		glPushMatrix();
		glColor3f(0.0,1.0,0.0);
		glTranslated((GLfloat)tX, (GLfloat)tY,(GLfloat)tZ);
		glRotated(-90, 1,0,0);
		glScalef(scalesize,scalesize,scalesize);
		glTranslated(0,0,1);
		gluCylinder(quadratic,0.05,0.0001,0.1,20,20);    //Head of Arrow
		glTranslated(0,0,-1);
		gluCylinder(quadratic,0.03,0.03,1.0,20,20);  //Axis
		glPopMatrix();

		//Y-Axis
		glPushMatrix();
		glColor3f(0.0,0.0,1.0);
		glTranslated((GLfloat)tX, (GLfloat)tY,(GLfloat)tZ);
		glRotated(90,0 ,0,1);
		glScalef(scalesize,scalesize,scalesize);
		glTranslated(0,0,1);
		gluCylinder(quadratic,0.05,0.0001,0.1,20,20);  //Axis
		glTranslated(0,0,-1);
		gluCylinder(quadratic,0.03,0.03,1.0,20,20);  //Head of Arrow
		glPopMatrix();

		//Z-Axis
		glPushMatrix();
		glTranslated((GLfloat)tX, (GLfloat)tY,(GLfloat)tZ);
		glColor3f(1.0,0.0,0.0);
		glRotated(90, 0,1,0);
		glScalef(scalesize,scalesize,scalesize);
		glTranslated(0,0,1);
		gluCylinder(quadratic,0.05,0.0001,0.1,20,20);  //Axis
		glTranslated(0,0,-1);
		gluCylinder(quadratic,0.03,0.03,1.0,20,20);  //Head of Arrow
		glPopMatrix();
		glDepthMask(GL_FALSE);


	};
	glPopMatrix();


	//Draw interaction
	glPushMatrix();
	glTranslatef((GLfloat)X, (GLfloat)Y,(GLfloat)Z);
	
	//In order to draw when bigger than the atom
	float show=0.0033*((Tensor^)MainICollection[ID])->A->isotope->Radius; 
	if(factor*abs(eigen[0])>show ||  factor*abs(eigen[1])>show ||  factor*abs(eigen[2])>show)
	  {
      //Rotation according to eigen vector
	  TransformMatrix[0]=(GLfloat)eigenvec[0];TransformMatrix[1]=(GLfloat)eigenvec[1];
	  TransformMatrix[2]=(GLfloat)eigenvec[2]; TransformMatrix[3]=(GLfloat)0.0;

	  TransformMatrix[4]=(GLfloat)eigenvec[3];TransformMatrix[5]=(GLfloat)eigenvec[4];
	  TransformMatrix[6]=(GLfloat)eigenvec[5];TransformMatrix[7]=(GLfloat)0.0;

	  TransformMatrix[8]=(GLfloat)eigenvec[6];TransformMatrix[9]=(GLfloat)eigenvec[7];
	  TransformMatrix[10]=(GLfloat)eigenvec[8];TransformMatrix[11]=(GLfloat)0.0;
	
	  TransformMatrix[12]=TransformMatrix[13]=TransformMatrix[14]=0.0;TransformMatrix[15]=(GLfloat)1.0;
	  glMultMatrixf(TransformMatrix);
	  
	  //Drawing of eigen values lines
	  if(EigenAxesFlag)
	  {
	  //Color red for positive eigenvalue and blue for negative eigenvalue
	  if(eigen[0]<=0)  glColor3f(0.0f,0.0f,1.0f); else glColor3f(1.0f,0.0f,0.0f);
	  glPushMatrix();
	  glRotated(90,0,1,0);
	  glTranslatef(0.0f,0.0f,-factor*abs(eigen[0]));
	  glScalef(1.0f,1.0f,2.0*factor*abs(eigen[0]));
	  glCallList(CylinderHighRes);
	  glPopMatrix();

	  //Color red for positive eigenvalue and blue for negative eigenvalue
	  if(eigen[1]<=0)  glColor3f(0.0f,0.0f,1.0f); else glColor3f(1.0f,0.0f,0.0f);
	  glPushMatrix();
	  glRotated(-90,1,0,0);
	  glTranslatef(0.0f,0.0f,-factor*abs(eigen[1]));
	  glScalef(1.0f,1.0f,2.0*factor*abs(eigen[1]));
	  glCallList(CylinderHighRes);
	  glPopMatrix();

	  //Color red for positive eigenvalue and blue for negative eigenvalue
	  if(eigen[2]<=0)  glColor3f(0.0f,0.0f,1.0f); else glColor3f(1.0f,0.0f,0.0f);
	  glPushMatrix();
	  glRotated(90,0,0,1);
	  glTranslatef(0.0f,0.0f,-factor*abs(eigen[2]));
	  glScalef(1.0f,1.0f,2.0*factor*abs(eigen[2]));
	  glCallList(CylinderHighRes);
	  glPopMatrix();
	  }


	  //Drawing of eigen values ellispoid
	  if(EigenEllipsoidsFlag)
	  {
	  glPushMatrix();
	  glScaled(factor*abs(eigen[0]),  factor*abs(eigen[1]),	  factor*abs(eigen[2])); 

	  //Setting up the material and coloring of the ellipsoids
	  GLfloat material[4];
      material[0] =(GLfloat)color[0]; material[1] =(GLfloat)color[1];
	  material[2] =(GLfloat)color[2]; material[3]=  0.7f;
	  GLfloat white[4]={0.8f,0.8f,0.8f,0.0f};
	  glMaterialfv(GL_FRONT, GL_SPECULAR, white);
      glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material);

	  //Drawing of wire
	  glColor4f((GLfloat)color[0],(GLfloat)color[1],(GLfloat)color[2],(GLfloat)(alpha+0.1f));
	  glCallList(WiredSphereHighRes);

      //Drawing of solid ellispoid
	  glColor4f((GLfloat)color[0],(GLfloat)color[1],(GLfloat)(color[2]+0.1f),(GLfloat)alpha);
	  glCallList(SelectedSphereHighRes);
	  glPopMatrix();
	  };

	  };

	glPopMatrix();

	
}

/**
* @brief Render bonds and j-couplings
*/
System::Void OpenGLForm::COpenGL::Render_Bonds(System::Void)
{
 double RADDEG = 57.29578;

  //Factor for bonds
 CBfinfactor=2.0*(CBfactor-0.5)*(MainICollection->InterKindMax[InteractionKind::CBond]-1.6f)+1.6f;

 //Factor for J-couplings
 double absMax=MainICollection->InterKindMax[InteractionKind::Jcoupling];
 double absMin=MainICollection->InterKindMin[InteractionKind::Jcoupling];

 JCfinfactor= JCfactor*absMax;

 //Loop over all the interactions to draw bond like interactions
 for each(int i in MainICollection->Keys)
 {
	 //Chemical bonds represented as cylinder without bases
	 if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::CBond && CBFlag)
	 {
		 //Calculate the distance
		 double distance=((Tensor^)MainICollection[i])->matrix3x3[0];
		 
		 //For bonds only smaller from the factor of the trackbar
		 if((distance-epsilo)<=CBfinfactor)
		 {
			 glPushMatrix();
			 Atom^ A=((Tensor^)MainICollection[i])->A;
			 Atom^ B=((Tensor^)MainICollection[i])->B;

			 //Origin of the bond
			 glTranslatef((GLfloat)A->X,(GLfloat)A->Y,(GLfloat)A->Z);

			 //Orientation
			 double az = atan2 ((B->Y-A->Y), (B->X-A->X)) * RADDEG;
			 double el = asin  ((B->Z-A->Z) / distance) * RADDEG;
			 glRotated ( az, 0,0,1);
			 glRotated (-el, 0,1,0);
			 glRotatef ( 90, 0,1,0);

			 //Coloring
			 glColor3f((GLfloat)CBcolor[0],(GLfloat)CBcolor[1], (GLfloat)CBcolor[2]);

			 //If bond is hydrogen bond(between this elements) draw from small spheres
			 if(((A->isotope->Element=="O" || A->isotope->Element=="F" || A->isotope->Element=="N") &&
				 B->isotope->Element=="H") ||
				 ((B->isotope->Element=="O" || B->isotope->Element=="F" || B->isotope->Element=="N") &&
				 A->isotope->Element=="H") )
			 {
				 //Spheres only for distances bigger than 1.2f
				 if((distance-epsilo)>=1.2f && (distance-epsilo)<=2.5f)
					 for(float j=0.0;j<distance;j+=0.075f)
					 {
						 glPushMatrix();
						 glTranslatef((GLfloat)0.0f,(GLfloat)0.0f, (GLfloat)j);
						 glScalef((GLfloat)0.035,(GLfloat)0.035, (GLfloat)0.035);
						 if(LowRes) glCallList(SphereLowRes); else glCallList(SphereHighRes);
						 glPopMatrix();
					 }
				 else if((distance-epsilo)<1.2f)
				 {
					 glScalef(1.0,1.0, distance);
					 if(LowRes) glCallList(CylinderLowRes); else glCallList(CylinderHighRes);
				 }
			 }
			 else
			 {
				 glScalef(1.0,1.0, distance);
				 if(LowRes) glCallList(CylinderLowRes); else glCallList(CylinderHighRes);
			 }

			 glPopMatrix();
		 };
	 };

   //J-couplings are represented as dashed cylinder with bases
	 if(((Interaction^)MainICollection[i])->IntKind==InteractionKind::Jcoupling && JCFlag)
	 {
		 //For j-couplings only larger from the factor of the trackbar
		 if(abs(((Tensor^)MainICollection[i])->matrix3x3[0])>JCfinfactor) 
		 {
			 glPushMatrix();
			 Atom^ A=((Tensor^)MainICollection[i])->A;
			 Atom^ B=((Tensor^)MainICollection[i])->B;
			 double distance=sqrt((A->X-B->X)*(A->X-B->X)+(A->Y-B->Y)*(A->Y-B->Y)+(A->Z-B->Z)*(A->Z-B->Z));

			 //Origin of the bond
			 glTranslatef((GLfloat)A->X,(GLfloat)A->Y,(GLfloat)A->Z);

			 //Orientation
			 double az = atan2 ((B->Y-A->Y), (B->X-A->X)) * RADDEG;
			 double el = asin  ((B->Z-A->Z) / distance) * RADDEG;
			 glRotated ( az, 0,0,1);
			 glRotated (-el, 0,1,0);
			 glRotatef ( 90, 0,1,0);

			 //Factor for Coloring
			 double magn=absMax>absMin ? (abs(((Tensor^)MainICollection[i])->matrix3x3[0])-absMin)/(absMax-absMin) : 0.0;
			 glColor3f((GLfloat)(magn*JCcolor[3]+(1.0-magn)*JCcolor[0]),
					   (GLfloat)(magn*JCcolor[4]+(1.0-magn)*JCcolor[1]),
					   (GLfloat)(magn*JCcolor[5]+(1.0-magn)*JCcolor[2]));

			 //Draw the dashed cylinders
			 for(float j=0.0;j<distance;j+=0.2f)
			 {
				 glPushMatrix();
				 glTranslatef((GLfloat)0.0f,(GLfloat)0.0f, (GLfloat)j);
				 if(LowRes) glCallList(FullCylinderLowRes); else glCallList(FullCylinderHighRes);
				 glPopMatrix();

			 };

			 glPopMatrix();
		 };
	 };
 };

}
#pragma endregion

/////Render electron cloud plane//////////
#pragma region Render electron cloud plane
/**
* @brief Render electron cloud plane. Usually a verysmall number of electron are drawn,
* so not so much optimization is needed
*/
System::Void OpenGLForm::COpenGL::Render_Electrons()
{	
	 double RADDEG = 57.29578;

 	//Render Floor 
    Render_Floor(-15, 1);

	//Set Clip Pane
	glColor4f(1, 1, 1, 0.05f);
	glEnable(GL_STENCIL_TEST);	//Enable using the stencil buffer
	glColorMask(0, 0, 0, 0);	//Disable drawing colors to the screen
	glDepthMask(false); 
	glDisable(GL_DEPTH_TEST);	//Disable depth testing


	glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE);
	glStencilFunc(GL_ALWAYS, 1, 0xffffffff);

	glDepthMask(true); 
	glColorMask(1, 1, 1, 1); //Enable drawing colors to the screen
	glEnable(GL_DEPTH_TEST); //Enable depth testing
	//Make the stencil test pass only when the pixel is 1 in the stencil buffer
	glStencilFunc(GL_EQUAL, 1, 1);
	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP); //Make the stencil buffer not change

	double eqr[] = {0.0f,-1.0f, 0.0f, -18.0f};  
	glEnable(GL_CLIP_PLANE0);                       // Enable Clip Plane For Removing Artifacts
	// (When The Object Crosses The Floor)
	glClipPlane(GL_CLIP_PLANE0, eqr);               // Equation For Reflected Objects

	
	//Matrices for Electrons
	ArrayList^ ElectronIndex=gcnew ArrayList();
	ArrayList^ Xcoord=gcnew ArrayList();
	ArrayList^ Ycoord=gcnew ArrayList();
	ArrayList^ Zcoord=gcnew ArrayList();
	ArrayList^ index=gcnew ArrayList();

	//Find electron and IDs of electron and seve in a array
	for each(int ii in MainCollection->Keys)
	{ 
		if(((Atom^)SystemModel->AtomCollection[ii])->isotope->Element=="e") 
		{
			ElectronIndex->Add(((Atom^)SystemModel->AtomCollection[ii]));
			index->Add(((Atom^)SystemModel->AtomCollection[ii])->ID);
		}
	};

	//Calculate coordinates of electron. These are calculated according to: 2-in line, 
	// 3-in triangle, >4-in quad
	for(int ii=0;ii<ElectronIndex->Count;ii++)
	{

		Xcoord->Add(-7.0f*!((ElectronIndex->Count)%2 &&(ElectronIndex->Count-1==ii))+14.0f*(ii%2));
		Ycoord->Add(-25.0f);
		Zcoord->Add(-80.0f-10.0f*((ElectronIndex->Count-1)/2)+20.0f*(ii/2));
	};
	
	//Loop over all electrons
	for(int ii=0;ii<ElectronIndex->Count;ii++)
	{
		
		glPushMatrix();
		glTranslated(Convert::ToDouble(Xcoord[ii]),
					 Convert::ToDouble(Ycoord[ii]),
					 Convert::ToDouble(Zcoord[ii]));

		///////////////////////////////////////
		//Draw the exchange couplings without any rotation on them

		glDisable(GL_LIGHTING);

		//Calcalate max for exchange coupling
		double absMax=MainICollection->InterKindMax[InteractionKind::Exchange];
		ECfinfactor= ECfactor*absMax;

		//Loop over all interaction to draw Exchange couplings
		for each(int jj in MainICollection->Keys)
		{
			//Draw exchange coupling
			if(((Interaction^)MainICollection[jj])->IntKind==InteractionKind::Exchange && ECFlag &&
				((Interaction^)MainICollection[jj])->A==((Atom^)ElectronIndex[ii]))
			{
				//Coordinates of the first electron
				double AX=Convert::ToDouble(Xcoord[ii]);
			    double AY=Convert::ToDouble(Ycoord[ii]);
				double AZ=Convert::ToDouble(Zcoord[ii]);

				//Find the coordinates of the second electron 
				int kk;
				int test=((Interaction^)MainICollection[jj])->B->ID;
				for(int k=0;k<ElectronIndex->Count;k++) 
				{
					if(Convert::ToInt32(index[k])==((Interaction^)MainICollection[jj])->B->ID)
					{
						kk=k;
					}
				}
				double BX=Convert::ToDouble(Xcoord[kk]);
				double BY=Convert::ToDouble(Ycoord[kk]);
				double BZ=Convert::ToDouble(Zcoord[kk]);


				glPushMatrix();
				//Calculate distance 
				double distance=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY)+(AZ-BZ)*(AZ-BZ));
				if (distance<0.000001) distance=1.0;

				//The enclosing matrix is already translated to the first electron.

				//Orientation
				double az = atan2 ((BY-AY), (BX-AX)) * RADDEG;
				double el = asin  ((BZ-AZ) / distance) * RADDEG;
				glRotated ( az, 0,0,1);
				glRotated (-el, 0,1,0);
				glRotatef ( 90, 0,1,0);

				//Factor for Coloring
				double minimum=MainICollection->InterKindMin[InteractionKind::Exchange];
				double maximum=MainICollection->InterKindMax[InteractionKind::Exchange];
				double magn=maximum>minimum ? (abs(((Tensor^)MainICollection[jj])->matrix3x3[0])-minimum)/(maximum-minimum) : 0.0;
				glColor3f((GLfloat)(magn*ECcolor[3]+(1.0-magn)*ECcolor[0]),
						  (GLfloat)(magn*ECcolor[4]+(1.0-magn)*ECcolor[1]),
						  (GLfloat)(magn*ECcolor[5]+(1.0-magn)*ECcolor[2]));

				//Draw the spring if bigger than the exchange coupling factor of trackbar
				if(abs(((Tensor^)MainICollection[jj])->matrix3x3[0])>ECfinfactor)
				{
					glLineWidth((GLfloat)4.0);

					//Draw a spring
					glBegin(GL_LINE_STRIP);
					int steps=100;
					for(int i=0;i<steps;i++)
					glVertex3f((GLfloat)(0.8*cos(i*distance/(steps/5))),
							   (GLfloat)(0.8*sin(i*distance/(steps/5))),
							   (GLfloat)(i*distance/steps));

					
					glEnd();

					// Back to default
					glLineWidth((GLfloat)1.0);
				};

				glPopMatrix();
			}

		}

		glEnable(GL_LIGHTING);

		/////////////////////////////////////
		//Draw the rest thing rotating them
		glMultMatrixf(TransformMat);
		glColor3f(0.7f,0.5f,0.7f);

		//Draw metaballs/electron
		glPushMatrix();
		glTranslated(0.40f,-0.6f,-0.3f); //For centering Metaball
		glScaled(0.4,0.4,0.4);
		Render_MetaBalls();	
		glPopMatrix();

		//Render electron ID
		if(AtomIDFlag) 	Render_IDText(0.0, 0.0, 0.0,2.0f, ((Atom^)ElectronIndex[ii])->getID());

		//Draw interactions elllipsoids of electrons
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable (GL_BLEND);
		glDepthMask(GL_FALSE);
		
		//Calculate the factors for g-tensor and zero-field
		double validminCB=3.0f; //standard
		GTfinfactor= GTfactor*validminCB/(abs(MainICollection->InterKindMax[InteractionKind::GTensor])>
			epsilo ?
			(MainICollection->InterKindMax[InteractionKind::GTensor]) 
			: 1.0f);
		ZFfinfactor= ZFfactor*validminCB/(abs(MainICollection->InterKindMax[InteractionKind::ZFS])>
			epsilo ?
			(MainICollection->InterKindMax[InteractionKind::ZFS]) 
			: 1.0f);

		//Loop over interactions to find g-tensors and zero-fields
		for each(int jj in MainICollection->Keys)
		{
			if (((Interaction^)MainICollection[jj])->IntKind==InteractionKind::GTensor && GTFlag &&
				((Interaction^)MainICollection[jj])->A==((Atom^)ElectronIndex[ii]))
				Render_Atom_Properties(0.0, 0.0, 0.0,
				((Tensor^)MainICollection[jj])->dcm, 
				((Tensor^)MainICollection[jj])->eigen, 
				GTfinfactor, GTcolor, 0.6f, jj);

			if (((Interaction^)MainICollection[jj])->IntKind==InteractionKind::ZFS && ZFFlag &&
				((Interaction^)MainICollection[jj])->A==((Atom^)ElectronIndex[ii]))
				Render_Atom_Properties(0.0, 0.0, 0.0,
				((Tensor^)MainICollection[jj])->dcm, 
				((Tensor^)MainICollection[jj])->eigen, 
				ZFfinfactor, ZFcolor, 0.6f, jj);
		}

		glPopMatrix();

		glDepthMask(GL_TRUE);
	};
	
	//Deactivate finally the clip pane
	GLfloat white[4]={0.8f,0.8f,0.8f,0.0f};
	glMaterialfv(GL_FRONT, GL_SPECULAR, white);
	glDisable(GL_STENCIL_TEST); //Disable using the stencil buffer
	glDisable(GL_CLIP_PLANE0);  

	

	

}

/**
* @brief Render a plane with name "electrons" in z and y coordinates
* @param y Coordinate y.
* @param z Coordinate z.
*/
System::Void OpenGLForm::COpenGL::Render_Floor(int z ,int y)
{

	//Blend the floor onto the screen
	glDisable(GL_DEPTH_TEST);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4f(0, 0, 0, 0.85f);

	glPushMatrix();
	int FLOOR_SIZE=200; //size of the plane
	glBegin(GL_QUADS);
	glNormal3d(0, 0, 1);
	glVertex3d(-y*FLOOR_SIZE / 2, z, y*FLOOR_SIZE / 2);
	glVertex3d(-y*FLOOR_SIZE / 2, z, -y*FLOOR_SIZE / 2);
	glVertex3d(y*FLOOR_SIZE / 2, z, -y*FLOOR_SIZE / 2);
	glVertex3d(y*FLOOR_SIZE / 2, z, y*FLOOR_SIZE / 2);
	glEnd();
	glPopMatrix();

	glDisable (GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	//Render "electrons" name///////////////////////////////////////////////////
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, (GLfloat)mWidth, 0.0,(GLfloat)mHeight);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	
	glListBase(nFontList);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	//Write the name in a string before givng it to OpenGL machinery
	char buf[10];
	glColor3f(1.0f, 1.0f, 1.0f);
	glRasterPos2i((GLint)10, (GLint)12);
	sprintf_s(buf,"ELECTRONS");
	glCallLists (9, GL_UNSIGNED_BYTE, &buf);

	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);

	glMatrixMode(GL_PROJECTION); //swapped this with...
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW); //...this
	glPopMatrix();

}

/**
* @brief Initiallizing and creating a metaball of spesific shape
*/
System::Void OpenGLForm::COpenGL::Init_Metaballs()
{
	const int minGridSize=10;
	int gridSize=20; //Was 40, changed for speed
	const int numMetaballs=8;
	METABALL metaballs[numMetaballs];

	//set up grid
	if(!cubeGrid->CreateMemory()){};
    //		return false;
	if(!cubeGrid->Init(gridSize)){};
	//		return false;

	for(int i=0; i<numMetaballs; i++)
		metaballs[i].Init(VECTOR3D(0.0f, 0.0f, 0.0f), 5.0f+float(i));

	//Using 8 balls putting them "carefully" randomly  
	metaballs[0].position.x=-6.0f;
	metaballs[0].position.y=0.0f;
	metaballs[0].position.z=0.0f;

	metaballs[1].position.x=6.0f;
	metaballs[1].position.y=0.0f;
	metaballs[1].position.z=0.0f;

	metaballs[2].position.x=0.0f;
	metaballs[2].position.y=5.0f;
	metaballs[2].position.z=0.0f;

	metaballs[3].position.x=0.0f;
	metaballs[3].position.y=-6.5f;
	metaballs[3].position.z=0.0f;

	metaballs[4].position.x=0.0f;
	metaballs[4].position.y=0.0f;
	metaballs[4].position.z=5.0f;

	metaballs[5].position.x=0.0f;
	metaballs[5].position.y=0.0f;
	metaballs[5].position.z=-6.0f;

	metaballs[6].position.x=6.0f;
	metaballs[6].position.y=6.0f;
	metaballs[6].position.z=6.0f;

	metaballs[7].position.x=-5.0f;
	metaballs[7].position.y=5.0f;
	metaballs[7].position.z=5.0f;

		//clear the field
	for(int i=0; i<cubeGrid->numVertices; i++)
	{
		cubeGrid->vertices[i].value=0.0f;
		cubeGrid->vertices[i].normal.LoadZero();
	}

	//evaluate the scalar field at each point
	VECTOR3D ballToPoint;
	float squaredRadius;
	VECTOR3D ballPosition;
	float normalScale;

	for(int i=0; i<numMetaballs; i++)
	{
		squaredRadius=metaballs[i].squaredRadius;
		ballPosition=metaballs[i].position;

		//VC++6 standard does not inline functions
		//by inlining these maually, in this performance-critical area,
		//almost a 100% increase in speed is found
		for(int j=0; j<cubeGrid->numVertices; j++)
		{
			//ballToPoint=cubeGrid.vertices[j].position-ballPosition;
			ballToPoint.x=cubeGrid->vertices[j].position.x-ballPosition.x;
			ballToPoint.y=cubeGrid->vertices[j].position.y-ballPosition.y;
			ballToPoint.z=cubeGrid->vertices[j].position.z-ballPosition.z;

			//get squared distance from ball to point
			//float squaredDistance=ballToPoint.GetSquaredLength();
			float squaredDistance=	ballToPoint.x*ballToPoint.x +
				ballToPoint.y*ballToPoint.y +
				ballToPoint.z*ballToPoint.z;
			if(squaredDistance==0.0f)
				squaredDistance=0.0001f;

			//value = r^2/d^2
			cubeGrid->vertices[j].value+=squaredRadius/squaredDistance;

			//normal = (r^2 * v)/d^4
			normalScale=squaredRadius/(squaredDistance*squaredDistance);
			//cubeGrid.vertices[j].normal+=ballToPoint*normalScale;
			cubeGrid->vertices[j].normal.x+=ballToPoint.x*normalScale;
			cubeGrid->vertices[j].normal.y+=ballToPoint.y*normalScale;
			cubeGrid->vertices[j].normal.z+=ballToPoint.z*normalScale;
		}
	}

}

/**
* @brief Final Render and orientation of a metaball 
*/
System::Void OpenGLForm::COpenGL::Render_MetaBalls()
{
	float threshold=1.1f;
	float TransformMatrix[16];


	//Set orientation "randomly nice"
	glPushMatrix();
	TransformMatrix[0]=(GLfloat)0.5159281;TransformMatrix[1]=(GLfloat)0.8540189;
	TransformMatrix[2]=(GLfloat)0.06682559; TransformMatrix[3]=(GLfloat)0.0;

	TransformMatrix[4]=(GLfloat)-0.6884727;TransformMatrix[5]=(GLfloat)0.4598062;
	TransformMatrix[6]=(GLfloat)-0.5608736;TransformMatrix[7]=(GLfloat)0.0;

	TransformMatrix[8]=(GLfloat)-0.5097246;TransformMatrix[9]=(GLfloat)0.2433634;
	TransformMatrix[10]=(GLfloat)0.8251974;TransformMatrix[11]=(GLfloat)0.0;

	TransformMatrix[12]=TransformMatrix[13]=TransformMatrix[14]=0.0;TransformMatrix[15]=(GLfloat)1.0;
	glMultMatrixf(TransformMatrix);

	glScalef(0.5f, 0.5f, 0.5f);
	cubeGrid->DrawSurface(threshold);

	glPopMatrix();


}
#pragma endregion

/////Render Navigation axes///////////////
#pragma region Render Navigation axes
/**
* @brief Render of navigation axes in the edge fo the window
*/
System::Void OpenGLForm::COpenGL::Render_Navigation()
{
	GLint viewport[4]; //var to hold the viewport info
	GLdouble modelview[16]; //var to hold the modelview info
	GLdouble projection[16]; //var to hold the projection matrix info
	GLfloat winX, winY, winZ; //variables to hold screen x,y,z coordinates
	GLdouble worldX, worldY, worldZ; //variables to hold world x,y,z coordinates
	float scalesize=15*0.75;  //size of the Navigation Axes

	glDisable(GL_FOG);
	
	//Extract coordinates form OpenGL
	glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info
	glGetDoublev( GL_PROJECTION_MATRIX, projection ); //get the projection matrix info
	glGetIntegerv( GL_VIEWPORT, viewport ); //get the viewport info

	//Setting the postion of navigation axes
	winX = (GLfloat) 0.05*mWidth;
	winY = (GLfloat) 0.10*mHeight;
	winZ = (GLfloat) 0.995;

	//get the world coordinates from the screen coordinates
	gluUnProject( winX, winY, winZ, modelview, projection, viewport, &worldX, &worldY, &worldZ);

	// DRAW WORLD ORIGIN AXIS
	glPushMatrix();
	glShadeModel(GL_SMOOTH);
	glTranslated(worldX, worldY, worldZ);
	glScalef(scalesize,scalesize,scalesize);
	gluQuadricDrawStyle(quadratic, GLU_FILL);
    gluQuadricNormals(quadratic, GLU_SMOOTH );

	//X-Axis
	glPushMatrix();
	glColor3f(0.0,1.0,0.0);
	glRotated(-90, 1,0,0);
	glTranslated(0,0,1);
	gluCylinder(quadratic,0.15,0.0001,0.1,20,20);  //Axis
	glTranslated(0,0,-1);
	gluCylinder(quadratic,0.1,0.1,1.0,20,20);  //Head of Arrow
	glPopMatrix();

	//Y-Axis
	glPushMatrix();
	glColor3f(0.0,0.0,1.0);
	glRotated(90,0 ,0,1);
	glTranslated(0,0,1);
	gluCylinder(quadratic,0.15,0.0001,0.1,20,20);  //Axis
	glTranslated(0,0,-1);
	gluCylinder(quadratic,0.1,0.1,1.0,20,20);  //Head of Arrow
	glPopMatrix();

	//Z-Axis
	glPushMatrix();
	glColor3f(1.0,0.0,0.0);
	glRotated(90, 0,1,0);
	glTranslated(0,0,1);
	gluCylinder(quadratic,0.15,0.0001,0.1,20,20);  //Axis
	glTranslated(0,0,-1);
	gluCylinder(quadratic,0.1,0.1,1.0,20,20);  //Head of Arrow
	glPopMatrix();

	glPopMatrix();

	glEnable(GL_FOG);
}
#pragma endregion

///// Selection Buffer///////////////////
#pragma region Selection Buffer
/**
* @brief Render a selection buffer according to the position of the mouse.
* @param MouseX Mouse pointer x Coordinate.
* @param MouseY Mouse pointer y Coordinate.
*/
System::Void OpenGLForm::COpenGL::SelectionBuffer(GLfloat MouseX,GLfloat MouseY)
{   
	   
	long len=1024;
    GLuint selectBuf[1024];
	GLint viewport[4];
    GLint mHits;


	//Set viewport and selection buffer
	glGetIntegerv(GL_VIEWPORT, viewport);
	glSelectBuffer(len, selectBuf);

    //Set render mode and initialize names
    glRenderMode(GL_SELECT);
	glInitNames();
    glPushName(-1);

	
	//For setting perspective
    glMatrixMode (GL_PROJECTION);
    //glPushMatrix ();
    glLoadIdentity ();
   

	//Create 5x5 pixel picking region near cursor location and perspective
    gluPickMatrix ((GLdouble) MouseX, (GLdouble) (viewport[3] - MouseY), 5.0, 5.0, viewport);
 

	//If using stereo, use right eye
	if(stereoview) Calc_eye(2, false);
	  else Calc_eye(0, false);
	
	//Zoom & Translation, Rotation and Centering
	glTranslatef((GLfloat)(TranslX/abs(Zoom)), (GLfloat)(TranslY/abs(Zoom)),(GLfloat) Zoom);
	glMultMatrixf(TransformMat);
	glTranslatef((GLfloat)(-(MainCollection->MaxX+MainCollection->MinX)/2.0f),
		         (GLfloat)(-(MainCollection->MaxY+MainCollection->MinY)/2.0f),
		         (GLfloat)(-(MainCollection->MaxZ+MainCollection->MinZ)/2.0f));

	glCallList(AllAtoms);

    //Finish Rendering but not swapp buffer(not needed to see it!)
    glFlush();

	//Find number of 3D objects hitted
    mHits=glRenderMode (GL_RENDER);
	 

	if (mHits) //For number of selected atoms more than 0
	{
		int n=0;
		double minz=selectBuf[1];

		//Find the minimum of the second element of each group of 4 elements
		for(int i=1;i<mHits;i++) if (selectBuf[1+i*4]<minz) 
		{
			n=i;
			minz=selectBuf[1+i*4];
		};
		
		//Select the minimum and take the ID of the 3d object saved in the fourth element
		if((int)selectBuf[3+n*4]>-1) ToBeSelectedAtom=selectBuf[3+n*4];
	}
	else ToBeSelectedAtom=-1; //Nothing selected
	
}
#pragma endregion

///// Mouse Handlers/////////////////////
#pragma region Mouse Handlers
/**
* @brief Delete the selected atom in OpneGL panel.
* @param sender The usual sender object for Event Handlers.
* @param e The usual EventArguments object for Event Handlers.
*/
void  OpenGLForm::COpenGL::Del_KeyDown(System::Object^  sender, System::Windows::Forms::KeyEventArgs^  e) 
{
	// For coordinates columns, determine whether the keystroke is the delete button.
	if (e->KeyCode == Keys::Delete && SelectedAtom>-1)
	{
		((SpinachGUI::MainForm^)mainForm)->dataGridView1->Rows->RemoveAt(SelectedAtom);
	};

 }

/**
* @brief Used for focus of the OpenGL panel when context menu closes.
* @param sender The usual sender object for Event Handlers.
* @param e The usual EventArguments object for Event Handlers.
*/
void OpenGLForm::COpenGL::OnMouseEnter(System::Object^ sender, System::EventArgs^  e)
{
	this->Focus();
	focused=true;
	
}

/**
* @brief Used for unfocus of the OpenGL panel when context menu opens.
* @param sender The usual sender object for Event Handlers.
* @param e The usual EventArguments object for Event Handlers.
*/
void OpenGLForm::COpenGL::OnMouseLeave(System::Object^ sender, System::EventArgs^  e)
{
	focused=false;
}

/**
* @brief For handling the wheel rotation of the mouse and used for zooming of the model.
* @param sender The usual sender object for Event Handlers.
* @param e The usual EventArguments object for Event Handlers.
*/
void OpenGLForm::COpenGL::OnMouseWheel(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	//Zooming using the wheel of the mouse
	Zoom += e->Delta/(WHEEL_DELTA/1);
	Render();	
}

/**
* @brief When the mouse move inside the OpenGL panel.
* @param sender The usual sender object for Event Handlers.
* @param e The usual EventArguments object for Event Handlers.
*/
void OpenGLForm::COpenGL::OnMouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{  
  //In order to not confused with the travelling
  if(!((SpinachGUI::MainForm^)mainForm)->OpenGLtimer->Enabled)
	{
	//For creating selecting buffer and finding the 3D object ID
	if(SelectFlag) SelectionBuffer((GLfloat) e->X, (GLfloat) e->Y);

	//Left mouse Button Rotation
	if( e->Button == ::MouseButtons::Left) 
	{
		//Use of ArcBall drag
		Quat4fT^ ThisQuat=gcnew Quat4fT();
		MousePt->X=(GLfloat) e->X;
		MousePt->Y=(GLfloat) e->Y;
		ArcBall->drag(ArcMargin+MousePt, ThisQuat);

		// Update End Vector And Get Rotation As Quaternion
		Matrix3fSetRotationFromQuat4f(ThisRot, ThisQuat);		// Convert Quaternion Into Matrix3fT
		Matrix3fMulMatrix3f(ThisRot, LastRot);				   // Accumulate Last Rotation Into This One
		Matrix4fSetRotationFromMatrix3f(Transform, ThisRot);
		ConMat();
	}	
	else if(e->Button == ::MouseButtons::Right) //Right mouse Button Translation
	{
		ThisRot=LastRot;
		TranslX+=(-MousePt->X+e->X);
		TranslY+=(MousePt->Y-e->Y);
		MousePt->X=(GLfloat) e->X;
		MousePt->Y=(GLfloat) e->Y;
	};

	//Render normally without change or redraw ellispoids if the 
	//selected atom has deactivated interaction showing.
	if(SelectedAtomsInteractionsFlag) Render();	
	else
	{
		UpdateDisplay(UpdateType::BONDS_LIKE);
		UpdateDisplay(UpdateType::ELLIPSOIDS);
	};

  };
}

/**
* @brief When a mouse button is down but not yet up.
* @param sender The usual sender object for Event Handlers.
* @param e The usual EventArguments object for Event Handlers.
*/
void OpenGLForm::COpenGL::OnMouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	//Prepare Context menu empting it
	this->ContextMenuStrip=nullptr;

	//Update the selection in the Atom's Table
	if(SelectFlag)
	{
		SelectedAtom=ToBeSelectedAtom;
		((SpinachGUI::MainForm^)mainForm)->dataGridView1->ClearSelection();
		if(ToBeSelectedAtom>-1) 
			((SpinachGUI::MainForm^)mainForm)->dataGridView1->Rows[SelectedAtom]->Selected=true;		 

	};

	//Save the initial point when the mouse clicks
	MousePt->X  = (GLfloat) e->X;
	MousePt->Y  = (GLfloat) e->Y;


	if (e->Button == ::MouseButtons::Left || e->Button == ::MouseButtons::Right) 
	{
		//Create context menu
		if((ToBeSelectedAtom>-1) &&  e->Button == ::MouseButtons::Right)
		{
			//Cleaning the contect menu
			this->ContextMenuStrip=((SpinachGUI::MainForm^)mainForm)->OpenGLcontextMenu;
			this->ContextMenuStrip->Items->Clear();

			//Creating the menu items of the context menu
			for each(int i in MainICollection->Keys)
			{
				//Create items only for the selected atom and no bond ones
				if(((Interaction^)MainICollection[i])->IntKind!=InteractionKind::CBond &&
					(((Interaction^)MainICollection[i])->A->getID()==ToBeSelectedAtom+1 ||
					((Interaction^)MainICollection[i])->B->getID()==ToBeSelectedAtom+1))
				{
					this->ContextMenuStrip->Items->Add("Edit " + ((Interaction^)MainICollection[i])->IntKind.ToString()+
						" (ID: "+((Interaction^)MainICollection[i])->getID()+")");
					this->ContextMenuStrip->Items[this->ContextMenuStrip->Items->Count-1]->Tag=
						((Interaction^)MainICollection[i])->getID();
				};
			};
		};

		//For faster rendering
		LowRes=true; 
		//Arcball initial position
		ArcBall->click(ArcMargin+MousePt);
	};

	//Update roation matrix
	LastRot = ThisRot;	

	//Good for checking if mouse moved during the mouse down
	OnMouseMove(sender,e);

	//For updating the mouse icon
	((SpinachGUI::MainForm^)mainForm)->OnMouseDownOpenGL(sender,e);
}

/**
* @brief When a mouse button is up after down.
* @param sender The usual sender object for Event Handlers.
* @param e The usual EventArguments object for Event Handlers.
*/
void OpenGLForm::COpenGL::OnMouseUp(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	//Back to high resolution
	LowRes=false;

	//For updating the mouse icon
	((SpinachGUI::MainForm^)mainForm)->OnMouseUpOpenGL(sender,e);

	//Dangerous crashes the program, but good for mesuring full speed of FPS when you click
	//while(true) Render();
}

#pragma endregion

///// Setting OpenGL System/////////////
#pragma region Setting OpenGL System
/**
* @brief Setting Pixel Format for OpenGL.
* @param hdc Windows machinery object.
*/
GLint OpenGLForm::COpenGL::MySetPixelFormat(HDC hdc)
{
	static	PIXELFORMATDESCRIPTOR pfd=				// pfd Tells Windows How We Want Things To Be
	{
		sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
		1,											// Version Number
		PFD_DRAW_TO_WINDOW |						// Format Must Support Window
		PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
		PFD_DOUBLEBUFFER,							// Must Support Double Buffering
		PFD_TYPE_RGBA,								// Request An RGBA Format
		32,										    // Select Our Color Depth
		0, 0, 0, 0, 0, 0,							// Color Bits Ignored
		0,											// No Alpha Buffer
		0,											// Shift Bit Ignored
		0,											// No Accumulation Buffer
		0, 0, 0, 0,									// Accumulation Bits Ignored
		32,											// 16Bit Z-Buffer (Depth Buffer)  
		0,											// No Stencil Buffer
		0,											// No Auxiliary Buffer
		PFD_MAIN_PLANE,								// Main Drawing Layer
		0,											// Reserved
		0, 0, 0										// Layer Masks Ignored
	};

	GLboolean StereoSupport=GL_FALSE;
	glGetBooleanv(GL_STEREO,&StereoSupport);
	if(StereoSupport) 	pfd.dwFlags |= PFD_STEREO;

	GLint  iPixelFormat; 

	// get the device context's best, available pixel format match 
	if((iPixelFormat = ChoosePixelFormat(hdc, &pfd)) == 0)
	{
		MessageBox::Show("ChoosePixelFormat Failed");
		return 0;
	}

	// make that match the device context's current pixel format 
	if(SetPixelFormat(hdc, iPixelFormat, &pfd) == FALSE)
	{
		MessageBox::Show("SetPixelFormat Failed");
		return 0;
	}

	if((m_hglrc = wglCreateContext(m_hDC)) == NULL)
	{
		MessageBox::Show("wglCreateContext Failed");
		return 0;
	}

	if((wglMakeCurrent(m_hDC, m_hglrc)) == NULL)
	{
		MessageBox::Show("wglMakeCurrent Failed");
		return 0;
	}
	//wglswapI
	
	return 1;
}

/**
* @brief Setting Forn Format for OpenGL.
* @param hdc Windows machinery object.
*/
GLint OpenGLForm::COpenGL::MySetFontFormat(HDC hdc)
{
 
 // Set up the Font characteristics
 HFONT hFont;
 LOGFONT logfont={};

 logfont.lfHeight = -20;
 logfont.lfWidth = 0;
 logfont.lfEscapement = 0;
 logfont.lfOrientation = 0;
 //logfont.lfWeight = FW_BOLD;
 logfont.lfItalic = FALSE;
 logfont.lfUnderline = FALSE;
 logfont.lfStrikeOut = FALSE;
 logfont.lfCharSet = ANSI_CHARSET;
 logfont.lfOutPrecision = OUT_DEFAULT_PRECIS;
 logfont.lfClipPrecision = CLIP_DEFAULT_PRECIS;
 logfont.lfQuality = DEFAULT_QUALITY;
 logfont.lfPitchAndFamily = DEFAULT_PITCH;
 wcscpy_s(logfont.lfFaceName,L"Arial");
 //_tcsncpy_s(logfont.lfFaceName, LF_FACESIZE, _T("Arial"), 7);  
 // Create the font and display list
 hFont = CreateFontIndirect(&logfont);
 
 SelectObject (hdc, hFont);
 // Create display lists for glyphs 0 through 128
 nFontList = glGenLists(128);

 wglUseFontBitmaps(hdc, 0, 128, nFontList);
 DeleteObject(hFont); // Don’t need original font anymore
 // Black Background
 //glClearColor(0.0f, 0.0f, 0.0f, 1.0f );

 return 1;
}

/**
*  @brief Setting Lighting(artistic) Format for OpenGL.
*/
System::Void  OpenGLForm::COpenGL::MySetLightingFormat(System::Void)
{
  GLfloat  ambient               = (GLfloat) 0.20;
  GLfloat  diffuse               = (GLfloat) 0.40;
  GLfloat  specular              = (GLfloat) 0.30;
  GLfloat ambientLight[] = {ambient/10.0f,ambient/10.0f,ambient/10.f,1.0};
  GLfloat diffuseLight[] = {diffuse,diffuse,diffuse,1.0};
  GLfloat specularLight[] = {specular/4.0f,specular/4.0f,specular/4.0f,1.0};

  // setting the lighting parameters
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT2, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT2, GL_SPECULAR, specularLight);
  GLfloat position[] = { 0.0f, 0.0f, -1.0f, 1.0f};
  GLfloat pos[]      = { 20.0f, 0.0f, 200.0f, 1.0f};

  glLightfv(GL_LIGHT0, GL_POSITION, position);
  glLightfv(GL_LIGHT2, GL_POSITION, pos);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT2);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularLight);

}

/**
* @brief  Initialization of parameters of OpenGL.
*/
bool OpenGLForm::COpenGL::InitGL(GLvoid)					// All setup for opengl goes here
{

	glShadeModel(GL_SMOOTH);							// Enable smooth shading
	//glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor((GLfloat)BackgroundColor[0], (GLfloat) BackgroundColor[1],
		         (GLfloat) BackgroundColor[2],1.0f);    // White background
	
	glEnable(GL_DEPTH_TEST);							// Enables depth testing
	glDepthFunc(GL_LEQUAL);								// The type of depth testing to do
	glClearDepth(1.0);									// Depth buffer setup
	glClearStencil(0); 
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really nice perspective calculations
	
	//More anti-aliasing, not sure if it makes more difference
	glEnable (GL_LINE_SMOOTH);
	glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable (GL_POLYGON_SMOOTH);
	glHint (GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	glDepthMask(GL_TRUE);
	glMaterialf(GL_FRONT, GL_SHININESS, 10.0f);         //Create a dictionary of colours
	glEnable(GL_NORMALIZE);

	glEnable(GL_LIGHTING);											// Enable Lighting
	glEnable(GL_LIGHT2);
	glEnable(GL_LIGHT0);	
	//glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);  // Enable Default Light
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_FOG);
	glFogi(GL_FOG_MODE, GL_EXP);
    glFogf(GL_FOG_DENSITY, 0.005f);
	//glFogf(GL_FOG_START,   30.0f);
    //glFogf(GL_FOG_END,     200.0f);
	

	quadratic=gluNewQuadric();							// Create A Pointer To The Quadric Object
	gluQuadricDrawStyle(quadratic,GLU_FILL);
	
	//Build the primary display lists
	BuildDisplayLists();
	return TRUE;										// Initialisation went ok
}

/**
* @brief  Swap Buffer in OpenGL panel, usually after finishing rendering.
*/
System::Void OpenGLForm::COpenGL::SwapOpenGLBuffers(System::Void)
{
	SwapBuffers(m_hDC) ;
}

/**
* @brief Resize of the windows containing the OpenGL panel.
* @param parentForm The parent Windows form.
* @param origin Point of origin.
* @param iWidth Width of window.
* @param iHeight Height of window.
*/
System::Void OpenGLForm::COpenGL::ReSizeWindow(Form ^ parentForm, System::Drawing::Point origin, GLsizei iWidth, GLsizei iHeight)
{
	//Setting location and size of the panel 
	this->Location= origin;
	this->Height = iHeight;
	this->Width = iWidth;
	mWidth=iWidth;mHeight=iHeight;

	//Setting the boundaries of ArcBall
	ArcMargin->X=(GLfloat)(-0.05*mWidth);
	ArcMargin->Y=(GLfloat)(-0.05*mHeight);
	ArcBall->setBounds((GLfloat)iWidth+2*ArcMargin->X,(GLfloat) iHeight+2*ArcMargin->Y);
	
	// Setting the size to the final OpenGl panel
	ReSizeGLScene(iWidth, iHeight);
	Render();
}

/**
* @brief Resize the OpenGL panel.
* @param width Width of OpenGL panel.
* @param height Height of OpenGL panel.
*/
GLvoid OpenGLForm::COpenGL::ReSizeGLScene(GLsizei width, GLsizei height)		
{
	if (height==0)	height=1;							// Prevent A Divide By Zero By Making Height Equal One
	width=System::Math::Max(1,width); height=System::Math::Max(1,height);
	glViewport(0,0,width,height);						// Reset The Current Viewport
	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glLoadIdentity();									// Reset The Projection Matrix
	
	// Calculate The Aspect Ratio Of The Window
	gluPerspective(45.0f,(GLfloat)width/(GLfloat)height,0.1f,300.0f);
	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Modelview Matrix
}
#pragma endregion

/////Bitmap export////////////////////
/**
* @brief Export the contens of the  OpenGL panel in a bitmap file.
*/
System::Drawing::Bitmap^  OpenGLForm::COpenGL::Save_BMP()
{
	/* //Old Method not working so good
		static int id = 0;
	int xs = mWidth;
	int ys = mHeight-1;
	System::Drawing::Bitmap^ nbitmap=gcnew System::Drawing::Bitmap(xs+1,ys);
	GLubyte *buf = new GLubyte [3*(xs+2)*(ys+2)];

	if (buf==NULL) {return nbitmap;}

	bool EyeFlip = false;
	Render();
	//glReadBuffer((!EyeFlip)?GL_BACK:GL_BACK_RIGHT);
	glReadBuffer((!EyeFlip)?GL_FRONT:GL_FRONT_RIGHT);
	glReadPixels(0,2,xs, ys,GL_RGB,GL_UNSIGNED_BYTE, buf);//76

	int ti;
	for(int i=0;i<5;i++)//Not Perfect but not bad solution
	{
	int temp=(xs+ys/2*(xs+1))*3-i*ys/2;
	if(buf[temp]!=255 || buf[temp+1]!=255 || buf[temp+2]!=255) ti=i;
	}

	for(int j=0;j<ys;j++){
		for(int i=0;i<(xs+1);i++){
			nbitmap->SetPixel(i,ys-1-j,System::Drawing::Color::FromArgb(255,buf[(i+j*(xs+1))*3-ti*j],
																			buf[(i+j*(xs+1))*3+1-ti*j],
																			buf[(i+j*(xs+1))*3+2-ti*j]));
		}
	};

	return nbitmap;
	free(buf);
	buf=NULL;
	*/

	//Create the bitmap file
	System::Drawing::Bitmap^ bmp = gcnew System::Drawing::Bitmap(((SpinachGUI::MainForm^)mainForm)->splitContainer2->Panel1->ClientSize.Width, 
		                                                         ((SpinachGUI::MainForm^)mainForm)->splitContainer2->Panel1->ClientSize.Height);
	//Create the datat file
    System::Drawing::Imaging::BitmapData^ data =
                bmp->LockBits(((SpinachGUI::MainForm^)mainForm)->splitContainer2->Panel1->ClientRectangle, 
				              System::Drawing::Imaging::ImageLockMode::WriteOnly, 
							  System::Drawing::Imaging::PixelFormat::Format24bppRgb);

    //Export OpenGL pixel to the data file        
	glReadPixels(0, 0, mWidth , mHeight, GL_BGR_EXT, GL_UNSIGNED_BYTE, data->Scan0.ToPointer());
 
    //Move data to bmp and rotate it
	bmp->UnlockBits(data);
    bmp->RotateFlip(RotateFlipType::RotateNoneFlipY);
	
    return bmp;
}

/**
* @brief Draws the molecular magnetic susceptibility tensor as an ellipsoid at the origin of the
* atom coordinates. The largest semi-axis is SUSfactor*validminCB, as for the electron tensors.
*/
System::Void OpenGLForm::COpenGL::Render_Susceptibility()
{
	if (!SUSFlag || SystemModel->Susceptibility==nullptr) return;
	Tensor^ chi=SystemModel->Susceptibility;
	double maxeigen=0.0;
	for (int i=0; i<3; i++) maxeigen=System::Math::Max(maxeigen, System::Math::Abs(chi->eigen[i]));
	double validminCB=3.0f; //standard, as for the electron tensors
	SUSfinfactor=SUSfactor*validminCB/(maxeigen>epsilo ? maxeigen : 1.0f);
	Render_Atom_Properties(0.0, 0.0, 0.0, chi->dcm, chi->eigen, SUSfinfactor, SUScolor, 0.6f, -1);
}
