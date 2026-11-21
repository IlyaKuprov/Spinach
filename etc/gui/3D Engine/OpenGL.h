#pragma once
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <vector>
#include "MetaBall.h"
#include "ArcBall.h"
#include "../System Core/Model.h"


using namespace System::Windows::Forms;
using namespace std;
using namespace System::IO;
using namespace System::Collections;

class GLUquadric{};//for stop linking warning 

/** \enum UpdateType
* Used to separate different kind of rebuilding call lists
*/
public enum class UpdateType {
	ALL,
	BONDS_LIKE,
	ELLIPSOIDS,
	ATOMS
};

namespace OpenGLForm 
{

	/** \class COpenGL
	* Main COpenGL class for the rendering OpenGL windows, inherited from  windows panel
	*/
	public ref class COpenGL: public System::Windows::Forms::Panel
	{
	public:
		COpenGL(Form ^ ,System::Drawing::Point , GLsizei , GLsizei);     
		~COpenGL(System::Void);											
		System::Void Home(System::Void);
		System::Void Reset(System::Void);
		System::Void ZoomtoFillScreen(System::Void);
		System::Void Render(System::Void);
		System::Void UpdateDisplay(UpdateType type);
		System::Void TotalRender(System::Void);
		System::Void Calc_eye(int, bool);
		System::Void ReSizeWindow(Form ^,System::Drawing::Point , GLsizei , GLsizei);
		System::Void SwapOpenGLBuffers(System::Void);
		System::Drawing::Bitmap^ Save_BMP();
		void OnMouseWheel(System::Object^ , System::Windows::Forms::MouseEventArgs^  );
	    System::Void RotateToPlane(char , char );
		System::Void Travelling(System::Void);

	private:
		System::Void Render_Atoms(System::Void);
		System::Void Render_SelectedAtoms(System::Void);
		System::Void BuildDisplayLists(System::Void);
		System::Void Render_Bonds(System::Void);
		System::Void Render_Atom_Properties(double, double, double, 
			                                ObservableCollection<double>^ , 
											ObservableCollection<double>^, 
											double, double *, double, int );
		System::Void Render_Navigation();
        System::Void Render_Floor(int, int);
		System::Void Render_Electrons();
		System::Void Render_Susceptibility();
		System::Void Render_MetaBalls();
		System::Void Render_IDText(double x, double y, double z, double radius, int ID);
		System::Void Render_Transparent_Atom_Properties(int);
		System::Void Render_Transparency(double, bool);
		System::Void ConMat(System::Void);
		System::Void Init_Metaballs(System::Void);
		GLvoid ReSizeGLScene(GLsizei , GLsizei );	
		System::Void  MySetLightingFormat(System::Void);
		System::Void SelectionBuffer(GLfloat ,GLfloat );
		void OnMouseEnter(System::Object^ , System::EventArgs^  );
		void OnMouseLeave(System::Object^ , System::EventArgs^  );
		void OnMouseMove(System::Object^  , System::Windows::Forms::MouseEventArgs^  );
		void OnMouseDown(System::Object^  , System::Windows::Forms::MouseEventArgs^  );
		void OnMouseUp(System::Object^  , System::Windows::Forms::MouseEventArgs^ );
		void Del_KeyDown(System::Object^ , System::Windows::Forms::KeyEventArgs^  ) ;
		GLint MySetPixelFormat(HDC );
		GLint MySetFontFormat(HDC );
		bool InitGL(GLvoid);
		System::ComponentModel::BackgroundWorker^  OpenGLbackgroundWorker;

		

     public:
		double CSfactor, HCfactor, GTfactor, CBfactor, SUSfactor;				///< Multiplication factor for the size of ellipsoids
		double JCfactor, QCfactor, ECfactor, ZFfactor;				///< Multiplication factor for the size of ellipsoids
		double CSfinfactor, HCfinfactor, GTfinfactor, CBfinfactor, SUSfinfactor;	///< Final multiplication factor for the size of ellipsoids
		double JCfinfactor, QCfinfactor, ECfinfactor, ZFfinfactor;	///< Final multiplication factor for the size of ellipsoids
		double *CScolor, *HCcolor, *GTcolor, *JCcolor, *SUScolor;				///< RGB Matrices for the color of interaction kind
		double *CBcolor, *QCcolor, *ECcolor, *ZFcolor;				///< RGB Matrices for the color of interaction kind
		bool CSFlag, HCFlag, GTFlag, JCFlag, CBFlag, SUSFlag;				///< Flags for showing kind of interactions
		bool DCFlag, QCFlag, ECFlag, ZFFlag;						///< Flags for showing kind of interactions
		double *BackgroundColor;									///< RGB Matrix for the color
		bool AxesFlag;												///< Flag for showing navigation axes
		bool AtomIDFlag;											///< Flag for showing ID of all atoms
		bool EigenAxesFlag;											///< Flag for showing interaction axes
		bool EigenEllipsoidsFlag;									///< Flag for showing interaction ellipsoids
		bool SelectedAtomsInteractionsFlag;							///< Flag for showing the interaction of the selected atom
		bool ReferenceFrameAxesFlag;								///< Flag for showing current reference frame axis of interaction
		bool focused;												///< Flag used for context menu showing
		bool stereoview;											///< Flag for stero viewing
		Model^ SystemModel;											///< Current system model
		int ToBeSelectedAtom;										///< Atom to be selected

	private:
		double last_update;					///< Time of last redraw
		int loop_count;						///< Number of times through display loop since redraw
		double FPSrate;						///< FPS Counter
		GLuint MySphere;					///< Non Default sphere for Call list
		GLuint SphereHighRes;				///< Sphere with high Resolution for Call list
		GLuint SphereLowRes;				///< Sphere with low Resolution for Call list
		GLuint SelectedSphereHighRes;		///< Selected Sphere with high Resolution for Call list
		GLuint SelectedSphereLowRes;		///< Selected Sphere with low Resolution for Call list
		GLuint CylinderHighRes;				///< Cylinder with high Resolution for Call list
		GLuint CylinderLowRes;				///< Cylinder with low Resolution for Call list
		GLuint FullCylinderHighRes;			///< Cylinder with bases with high Resolution for Call list
		GLuint FullCylinderLowRes;			///< Cylinder with bases with low Resolution for Call list
		GLuint AllAtoms;					///< Rendering All atoms for Call list
		GLuint AllBonds_Like;				///< Rendering All bonds and j-couplings for Call list
		GLuint All_Ellipsoids;				///< Rendering All ellipsoidal interactions for Call list
		GLuint Navigation;					///< Rendering navigation axes for Call list
		GLuint WiredSphereLowRes;			///< Wired Sphere with high Resolution for Call list
		GLuint WiredSphereHighRes;			///< Wired Sphere with low Resolution for Call list
		HDC m_hDC;							///< For OpenGL window panel
		HGLRC m_hglrc;						///< For OpenGL window panel
		GLUquadricObj *quadratic;			///< Quadratic object ofr drawing of spheres
		GLuint nFontList;					///< Fonts for OpenGL text
		Point2fT^ MousePt;					///< Current Mouse pointer coordinates
		Point2fT^ ArcMargin;				///< Margins for the Arcball algorithm
        Matrix4fT^ Transform;				///< Current matrix for transformation
		GLfloat* TransformMat;				///< Current matrix for transformation used for migration to OpenGL format
		GLfloat* SaveMat;					///< Matrix for saving current frame in travelling
        Matrix3fT ^LastRot, ^ThisRot;		///< Matrix with last rotation and current rotation
		ArcBallT^ ArcBall;					///< ArdBall object for 3D roation using the mouse
		int mWidth,mHeight;					///< Width and Height of the OpenGL Window
		bool SelectFlag;					///< Flag for selecting atoms
		int SelectedAtom;					///< ID of Selected Atom- zero indexed as coming from OpenGL machinery
		double Zoom;						///< Value of Zooming -distance for the center, focal pont
		double TranslX,TranslY;				///< Translation X, Y of model 
		Form^ mainForm;						///< Related to MainForm of Application
		CUBE_GRID *cubeGrid;				///< Cube grid for use from metaballs object
		bool LowRes;						///< Flag for using low resolution during rotation
		int Spherediv;						///< Number of quad used for spheres
	};
}

