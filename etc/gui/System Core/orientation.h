#pragma once
#include <Dense>
#include <Geometry>

using namespace Eigen;


namespace SpinachGUI {

#pragma region "Orientation related classes"

	/** \class  EulerAngles
	*  Class for storing a 3 dimentional rotation using Euler angles for parameterization
	*/
	public ref class EulerAngles {

	public:
		/**
		* @brief Contrustor setting all angles to zero.
		*/
		EulerAngles() : alpha(0.0),beta(0.0),gamma(0.0) { }

		/**
		* @brief Contrustor setting all angles.
		* @param _alpha Value of alpha angle.
		* @param _beta Value of beta angle.
		* @param _gamma Value of gamma angle.
		*/
		EulerAngles(double _alpha,double _beta, double _gamma): alpha(_alpha),beta(_beta),gamma(_gamma) {}

		/**
		* @brief Setter for all angles
		* @param _alpha Value of alpha angle.
		* @param _beta Value of beta angle.
		* @param _gamma Value of gamma angle.
		*/
		void Set(double _alpha,double _beta, double _gamma) {alpha=_alpha;  beta= _beta; gamma=_gamma; }
		void Normalize() {};

		/**
		* @brief Getter of normalized eulerangles object.
		* @result Returns normalized eulerangles object.
		*/
		EulerAngles^ Normalized() { return this;}

		/**
		* @brief Operator for checking if two object are equal
		*/
		bool operator==(EulerAngles o)  { return alpha == o.alpha &&  beta == o.beta && gamma == o.gamma;}

		/**
		* @brief Operator for checking if two object are not equal
		*/
		bool operator!=(EulerAngles o) { return alpha != o.alpha ||  beta != o.beta || gamma != o.gamma;}

		double alpha,beta,gamma; ///< Alpha beta gamma angles
	};

	/** \class  Matrix3x3
	*  Class for storing a 3 dimentional rotation using a 3x3 matrix, DCM for parameterization
	*/
	public ref class Matrix3x3 {

	public:

		/**
		* @brief Contrustor setting all elements to zero and initialize Eigen lib matrix.
		*/
		Matrix3x3() : xx(0),xy(0),xz(0),yx(0),yy(0),yz(0),zx(0),zy(0),zz(0) { matrix3d=new  Eigen::Matrix3d();}


		/**
		* @brief Getter of normalized  Matrix3x3 object.
		* @result Returns normalized  Matrix3x3 object.
		*/
		Matrix3x3^ Normalized() { return this;}

		void Normalize() {};
		void toEigen();
		void fromEigen();

		double xx,xy,xz,yx,yy,yz,zx,zy,zz; ///< Elements of the matrix
		Eigen::Matrix3d* matrix3d;				   ///< Eigen lib matrix
	};

	/** \class Quaternion
	*  Class for storing a 3 dimentional rotation using quaternion for parameterization
	*/
	public ref class Quaternion {

	public:

		/**
		* @brief Contrustor setting all elements to zero and initialize Eigen lib quaternion.
		*/
		Quaternion() : x(0),y(0),z(0),w(0) {quaterniond= new Eigen::Quaterniond();}

		/**
		* @brief Contrustor setting all elements to zero and initialize Eigen lib quaternion.
		* @param _w Value of quaternion w parameter.
		* @param _x Value of quaternion x parameter.
		* @param _y Value of quaternion y parameter.
		* @param _z Value of quaternion z parameter.
		*/
		Quaternion(double _w, double _x, double _y, double _z) : x(_x),y(_y),z(_z),w(_w) {quaterniond= new Eigen::Quaterniond();}

		/**
		* @brief Copy Contrustor.
		*/
		Quaternion(Quaternion ^quat):x(quat->x),y(quat->y),z(quat->z),w(quat->w){quaterniond=new Eigen::Quaterniond();} 


		Quaternion^ Normalized();
		void Normalize() {};
		void toEigen();
		void fromEigen();

		double x,y,z,w;					 ///< Elements of the quaternion
		Eigen::Quaterniond* quaterniond; ///< Eigen lib quaternion

	};

	/** \class  AngleAxis
	*  Class for storing a 3 dimentional rotation using angle axis for parameterization
	*/
	public ref class AngleAxis {

	public:

		/**
		* @brief Contrustor setting all elements to zero and initialize Eigen lib angleaxis.
		*/
		AngleAxis() : x(0),y(0),z(0),angle(0) {angleaxisd=new Eigen::AngleAxisd(); }

		/**
		* @brief Contrustor setting all elements to zero and initialize Eigen lib angleaxis.
		* @param _x Value of angleaxis x parameter.
		* @param _y Value of angleaxis y parameter.
		* @param _z Value of angleaxis z parameter.
		* @param _angle Value of angleaxis angle parameter.
		*/
		AngleAxis(double _x, double _y, double _z, double _angle): x(_x),y(_y),z(_z),angle(_angle) {angleaxisd=new Eigen::AngleAxisd();}

		/**
		* @brief Copy Contrustor.
		*/
		AngleAxis(AngleAxis ^angleax):x(angleax->x),y(angleax->y),z(angleax->z),angle(angleax->angle){angleaxisd=new Eigen::AngleAxisd();} 


		void Normalize(){};
		AngleAxis^ Normalized();
		void toEigen();
		void fromEigen();

		double x,y,z,angle;				///< Elements of the angleaxis
		Eigen::AngleAxisd* angleaxisd;  ///< Eigen lib angleaxis

	};

	/** \class  Type
	* Enumeration of the four conventions on storing rotations
	*/
	public enum  Type {
		EULER,
		ANGLE_AXIS,
		QUATERNION,
		DCM
	};

#pragma endregion

	/** \class  Orientation
	*  Class for storing a 3 dimentional rotation using different conventions
	*/
	public ref  class Orientation {
	public:
		//        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

#pragma region "Constructors, Destructor, getters"

		/**
		* @brief Constructs the identity rotation.
		*/
		Orientation(): mType(QUATERNION) {myQuaternion=gcnew Quaternion(1,0,0,0);Invariant();}

		/**
		* @brief Constructs an orientation from euler angles.
		* @param ea Euler angles object.
		*/
		Orientation(EulerAngles^ ea) : mType(EULER) {mEuler=ea;Invariant();};


		/**
		* @brief Constructs an orientation from a matrix.
		* @param m DCM matrix object.
		*/
		Orientation(Matrix3x3^ m): mType(DCM) {myMatrix=m;Invariant();};

		/**
		* @brief Constructs an orientation from a quaternion.
		* @param q Quaternion object.
		*/
		Orientation(Quaternion^ q): mType(QUATERNION)  {myQuaternion=q;Invariant();};

		/**
		* @brief Constructs an orientation from an AngleAxis.
		* @param aa Angle axis object.
		*/
		Orientation(AngleAxis^ aa): mType(ANGLE_AXIS)  {myAngleAxis=aa;Invariant();};

		/**
		* @brief Destructor
		*/
		~Orientation(){};


		///Normalise the contained rotation.
		Orientation^ Normalized() ;

		/**
		* @brief Getter of the type of orietation
		* @result Returns the convention that is being used to store the rotation.
		*/
		SpinachGUI::Type GetType() {return mType;}

		//Getter to different forms
		EulerAngles^ GetAsEuler();
		Matrix3x3^ GetAsDCM();
		AngleAxis^  GetAsAngleAxis();
		Quaternion^ GetAsQuaternion();
#pragma endregion


#pragma region "Rotation Converions"

		/**
		* @brief Convertion function from Euler angles to Euler angles object.
		* @param rot Euler angles object.
		* @result Returns the converted to Euler angles object.
		*/
		EulerAngles^ ConvertToEuler(EulerAngles^ rot){return rot;};

		/**
		* @brief Convertion function from DCM matrix to Euler angles object.
		* @param rot DCM matrix  object.
		* @result Returns the converted to Euler angles object.
		*/
		EulerAngles^ ConvertToEuler(Matrix3x3^   rot){return ConvertToEuler(ConvertToQuaternion(rot));};

		EulerAngles^ ConvertToEuler(Quaternion^ rot);

		/**
		* @brief Convertion function from AngleAxis to Euler angles object.
		* @param rot AngleAxis object.
		* @result Returns the converted to Euler angles object.
		*/
		EulerAngles^ ConvertToEuler(AngleAxis^   rot){return ConvertToEuler(ConvertToQuaternion(rot));};

		/**
		* @brief Convertion function from Euler angles to DCM matrix object.
		* @param rot Euler angles object.
		* @result Returns the converted to DCM matrix object.
		*/
		Matrix3x3^    ConvertToDCM(EulerAngles^ rot){return ConvertToDCM(ConvertToQuaternion(rot));};

		/**
		* @brief Convertion function from DCM matrix to DCM matrix object.
		* @param rot  DCM matrix object.
		* @result Returns the converted to DCM matrix object.
		*/
		Matrix3x3^    ConvertToDCM(Matrix3x3^   rot){return rot;};

		Matrix3x3^    ConvertToDCM(Quaternion^ );

		/**
		* @brief Convertion function from AngleAxis to DCM matrix object.
		* @param rot  AngleAxis object.
		* @result Returns the converted to DCM matrix object.
		*/
		Matrix3x3^    ConvertToDCM(AngleAxis^   rot){return ConvertToDCM(ConvertToQuaternion(rot));};

		Quaternion^ ConvertToQuaternion(EulerAngles^ );
		Quaternion^ ConvertToQuaternion(Matrix3x3^   );

		/**
		* @brief Convertion function from Quaternion to Quaternion object.
		* @param rot Quaternion object.
		* @result Returns the converted to Quaternion object.
		*/
		Quaternion^ ConvertToQuaternion(Quaternion^ rot){return rot;};

		Quaternion^ ConvertToQuaternion(AngleAxis^   );


		/**
		* @brief Convertion function from EulerAngles to AngleAxis object.
		* @param rot DCM matrix object.
		* @result Returns the converted to AngleAxis object.
		*/
		AngleAxis^  ConvertToAngleAxis(EulerAngles^ rot){return ConvertToAngleAxis(ConvertToQuaternion(rot));}

		/**
		* @brief Convertion function from DCM matrix to AngleAxis object.
		* @param rot EulerAngles object.
		* @result Returns the converted to AngleAxis object.
		*/
		AngleAxis^  ConvertToAngleAxis(Matrix3x3^ rot){return ConvertToAngleAxis(ConvertToQuaternion(rot));};

		AngleAxis^  ConvertToAngleAxis(Quaternion^ );

		/**
		* @brief Convertion function from AngleAxis to AngleAxis object.
		* @param rot AngleAxis object.
		* @result Returns the converted to AngleAxis object.
		*/
		AngleAxis^  ConvertToAngleAxis(AngleAxis^ rot){return rot;};

#pragma endregion

	private:

		//Main variables of each orientation
		SpinachGUI::Type mType;
		EulerAngles^ mEuler;
		AngleAxis^ myAngleAxis;
		Quaternion^ myQuaternion;
		Matrix3x3^ myMatrix;

		//Eigen-based Variables
		Eigen::AngleAxisd* mAngleAxis;
		Eigen::Quaterniond* mQuaternion;
		Eigen::Matrix3d* mMatrix;

		void Invariant();
	};
}

