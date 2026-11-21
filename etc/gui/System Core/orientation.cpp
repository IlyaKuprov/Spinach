#include "../StdAfx.h"
#include "orientation.h"

using namespace SpinachGUI;

/**
*\file orientation.cpp
*/

#pragma region "Matrix3x3, Quaternion, Angleaxis functions"

/**
* @brief Function for setting the Eigen lib object values to the elements of the matrix.
*/
void Matrix3x3::toEigen()
{
	(*matrix3d)(0)=xx;
	(*matrix3d)(1)=xy;
	(*matrix3d)(2)=xz;
	(*matrix3d)(3)=yx;
	(*matrix3d)(4)=yy;
	(*matrix3d)(5)=yz;
	(*matrix3d)(6)=zx;
	(*matrix3d)(7)=zy;
	(*matrix3d)(8)=zz;
};

/**
* @brief Function for setting the elements of the matrix to the Eigen lib object values. 
*/
void Matrix3x3::fromEigen()
{
	xx=(*matrix3d)(0);
	xy=(*matrix3d)(1);
	xz=(*matrix3d)(2);
	yx=(*matrix3d)(3);
	yy=(*matrix3d)(4);
	yz=(*matrix3d)(5);
	zx=(*matrix3d)(6);
	zy=(*matrix3d)(7);
	zz=(*matrix3d)(8);
};

/**
* @brief Function for setting the Eigen lib object values to the elements of the quaternion.
*/
void SpinachGUI::Quaternion::toEigen()
{
	quaterniond->x()=x;
	quaterniond->y()=y;
	quaterniond->z()=z;
	quaterniond->w()=w;
};

/**
* @brief Function for setting the elements of the quaternion to the Eigen lib object values. 
*/
void SpinachGUI::Quaternion::fromEigen()
{
	x=quaterniond->x();
	y=quaterniond->y();
	z=quaterniond->z();
	w=quaterniond->w();
};

/**
* @brief Function for setting the Eigen lib object values to the elements of the Angle axis.
*/
void SpinachGUI::AngleAxis::toEigen()
{
	angleaxisd->axis().x()=x;
	angleaxisd->axis().y()=y;
	angleaxisd->axis().z()=z;
	angleaxisd->angle()=angle;
};

/**
* @brief Function for setting the elements of the Angle axis to the Eigen lib object values.
*/
void SpinachGUI::AngleAxis::fromEigen()
{
	x=angleaxisd->axis().x();
	y=angleaxisd->axis().y();
	z=angleaxisd->axis().z();
	angle=angleaxisd->angle();
};

/**
* @brief Getter of normalized quaternion object.
* @result Returns normalized quaternion object.
*/
SpinachGUI::Quaternion^ SpinachGUI::Quaternion::Normalized() 
{
	//Pass to the Eigen library object
	this->toEigen();

	//Normalize
	quaterniond->normalize();

	//Take it from Eigen library object and return
	this->fromEigen();
	return this;
}

/**
* @brief Getter of normalized AngleAxis object.
* @result Returns normalized AngleAxis object.
*/
SpinachGUI::AngleAxis^ SpinachGUI::AngleAxis::Normalized()
{
	//Pass to the Eigen library object
	this->toEigen();

	//Normalize
	angleaxisd->axis().normalize();

	//Take it from Eigen library object and return
	this->fromEigen();
	return this;
}

#pragma endregion

#pragma region "Orientation various functions"

/**
* @brief Function for giving the normalized orientation.
* @result Returns a new orientation which is the normalized current orrientation.
*/
Orientation^ Orientation::Normalized()  
{

	//Check what type of orientation and then create a new one which is normalized
	if(mType == EULER) 	return gcnew Orientation(mEuler);
	if(mType == DCM) 	return gcnew Orientation(myMatrix);
	if(mType == QUATERNION) return gcnew Orientation(myQuaternion->Normalized());	
	if(mType == ANGLE_AXIS) return gcnew Orientation(myAngleAxis->Normalized());

	//In case of unknown type, shows message and returns a new empty quaternion orientation
	System::Windows::Forms::MessageBox::Show("Orientation with Unknown Type, returning Quaternion Orientation",
		"Error",  System::Windows::Forms::MessageBoxButtons::OK, 
		System::Windows::Forms::MessageBoxIcon::Error) ;
	return gcnew Orientation();
}


/**
* @brief Function for checking if elements are zero.
*/
void Orientation::Invariant() 
{
	//In case of unknown type, Shows message
	if(!(mType == ANGLE_AXIS || mType == DCM || mType == QUATERNION || mType == EULER)) 
	{
		System::Windows::Forms::MessageBox::Show("mType set to invalid value.",
			"Error",  System::Windows::Forms::MessageBoxButtons::OK, 
			System::Windows::Forms::MessageBoxIcon::Error) ;
	};

	if(mType == ANGLE_AXIS) 
	{
		//Pass to Eigen lib object and set a new vector form the axis
		myAngleAxis->toEigen();
		Eigen::Vector3d axis = myAngleAxis->angleaxisd->axis();

		//Check if the axis is zero give message
		if(axis.x() == 0 && axis.y() == 0 && axis.z() == 0) 
		{
			System::Windows::Forms::MessageBox::Show(
				"An orientation in angle-axis notation has (0,0,0) as it's axis vector",
				"Error",  System::Windows::Forms::MessageBoxButtons::OK, 
				System::Windows::Forms::MessageBoxIcon::Error) ;
		};
	} 
	else if(mType == QUATERNION) 
	{
		//Pass to Eigen lib object
		myQuaternion->toEigen();

		//Check if parameters are zero give message
		if(myQuaternion->quaterniond->x() == 0 && myQuaternion->quaterniond->y() == 0 && 
			myQuaternion->quaterniond->z() == 0 && myQuaternion->quaterniond->w() == 0)
		{
			System::Windows::Forms::MessageBox::Show(
				"An orientation in quaternion notation is (0,0,0,0)",
				"Error",  System::Windows::Forms::MessageBoxButtons::OK, 
				System::Windows::Forms::MessageBoxIcon::Error) ;
		};
	};
}


/**
* @brief Function for returning orientation as Euler angles.
* @result Returns orientation object to Euler angles object.
*/
SpinachGUI::EulerAngles^ Orientation::GetAsEuler() 
{
	//Set new orientation as the current one normalized
	Orientation^ normalized=this->Normalized();

	//Depending type, convert the normalized orientation to Euler angles form
	if(mType == EULER) return ConvertToEuler(normalized->mEuler);				
	if(mType == DCM)   return ConvertToEuler(normalized->myMatrix);				
	if(mType == ANGLE_AXIS) return ConvertToEuler(normalized->myAngleAxis);		
	if(mType == QUATERNION) return ConvertToEuler(normalized->myQuaternion);

	//In case of unknown type, shows message and returns a new empty euler angles orientation
	System::Windows::Forms::MessageBox::Show("Requesting Unknown Orientation Type",
		"Error",  System::Windows::Forms::MessageBoxButtons::OK, 
		System::Windows::Forms::MessageBoxIcon::Error) ;
	return gcnew EulerAngles();
}


/**
* @brief Function for returning orientation as DCM matrix.
* @result Returns orientation object to DCM matrix object.
*/
SpinachGUI::Matrix3x3^ Orientation::GetAsDCM() 
{
	//Set new orientation as the current one normalized
	Orientation^ normalized=Normalized();

	//Depending type, convert the normalized orientation to DCM matrix form
	if(mType == EULER) return ConvertToDCM(normalized->mEuler);				
	if(mType == DCM)   return ConvertToDCM(normalized->myMatrix);				
	if(mType == ANGLE_AXIS) return ConvertToDCM(normalized->myAngleAxis);		
	if(mType == QUATERNION)  return ConvertToDCM(normalized->myQuaternion);

	//In case of unknown type, returns a new empty DCM matrix orientation
	return gcnew Matrix3x3();
}


/**
* @brief Function for returning orientation as Angle axis.
* @result Returns orientation object to angle axis object.
*/
SpinachGUI::AngleAxis^ Orientation::GetAsAngleAxis() 
{
	//Set new orientation as the current one normalized
	Orientation^ normalized=Normalized();

	//Depending type, convert the normalized orientation to angle axis form
	if(mType == EULER) return ConvertToAngleAxis(normalized->mEuler);				
	if(mType == DCM)   return ConvertToAngleAxis(normalized->myMatrix);				
	if(mType == ANGLE_AXIS) return ConvertToAngleAxis(normalized->myAngleAxis);		
	if(mType == QUATERNION) return ConvertToAngleAxis(normalized->myQuaternion);

	//In case of unknown type, shows message and returns a new empty angle axis orientation
	System::Windows::Forms::MessageBox::Show("Requesting Unknown Orientation Type",
		"Error",  System::Windows::Forms::MessageBoxButtons::OK, 
		System::Windows::Forms::MessageBoxIcon::Error) ;
	return gcnew SpinachGUI::AngleAxis();		
}



/**
* @brief Function for returning orientation as quaternion.
* @result Returns orientation object to quaternion object.
*/
SpinachGUI::Quaternion^ Orientation::GetAsQuaternion() 
{
	//Set new orientation as the current one normalized
	Orientation^ normalized=Normalized();

	//Depending type, convert the normalized orientation to quaternion form
	if(mType == EULER) return ConvertToQuaternion(normalized->mEuler);				
	if(mType == DCM)   return ConvertToQuaternion(normalized->myMatrix);				
	if(mType == ANGLE_AXIS) return ConvertToQuaternion(normalized->myAngleAxis);		
	if(mType == QUATERNION) return ConvertToQuaternion(normalized->myQuaternion);

	//In case of unknown type, shows message and returns a new empty quaternion orientation
	System::Windows::Forms::MessageBox::Show("Requesting Unknown Orientation Type",
		"Error",  System::Windows::Forms::MessageBoxButtons::OK, 
		System::Windows::Forms::MessageBoxIcon::Error) ;
	return gcnew SpinachGUI::Quaternion();
}

#pragma endregion

#pragma region "Orientation Convertion functions"

/**
* @brief Convertion function from  Quaternion to Euler angles  object.
* @param rot Quaternion object.
* @result Returns the converted to Euler angles object.
*/
SpinachGUI::EulerAngles^ Orientation::ConvertToEuler(SpinachGUI::Quaternion^ rot) 
{
	double PI=3.141592653589793238462643;
	Vector3d z_axis=Vector3d(0,0,1);
	Vector3d x_axis=Vector3d(1,0,0);

	//Pass to Eigen lib object
	rot->toEigen();

	//Set z,x axis form the quaternion
	z_axis=*(rot->quaterniond)*z_axis;
	x_axis=*(rot->quaterniond)*x_axis;

	//Find the gamma from the angle of z_axis
	double gamma=atan2(z_axis.y(),z_axis.x());

	//Find beta from the angle of z_axis using the length 
	double beta=atan2(sqrt(z_axis.x()*z_axis.x() + z_axis.y()*z_axis.y()),z_axis.z());

	//Use γ and β to rotate V2 back onto the point 0,0,1
	Quaterniond betaTwist(AngleAxisd(-beta, Vector3d(0,1,0)));
	Quaterniond gammaTwist (AngleAxisd(-gamma,Vector3d(0,0,1)));

	//Calculate x axis form the twists
	x_axis=(betaTwist*gammaTwist)*x_axis;

	//Find alpha from the angle of x_axis
	double alpha = atan2(x_axis.y(),x_axis.x());

	//Convention for alpha angle:Transcate to first quarter
	if(alpha < 0 || alpha >= 2*PI)  alpha = alpha-2*PI*floor(alpha/(2*PI));
	if(alpha >= 2*PI) alpha = 0;

	//Convention for beta angle:Transcate to first two quarters
	if(beta < 0 || beta > PI)     beta = beta-  PI*floor(beta/PI);

	//Convention for gamma angle:Transcate to first quarter
	if(gamma < 0 || gamma >= 2*PI)  gamma = gamma-2*PI*floor(gamma/(2*PI));
	if(gamma >= 2*PI) gamma = 0;

	return gcnew EulerAngles(gamma,beta,alpha);
}

/**
* @brief Convertion function from  Euler angles to Quaternion object.
* @param rot Euler angles object.
* @result Returns the converted to Quaternion object.
*/
SpinachGUI::Quaternion^ Orientation::ConvertToQuaternion(EulerAngles^ rot) 
{
	//Set three new quaternion using each angle axis in the x,y,z axes respectively
	const double alpha=rot->alpha, beta=rot->beta, gamma=rot->gamma;
	Quaterniond q1(AngleAxisd(alpha,Vector3d(0,0,1)));
	Quaterniond q2(AngleAxisd(beta, Vector3d(0,1,0)));
	Quaterniond q3(AngleAxisd(gamma,Vector3d(0,0,1)));

	//Create a new empty quaternion object
	Quaternion^ result=gcnew Quaternion();

	//Set Eigen lib object using the Eigen quaternion form  multipling the three previous quaternions
	*(result->quaterniond)=q3*q2*q1;

	//Pass to quaternion form from eigen lib and return
	result->fromEigen();
	return result;
}

/**
* @brief Convertion function from  Quaternion to DCM matrix  object.
* @param rot Quaternion object.
* @result Returns the converted to DCM matrix object.
*/
SpinachGUI::Matrix3x3^  Orientation::ConvertToDCM(SpinachGUI::Quaternion^ rot)
{
	//Pass to Eigen lib object
	rot->toEigen();

	//Create a new empty DCM matrix object
	Matrix3x3^ result=gcnew Matrix3x3();

	//Set Eigen lib object using the Eigen quaternion form of rot
	*result->matrix3d=rot->quaterniond->toRotationMatrix();

	//Pass to DCM matrix form from eigen lib and return
	result->fromEigen(); 
	return result;
};

/**
* @brief Convertion function from  DCM matrix to Quaternion  object.
* @param rot DCM matrix object.
* @result Returns the converted to Quaternion object.
*/
SpinachGUI::Quaternion^  Orientation::ConvertToQuaternion(SpinachGUI::Matrix3x3^    rot)
{
	//Pass to Eigen lib object
	rot->toEigen();

	//Create a new empty quaternion object
	SpinachGUI::Quaternion^ result=gcnew SpinachGUI::Quaternion();

	//Set Eigen lib object using the Eigen matrix3d form of rot
	result->quaterniond=new Quaterniond(*rot->matrix3d);

	//Pass to quaternion form from eigen lib and return
	result->fromEigen(); 
	return result;
};

/**
* @brief Convertion function from  AngleAxis to Quaternion  object.
* @param rot AngleAxis object.
* @result Returns the converted to Quaternion object.
*/
SpinachGUI::Quaternion^  Orientation::ConvertToQuaternion(SpinachGUI::AngleAxis^   rot)
{
	//Pass to Eigen lib object
	rot->toEigen();

	//Create a new empty quaternion object
	SpinachGUI::Quaternion^ result=gcnew SpinachGUI::Quaternion();

	//Set Eigen lib object using the Eigen angleaxis form of rot
	result->quaterniond=new Quaterniond(*rot->angleaxisd);

	//Pass to quaternion form from eigen lib and return
	result->fromEigen(); 
	return result;
};

/**
* @brief Convertion function from Quaternion to AngleAxis object.
* @param rot Quaternion object.
* @result Returns the converted to AngleAxis object.
*/
SpinachGUI::AngleAxis^ Orientation::ConvertToAngleAxis(SpinachGUI::Quaternion^ rot)
{
	//Pass to Eigen lib object
	rot->toEigen();

	//Create a new empty Angle axis object
	SpinachGUI::AngleAxis^ result=gcnew SpinachGUI::AngleAxis();

	//Set Eigen lib object using the Eigen quaternion form of rot
	result->angleaxisd=new AngleAxisd(*rot->quaterniond);

	//Pass to angleaxis form from eigen lib and return
	result->fromEigen(); 
	return result;	
};

#pragma endregion


