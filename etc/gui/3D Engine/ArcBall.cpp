/** KempoApi: The Turloc Toolkit *****************************/
/**  Arcball class for mouse manipulation.                  **/
/**                              (C) 1999-2003 Tatewake.com **/
/**   History:                                              **/
/**   08/17/2003 - (TJG) - Creation                         **/
/**   09/23/2003 - (TJG) - Bug fix and optimization         **/
/**   09/25/2003 - (TJG) - Version for NeHe Basecode users  **/
/**   15/11/2012 - (Andreas Biternas) - C++/CLI version     **/
/*************************************************************/
#include "../StdAfx.h"
#include <windows.h>											// Header File For Windows
#include <gl/gl.h>												// Header File For The OpenGL32 Library
#include <gl/glu.h>												// Header File For The GLu32 Library
#include "math.h"                                               // Needed for sqrtf
#include "ArcBall.h"                                            // ArcBall header

/**
*\file ArcBall.cpp
*/

///Arcball sphere constants:
///Diameter is       2.0f
///Radius is         1.0f
///Radius squared is 1.0f

/**
Map point ot sphere
*/
void ArcBall_t::_mapToSphere(Point2fT^ NewPt, Vector3fT^ NewVec)
{
    Point2fT^ TempPt=gcnew Point2fT() ;
    GLfloat length;

    //Copy paramter into temp point
	TempPt ->X = NewPt->X;
    TempPt ->Y = NewPt->Y;

    //Adjust point coords and scale down to range of [-1 ... 1]
    TempPt->X  =        (TempPt->X * this->AdjustWidth)  - 1.0f;
    TempPt->Y  = 1.0f - (TempPt->Y * this->AdjustHeight);

    //Compute the square of the length of the vector to the point from the center
    length      = (TempPt->X * TempPt->X) + (TempPt->Y * TempPt->Y);

    //If the point is mapped outside of the sphere... (length > radius squared)
    if (length > 1.0f)
    {
        GLfloat norm;

        //Compute a normalizing factor (radius / sqrt(length))
        norm    = 1.0f / Funcsqrt(length);

        //Return the "normalized" vector, a point on the sphere
        NewVec->X = TempPt->X * norm;
        NewVec->Y = TempPt->Y * norm;
        NewVec->Z = 0.0f;
    }
    else    //Else it's on the inside
    {
        //Return a vector to a point mapped inside the sphere sqrt(radius squared - length)
        NewVec->X = TempPt->X;
        NewVec->Y = TempPt->Y;
        NewVec->Z = Funcsqrt(1.0f - length);
    }
}

/**
Create/Destroy
*/
ArcBall_t::ArcBall_t(GLfloat NewWidth, GLfloat NewHeight)
{
     StVec=gcnew Vector3fT();          //Saved click vector
     EnVec=gcnew Vector3fT();  
	//Clear initial values
    this->StVec->X     =
    this->StVec->Y     = 
    this->StVec->Z     = 

    this->EnVec->X     =
    this->EnVec->Y     = 
    this->EnVec->Z     = 0.0f;

    //Set initial bounds
    this->setBounds(NewWidth, NewHeight);
}

/**
Mouse down
*/
void    ArcBall_t::click(Point2fT^ NewPt)
{
    //Map the point to the sphere
    this->_mapToSphere(NewPt, this->StVec);
}

/**
Mouse drag, calculate rotation
*/
void    ArcBall_t::drag(Point2fT^ NewPt, Quat4fT^ NewRot)
{
    //Map the point to the sphere
	this->_mapToSphere(NewPt, this->EnVec);
	
	

    //Return the quaternion equivalent to the rotation
    if (NewRot)
    {
        Vector3fT^  Perp=gcnew Vector3fT();

        //Compute the vector perpendicular to the begin and end vectors
        Vector3fCross(Perp, this->StVec, this->EnVec);

        //Compute the length of the perpendicular vector	    
        if (Vector3fLength(Perp) > Epsilon)    //if its non-zero
        {
            //We're ok, so return the perpendicular vector as the transform after all
            NewRot->X = Perp->X;
            NewRot->Y = Perp->Y;
            NewRot->Z = Perp->Z;
            //In the quaternion values, w is cosine (theta / 2), where theta is rotation angle
            NewRot->W= Vector3fDot(this->StVec, this->EnVec);
        }
        else                                    //if its zero
        {
            //The begin and end vectors coincide, so return an identity transform
            NewRot->X = 
            NewRot->Y = 
            NewRot->Z = 
            NewRot->W = 0.0f;
        }
    }
}



