/** KempoApi: The Turloc Toolkit *****************************/
/**  Arcball class for mouse manipulation.                  **/
/**                              (C) 1999-2003 Tatewake.com **/
/**   History:                                              **/
/**   08/17/2003 - (TJG) - Creation                         **/
/**   09/23/2003 - (TJG) - Bug fix and optimization         **/
/**   09/25/2003 - (TJG) - Version for NeHe Basecode users  **/
/**   15/11/2012 - (Andreas Biternas) - C++/CLI version     **/
/*************************************************************/
#pragma once
# include "math.h"
using namespace System;
using namespace System::Runtime::InteropServices;
# include "assert.h"


///Custom math, or speed overrides
#define Funcsqrt    sqrtf

///utility macros
///assuming IEEE-754(GLfloat), which i believe has max precision of 7 bits
# define Epsilon 1.0e-5

///Math types derived from the KempoApi tMath library

///A generic 2-element tuple that is represented by Single-precision floating point x,y coordinates. 
typedef public ref struct Tuple2f_t
{
	GLfloat X, Y;

	Tuple2f_t^ operator+(Tuple2f_t^ a) 
	{ //outside the class
		Tuple2f_t^ result=gcnew Tuple2f_t();
		result->X=a->X+X;
		result->Y=a->Y+Y;
		return result;
	};


} Tuple2fT;     

///A generic 3-element tuple that is represented by Single precision-floating point x,y,z coordinates. 
typedef public ref struct Tuple3f_t
{
	GLfloat X, Y, Z;
} Tuple3fT; 

///A 4-element tuple represented by Single-precision floating point x,y,z,w coordinates. 
typedef public ref struct Tuple4f_t
{
	GLfloat X, Y, Z, W;
} Tuple4fT;      

///A Single precision floating point 3 by 3 matrix. 
typedef public ref struct Matrix3f_t
{
	//column major
	GLfloat XX; GLfloat XY;  GLfloat XZ;
	GLfloat YX; GLfloat YY;  GLfloat YZ;
	GLfloat ZX; GLfloat ZY;  GLfloat ZZ;
} Matrix3fT;     

///A Single precision floating point 4 by 4 matrix.
typedef public ref struct Matrix4f_t
{
	//column major
	GLfloat XX;  GLfloat XY;  GLfloat XZ;   GLfloat XW; 
	GLfloat YX;  GLfloat YY;  GLfloat YZ;   GLfloat YW; 
	GLfloat ZX;  GLfloat ZY;  GLfloat ZZ;   GLfloat ZW; 
	GLfloat TX;  GLfloat TY;  GLfloat TZ;   GLfloat TW; 
} Matrix4fT;      

///"Inherited" types

#define Point2fT    Tuple2fT   ///A 2 element point that is represented by Single precision floating point x,y coordinates. 
#define Quat4fT     Tuple4fT   ///A 4 element unit quaternion represented by Single precision floating point x,y,z,w coordinates. 
#define Vector2fT   Tuple2fT   ///A 2-element vector that is represented by Single-precision floating point x,y coordinates. 
#define Vector3fT   Tuple3fT   ///A 3-element vector that is represented by Single-precision floating point x,y,z coordinates. 




//Math functions

/**
* Sets the value of this tuple to the vector sum of itself and tuple t1.
* @param NewObj This tuple point.
* @param t1  the other tuple
*/
inline static void Point2fAdd(Point2fT^ NewObj,  Tuple2fT^ t1)
{
	assert(NewObj && t1);

	NewObj->X += t1->X;
	NewObj->Y += t1->Y;
}

/**
* Sets the value of this tuple to the vector difference of itself and tuple t1 (this = this - t1).
* @param NewObj This tuple point.
* @param t1 The other tuple.
*/
inline static void Point2fSub(Point2fT^ NewObj,  Tuple2fT^ t1)
{
	assert(NewObj && t1);

	NewObj->X -= t1->X;
	NewObj->Y -= t1->Y;
}

/**
* Sets this vector to be the vector cross product of vectors v1 and v2.
* @param Result The Resulted vector.
* @param v1 the first vector.
* @param v2 the second vector.
*/
inline static void Vector3fCross(Vector3fT^ Result, Vector3fT^ v1, Vector3fT^ v2)
{
	//Vector3fT Result; //safe not to initialize

	assert(Result && v1 && v2);

	/// store on stack once for aliasing-safty
	/// i.e. safe when a.cross(a, b)

	Result->X = (v1->Y * v2->Z) - (v1->Z * v2->Y);
	Result->Y = (v1->Z * v2->X) - (v1->X * v2->Z);
	Result->Z = (v1->X * v2->Y) - (v1->Y * v2->X);

	//copy result back
	//*NewObj = Result;
}

/**
* Computes the dot product of the this vector and vector v1.
* @param NewObj This vector.
* @param  v1 the other vector.
*/
inline 	static GLfloat Vector3fDot(const Vector3fT^ NewObj, const Vector3fT^ v1)
{
	assert(NewObj && v1);

	return  (NewObj->X * v1->X) +
		(NewObj->Y * v1->Y) +
		(NewObj->Z * v1->Z);
}

/**
* Returns the squared length of this vector.
* @param NewObj This vector.
* @return the squared length of this vector.
*/
inline static GLfloat Vector3fLengthSquared( Vector3fT^ NewObj)
{
	assert(NewObj);

	return  (NewObj->X * NewObj->X) +
		(NewObj->Y * NewObj->Y) +
		(NewObj->Z * NewObj->Z);
}

/**
* Returns the length of this vector.
* @param NewObj This vector.
* @return the length of this vector
*/
inline static GLfloat Vector3fLength(Vector3fT^ NewObj)
{
	assert(NewObj);

	return Funcsqrt(Vector3fLengthSquared(NewObj));
}

/**
* Sets this Matrix3 to zero.
* @param NewObj This matrix.
*/
inline static void Matrix3fSetZero(Matrix3fT^ NewObj)
{
	NewObj->XX = NewObj->XY = NewObj->XZ = 0.0f;
	NewObj->YX = NewObj->YY = NewObj->YZ = 0.0f;
	NewObj->ZX = NewObj->ZY = NewObj->ZZ = 0.0f;
}

/**
* Sets this Matrix3 to identity.
* @param NewObj This matrix.
*/
inline static void Matrix3fSetIdentity(Matrix3fT^ NewObj)
{
	Matrix3fSetZero(NewObj);

	///then set diagonal as 1
	NewObj->XX = NewObj->YY = NewObj->ZZ = 1.0f;
}

/**
* Sets the value of this matrix to the matrix conversion of the
* quaternion argument. 
* @param NewObj This matrix.
* @param q1 the quaternion to be converted 
*/
//$hack this can be optimized some(if s == 0)
inline static void Matrix3fSetRotationFromQuat4f(Matrix3fT^% NewObj, Quat4fT^ q1)
{
	GLfloat n, s;
	GLfloat xs, ys, zs;
	GLfloat wx, wy, wz;
	GLfloat xx, xy, xz;
	GLfloat yy, yz, zz;

	assert(NewObj && q1);

	Matrix3fT^ tNewObj=gcnew Matrix3fT();


	n = (q1->X * q1->X) + (q1->Y * q1->Y) + (q1->Z * q1->Z) + (q1->W * q1->W);
	s = (n > 0.0f) ? (2.0f / n) : 0.0f;

	xs = q1->X * s;  ys = q1->Y * s;  zs = q1->Z * s;
	wx = q1->W * xs; wy = q1->W * ys; wz = q1->W * zs;
	xx = q1->X * xs; xy = q1->X * ys; xz = q1->X * zs;
	yy = q1->Y * ys; yz = q1->Y * zs; zz = q1->Z * zs;

	tNewObj->XX = 1.0f - (yy + zz); tNewObj->YX =         xy - wz;  tNewObj->ZX =         xz + wy;
	tNewObj->XY =         xy + wz;  tNewObj->YY = 1.0f - (xx + zz); tNewObj->ZY =         yz - wx;
	tNewObj->XZ =         xz - wy;  tNewObj->YZ =         yz + wx;  tNewObj->ZZ = 1.0f - (xx + yy);

	NewObj=tNewObj;
}

/**
* Sets the value of this matrix to the result of multiplying itself
* with matrix m1. 
* @param NewObj This matrix.
* @param m1 the other matrix 
*/
inline static void Matrix3fMulMatrix3f(Matrix3fT^% NewObj, Matrix3fT^ m1)
{
	Matrix3fT^ Result=gcnew Matrix3fT(); //safe not to initialize

	assert(NewObj && m1);

	// alias-safe way.
	Result->XX = (NewObj->XX * m1->XX) + (NewObj->YX * m1->XY) + (NewObj->ZX * m1->XZ);
	Result->YX = (NewObj->XX * m1->YX) + (NewObj->YX * m1->YY) + (NewObj->ZX * m1->YZ);
	Result->ZX = (NewObj->XX * m1->ZX) + (NewObj->YX * m1->ZY) + (NewObj->ZX * m1->ZZ);

	Result->XY = (NewObj->XY * m1->XX) + (NewObj->YY * m1->XY) + (NewObj->ZY * m1->XZ);
	Result->YY = (NewObj->XY * m1->YX) + (NewObj->YY * m1->YY) + (NewObj->ZY * m1->YZ);
	Result->ZY = (NewObj->XY * m1->ZX) + (NewObj->YY * m1->ZY) + (NewObj->ZY * m1->ZZ);

	Result->XZ = (NewObj->XZ * m1->XX) + (NewObj->YZ * m1->XY) + (NewObj->ZZ * m1->XZ);
	Result->YZ = (NewObj->XZ * m1->YX) + (NewObj->YZ * m1->YY) + (NewObj->ZZ * m1->YZ);
	Result->ZZ = (NewObj->XZ * m1->ZX) + (NewObj->YZ * m1->ZY) + (NewObj->ZZ * m1->ZZ);

	//copy result back to this
	NewObj = Result;
}

/**
* Sets the value of this matrix to the matrix m1.
* @param NewObj This matrix.
* @param m1 the other matrix 
*/
inline static void Matrix4fSetRotationScaleFromMatrix4f(Matrix4fT^ NewObj, Matrix4fT^ m1)
{
	assert(NewObj && m1);

	NewObj->XX = m1->XX; NewObj->YX = m1->YX; NewObj->ZX = m1->ZX;
	NewObj->XY = m1->XY; NewObj->YY = m1->YY; NewObj->ZY = m1->ZY;
	NewObj->XZ = m1->XZ; NewObj->YZ = m1->YZ; NewObj->ZZ = m1->ZZ;
}

/**
* Performs SVD on this matrix and gets scale and rotation.
* Rotation is placed into rot3, and rot4.
* @param NewObj This matrix.
* @param rot3 the rotation factor(Matrix3d). if null, ignored
* @param rot4 the rotation factor(Matrix4) only upper 3x3 elements are changed. if null, ignored
* @return scale factor
*/
inline static GLfloat Matrix4fSVD(Matrix4fT^ NewObj, Matrix3fT^ rot3, Matrix4fT^ rot4)
{
	GLfloat s, n;

	assert(NewObj);

	/// this is a simple svd.
	/// Not complete but fast and reasonable.
	/// See comment in Matrix3d.

	s = Funcsqrt(
		( (NewObj->XX * NewObj->XX) + (NewObj->XY * NewObj->XY) + (NewObj->XZ * NewObj->XZ) + 
		(NewObj->YX * NewObj->YX) + (NewObj->YY * NewObj->YY) + (NewObj->YZ * NewObj->YZ) +
		(NewObj->ZX * NewObj->ZX) + (NewObj->ZY * NewObj->ZY) + (NewObj->ZZ * NewObj->ZZ) ) / 3.0f );

	if (rot3)   //if pointer not null
	{
		//this->getRotationScale(rot3);
		rot3->XX = NewObj->XX; rot3->XY = NewObj->XY; rot3->XZ = NewObj->XZ;
		rot3->YX = NewObj->YX; rot3->YY = NewObj->YY; rot3->YZ = NewObj->YZ;
		rot3->ZX = NewObj->ZX; rot3->ZY = NewObj->ZY; rot3->ZZ = NewObj->ZZ;

		// zero-div may occur.

		n = 1.0f / Funcsqrt( (NewObj->XX * NewObj->XX) +
			(NewObj->XY * NewObj->XY) +
			(NewObj->XZ * NewObj->XZ) );
		rot3->XX *= n;
		rot3->XY *= n;
		rot3->XZ *= n;

		n = 1.0f / Funcsqrt( (NewObj->YX * NewObj->YX) +
			(NewObj->YY * NewObj->YY) +
			(NewObj->YZ * NewObj->YZ) );
		rot3->YX *= n;
		rot3->YY *= n;
		rot3->YZ *= n;

		n = 1.0f / Funcsqrt( (NewObj->ZX * NewObj->ZX) +
			(NewObj->ZY * NewObj->ZY) +
			(NewObj->ZZ * NewObj->ZZ) );
		rot3->ZX *= n;
		rot3->ZY *= n;
		rot3->ZZ *= n;
	}

	if (rot4)   //if pointer not null
	{
		if (rot4 != NewObj)
		{
			Matrix4fSetRotationScaleFromMatrix4f(rot4, NewObj);  // private method
		}

		// zero-div may occur.

		n = 1.0f / Funcsqrt( (NewObj->XX * NewObj->XX) +
			(NewObj->XY * NewObj->XY) +
			(NewObj->XZ * NewObj->XZ) );
		rot4->XX *= n;
		rot4->XY *= n;
		rot4->XZ *= n;

		n = 1.0f / Funcsqrt( (NewObj->YX * NewObj->YX) +
			(NewObj->YY * NewObj->YY) +
			(NewObj->YZ * NewObj->YZ) );
		rot4->YX *= n;
		rot4->YY *= n;
		rot4->YZ *= n;

		n = 1.0f / Funcsqrt( (NewObj->ZX * NewObj->ZX) +
			(NewObj->ZY * NewObj->ZY) +
			(NewObj->ZZ * NewObj->ZZ) );
		rot4->ZX *= n;
		rot4->ZY *= n;
		rot4->ZZ *= n;
	}

	return s;
}

/**
* Sets the value of this matrix to the matrix m1. 
* @param NewObj This matrix.
* @param m1 the other matrix 
*/
inline static void Matrix4fSetRotationScaleFromMatrix3f(Matrix4fT^ NewObj, Matrix3fT^ m1)
{
	assert(NewObj && m1);

	NewObj->XX = m1->XX; NewObj->YX = m1->YX; NewObj->ZX = m1->ZX;
	NewObj->XY = m1->XY; NewObj->YY = m1->YY; NewObj->ZY = m1->ZY;
	NewObj->XZ = m1->XZ; NewObj->YZ = m1->YZ; NewObj->ZZ = m1->ZZ;
}

/**
* Sets the value of this matrix to the result of multiplying itself
* with value scale. 
* @param NewObj This matrix.
* @param scale the value to scale matrix to 
*/
inline static void Matrix4fMulRotationScale(Matrix4fT^ NewObj, GLfloat scale)
{
	assert(NewObj);

	NewObj->XX *= scale; NewObj->YX *= scale; NewObj->ZX *= scale;
	NewObj->XY *= scale; NewObj->YY *= scale; NewObj->ZY *= scale;
	NewObj->XZ *= scale; NewObj->YZ *= scale; NewObj->ZZ *= scale;
}

/**
* Sets the rotational component (upper 3x3) of this matrix to the matrix
* values in the T precision Matrix3d argument; the other elements of
* this matrix are unchanged; a singular value decomposition is performed
* on this object's upper 3x3 matrix to factor out the scale, then this
* object's upper 3x3 matrix components are replaced by the passed rotation
* components, and then the scale is reapplied to the rotational
* component
* @param NewObj This matrix.
* @param m1 T precision 3x3 matrix
*/
inline static void Matrix4fSetRotationFromMatrix3f(Matrix4fT^ NewObj,  Matrix3fT^ m1)
{
	GLfloat scale;

	assert(NewObj && m1);

	scale = Matrix4fSVD(NewObj, nullptr, nullptr);

	Matrix4fSetRotationScaleFromMatrix3f(NewObj, m1);
	Matrix4fMulRotationScale(NewObj, scale); 
}

/** \class ArcBall
* Main Arcball class for the Shoemake algorithm 
*/
typedef public ref class ArcBall_t
{
protected:
	inline
		void _mapToSphere(Point2fT^ NewPt, Vector3fT^ NewVec) ;

public:
	///Create/Destroy
	ArcBall_t(GLfloat NewWidth, GLfloat NewHeight);
	~ArcBall_t() { /* nothing to do */ };

	///Set new bounds
	inline	void    setBounds(GLfloat NewWidth, GLfloat NewHeight)
	{
		//assert((NewWidth > 1.0f) && (NewHeight > 1.0f));

		//Set adjustment factor for width/height
		this->AdjustWidth  = NewWidth > 1.0f ? 2.0f / (NewWidth - 1.0f) : 1.0f;
		this->AdjustHeight = NewHeight > 1.0f ? 2.0f / (NewHeight - 1.0f) : 1.0f;
	}

	///Mouse down
	void    click(Point2fT^ NewPt);

	///Mouse drag, calculate rotation
	void    drag(Point2fT^ NewPt, Quat4fT^ NewRot);

protected:
	Vector3fT^   StVec;          ///Saved click vector
	Vector3fT^   EnVec;          ///Saved drag vector
	GLfloat     AdjustWidth;    ///Mouse bounds width
	GLfloat     AdjustHeight;   ///Mouse bounds height

} ArcBallT;