//////////////////////////////////////////////////////////////////////////////////////////
//	VECTOR3D.h
//	Class declaration for a 3d vector
//	Downloaded from: www.paulsprojects.net
//	Created:	20th July 2002
//	Modified:	8th November 2002	-	Changed Constructor layout
//									-	Some speed Improvements
//									-	Corrected Lerp
//				7th January 2003	-	Added QuadraticInterpolate
//
//	Copyright (c) 2006, Paul Baker
//	Distributed under the New BSD Licence. (See accompanying file License.txt or copy at
//	http://www.paulsprojects.net/NewBSDLicense.txt)
//////////////////////////////////////////////////////////////////////////////////////////	

#pragma once

// The marching-cubes renderer stores positions and surface normals as XYZ triples.
struct VECTOR3D
{
    float x, y, z;

    VECTOR3D() : x(0), y(0), z(0) {}
    VECTOR3D(float x, float y, float z) : x(x), y(y), z(z) {}

    void Set(float a, float b, float c) { x=a; y=b; z=c; }
    void LoadZero() { Set(0,0,0); }
    operator const float*() const { return &x; }
};
