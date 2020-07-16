// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __OPENMESHINCLUDES_H
#define __OPENMESHINCLUDES_H

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC system_header
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/McDecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>

#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Geometry/QuadricT.hh>

//#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>

#ifdef PI
#undef PI
#endif
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

//===========================================================
//===========================================================
//===========================================================
struct TriTraits : public OpenMesh::DefaultTraits
{
  /// Use double precision points
  typedef OpenMesh::Vec3d Point;
  /// Use double precision Normals
  typedef OpenMesh::Vec3d Normal;
  /// Use RGBA Color
  typedef OpenMesh::Vec4f Color;  
 
  // need status attributes for decimation
  VertexAttributes( OpenMesh::Attributes::Status );
  EdgeAttributes( OpenMesh::Attributes::Status );
  FaceAttributes( OpenMesh::Attributes::Status );  
};

typedef OpenMesh::TriMesh_ArrayKernelT<TriTraits>  TriMesh;


//Shortcuts for different Handles
typedef TriMesh::VertexHandle          VH;
typedef TriMesh::FaceHandle            FH;
typedef TriMesh::EdgeHandle            EH;
typedef TriMesh::HalfedgeHandle        HEH;
//Iterators
typedef TriMesh::VertexIter            VIter;
typedef TriMesh::FaceIter              FIter;
typedef TriMesh::EdgeIter              EIter;
typedef TriMesh::HalfedgeIter          HIter;
//Circulators
typedef TriMesh::VertexVertexIter      VVIter;
typedef TriMesh::VertexFaceIter        VFIter;
typedef TriMesh::VertexEdgeIter        VEIter;
typedef TriMesh::VertexOHalfedgeIter   VOHIter;
typedef TriMesh::VertexIHalfedgeIter   VIHIter;
typedef TriMesh::FaceVertexIter        FVIter;
typedef TriMesh::FaceFaceIter          FFIter;
typedef TriMesh::FaceEdgeIter          FEIter;
//3D-Point
typedef OpenMesh::Vec3d      Vec3d;
typedef OpenMesh::Vec2d      Vec2d;
typedef OpenMesh::Vec3i      Vec3i;
typedef OpenMesh::Vec2i      Vec2i;

#endif
