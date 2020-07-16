// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Various functions for IO using the VTK library
 * \author Sassen
 */

#ifndef VTKIO_H
#define VTKIO_H

#ifdef GOAST_WITH_VTK

#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkAOSDataArrayTemplate.h>
#include <vtkVertex.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkLine.h>

#include <goast/Core/Topology.h>
#include <goast/Core/LocalMeshGeometry.h>

/**
 * \brief Save a triangle mesh as VTP file
 * \tparam VectorType the type of vectors used
 * \param Topology the connectivity of the mesh
 * \param Geometry the nodal positions
 * \param path the name of the VTP file to write
 */
template<typename VectorType>
void saveAsVTP( const MeshTopologySaver &Topology, const VectorType &Geometry, const std::string &path ) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  for ( unsigned int i = 0; i < Topology.getNumVertices(); ++i ) {
    VectorType pt( 3 );
    getXYZCoord( Geometry, pt, i );
    points->InsertNextPoint( pt[0], pt[1], pt[2] );
  }

  vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
  for ( int t = 0; t < Topology.getNumFaces(); ++t ) {
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
    triangle->GetPointIds()->SetId( 0, Topology.getNodeOfTriangle( t, 0 ));
    triangle->GetPointIds()->SetId( 1, Topology.getNodeOfTriangle( t, 1 ));
    triangle->GetPointIds()->SetId( 2, Topology.getNodeOfTriangle( t, 2 ));
    triangles->InsertNextCell( triangle );
  }


  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetPolys( triangles );

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName( path.c_str());
  writer->SetInputData( polydata );
  writer->Write();
}

/**
 * \brief Save a triangle mesh as VTP file with vertex-, face-, or mesh-wise scalar values
 * \tparam VectorType the type of vectors used
 * \param Topology the connectivity of the mesh
 * \param Geometry the nodal positions
 * \param path the name of the VTP file to write
 * \param Scalars the scalar values, the size of the vector has to be either the number of vertices or the number of faces or one, based on this it will be determined whether it are vertex-, face-, or mesh-wise scalars
 * \param scalarsName the display name of the scalar values
 */
template<typename VectorType>
void saveAsVTP( const MeshTopologySaver &Topology, const VectorType &Geometry, const std::string &path,
                const VectorType &Scalars, const std::string &scalarsName = "Scalars" ) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  for ( unsigned int i = 0; i < Topology.getNumVertices(); ++i ) {
    VectorType pt( 3 );
    getXYZCoord( Geometry, pt, i );
    points->InsertNextPoint( pt[0], pt[1], pt[2] );
  }

  vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
  for ( int t = 0; t < Topology.getNumFaces(); ++t ) {
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
    triangle->GetPointIds()->SetId( 0, Topology.getNodeOfTriangle( t, 0 ));
    triangle->GetPointIds()->SetId( 1, Topology.getNodeOfTriangle( t, 1 ));
    triangle->GetPointIds()->SetId( 2, Topology.getNodeOfTriangle( t, 2 ));
    triangles->InsertNextCell( triangle );
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetPolys( triangles );

  vtkAOSDataArrayTemplate<typename VectorType::Scalar> *data = vtkAOSDataArrayTemplate<typename VectorType::Scalar>::New();
  data->SetName( scalarsName.c_str());
  data->SetNumberOfValues( Scalars.size());
  for ( int i = 0; i < Scalars.size(); i++ ) {
    data->SetValue( i, Scalars[i] );
  }

  if ( Scalars.size() == Topology.getNumVertices()) {
    polydata->GetPointData()->AddArray( data );
    polydata->GetPointData()->SetActiveAttribute( 0, 0 );
  }
  else if ( Scalars.size() == Topology.getNumFaces()) {
    polydata->GetCellData()->AddArray( data );
    polydata->GetCellData()->SetActiveAttribute( 0, 0 );
  }
  else if ( Scalars.size() == 1 ) {
    polydata->GetFieldData()->AddArray( data );
  }
  else
    throw std::length_error( "saveAsVTP: Scalar values'" + scalarsName +
                             "' neither have the same number as vertices or faces." );

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName( path.c_str());
  writer->SetInputData( polydata );
  writer->Write();
}

/**
 * \brief Save a triangle mesh as VTP file with multiple vertex-, face-, or mesh-wise scalar values
 * \tparam VectorType the type of vectors used
 * \param Topology the connectivity of the mesh
 * \param Geometry the nodal positions
 * \param path the name of the VTP file to write
 * \param Scalars a map containing the scalar values, the keys will be used as display name, the sizes of the vectors have to be either the number of vertices or the number of faces or one, based on this it will be determined whether it are vertex-, face-, or mesh-wise scalars
 * \param activeScalar the name of the scalar value, which will be shown as default when opening the VTP file
 */
template<typename VectorType>
void saveAsVTP( const MeshTopologySaver &Topology, const VectorType &Geometry, const std::string &path,
                const std::map<std::string, VectorType> &Scalars, const std::string &activeScalar = "" ) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  for ( unsigned int i = 0; i < Topology.getNumVertices(); ++i ) {
    VectorType pt( 3 );
    getXYZCoord( Geometry, pt, i );
    points->InsertNextPoint( pt[0], pt[1], pt[2] );
  }

  vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
  for ( int t = 0; t < Topology.getNumFaces(); ++t ) {
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
    triangle->GetPointIds()->SetId( 0, Topology.getNodeOfTriangle( t, 0 ));
    triangle->GetPointIds()->SetId( 1, Topology.getNodeOfTriangle( t, 1 ));
    triangle->GetPointIds()->SetId( 2, Topology.getNodeOfTriangle( t, 2 ));
    triangles->InsertNextCell( triangle );
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetPolys( triangles );

  for ( auto const &d : Scalars ) {
    if ( d.second.size() == 3 * Topology.getNumVertices()) {
      vtkSmartPointer<vtkAOSDataArrayTemplate<typename VectorType::Scalar>> data
              = vtkSmartPointer<vtkAOSDataArrayTemplate<typename VectorType::Scalar>>::New();
      data->SetName( d.first.c_str());
      data->SetNumberOfComponents( 3 );
      data->SetNumberOfTuples( d.second.size() / 3 );
      VectorType localVec( 3 );
      for ( int i = 0; i < Topology.getNumVertices(); i++ ) {
        getXYZCoord( d.second, localVec, i );
        data->SetTuple3( i, localVec[0], localVec[1], localVec[2] );
      }

      polydata->GetPointData()->AddArray( data );
      if ( d.first == activeScalar )
        polydata->GetPointData()->SetActiveAttribute( activeScalar.c_str(), 0 );
    }
    else {
      vtkSmartPointer<vtkAOSDataArrayTemplate<typename VectorType::Scalar>> data
              = vtkSmartPointer<vtkAOSDataArrayTemplate<typename VectorType::Scalar>>::New();
      data->SetName( d.first.c_str());
      data->SetNumberOfValues( d.second.size());
      for ( int i = 0; i < d.second.size(); i++ ) {
        data->SetValue( i, d.second[i] );
      }

      if ( d.second.size() == Topology.getNumVertices()) {
        polydata->GetPointData()->AddArray( data );
        if ( d.first == activeScalar )
          polydata->GetPointData()->SetActiveAttribute( activeScalar.c_str(), 0 );
      }
      else if ( d.second.size() == Topology.getNumFaces()) {
        polydata->GetCellData()->AddArray( data );
        if ( d.first == activeScalar )
          polydata->GetCellData()->SetActiveAttribute( activeScalar.c_str(), 0 );
      }
      else if ( d.second.size() == 1 ) {
        polydata->GetFieldData()->AddArray( data );
      }
      else
        throw std::length_error( "saveAsVTP: Scalar values'" + d.first +
                                 "' neither have the same number as vertices or faces." );
    }


  }

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName( path.c_str());
  writer->SetInputData( polydata );
  writer->Write();
}

/**
 * \brief Save the edges of a triangle mesh as VTP with edge-wise scalar values
 * \tparam VectorType the type of vectors used
 * \param Topology the connectivity of the mesh
 * \param Geometry the nodal positions
 * \param path the name of the VTP file to write
 * \param Scalars the scalar values, the size of the vector has to be the same as the number of edges
 * \param scalarsName the display name of the scalar values
 */
template<typename VectorType>
void saveWireframeAsVTP( const MeshTopologySaver &Topology, const VectorType &Geometry, const std::string &path,
                         const VectorType &Scalars, const std::string &scalarsName = "Scalars" ) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  if ( Scalars.size() != Topology.getNumEdges())
    throw std::length_error( "saveWireframeAsVTP: Scalar values '" + scalarsName +
                             "' must have same number as edges." );

  for ( unsigned int i = 0; i < Topology.getNumVertices(); ++i ) {
    VectorType pt( 3 );
    getXYZCoord( Geometry, pt, i );
    points->InsertNextPoint( pt[0], pt[1], pt[2] );
  }

  vtkSmartPointer<vtkCellArray> edges = vtkSmartPointer<vtkCellArray>::New();
  for ( int e = 0; e < Topology.getNumEdges(); ++e ) {
    vtkSmartPointer<vtkLine> edge = vtkSmartPointer<vtkLine>::New();
    edge->GetPointIds()->SetId( 0, Topology.getAdjacentNodeOfEdge( e, 0 ));
    edge->GetPointIds()->SetId( 1, Topology.getAdjacentNodeOfEdge( e, 1 ));
    edges->InsertNextCell( edge );
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetPolys( edges );

  vtkAOSDataArrayTemplate<typename VectorType::Scalar> *data = vtkAOSDataArrayTemplate<typename VectorType::Scalar>::New();
  data->SetName( scalarsName.c_str());
  data->SetNumberOfValues( Scalars.size());
  for ( int i = 0; i < Scalars.size(); i++ ) {
    data->SetValue( i, Scalars[i] );
  }

  polydata->GetCellData()->AddArray( data );
  polydata->GetCellData()->SetActiveAttribute( 0, 0 );

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName( path.c_str());
  writer->SetInputData( polydata );
  writer->Write();
}

/**
 * \brief Save the edges of a triangle mesh as VTP with edge-wise scalar values
 * \tparam VectorType the type of vectors used
 * \param Topology the connectivity of the mesh
 * \param Geometry the nodal positions
 * \param path the name of the VTP file to write
 * \param Scalars a map containing the scalar values, the keys will be used as display name, the size of the vector has to be the same as the number of edges
 */
template<typename VectorType>
void saveWireframeAsVTP( const MeshTopologySaver &Topology, const VectorType &Geometry, const std::string &path,
                         const std::map<std::string, VectorType> &Scalars ) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();


  for ( unsigned int i = 0; i < Topology.getNumVertices(); ++i ) {
    VectorType pt( 3 );
    getXYZCoord( Geometry, pt, i );
    points->InsertNextPoint( pt[0], pt[1], pt[2] );
  }

  vtkSmartPointer<vtkCellArray> edges = vtkSmartPointer<vtkCellArray>::New();
  for ( int e = 0; e < Topology.getNumEdges(); ++e ) {
    vtkSmartPointer<vtkLine> edge = vtkSmartPointer<vtkLine>::New();
    edge->GetPointIds()->SetId( 0, Topology.getAdjacentNodeOfEdge( e, 0 ));
    edge->GetPointIds()->SetId( 1, Topology.getAdjacentNodeOfEdge( e, 1 ));
    edges->InsertNextCell( edge );
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetPolys( edges );

  for ( auto const &d : Scalars ) {
    if ( d.second.size() != Topology.getNumEdges())
      throw std::length_error( "saveAsVTP: Scalar values'" + d.first +
                               "' neither have the same number as vertices or faces." );

    vtkSmartPointer<vtkAOSDataArrayTemplate<typename VectorType::Scalar>> data
            = vtkSmartPointer<vtkAOSDataArrayTemplate<typename VectorType::Scalar>>::New();
    data->SetName( d.first.c_str());
    data->SetNumberOfValues( d.second.size());
    for ( int i = 0; i < d.second.size(); i++ ) {
      data->SetValue( i, d.second[i] );
    }

    polydata->GetCellData()->AddArray( data );
  }


  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName( path.c_str());
  writer->SetInputData( polydata );
  writer->Write();
}

/**
 * \brief Save a selection of the edges of a triangle mesh as VTP
 * \tparam VectorType the type of vectors used
 * \param Topology the connectivity of the mesh
 * \param Geometry the nodal positions
 * \param path the name of the VTP file to write
 * \param Edges a vector containing the indices of the edges to be written
 */
template<typename VectorType>
void savePartialWireframeAsVTP( const MeshTopologySaver &Topology, const VectorType &Geometry, const std::string &path,
                                const std::vector<int> &Edges ) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();


  for ( unsigned int i = 0; i < Topology.getNumVertices(); ++i ) {
    VectorType pt( 3 );
    getXYZCoord( Geometry, pt, i );
    points->InsertNextPoint( pt[0], pt[1], pt[2] );
  }

  vtkSmartPointer<vtkCellArray> edges = vtkSmartPointer<vtkCellArray>::New();
  for ( int e = 0; e < Topology.getNumEdges(); ++e ) {
    if ( std::find( Edges.begin(), Edges.end(), e ) == Edges.end())
      continue;

    vtkSmartPointer<vtkLine> edge = vtkSmartPointer<vtkLine>::New();
    edge->GetPointIds()->SetId( 0, Topology.getAdjacentNodeOfEdge( e, 0 ));
    edge->GetPointIds()->SetId( 1, Topology.getAdjacentNodeOfEdge( e, 1 ));
    edges->InsertNextCell( edge );
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetPolys( edges );



  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName( path.c_str());
  writer->SetInputData( polydata );
  writer->Write();
}

#endif

#endif //VTKIO_H
