// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include <goast/Core.h>
#include <goast/DiscreteShells/DeformationTransfer.h>


/** 
 * \brief Deformation transfer based on deformation gradients
 * \author Heeren
 * 
 * Implementation of the method proposed in \cite BoSuPa06 (see all \cite SuPo04 )
 * 
 * Consider a reference shape \f$ r \f$ and a corresponding deformed shape \f$ \tilde r = \phi(r) \f$,
 * where \f$ \phi \f$ is some unknown deformation, and a target shape \f$ t \f$. 
 * The algorithm computes the corresponding deformed shape \f$ \tilde t = \phi(t) \f$.
 * 
 * All shapes are supposed to be given as nodal positions of meshes that are in dense correspondence. 
 * 
 */
int main(int argc, char *argv[]){

  try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE DEFORMATION TRANSFER" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
        
    // load two end meshes and variational mesh
    TriMesh refMesh, targetMesh, deformedRefMesh;
    OpenMesh::IO::read_mesh(refMesh, "../../data/faces/neutral.ply"); 
    OpenMesh::IO::read_mesh(targetMesh, "../../data/faces/smile.ply");
    OpenMesh::IO::read_mesh(deformedRefMesh, "../../data/faces/disgust.ply"); 
    
    //! get geometries
    typename DefaultConfigurator::VectorType refGeom, targetGeom, deformedRefGeom, deformedTargetGeom;
    getGeometry( refMesh, refGeom );
    getGeometry( targetMesh, targetGeom );
    getGeometry( deformedRefMesh, deformedRefGeom );
    
    // set up topology and define boundary mask (of fixed vertices)
    MeshTopologySaver Topol( refMesh );
    std::vector<int> mask;  
    Topol.fillPartialBoundaryMask( mask, 6137 );
    
    // perform deformation transfer
    DeformationTransferOperator<DefaultConfigurator>( Topol, mask ).execute( refGeom, deformedRefGeom, targetGeom, deformedTargetGeom );
    
    // save transfered variation
    setGeometry( deformedRefMesh, deformedTargetGeom );
    OpenMesh::IO::write_mesh(deformedRefMesh, "deformationTransfer_faces.ply");
    
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}