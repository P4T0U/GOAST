// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __MULTIRESOPTIMIZATION_H
#define __MULTIRESOPTIMIZATION_H

#include "MultiResOperator.h"
#include "MeshDecimation.h"

//==========================================================================================================
//  RECURSIVE MULTILEVEL OPERATORS
//==========================================================================================================

/**
 * \brief Abstract base class to prescribe structure of optimization operators that are supposed to be used with the recursive multilevel optimization.
 * \author Heeren
 *
 * Recursive multilevel operators are designed to have more than two hierarchical levels.
 *
 * The hierarchical scheme is controlled by the variables "numOfRecursiveLevels" and "theta" in the constructor,
 * where #(nodes on level k-1) = theta * #(nodes on level k),  0 < theta < 1.
 *
 * To make use of this scheme proceed as follows:
 *
 * (1) Write a class that inherits from MultiResOptimInterface<>.
 *
 * (2) Overwrite the pure virtual functions fillInputMeshes() and optimize(). In particular, have a look at documentation of MultiResOptimInterface<>::optimize()!
 *
 * (3) Construct an object of your class as well as one of RecursiveMultilevelOptimization<>, by parsing your class to the latter one in the constructor.
 *
 * (4) call apply( ) on the RecursiveMultilevelOptimization<> object.
 *
 * NOTE In particular, RecursiveMultilevelOptimization<> can also handle the case without any hierarchical scheme!
 */
template <typename ConfiguratorType>
class MultiResOptimInterface{

  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  int _numOfRecursiveLevels;
  RealType _theta;
  
public:
  MultiResOptimInterface( int numOfRecursiveLevels, RealType Theta ) : _numOfRecursiveLevels( numOfRecursiveLevels ), _theta(Theta) {}
  
  virtual ~MultiResOptimInterface(){}
  
  // All important optimization function that is called on every level.
  // The independent variables are given as meshes, the dependent variables as geometries.
  // Except for the lowest level (where the dep. variables have to be initialized) the dependent variables
  // are initializations for the optimization on this level and are supposed to be overwritten by the solution on this level.
  virtual void optimize( const MeshTopologySaver& /*Topology*/, 
			 const std::vector<TriMesh>& /*IndependentVariablesOnThisLevel*/,  //TODO change to passing inputGeometries, not meshes!?
			 std::vector<VectorType>& /*DependentVariablesOnThisLevel*/, 
			 int /*level*/ ) const = 0;
  
  // load independent variables, i.e. initial meshes, from file (or do nothing)
  virtual void fillInputMeshes( std::vector<TriMesh>& ) const = 0;
  
  int getNumOfRecursiveLevels() const {
    return _numOfRecursiveLevels;
  }
  
  RealType getTheta() const {
    return _theta;
  }
  
};

/**
 * \brief Recursive multilevel operator
 * \author Heeren
 *
 * Needs an operator derived from abstract class MultiResOptimInterface (see above) as argument in the constructor.
 *
 * Usage: Provide this operator, call apply().
 */
template <typename ConfiguratorType>
class RecursiveMultilevelOptimization {

protected:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  const MultiResOptimInterface<ConfiguratorType>& _optInterface;  
  bool _quiet;

public:
  RecursiveMultilevelOptimization( const MultiResOptimInterface<ConfiguratorType>& optInterface, bool Quiet = false ) : _optInterface(optInterface) , _quiet( Quiet ){}

  // recursive optimization
  void apply( std::vector<VectorType>& Solutions, const BitVector* bdryMask = NULL ) const {
        // write input meshes to inputMeshes
        std::vector<TriMesh> inputMeshes;
        _optInterface.fillInputMeshes( inputMeshes );

        // reference mesh
        TriMesh refMesh( inputMeshes[0] );
        MeshTopologySaver Topology( refMesh );

        //take care of boundary conditions
        if( bdryMask )
            Topology.setFixedNodes( *bdryMask );

        apply( Topology, inputMeshes, Solutions );
  }

  // recursive optimization
  void apply( const BitVector* bdryMask = NULL ) const {
        std::vector<VectorType> Solutions;
        apply( Solutions, bdryMask );
  }

  // recursive optimization
  //TODO change to passing inputGeometries, not meshes!?
  void apply( const MeshTopologySaver& Topology, const std::vector<TriMesh>& inputMeshes, std::vector<VectorType>& Solutions ) const {    
    optmizeOrDecimate( Topology, inputMeshes, Solutions, _optInterface.getNumOfRecursiveLevels() );
  }
  
protected:
    //TODO change to passing inputGeometries, not meshes!?
  void optmizeOrDecimate( const MeshTopologySaver& Topology, const std::vector<TriMesh>& inputMeshes, std::vector<VectorType>& Solutions, int k ) const {

    
    // if k < K compute initialization on coarser level
    if( k > 0 ){          
      
      // decimation
      std::vector<TriMesh> coarseMeshes;
      MultiResolutionOperator<ConfiguratorType> multiresOp( Topology, _quiet );
      multiresOp.computeDecimation( inputMeshes, _optInterface.getTheta(), coarseMeshes );
      
      std::vector<VectorType> coarseSolutions;  
      TriMesh coarseRefMesh( coarseMeshes[0] );
      setGeometry( coarseRefMesh,  multiresOp.getCoarseRefGeom() );
      MeshTopologySaver coarseTopology( coarseRefMesh );    

      // adapt Dirichlet boundary mask
      if( Topology.getNumFixedNodes() > 0 ){
        if( !_quiet ) std::cerr << "Adapt Dirichlet boundary conditions..." << std::endl;    
        BitVector coarseMask( coarseTopology.getNumVertices() );
        multiresOp.computeDecimatedBoundaryMask( Topology.getFixedNodes(), coarseMask );
        coarseTopology.setFixedNodes( coarseMask );
      }
    
      // recursive optimization
      optmizeOrDecimate( coarseTopology, coarseMeshes, coarseSolutions, k-1 );    
      
      // detail transfer
      Solutions.resize( coarseSolutions.size() );
      multiresOp.prolongate( coarseSolutions, Solutions );
    
      //TODO add saving methods for base mesh or anchors?
    }
        
    if( !_quiet ) std::cerr << std::endl << "-------------------------------------------------" << std::endl;
    if( !_quiet ) std::cerr << "------- START RECURSIVE LEVEL " << k << " -----------------" << std::endl;       
    if( !_quiet ) std::cerr << "-------------------------------------------------" << std::endl << std::endl;
    
    // optimization on this level
    _optInterface.optimize( Topology, inputMeshes, Solutions, k );             
        
    if( !_quiet ) std::cerr << std::endl <<  "-------------------------------------------------" << std::endl;
    if( !_quiet ) std::cerr << "Time for optimization on level  " << k << ": [UNKNOWN] seconds." << std::endl; 
    if( !_quiet ) std::cerr << "------ END RECURSIVE LEVEL " << k << " --------------------" << std::endl;
    if( !_quiet ) std::cerr <<  "-------------------------------------------------" << std::endl<< std::endl;
    
  }
  

};

#endif