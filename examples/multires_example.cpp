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
#include <goast/DiscreteShells.h>
#include <goast/GeodesicCalculus.h>
#include <goast/Multiresolution/MultiResOptimization.h>


//! \brief Test class for recursive multilevel computation.
//! \author Heeren
template < typename ConfiguratorType, typename ShellDeformationType>
class RecursiveMultilevelGeodesicOp : public MultiResOptimInterface<ConfiguratorType>{

protected:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const TriMesh& _startMesh, _endMesh;
  RealType _bendWeight;
  int _lengthOfGeodesic;
  int _initModeOnCoarseLevel;
  std::string _saveNameStem;
  
  const OptimizationParameters<ConfiguratorType>& _optPars;
  
public:  
  RecursiveMultilevelGeodesicOp( const TriMesh& StartMesh, 
                                 const TriMesh& EndMesh,
                                 RealType bendWeight,
                                 int lengthOfGeodesic, 
                                 const OptimizationParameters<ConfiguratorType>& optPars,
                                 RealType Theta,
                                 int initModeOnCoarseLevel = 2,
                                 std::string saveNameStem = "" )
  //! CAUTION Assume that we have only two levels since the OptimizationParameters class can only handle two levels so far!
  : MultiResOptimInterface<ConfiguratorType>( 1, Theta ), 
  _startMesh(StartMesh), _endMesh(EndMesh), _bendWeight(bendWeight), _lengthOfGeodesic(lengthOfGeodesic), _initModeOnCoarseLevel(initModeOnCoarseLevel), _saveNameStem(saveNameStem), _optPars(optPars) { }
  
  virtual ~RecursiveMultilevelGeodesicOp(){}
  
  //
  virtual void optimize ( const MeshTopologySaver& Topology, const std::vector<TriMesh>& inputMeshes, std::vector<VectorType>& Path, int level ) const {     
       
    std::cerr << "#nodes on this level = " << Topology.getNumVertices() << std::endl; 
    
    if( level == 0 )
        Path.resize( _lengthOfGeodesic-2 );
    
    if( Path.size() != _lengthOfGeodesic-2 )
        throw BasicException( "RecursiveMultilevelGeodesicOp::optimize(): input vector has wrong size!" );
     
    VectorType startShape, endShape;
    getGeometry( inputMeshes[0], startShape );
    getGeometry( inputMeshes[1], endShape );
      
    // cast input into single vector
    int numLocalDofs = 3 * Topology.getNumVertices();
    int numTotalDofs = (_lengthOfGeodesic-2) * numLocalDofs;
    VectorType concPath( numTotalDofs );   
    if( level > 0 ){
        for( int k = 0; k < _lengthOfGeodesic-2; k++ ){
            if( Path[k].size() != numLocalDofs )
              throw BasicException( "RecursiveMultilevelGeodesicOp::optimize(): input vector component has wrong size!" );
            for( int i = 0; i < numLocalDofs; i++ )
                concPath[k*numLocalDofs+i] = Path[k][i];
        }
    }    
    
    // save initialization?
    if( (_saveNameStem.size() > 0) && (level > 0) ){
      TriMesh temp( Topology.getGrid() );
      for( int k = 0; k < _lengthOfGeodesic-2; k++ ){
        setGeometry( temp, concPath.segment(k*numLocalDofs, numLocalDofs) );
        std::ostringstream filename;
        filename << _saveNameStem << "_init_level" << level << "_" << k+1 << ".ply";
        OpenMesh::IO::write_mesh( temp, filename.str() );  
      }
    }
    
    _optPars.setLevel( level );
    ShellDeformationType W(Topology, _bendWeight );
    GeodesicInterpolation<DefaultConfigurator> interpolationOp( Topology, startShape, endShape, W, _lengthOfGeodesic, _optPars, false );
    
    // take boundary condition into account
    std::vector<int> bdryMask;
    if( Topology.getNumFixedNodes() > 0 ){        
        bdryMask.reserve( Topology.getNumFixedNodes() );
        for( int i = 0; i < Topology.getNumVertices(); i++ )
            if( Topology.isFixedNode(i) )
                bdryMask.push_back( i );
        extendBoundaryMask( Topology.getNumVertices(), bdryMask ); 
        interpolationOp.setBoundaryMask( bdryMask );
    }
    
    // initialize with linear interpoaltion 
    int initMode = ( level == 0 ) ? _initModeOnCoarseLevel : -1;
    interpolationOp.execute( concPath, initMode );
    
    // cast back
    for( int k = 0; k < _lengthOfGeodesic-2; k++ ){
      if( Path[k].size() != numLocalDofs )
        Path[k].resize( numLocalDofs );
      for( int i = 0; i < numLocalDofs; i++ )
        Path[k][i] = concPath[k*numLocalDofs+i];
    }
    
    // save final result?
    if( _saveNameStem.size() > 0 ){
      TriMesh temp( Topology.getGrid() );
      for( int k = 0; k < _lengthOfGeodesic-2; k++ ){
        setGeometry( temp, concPath.segment(k*numLocalDofs, numLocalDofs) );
        std::ostringstream filename;
        filename << _saveNameStem << "_final_level" << level << "_" << k+1 << ".ply";
        OpenMesh::IO::write_mesh( temp, filename.str() );  
      }
    }
    
  }  
  
  //
  virtual void fillInputMeshes( std::vector<TriMesh>& inputMeshes ) const {
    inputMeshes.clear();
    inputMeshes.push_back( _startMesh );
    inputMeshes.push_back( _endMesh );
  }

  
};

/** 
 * \brief Computation of discrete geodesic interpolation in shell space via multiresolution scheme
 * \author Heeren
 * 
 * Geodesic interpolation as described in geodesic_interpolation_shells.cpp but using
 * the multi resolution scheme (with two levels, i.e. coarse level and original level) 
 * proposed in Botsch et al. \cite BoSuPa06
 * 
 */
int main(int argc, char *argv[])
{

try{

  TriMesh mesh1, mesh2;
  OpenMesh::IO::read_mesh(mesh1, "../../data/finger/finger0.ply");
  OpenMesh::IO::read_mesh(mesh2, "../../data/finger/finger1.ply");
    
  MeshTopologySaver Topology( mesh1 );
  double bendWeight = 0.01;
  int lengthOfGeodesic = 5;
  double Theta = 0.05;
  int initModeOnCoarseLevel = 2; // linear interpolation
  
  // define optimization parameters
  OptimizationParameters<DefaultConfigurator> optPars;
  // coarse level
  optPars.setGradientIterations( 200, 0 );
  //optPars.setBFGSIterations( 100, 0 );
  optPars.setNewtonIterations( 50, 0 );
  // fine level
  //optPars.setGradientIterations( 50, 1 );
  //optPars.setBFGSIterations( 25, 1 );
  optPars.setNewtonIterations( 10, 1 );
  //optPars.setQuietMode( SHOW_ALL );
  
  // Dirichlet boundary nodes
  std::vector<int> indices;
  Topology.getVertexNRingVertices( 10, 1302,  indices );
  std::cerr << "#bdry nodes = " << indices.size() << std::endl;
  BitVector mask( Topology.getNumVertices() );
  mask.setAll( false );
  for( int i = 0; i < indices.size(); i++ )
      mask.set( indices[i], true );


  std::cerr << "Start multilevel gedoesic test."  << std::endl;
  typedef ShellDeformation<DefaultConfigurator, NonlinearMembraneDeformation<DefaultConfigurator>, SimpleBendingDeformation<DefaultConfigurator> > ShellDeformationType;                                 
  RecursiveMultilevelGeodesicOp<DefaultConfigurator, ShellDeformationType> op( mesh1, mesh2, bendWeight, lengthOfGeodesic, optPars, Theta, initModeOnCoarseLevel, "geodMultiresTest_fixed" );    
  RecursiveMultilevelOptimization<DefaultConfigurator>( op ).apply( &mask );
 
  
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}