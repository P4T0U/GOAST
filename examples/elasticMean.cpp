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


//==============================================================================================================
//==============================================================================================================

//! \brief Class to compute elastic average by multiresolution optimization.
//! \author Heeren
template < typename ConfiguratorType, typename ShellDeformationType>
class RecursiveMultilevelElasticAverageOp : public MultiResOptimInterface<ConfiguratorType>{

protected:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const std::vector<TriMesh>& _inputMeshes;
  RealType _bendWeight;
  int _numOfShapes;
  std::string _saveNameStem;  
  const OptimizationParameters<ConfiguratorType>& _optPars;
  
public:
  RecursiveMultilevelElasticAverageOp( const std::vector<TriMesh>& InputMeshes, RealType bendWeight, const OptimizationParameters<DefaultConfigurator>& optPars, RealType Theta, std::string saveNameStem = ""  )  
  //! CAUTION Assume that we have only two levels since the OptimizationParameters class can only handle two levels so far!
  : MultiResOptimInterface<ConfiguratorType>( 1, Theta ), _inputMeshes(InputMeshes), _bendWeight(bendWeight), _numOfShapes(InputMeshes.size()), _saveNameStem(saveNameStem), _optPars(optPars) { }
  
  virtual ~RecursiveMultilevelElasticAverageOp(){}
  
  //
  virtual void optimize ( const MeshTopologySaver& Topology, const std::vector<TriMesh>& inputMeshes, std::vector<VectorType>& Average, int level ) const {     
       
    std::cerr << "#nodes on this level = " << Topology.getNumVertices() << std::endl; 
    
    if( level == 0 )
        Average.resize( 1 );
    
    if( Average.size() != 1 )
        throw BasicException( "RecursiveMultilevelElasticAverageOp::optimize(): input vector has wrong size!" );
     
    int numDOFs = 3 * Topology.getNumVertices();
    // Extract geometry
    VectorType shapes(_numOfShapes * numDOFs);   
    for (auto it = inputMeshes.begin(); it != inputMeshes.end(); ++it) {
      VectorType Geometry;
      getGeometry(*it, Geometry);
      shapes.block((it - inputMeshes.begin()) * numDOFs, 0, numDOFs, 1) = Geometry;
    }
      
    // cast input into single vector
    VectorType mean = Average[0];     
    // initialization on coarse level
    if( level == 0 )
        mean = shapes.segment( 0, numDOFs );

    // save initialization?
    if( _saveNameStem.size() > 0 ){
      TriMesh temp( Topology.getGrid() );
      setGeometry( temp, mean );
      std::ostringstream filename;
      filename << _saveNameStem << "_init_level" << level << ".ply";
      OpenMesh::IO::write_mesh( temp, filename.str() );  
    }
    
    _optPars.setLevel( level );
    ShellDeformationType W(Topology, _bendWeight );
    VectorType alpha = VectorType::Constant(_numOfShapes, 1.);
    ElasticMean<ConfiguratorType> elasticMeanOp(Topology, W, shapes, alpha, _numOfShapes, _optPars );
    
    // take boundary condition into account
    std::vector<int> bdryMask;
    if( Topology.getNumFixedNodes() > 0 ){        
        bdryMask.reserve( Topology.getNumFixedNodes() );
        for( int i = 0; i < Topology.getNumVertices(); i++ )
            if( Topology.isFixedNode(i) )
                bdryMask.push_back( i );
        std::cerr << "\n#bdry nodes on this level = " << bdryMask.size() << std::endl; 
        extendBoundaryMask( Topology.getNumVertices(), bdryMask ); 
        elasticMeanOp.setBoundaryMask( bdryMask );
    }
    
    // perform optimization
    elasticMeanOp.execute( mean );
    
    // cast back
    Average[0] = mean;
    
    // save final result?
    if( _saveNameStem.size() > 0 ){
      TriMesh temp( Topology.getGrid() );
      setGeometry( temp, mean );
      std::ostringstream filename;
      filename << _saveNameStem << "_final_level" << level << ".ply";
      OpenMesh::IO::write_mesh( temp, filename.str() );  
    }
    
  }  
  
  //
  virtual void fillInputMeshes( std::vector<TriMesh>& inputMeshes ) const {
    inputMeshes.clear();
    for( int i = 0; i < _numOfShapes; i++ )
        inputMeshes.push_back( _inputMeshes[i] );
  }

  
};


//==============================================================================================================
typedef ShellDeformation<DefaultConfigurator, NonlinearMembraneDeformation<DefaultConfigurator>, SimpleBendingDeformation<DefaultConfigurator> > ShellDeformationType;
typedef typename DefaultConfigurator::VectorType VectorType;



/** 
 * \brief Computation of elastic mean (either as single or muli resolution method)
 * \author Heeren
 * 
 * For \f$ n = 10 \f$ input shapes \f$ s_1, \ldots, s_n \f$, 
 * compute minimizer of \f$ F[s] = \sum_{i=1}^n W[s_i, s] \f$,
 * where W is an elastic shell deformation energy.
 * 
 * Optimization can be performed on a single level or using a multi resolution scheme
 * which first computes a solution on a coarser level and prolongates the coarse solution to 
 * the original level to obtain a good initialization.
 * 
 * All shapes are supposed to be given as nodal positions of meshes that are in dense correspondence. 
 * 
 */
int main(int argc, char *argv[])
{

try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE ELASTIC MEAN (SINGLE VS. MULTI RESOLUTION)" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
    
    // Read meshes from disc
    std::vector<TriMesh> meshes;
    int numOfShapes = 10;    
    bool multiresOptim = false;

    for (int k = 0; k < numOfShapes; k++) {
        TriMesh mesh;
        if (!OpenMesh::IO::read_mesh(mesh, "../../data/cactus/cactus" + std::to_string(k) + ".ply")) {
                    std::cerr << "read error\n";
                    exit(1);
        }
        meshes.push_back(mesh);
    }

    // Extract topology
    MeshTopologySaver Topol(meshes[0]);
    int numDOFs = 3 * Topol.getNumVertices();
    std::cerr << "#DOFS = " << numDOFs << std::endl;
    
    // Extract geometry
    VectorType shapes(numOfShapes * numDOFs);
    VectorType Geometry;
    for (auto it = meshes.begin(); it != meshes.end(); ++it) {
                getGeometry(*it, Geometry);
                shapes.block((it - meshes.begin()) * numDOFs, 0, numDOFs, 1) = Geometry;
    }
    
    // Create mask for bottom of cactus
    std::vector<int> mask;
    Topol.getVertexNRingVertices(10, 619, mask);
    std::cerr << "#bdry nodes = " << mask.size() << std::endl;    

    double bendWeight = 0.001;
    ShellDeformationType W(Topol, bendWeight);
    
    // optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setGradientIterations( 200 );
    optPars.setBFGSIterations( 100 );
    optPars.setNewtonIterations( 50 );
    optPars.setQuietMode(SHOW_ALL);

    // weights of elastic mean are all set to 1
    VectorType alpha = VectorType::Constant(numOfShapes, 1.);

    // initialize mean with first shape
    VectorType mean;
    getGeometry(meshes[0], mean);
            
    // perform single level or multi level optimization
    if( !multiresOptim ){
      // define elastic mean operato and set boundary mask
      ElasticMean<DefaultConfigurator> elasticMeanOp(Topol, W, shapes, alpha, numOfShapes, optPars);    
      
      extendBoundaryMask(Topol.getNumVertices(), mask);
      elasticMeanOp.setBoundaryMask(mask);
                
      // run optimization
      auto t_start = std::chrono::high_resolution_clock::now();
      std::cerr << "Start SINGLE level optimization for elastic mean computation." << std::endl;
      elasticMeanOp.execute(mean);    
      auto t_end = std::chrono::high_resolution_clock::now();
            std::cout << std::fixed << "Done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count()
            << "seconds." << std::endl;
            
      // save output
      TriMesh output(meshes[0]);
      setGeometry(output, mean);
      if (!OpenMesh::IO::write_mesh(output, "cactus_elastic_mean.ply") ){
                std::cerr << "write error\n";
                exit(1);
    }      
    }
    else{
      // run multiresolution optimization        
      double Theta = 0.1;
      RecursiveMultilevelElasticAverageOp<DefaultConfigurator, ShellDeformationType> multiresElasticMeanOp( meshes, bendWeight, optPars, Theta, "cactus_elastic_mean" );
      RecursiveMultilevelOptimization<DefaultConfigurator> multiResOptim( multiresElasticMeanOp );
      
      // cast boundary mash into BitVector
      BitVector binMask( Topol.getNumVertices() );
      binMask.setAll( false );
      for( int i = 0; i < mask.size(); i++ )
        binMask.set( mask[i], true );  
      
      // reset optim params for two different levels 
      optPars.setGradientIterations( 200, 0 );
      optPars.setBFGSIterations( 100, 0 );
      optPars.setNewtonIterations( 50, 0 );
      optPars.setQuietMode(SHOW_TERMINATION_INFO, 0 );
      
      optPars.setGradientIterations( 25, 1 );
      optPars.setBFGSIterations( 0, 1 );
      optPars.setNewtonIterations( 10, 1 );
      optPars.setQuietMode(SHOW_ALL, 1 );
      
      auto t_start = std::chrono::high_resolution_clock::now();
      std::cerr << "Start MULTI level optimization for elastic mean computation." << std::endl;
      multiResOptim.apply( &binMask );
      auto t_end = std::chrono::high_resolution_clock::now();
            std::cout << std::fixed << "Done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count()
            << "seconds." << std::endl;        
    } 
  
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}