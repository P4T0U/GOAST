// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include <iostream>
#include <string>

#include <goast/Core.h>
#include <goast/Optimization/Objectives.h>
#include <goast/Optimization/Constraints.h>

#include <goast/external/ipoptNonlinearConstraintSolver.h>
#include <goast/external/vtkIO.h>


//================================================================================
// Global typdefs for easier usability
typedef DefaultConfigurator ConfiguratorType;

typedef typename ConfiguratorType::RealType RealType;
typedef typename ConfiguratorType::VectorType VectorType;
typedef typename ConfiguratorType::SparseMatrixType MatrixType;
typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
typedef typename ConfiguratorType::TripletType TripletType;
typedef typename ConfiguratorType::VecType VecType;
typedef typename ConfiguratorType::MatType MatType;
typedef typename ConfiguratorType::TensorType TensorType;
typedef std::vector<TripletType> TripletListType;
//================================================================================

//================================================================================
template<typename T>
int sgn( T val ) {
  return ( T( 0 ) < val ) - ( val < T( 0 ));
}

template<typename ConfiguratorType>
class VolumeFunctional : public ConstraintsOp<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const MeshTopologySaver &_Topology;
  VectorType _nodalAreas;

public:
  VolumeFunctional( const MeshTopologySaver &Topology, const VectorType &Geometry ) : _Topology( Topology ) {
    computeNodalAreas<ConfiguratorType>( _Topology, Geometry, _nodalAreas );

  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    assert( Arg.size() == _Topology.getNumVertices());

    VectorType indicatorFunction( Arg );
    indicatorFunction.array() += 1.;
    indicatorFunction.array() = indicatorFunction.array().pow( 2 );
    indicatorFunction.array() /= 4.;

    Dest.resize( 1 );
    Dest.setZero();
    Dest[0] = _nodalAreas.dot( indicatorFunction );
  }

  int getTargetDimension() const override {
    return 1;
  }
};

template<typename ConfiguratorType>
class VolumeGradient : public ConstraintsGradient<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_Topology;
  VectorType _nodalAreas;

public:
  VolumeGradient( const MeshTopologySaver &Topology, const VectorType &Geometry ) : _Topology( Topology ) {
    computeNodalAreas<ConfiguratorType>( _Topology, Geometry, _nodalAreas );

  }

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    assert( Arg.size() == _Topology.getNumVertices());

    if ( Dest.rows() != 1 || Dest.cols() != _Topology.getNumVertices())
      Dest.resize( 1, _Topology.getNumVertices());
    Dest.setZero();

    TripletListType tripletList;
    tripletList.reserve( _Topology.getNumVertices());

    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.begin(), tripletList.end());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &tripletList ) const override {
    for ( int i = 0; i < _Topology.getNumVertices(); i++ )
      tripletList.emplace_back( 0, i, _nodalAreas[i] * ( Arg[i] + 1. ) / 2. );
  }

  int getTargetDimension() const override {
    return 1;
  }

  int getNNZ() const override {
    return _Topology.getNumVertices();
  }
};

template<typename ConfiguratorType>
class VolumeHessian : public ConstraintsHessian<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::TensorType TensorType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_Topology;
  VectorType _nodalAreas;

public:
  explicit VolumeHessian( const MeshTopologySaver &Topology, const VectorType &Geometry ) : _Topology( Topology ) {
    computeNodalAreas<ConfiguratorType>( _Topology, Geometry, _nodalAreas );
  }


  void apply( const VectorType &Arg, TensorType &Dest ) const override {
    Dest.resize( 1, Arg.size(), Arg.size());
    Dest.setZero();

    std::vector<TripletListType> tripletLists;
    pushTriplets( Arg, tripletLists );

    Dest.setFromTriplets( tripletLists );
  }

  void apply( const VectorType &Arg, std::vector<SparseMatrixType> &Dest ) const {
    Dest.resize( 1 );

    std::vector<TripletListType> tripletLists;
    pushTriplets( Arg, tripletLists );

    Dest[0].resize( Arg.size(), Arg.size());
    Dest[0].setFromTriplets( tripletLists.begin(), tripletLists.end());
  }

  void pushTriplets( const VectorType &Arg, std::vector<TripletListType> &tripletLists ) const {
    tripletLists.resize( 1 );

    for ( int i = 0; i < _Topology.getNumVertices(); i++ )
      tripletLists[0].emplace_back( i, i, _nodalAreas[i] / 2. );
  }

  void pushTriplets( const VectorType &Arg, TripletListType tripletList, const VectorType &Lambda ) const {
    assert( Lambda.size() == 1 );

    for ( int i = 0; i < _Topology.getNumVertices(); i++ )
      tripletList.emplace_back( i, i, Lambda[0] * _nodalAreas[i] / 2. );
  }


  int getTargetDimension() const override {
    return 1;
  }


};

template<typename ConfiguratorType>
class ModicaMortolaFunctional : public ObjectiveOp<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  const MeshTopologySaver &_Topology;
  const MatrixType &_stiffnessMatrix;
  VectorType _nodalAreas;
  const RealType _epsilon;

public:
  ModicaMortolaFunctional( const MeshTopologySaver &Topology, const VectorType &Geometry,
                           const MatrixType &stiffnessMatrix, const RealType epsilon ) : _Topology( Topology ),
                                                                                         _stiffnessMatrix(
                                                                                                 stiffnessMatrix ),
                                                                                         _epsilon( epsilon ) {
    computeNodalAreas<ConfiguratorType>( _Topology, Geometry, _nodalAreas );

  }

  void apply( const VectorType &Arg, RealType &Dest ) const override {
    assert( Arg.size() == _Topology.getNumVertices());


    VectorType nonlinearPart( _Topology.getNumVertices());
    for ( int i = 0; i < _Topology.getNumVertices(); i++ )
      nonlinearPart[i] = ( 9. / 16. ) * ( Arg[i] * Arg[i] - 1 ) * ( Arg[i] * Arg[i] - 1 );

    RealType stiffnessPart = ( Arg.transpose() * _stiffnessMatrix * Arg );

    Dest = ( 1 / _epsilon ) * _nodalAreas.dot( nonlinearPart ) + _epsilon * stiffnessPart;
    Dest *= 1. / 2.;
  }

  int getTargetDimension() const override {
    return 1;
  }
};

template<typename ConfiguratorType>
class ModicaMortolaGradient : public ObjectiveGradient<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  const MeshTopologySaver &_Topology;
  const MatrixType &_stiffnessMatrix;
  VectorType _nodalAreas;
  const RealType _epsilon;

public:
  ModicaMortolaGradient( const MeshTopologySaver &Topology, const VectorType &Geometry,
                         const MatrixType &stiffnessMatrix, const RealType epsilon ) : _Topology( Topology ),
                                                                                       _stiffnessMatrix(
                                                                                               stiffnessMatrix ),
                                                                                       _epsilon( epsilon ) {
    computeNodalAreas<ConfiguratorType>( _Topology, Geometry, _nodalAreas );

  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    assert( Arg.size() == _Topology.getNumVertices());

    if ( Dest.size() != _Topology.getNumVertices())
      Dest.resize( _Topology.getNumVertices());
    Dest.setZero();

    // nonlinear part
    for ( int i = 0; i < _Topology.getNumVertices(); i++ )
      Dest[i] = _nodalAreas[i] * ( 9. / 16. ) * 4. * Arg[i] * ( Arg[i] * Arg[i] - 1 );

    Dest *= 1 / ( 2 * _epsilon );

    // stiffness part
    Dest += _epsilon * _stiffnessMatrix * Arg;
  }


};

template<typename ConfiguratorType>
class ModicaMortolaHessian : public ObjectiveHessian<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_Topology;
  const MatrixType &_stiffnessMatrix;
  VectorType _nodalAreas;
  const RealType _epsilon;

public:
  ModicaMortolaHessian( const MeshTopologySaver &Topology, const VectorType &Geometry,
                        const MatrixType &stiffnessMatrix, const RealType epsilon ) : _Topology( Topology ),
                                                                                      _stiffnessMatrix(
                                                                                              stiffnessMatrix ),
                                                                                      _epsilon( epsilon ) {
    computeNodalAreas<ConfiguratorType>( _Topology, Geometry, _nodalAreas );

  }

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    assert( Arg.size() == _Topology.getNumVertices());

    if ( Dest.rows() != _Topology.getNumVertices() || Dest.cols() != _Topology.getNumVertices())
      Dest.resize( _Topology.getNumVertices(), _Topology.getNumVertices());
    Dest.setZero();

    TripletListType tripletList;
    tripletList.reserve( _Topology.getNumVertices());

    pushDiagonalTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.begin(), tripletList.end());
    Dest += _epsilon * _stiffnessMatrix;
  }


protected:
  void pushDiagonalTriplets( const VectorType &Arg, TripletListType &tripletList ) const {
    for ( int i = 0; i < _Topology.getNumVertices(); i++ )
      tripletList.emplace_back( i, i, _nodalAreas[i] / ( 2 * _epsilon ) * ( 9. / 16. ) * ( 12 * Arg[i] * Arg[i] - 4 ));
  }
};

//================================================================================



//#################################
//#################################
int main( int argc, char *argv[] ) {
  try {
    Eigen::IOFormat CleanFmt( Eigen::StreamPrecision, 0, ", ", ",", "", "", "(", ")" );

    // Default parameters
    std::string file = "../../data/cactus/cactus0.ply";

    std::string outputPrefix = "squares_";

    RealType Volume = .1;

    // ===== reading and preparing data ===== //
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Input data" << std::endl;
    std::cout << "=====================================================" << std::endl;
    if ( argc >= 2 ) {
      file = argv[1];
    }
    if ( argc >= 3 ) {
      outputPrefix = argv[2];
    }

    // Read mesh
    TriMesh mesh;
    if ( !OpenMesh::IO::read_mesh( mesh, file ))
      throw std::runtime_error( "Failed to read file: " + file );

    // Topology of the mesh
    MeshTopologySaver Topology( mesh );
    const int numVertices = Topology.getNumVertices();

    std::cout << " - Mesh size:" << std::endl;
    std::cout << " -- Number of vertices: " << Topology.getNumVertices() << std::endl;
    std::cout << " -- Number of edges: " << Topology.getNumEdges() << std::endl;
    std::cout << " -- Number of faces: " << Topology.getNumFaces() << std::endl;
    std::cout << std::endl;

    // Extract geometry
    VectorType Geometry;
    getGeometry( mesh, Geometry );

    // 1-ring about vertex as initialization
    std::vector<int> ring;
    Topology.getVertexNRingVertices( 2, 0, ring );

    // Color mesh according to stuff:
    std::map<std::string, VectorType> meshColorings;

    // Determine edge lengths
    VectorType edgeLengths;
    getEdgeLengths<DefaultConfigurator>( Topology, Geometry, edgeLengths );

    // Compute (approx.) smallest, biggest, and average diameter of triangle in mesh. This is necessary to properly
    // determine the epsilon in the Modica-Mortola functional.
    RealType minDiameter = std::numeric_limits<RealType>::infinity();
    RealType maxDiameter = -std::numeric_limits<RealType>::infinity();
    RealType avgDiameter = 0.;
    int diameterFace = -1;
    for ( int faceIdx = 0; faceIdx < Topology.getNumFaces(); ++faceIdx ) {
      std::array<int, 3> vertices{ Topology.getNodeOfTriangle( faceIdx, 0 ),
                                   Topology.getNodeOfTriangle( faceIdx, 1 ),
                                   Topology.getNodeOfTriangle( faceIdx, 2 ) };

      std::array<int, 3> edges{ Topology.getEdgeOfTriangle( faceIdx, 0 ),
                                Topology.getEdgeOfTriangle( faceIdx, 1 ),
                                Topology.getEdgeOfTriangle( faceIdx, 2 ) };

      RealType localDiameter = -std::numeric_limits<RealType>::infinity();

      // local diameter = length of longest edge
      for ( int i : { 0, 1, 2 } )
        if ( edgeLengths[edges[i]] > localDiameter )
          localDiameter = edgeLengths[edges[i]];

      if ( localDiameter < minDiameter ) {
        diameterFace = faceIdx;
        minDiameter = localDiameter;
      }
      if ( localDiameter > maxDiameter )
        maxDiameter = localDiameter;
      avgDiameter += localDiameter;
    }
    avgDiameter /= numVertices;
    std::cout << " - Diameter (min): " << std::scientific << minDiameter << std::endl;
    std::cout << " - Diameter (avg): " << std::scientific << avgDiameter << std::endl;
    std::cout << " - Diameter (max): " << std::scientific << maxDiameter << std::endl;
    std::cout << " - Diameter FaceIdx: " << diameterFace << std::endl;


    std::cout << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Functionals and initialization" << std::endl;
    std::cout << "=====================================================" << std::endl;
    MatrixType stiffnessMatrix, massMatrix;
    std::cout << " - Computing stiffness and mass matrix... " << std::flush;
    auto t_start = std::chrono::high_resolution_clock::now();

    computeStiffnessMatrix<ConfiguratorType>( Topology, Geometry, stiffnessMatrix, false );
    computeMassMatrix<ConfiguratorType>( Topology, Geometry, massMatrix, false );

    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done in " << std::chrono::duration<double, std::milli>( t_end - t_start ).count() << "ms."
              << std::endl;

    // Initialization, 1 on a neighborhood of the vertex, -1 everywhere else
    VectorType v = VectorType::Constant( numVertices, -1. );
    for ( int vIdx : ring )
      v[vIdx] = 1.;

    meshColorings["Initialization"] = v;

    // Volume functional
    VolumeFunctional<ConfiguratorType> volFct( Topology, Geometry );
    VolumeGradient<ConfiguratorType> volGrad( Topology, Geometry );
    VolumeHessian<ConfiguratorType> volHess( Topology, Geometry );

    std::cout << " - Initial volume: " << volFct( v )[0] << std::endl;

    // Perimeter functional
    ModicaMortolaFunctional<ConfiguratorType> perFct( Topology, Geometry, stiffnessMatrix, minDiameter );
    ModicaMortolaGradient<ConfiguratorType> perGrad( Topology, Geometry, stiffnessMatrix, minDiameter );
    ModicaMortolaHessian<ConfiguratorType> perHess( Topology, Geometry, stiffnessMatrix, minDiameter );

    std::cout << " - Intial perimeter: " << perFct( v ) << std::endl;

    std::cout << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Derivative tests" << std::endl;
    std::cout << "=====================================================" << std::endl;
//    std::cout << " - Testing derivative of volume functional... " << std::flush;
//    VectorValuedDerivativeTester<DefaultConfigurator>( volFct, volGrad, 1e-8, 1 ).plotRandomDirections( v, 50,
//                                                                                                        "test_volFct" );
//    std::cout << "Done." << std::endl;
//
//    std::cout << " - Testing second derivative of volume functional... " << std::flush;
//    TensorValuedDerivativeTester<DefaultConfigurator>( volGrad, volHess, 1e-8, 1 ).plotRandomDirections( v, 50,
//                                                                                                         "test_volGrad" );
//    std::cout << "Done." << std::endl;
//
//    std::cout << " - Testing derivative of perimeter functional... " << std::flush;
//    ScalarValuedDerivativeTester<DefaultConfigurator>( perFct, perGrad, 1e-8 ).plotRandomDirections( v, 50,
//                                                                                                     "test_perFct" );
//    std::cout << "Done." << std::endl;
//
//    std::cout << " - Testing second derivative of perimeter functional... " << std::flush;
//    VectorValuedDerivativeTester<DefaultConfigurator>( perGrad, perHess, 1e-8 ).plotRandomDirections( v, 50,
//                                                                                                      "test_perGrad" );
//    std::cout << "Done." << std::endl;

    std::cout << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Perimeter optimization" << std::endl;
    std::cout << "=====================================================" << std::endl;

    std::vector<RealType> nonlinearLowerBounds( 1, Volume );
    std::vector<RealType> nonlinearUpperBounds( 1, Volume );

    std::vector<RealType> linearLowerBounds( numVertices, -1. );
    std::vector<RealType> linearUpperBounds( numVertices, 1. );

    // the zeroth-vertex should be in the set
    linearLowerBounds[0] = 1.;

//    IpoptFirstOrderSolver<DefaultConfigurator> Solver( perFct, perGrad, volFct, volGrad, 1000, 1e-8, linearLowerBounds,
//                                                       linearUpperBounds, nonlinearLowerBounds, nonlinearUpperBounds );

    IpoptSecondOrderSolver<DefaultConfigurator, VolumeHessian> Solver( perFct, perGrad, perHess,
                                                                       volFct, volGrad, volHess,
                                                                       5000, 1e-8,
                                                                       linearLowerBounds, linearUpperBounds,
                                                                       nonlinearLowerBounds,
                                                                       nonlinearUpperBounds, 0, 5 );
    t_start = std::chrono::high_resolution_clock::now();
    Solver.solve( v, v );
    t_end = std::chrono::high_resolution_clock::now();

    std::cout << "It took " << std::fixed << std::chrono::duration<double, std::milli>( t_end - t_start ).count()
              << "ms." << std::endl;

    meshColorings["OptimizedPhaseField"] = v;

    for ( int vertexIdx = 0; vertexIdx < numVertices; vertexIdx++ )
      v[vertexIdx] = sgn( v[vertexIdx] );

    meshColorings["ThresholdedPhaseField"] = v;

#ifdef GOAST_WITH_VTK
    saveAsVTP( Topology, Geometry, outputPrefix + "phasefield.vtp", meshColorings );
#else
    std::vector<VectorType> VTKData;
    for ( const auto &mc : meshColorings )
      VTKData.push_back( mc.second );
    saveMeshAsLegacyVTK( Topology, Geometry, VTKData, 1, outputPrefix + "phasefield.vtk" );
#endif
  }
  catch ( std::exception &el ) {
    std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXCEPTION: " << std::endl << el.what() << std::endl
              << std::flush;
  }

  return 0;
}