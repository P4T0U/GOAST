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

//================================================================================

typedef DefaultConfigurator ConfiguratorType;

typedef typename ConfiguratorType::RealType   RealType;
typedef typename ConfiguratorType::VectorType VectorType;
typedef typename ConfiguratorType::SparseMatrixType MatrixType;
typedef typename ConfiguratorType::TripletType TripletType;
typedef typename ConfiguratorType::VecType    VecType;
typedef typename ConfiguratorType::MatType    MatType;
typedef std::vector<TripletType> TripletListType;

static double SCALINGFACTOR = 1.;

//================================================================================
// f(x,y) = 100 * (y - x^2)^2 + (1-x)^2
class TestEnergy : public BaseOp<VectorType, RealType>{
    
public:
  void apply( const VectorType& Arg, RealType & Dest ) const {
      RealType aux1 = Arg[1] - Arg[0]*Arg[0];
      RealType aux2 = 1 - Arg[0];
      Dest = SCALINGFACTOR * ( 100*aux1*aux1 + aux2*aux2 );
  }    
};

class TestGradient : public BaseOp<VectorType, VectorType >{ 
    
public:
  void apply( const VectorType& Arg, VectorType& Dest ) const {
      if( Dest.size() != Arg.size() )
          Dest.resize( Arg.size() );
      Dest.setZero();
      RealType aux1 = Arg[1] - Arg[0]*Arg[0];
      RealType aux2 = 1 - Arg[0];
      Dest[0] -= SCALINGFACTOR * ( 400*aux1*Arg[0] + 2*aux2 );
      Dest[1] += SCALINGFACTOR * ( 200*aux1 );
  }    
};

class TestHessian : public BaseOp<VectorType, MatrixType >{ 

public:
  void apply( const VectorType& Arg, MatrixType& Dest ) const {
      
    if( (Dest.rows() != 2) || (Dest.cols() != 2) )
        Dest.resize( 2, 2 );
    Dest.setZero();  
    
    RealType aux1 = Arg[1] - Arg[0]*Arg[0];
    Dest.coeffRef(0,0) += SCALINGFACTOR * ( 800*Arg[0]*Arg[0] - 400*aux1 + 2 );
    Dest.coeffRef(0,1) -= SCALINGFACTOR * ( 400*Arg[0] ); 
    Dest.coeffRef(1,0) -= SCALINGFACTOR * ( 400*Arg[0] );
    Dest.coeffRef(1,1) += SCALINGFACTOR * ( 200. );
  }    
};
//================================================================================




//#################################
//#################################
int main(int argc, char *argv[])
{

try{

  std::cout << "Hello GOAST! :-)" << std::endl;
  std::cout << "This main runs a couple of basic tests to make sure the three mandatory dependencies are included correctly." << std::endl<< std::endl;

  std::cout << "=====================================================" << std::endl;
  std::cout << "Start mesh test (check OpenMesh)..." << std::endl;
  TriMesh mesh;
  
  if ( !OpenMesh::IO::read_mesh( mesh, "../../data/finger/finger0.ply" )) {
    std::cout << "read error" << std::endl;
    exit( 1 );
  }
  MeshTopologySaver topology( mesh );
  std::cout << "Mesh has " << topology.getNumVertices() << " vertices, " << topology.getNumFaces() << " faces and " << topology.getNumEdges() << " edges!" << std::endl;

  std::cout << "=====================================================" << std::endl;
  std::cout << "Start vector test (check Eigen)..." << std::endl;
  VectorType geometry;
  getGeometry( mesh, geometry );
  VecType node;
  getXYZCoord<VectorType, VecType>( geometry, node, 0 );
  std::cout << "The first node is given by  " << node << std::endl;


  std::cout << "=====================================================" << std::endl;
  std::cout << "Start linear solver test (check SuiteSparse)..." << std::endl;
  int num = 7;
  MatrixType mat(num,num);
  VectorType vec(num), sol(num);

  std::vector<TripletType> tripletList;
  for( int i = 0; i < num; i++ ){
        tripletList.push_back( TripletType(i,i,i+1) );
        vec[i] = 1.;
  }
  mat.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  //std::cout << mat << std::endl;

  LinearSolver<ConfiguratorType> LUSolver;
  LUSolver.prepareSolver( mat );
  LUSolver.backSubstitute( vec, sol );
  std::cout << "Solution of linear system is:\n" << sol << std::endl;


  std::cout << "=====================================================" << std::endl;
  std::cout << "Start energy test..." << std::endl;
  TestEnergy E;  
  TestGradient DE;
  TestHessian D2E;  
  
  VectorType start(2), solution(2);
  start[0] = 1.2; start[1] = 1.2;
  
  std::cout << "Start first derivative test..." << std::endl;
  ScalarValuedDerivativeTester<ConfiguratorType> ( E, DE, 1e-8 ).plotAllDirections ( start, "gradTest");
  
  std::cout << "Start second derivative test..." << std::endl;
  VectorValuedDerivativeTester<ConfiguratorType>( DE, D2E, 1e-5 ).plotAllDirections ( start, "hessTest" );

  
  std::cout << "=====================================================" << std::endl;
  std::cout << "Start optimization tests..." << std::endl;
    
  std::cout << "Start gradient descent with Armijo... " << std::endl;
  GradientDescent<ConfiguratorType> ( E, DE, 1000, 1e-8, ARMIJO, SHOW_TERMINATION_INFO ).solve( start, solution );
  std::cout << solution << std::endl << std::endl;

  std::cout << "Start BFGS with Armijo... " << std::endl;
  int Reset = 50;
  QuasiNewtonBFGS<ConfiguratorType> ( E, DE, 1000, 1e-8, ARMIJO, Reset, SHOW_TERMINATION_INFO ).solve( start, solution );
  std::cout << solution << std::endl << std::endl;

  std::cout << "Start Newton with Armijo..." << std::endl;
  start[0] = 1.2; start[1] = 1.2;
  NewtonMethod<ConfiguratorType> ( DE, D2E, 100, 1e-12, ARMIJO, SHOW_TERMINATION_INFO ).solve( start, solution );
  std::cout << solution << std::endl << std::endl;

  std::cout << "Start Newton with optimal..." << std::endl;
  start[0] = 1.2; start[1] = 1.2;
  NewtonMethod<ConfiguratorType> ( DE, D2E, 100, 1e-12, NEWTON_OPTIMAL, SHOW_TERMINATION_INFO ).solve( start, solution );
  std::cout << solution << std::endl << std::endl;


  } 
  catch ( BasicException &el ){
        std::cout << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}