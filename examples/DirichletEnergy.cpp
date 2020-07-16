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

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Minimization of Dirichlet's energy
 * \author Heeren
 *
 * We optimize Dirichlet's energy by:
 * a) solving the Euler-Lagrange-equation, i.e. a Poisson problem with Dirichlet bc.,
 * b) direct optimization via gradient descent or BFGS.
 *
 * Let \f$ \bar x \f$ be the geometry of a reference mesh \$f M \f$ and \f$ L = L[\bar x] \f$ the corresponding stiffness matrix.
 * Let \f$ x \f$ be the geometry of an arbitrary mesh, and consider \f$ x = \phi(\bar x)\f$ for the (unique!) pw. affine deformation \f$ \phi \f$.
 * Then Dirichlet's energy is given by
 * \f[ E[x] = \frac12 x^T L x = \frac12 \int_M |\nabla \phi \|^2 da \, . \f]
 * 
 */
int main(int argc, char *argv[])
{

try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE ON DIRICHLET ENERGY" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    // load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/plate4SD.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );

    // determine boundary mask and deform part of boundary
    std::vector<int> bdryMask;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        if( coords[0] < 0.05 || coords[0] > 0.95 )
            bdryMask.push_back( i );
        // deform part of boundary
        if( coords[0] > 0.95 ) {
            coords[2] += 0.2;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        }
    }
    int numBdryNodes = bdryMask.size();
    std::cerr << "num of bdry nodes = " << numBdryNodes << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMask );

    // save initialization
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initPlate.ply");

    // FIRST CONSIDER SOLUTION OF EULER-LAGRANGE-EQUATION
    std::cerr << "\n\na) OPTIMIZATION BY SOLVING EULER-LAGRANGE EQUATION" << std::endl;
    // assemble and mask stiffness matrix
    std::cerr << "Set up system matrix" << std::endl;
    typename DefaultConfigurator::SparseMatrixType StiffnessMatrix;
    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomRef, StiffnessMatrix );
    applyMaskToMajor<typename DefaultConfigurator::SparseMatrixType>( bdryMask, StiffnessMatrix );

    // set up right hand side and mask
    VectorType rhs( plateGeomDef );
    rhs.setZero();
    for( int i = 0; i < bdryMask.size(); i++ )
        rhs[bdryMask[i]] = plateGeomDef[bdryMask[i]];

    // set up linear system and solve
    std::cerr << "Set up linear system and solve" << std::endl;
    LinearSolver<DefaultConfigurator> directSolver;
    VectorType solution;
    directSolver.prepareSolver( StiffnessMatrix );
    directSolver.backSubstitute( rhs, solution );

    // get final energy value
    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomRef, StiffnessMatrix );
    VectorType temp = StiffnessMatrix * solution;
    std::cerr << "Final Dirichlet energy = " << 0.5 * temp.dot( solution ) << std::endl;

    // saving
    setGeometry( plate, solution );
    OpenMesh::IO::write_mesh(plate, "solutionDirichlet_EulerLangrange.ply");

    // NOW CONSIDER SOLUTION BY OPTIMIZING DIRICHLET ENERGY
    std::cerr << "\n\nb) OPTIMIZATION BY DIRECT MINIMIZATION" << std::endl;

    // set up energies
    DirichletEnergy<DefaultConfigurator> E( plateTopol, plateGeomRef, true );
    DirichletGradientDef<DefaultConfigurator> DE( plateTopol, plateGeomRef );
    DirichletHessianDef<DefaultConfigurator> D2E( plateTopol, plateGeomRef );
    typename DefaultConfigurator::RealType energy;

    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setGradientIterations( 100 );
    optPars.setBFGSIterations( 50 );
    optPars.setNewtonIterations( 10 );
    optPars.setQuietMode( SHOW_TERMINATION_INFO );
    
    // optmization with gradient descent
    std::cerr << "Start gradient descent... " << std::endl;
    VectorType initialization = plateGeomDef;
    GradientDescent<DefaultConfigurator> GD( E, DE, optPars );
    GD.setBoundaryMask( bdryMask );
    GD.solve( initialization, plateGeomDef );

    // optmization with BFGS
    std::cerr << "Start Quasi-Newton... " << std::endl;
    initialization = plateGeomDef;
    QuasiNewtonBFGS<DefaultConfigurator> QNM( E, DE, optPars );
    QNM.setBoundaryMask( bdryMask );
    QNM.solve( initialization, plateGeomDef );

    // saving
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "solutionDirichlet_DirectMinimization.ply");

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}