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


//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Parametric mean curvature motion
 * \author Heeren
 *
 * Let \f$ x = x(t) \f$ be the time-dependent parametrization of a smooth surface,
 * then the mean curvature formulation reads
 * \f[ \frac{d}{dt} x = - h(x) n(x) \, ,\f]
 * where \f$ h \in \R \f$ and \f$ n \in S^2 \f$ denote the mean curvature and the unit normal vector, respectively.
 *
 * Now let \f$ x \f$ be the geometry of a triangle mesh having n vertices,
 * let \f$ L = L[x] \f$ and \f$ M = M[x] \f$ be the mass and stiffness matrix, respectively.
 * Since we have that the Laplace-Beltrami operator of the embedding is precisely
 * the mean curvature vector, the vertex-wise (integrated) mean curvature vector is given as \f$ h = Lx \f$.
 *
 * Hence an implicit time-stepping scheme of the mean curvature motion reads:
 * Let \f$ x_0 \f$ and \f$ \tau > 0 \f$ be given, then we compute for \f$ k = 0,1,2, \ldots \f$:
 * \f[ M_k x_{k+1} - M_k x_k = - \tau L_k x_{k+1} \, . \f]
 * with \f$ M_k = M[x_k]\f$ and \f$ L_k = L[x_k] \f$.
 * Hence we have to solve the following linear system in every iteration step
 * \f[ ( M_k - \tau L_k ) x_{k+1} = M_k x_k \, . \f]
 * 
 */
int main(int argc, char *argv[])
{

try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE ON MEAN CURVATURE MOTION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    // load noisy Stanford bunny
    TriMesh bunny;
    OpenMesh::IO::read_mesh(bunny, "../../data/noisyStanfordBunny.ply");
    MeshTopologySaver Topol( bunny );
    std::cerr << "num of nodes = " << Topol.getNumVertices() << std::endl;
    VectorType bunnyGeom;
    getGeometry( bunny, bunnyGeom );

    // no boundary considered here
    std::vector<bool> mask( Topol.getNumVertices() );
    for( int i = 0; i < Topol.getNumVertices(); i++ )
        mask[i] = false;
    
    // determine parameters
    int steps = 5;
    RealType tau = 1.e-5;
    
    // implicit time-stepping, i.e. solve (M + tau*L) u^k = Mu^{k-1} for u^k
    for( int i = 0; i < steps; i++ ){

        std::cerr << "Start " << i+1 << "th MCM iteration of " << steps << std::endl;

        typename DefaultConfigurator::SparseMatrixType MassMatrix, StiffnessMatrix;
        // stiffness matrix
        computeStiffnessMatrix<DefaultConfigurator>( Topol, bunnyGeom, StiffnessMatrix, false );
        // full mass matrix
        computeMassMatrix<DefaultConfigurator>( Topol, bunnyGeom, MassMatrix, false );
        // alternatively: lumped mass matrix
        //computeLumpedMassMatrix<DefaultConfigurator>( Topol, bunnyGeom, MassMatrix, false );

        // set up solver to compute (M + tau*L)^{-1}
        LinearSolver<DefaultConfigurator> solver;
        typename DefaultConfigurator::SparseMatrixType SystemMatrix = MassMatrix + tau * StiffnessMatrix;
        solver.prepareSolver( SystemMatrix );

        // apply to (x,y,z)-components separately
        for(int j = 0; j < 3; j++){
            // get jth component
            VectorType geomComp( Topol.getNumVertices() );
            for( int k = 0; k < Topol.getNumVertices(); k++ )
                geomComp[k] = bunnyGeom[ j*Topol.getNumVertices() + k ];
            VectorType rhs = MassMatrix * geomComp;
            solver.backSubstitute( rhs, geomComp );

            // write back component
            for( int k = 0; k < Topol.getNumVertices(); k++ )
                if( !mask[k] )
                    bunnyGeom[j*Topol.getNumVertices() + k] = geomComp[k];
        }
    }

    // saving
    setGeometry( bunny, bunnyGeom );
    OpenMesh::IO::write_mesh(bunny, "bunny_smoothedByMCM.ply");

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}