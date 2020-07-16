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

//==============================================================================================================
typedef ShellDeformation<DefaultConfigurator, NonlinearMembraneDeformation<DefaultConfigurator>, SimpleBendingDeformation<DefaultConfigurator> > ShellDeformationType;
typedef typename DefaultConfigurator::VectorType VectorType;

/** 
 * \brief Compute sectional curvature in shell space.
 * \author Heeren
 * 
 */
int main(int argc, char *argv[])
{

try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE SECTIONAL CURVATURE" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
    
    TriMesh posMesh, UMesh, VMesh;
    OpenMesh::IO::read_mesh(posMesh, "../../data/sphere/unitSphere3.ply" );
    OpenMesh::IO::read_mesh(UMesh, "../../data/sphere/unitSphere3_defX.ply" );
    OpenMesh::IO::read_mesh(VMesh, "../../data/sphere/unitSphere3_defY.ply" );    
    std::cerr << "Num of vertices = " << VMesh.n_vertices() << std::endl;
        
    VectorType posGeom, UGeom, VGeom;
    getGeometry( posMesh, posGeom );  
    getGeometry( UMesh, UGeom );
    getGeometry( VMesh, VGeom );
    MeshTopologySaver Topol( posMesh );        

    // set parameters
    double Tau = 0.01; 
    std::vector<int> Mask;   
    double bendingWeight = 0.01;
    ShellDeformationType W( Topol, bendingWeight );    
        
    // set optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setGradientIterations( 25 );
    optPars.setNewtonIterations( 25 );
    //optPars.setVerbose();  
    optPars.setQuiet();
    bool quiet = true;
    optPars.setInitializationScheme( 2 ); 
        
    // compute initial deformation 
    std::cerr << "\nCompute deformation energies: " << std::endl;
    typename DefaultConfigurator::RealType energy;
    W.applyEnergy( posGeom, UGeom, energy );
    std::cerr << "1/2 g_p(u,u) ~ W[p,p+u] = " << energy << std::endl;
    W.applyEnergy( posGeom, VGeom, energy );
    std::cerr << "1/2 g_p(v,v) ~ W[p,p+v] = " << energy << std::endl << std::endl;
        
    // compute tangent vectors as displacements
    std::cerr << "Perform registration..." << std::endl;
    VectorType OptimalRBMdofs;
    MomentumRegistrationOperator<DefaultConfigurator>( Topol, posGeom ).execute( UGeom, OptimalRBMdofs );
    MomentumRegistrationOperator<DefaultConfigurator>( Topol, posGeom ).execute( VGeom, OptimalRBMdofs ); 
    UGeom -= posGeom;        
    VGeom -= posGeom;                    
        
    // start curvature computation
    std::cerr << "Start computation of sectional curvature for mu = " << bendingWeight << " and tau = " << Tau << std::endl;
    auto t_start = std::chrono::high_resolution_clock::now();
    double secCurv = sectionalCurvature<DefaultConfigurator>( Topol, posGeom, UGeom, VGeom, Tau, Mask, W, optPars, quiet );
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << std::fixed << "Done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << "seconds." << std::endl ;     
    std::cerr << "Sectional curvature for tau = " <<  Tau << " is " << secCurv << std::endl;     
  
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}