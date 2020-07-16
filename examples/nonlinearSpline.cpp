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
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Computation of nonlinear splines in shell space 
 * \author Heeren
 * 
 * Corresponds to Fig. 3 in \cite HeRuSc16
 * 
 */
int main(int argc, char *argv[])
{

try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE NONLINEAR SPLINES" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
    
    //! load keyframe  shapes 
    TriMesh mesh1, mesh2, mesh3;
    OpenMesh::IO::read_mesh(mesh1, "../../data/sphere/unitSphere1.ply"); 
    OpenMesh::IO::read_mesh(mesh2, "../../data/sphere/unitSphere1_defX.ply");
    OpenMesh::IO::read_mesh(mesh3, "../../data/sphere/unitSphere1_defY.ply");    
    
    VectorType geom1, geom2, geom3;
    getGeometry( mesh1, geom1 );
    getGeometry( mesh2, geom2 );
    getGeometry( mesh3, geom3 );
    
    MeshTopologySaver Topol( mesh1 );
    RealType bendingWeight = 0.01;
    ShellDeformationType W(Topol, bendingWeight );
    int numLocalDOFs = 3 * Topol.getNumVertices();
  
    // set inner optimization parameters
    OptimizationParameters<DefaultConfigurator> innerOptPars;
    innerOptPars.setGradientIterations( 50 );
    //innerOptPars.setBFGSIterations( 75 );
    innerOptPars.setNewtonIterations( 25 );
    //innerOptPars.setQuietMode( SHOW_ALL );  
    innerOptPars.setQuiet();
  
    // total length of discrete spline curve
    int numShapesAlongDiscreteSpline = 7;
    std::vector<int> FixTimePoints( numShapesAlongDiscreteSpline );
    for( int i = 0; i < numShapesAlongDiscreteSpline; i++ )
        FixTimePoints[i] = 0;
    
    // determine end shapes and middle shape as key frames, i.e. we have two segments
    int middleIdx = (numShapesAlongDiscreteSpline-1)/2;
    FixTimePoints[0] = 1;
    FixTimePoints[middleIdx] = 1;
    FixTimePoints[numShapesAlongDiscreteSpline-1] = 1;
    
    // initialize two segments as linear interpolation
    int lengthOfSegment = (numShapesAlongDiscreteSpline+1)/2;
    VectorType spline( numShapesAlongDiscreteSpline * numLocalDOFs );
    RealType tau = 1. / (lengthOfSegment-1);    
    // first segment
    for( int k = 0; k < lengthOfSegment; k++ )
      for( int i = 0; i < numLocalDOFs; i++ )
         spline[k*numLocalDOFs+i]  = geom1[i] + k * tau * (geom2[i] - geom1[i]);      
    // second segment
    for( int k = 0; k < lengthOfSegment; k++ )
      for( int i = 0; i < numLocalDOFs; i++ )
         spline[(lengthOfSegment+k-1)*numLocalDOFs+i]  = geom2[i] + k * tau * (geom3[i] - geom2[i]);
    
    // set up energy and gradient
    bool periodicBC = false;  
    SplineEnergy<DefaultConfigurator> E( Topol, W, FixTimePoints, periodicBC, innerOptPars );
    SplineGradient<DefaultConfigurator> DE( Topol, W, FixTimePoints, periodicBC, innerOptPars );
    DE.setQuietMode( true );
    
    // get initial energies
    typename DefaultConfigurator::RealType energy;
    E.apply(spline, energy);
    VectorType energies = E.getEnergies( );
    std::cerr << "Initial energies:" << std::endl;
    for(int i = 0; i < energies.size(); i++ )
        std::cerr << energies[i] << " ";
    std::cerr << std::endl << std::endl;
    
    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> outerOptPars;
    outerOptPars.setGradientIterations( 100 );
    outerOptPars.setBFGSIterations( 50 );
    outerOptPars.setQuietMode( SHOW_ALL ); 
    
    // optmization with gradient descent
    std::cerr << "Start gradient descent... " << std::endl;
    VectorType initialization = spline;
    GradientDescent<DefaultConfigurator> GD( E, DE, outerOptPars );  
    GD.solve( initialization, spline );   
      
    // optmization with BFGS
    std::cerr << "Start Quasi-Newton... " << std::endl;
    initialization = spline;
    QuasiNewtonBFGS<DefaultConfigurator> QNM( E, DE, outerOptPars );
    QNM.solve( initialization, spline );
    
    // get final energies
    E.apply(spline, energy);
    energies = E.getEnergies( );
    std::cerr << "Final energies:" << std::endl;
    for(int i = 0; i < energies.size(); i++ )
        std::cerr << energies[i] << " ";
    std::cerr << std::endl << std::endl;
  
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}