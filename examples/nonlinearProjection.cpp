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
#include <goast/Multiresolution.h>


typedef ShellDeformation<DefaultConfigurator, NonlinearMembraneDeformation<DefaultConfigurator>, SimpleBendingDeformation<DefaultConfigurator> > ShellDeformationType;

//#################################
//#################################
int main( int argc, char *argv[] ) {

  try {

    TriMesh mesh1, mesh2, mesh3;


    if ( !OpenMesh::IO::read_mesh( mesh1, "../../data/finger/finger0.ply" )) {
      std::cerr << "read error 1\n";
      exit( 1 );
    }

    if ( !OpenMesh::IO::read_mesh( mesh2, "../../data/finger/finger_var.ply" )) {
      std::cerr << "read error 2\n";
      exit( 1 );
    }

    if ( !OpenMesh::IO::read_mesh( mesh3, "../../data/finger/finger1.ply" )) {
      std::cerr << "read error 3\n";
      exit( 1 );
    }

    std::vector<TriMesh> meshes;
    meshes.push_back( mesh1 );
    meshes.push_back( mesh2 );
    meshes.push_back( mesh3 );

    std::cerr << "#nodes = " << meshes[0].n_vertices() << std::endl;

    std::cerr << "Start decimation test." << std::endl;
    SimultaneousDecimater<DefaultConfigurator> simDec( meshes );
    simDec.calcDecimation( 0.1 );
    std::cerr << "decimated" << std::endl;

    TriMesh coarseMesh1 = simDec.getCoarseMesh( 0 );
    TriMesh coarseMesh2 = simDec.getCoarseMesh( 1 );
    TriMesh coarseMesh3 = simDec.getCoarseMesh( 2 );

    double bendWeight = 0.001;

    OptimizationParameters<DefaultConfigurator> optParsProjInner, optParsProjOuter;
    optParsProjInner.setGradientIterations( 25 );
    optParsProjInner.setNewtonIterations( 10 );
    optParsProjInner.setQuiet();
    //optParsProjInner.setQuietMode( SHOW_ALL );

    optParsProjOuter.setGradientIterations( 10 );
    optParsProjOuter.setBFGSIterations( 0 );
    optParsProjOuter.setQuietMode( SHOW_ALL );


    MeshTopologySaver Topol( coarseMesh1 );
    typename DefaultConfigurator::VectorType shape, point, grad;
    typename DefaultConfigurator::RealType Energy;

    int numLocalDofs = 3 * Topol.getNumVertices();
    typename DefaultConfigurator::VectorType principalVariations( 3 * numLocalDofs );
    principalVariations.setZero();
    getGeometry( coarseMesh2, point );
    for ( int i = 0; i < numLocalDofs; i++ )
      principalVariations[i] = point[i];
    getGeometry( coarseMesh1, point );
    for ( int i = 0; i < numLocalDofs; i++ )
      principalVariations[numLocalDofs + i] = point[i];
    getGeometry( coarseMesh3, point );
    for ( int i = 0; i < numLocalDofs; i++ )
      principalVariations[2 * numLocalDofs + i] = point[i];

    getGeometry( coarseMesh2, shape );
    shape += point;
    shape /= 2.;

    point.resize( 2 );
    point[0] = 0.25;
    point[1] = 0.35;

    ShellDeformationType W( Topol, bendWeight );

    //std::cerr << "Apply energy..." << std::endl;
    NonlinearProjectionFunctional<DefaultConfigurator> F( Topol, W, principalVariations, shape, optParsProjInner );
    //F.apply( point, Energy );

    //std::cerr << "Apply gradient..." << std::endl;
    NonlinearProjectionGradient<DefaultConfigurator> DF( Topol, W, principalVariations, shape, optParsProjInner );
    //DF.apply( point, grad );


    std::cerr << "Start first derivative test..." << std::endl;
    auto t_start = std::chrono::high_resolution_clock::now();
//        ScalarValuedDerivativeTester<DefaultConfigurator>(F, DF, 1e-6).plotAllDirections(point, "gradTest");
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << std::fixed << "Done in " << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count()
              << "seconds." << std::endl;

    TriMesh output = Topol.getGrid();
    setGeometry( output, shape );
    OpenMesh::IO::write_mesh( output, "pT_shape_to_project.ply" );

    std::cerr << "Calculating coefficients of projection..." << std::endl;
    t_start = std::chrono::high_resolution_clock::now();

    NonlinearLocalProjection<DefaultConfigurator> projection( Topol, W, principalVariations, shape, optParsProjOuter,
                                                              optParsProjInner, false );
    typename DefaultConfigurator::VectorType projShape;
    projection.execute( point, projShape );

    t_end = std::chrono::high_resolution_clock::now();
    std::cout << std::fixed << "Done in " << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count()
              << "seconds." << std::endl;

    std::cerr << "Calculating projection..." << std::endl;
    t_start = std::chrono::high_resolution_clock::now();

    typename DefaultConfigurator::VectorType extendedAlpha;
    bool validExtension = isAlphaExtensionValid<DefaultConfigurator>( point, extendedAlpha );

    computeElasticAverage<DefaultConfigurator>( Topol, W, principalVariations, extendedAlpha, 3, optParsProjInner, NULL,
                                                shape );

    t_end = std::chrono::high_resolution_clock::now();
    std::cout << std::fixed << "Done in " << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count()
              << "seconds." << std::endl;

    setGeometry( output, shape );
    OpenMesh::IO::write_mesh( output, "pT_projected_shape.ply" );

  }
  catch ( BasicException &el ) {
    std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl
              << std::flush;
  }

  return 0;
}