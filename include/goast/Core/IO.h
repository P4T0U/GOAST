// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IO_H
#define IO_H

#include "Auxiliary.h"
#include "Topology.h"

//==========================================================================================================
// HELPER FUNCTIONS PRIMARILY FOR DERIVATIVE TESTS
//==========================================================================================================

bool runGnuplot(const std::string GnuplotCommandFileName) {
	std::string systemCommand = "gnuplot ";
	systemCommand += GnuplotCommandFileName;
	const bool failed = (system(systemCommand.c_str()) != EXIT_SUCCESS);
	if (failed)
		std::cerr << "runGnuplot: Calling gnuplot returned an error.\n";
	return !failed;
}


template<typename VectorType>
void generateDataFile ( const VectorType &x, const VectorType &y, const std::string &outputFileName ) {
        std::ofstream outDataFile ( outputFileName.c_str()  );
        for( int j=0; j<x.size(); ++j ) 
            outDataFile << std::fixed << std::setprecision( 12 ) << x[j] << " " << y[j] << std::endl;
        outDataFile.close();
}

template<typename VectorType>
void generatePNG( const VectorType &t, const VectorType &x, const VectorType &y, const std::string saveName ) {
        
        std::ostringstream plotname;
        plotname << "tempGnuPlotFile.dat";
        std::ofstream out ( std::string(plotname.str()).c_str()  );
        
        generateDataFile ( t, x, "tempEnergyFile.dat" );
        generateDataFile ( t, y, "tempGateauxFile.dat" );
        
        out << "set terminal png" << std::endl;
        out << "set xlabel \"\"" << std::endl;
        out << "set ylabel \"\"" << std::endl;
        out << "set output \"" << saveName <<"\"" << std::endl;
        out << "plot \"tempGateauxFile.dat\" title \"Gateaux derivative\" w l, \"tempEnergyFile.dat\" title \"Energy\" w l" << std::endl;        
        out.close();        
        runGnuplot ( plotname.str() );
}

//==========================================================================================================
// I/O FUNCTIONS
//==========================================================================================================

template<typename VectorType>
void printVector( const VectorType& Arg, int numEntriesPerLine ) {
    for( unsigned int i = 1; i <= Arg.size(); i++ ){
      std::cerr << std::scientific << std::setw(13) << Arg[i- 1] << " ";
      if( (i > 0) && (i%numEntriesPerLine == 0) ) std::cerr << std::endl;
    }        
    std::cerr << std::endl << std::endl;
}

// export to file
template<typename VectorType>
void writeToFile(const VectorType &array, const std::string filename )
{
    std::ofstream arrayData(filename, std::ios::out); 

    for(int i=0;i<array.size();i++)
        arrayData<<array(i)<<std::endl; //Outputs array to txtFile

    arrayData.close();
}

// read single string from file
std::string readSingleStringFromFile( const std::string filename ){
  std::ifstream readfile ( filename.c_str() );
  if (readfile.is_open()){
    while (! readfile.eof() ){
      std::string line;
      getline (readfile,line);
      if(! line.empty() )
	return line;
      else
        std::cerr << "ERROR in readSingleStringFromFile! Line is empty!\n";
    }
    readfile.close();
  }
  return "";
}

// read single double from file
double readSingleDoubleFromFile( const std::string filename ){
  std::ifstream readfile ( filename.c_str() );
  if (readfile.is_open()){
    while (! readfile.eof() ){
      std::string line;
      getline (readfile,line);
      if(! line.empty() )
	return std::atof( line.c_str() ) ;
      else
        std::cerr << "ERROR in readSingleDoubleFromFile! Line is empty!\n";
    }
    readfile.close();
  }
  return 0.;
}

// write integer list to file
template<typename ListType>
void writeListToFile( const std::string filename, const ListType& list ){    
  std::ofstream writefile(filename, std::ios::out); 
  for( uint i = 0; i < list.size(); i++ )
   writefile  <<  list[i]  <<  std::endl; 
  writefile.close();
}

// read integer list from file (assume there is one integer per line)
void readIntegerListFromFile( const std::string filename, std::vector<int>& list ){
  std::string line;
  std::ifstream readfile (filename.c_str());
  list.clear();

  if (readfile.is_open()){
    while (! readfile.eof() ){
      getline (readfile,line);      
      if(! line.empty() )
          list.push_back( std::atoi( line.c_str() ) );
    }
    readfile.close();
  }      
}

//! \todo Generic type of number (float/double/...)
//! \todo Generic number of dimensions (not just 3)
template<typename VecType, typename RealType=double>
std::vector<VecType> readPositionList(const std::string &filename) {
  std::string line;
  std::ifstream readfile( filename.c_str());

  std::vector<VecType> list;

  if ( readfile.is_open()) {
    while ( !readfile.eof()) {
      getline( readfile, line );
      if ( !line.empty()) {
        std::stringstream ss( line );
        std::string token;
        std::vector<RealType> v;
        while ( std::getline( ss, token, ' ' )) {
          v.push_back( std::stod( token ));
        }
        if (v.size() == 3) {
          list.emplace_back(v[0], v[1], v[2]);
        }
      }
    }
    readfile.close();
  }

  return list;
}

//=================================================================================
// SAVE MESH AS VTK (POSSIBLY WITH DATA)
//=================================================================================

// save mesh as vtk file
template<typename VectorType>
void saveMeshAsLegacyVTK ( const MeshTopologySaver & Topology, const VectorType& Geometry, std::string filename ) {
    std::vector<VectorType> Data;
    saveMeshAsLegacyVTK ( Topology, Geometry, Data, 0, filename );
}

// save mesh with (scalar) data as vtk file
// dataType = 1 is vertex data, dataType = 2 is face data
template<typename VectorType>
void saveMeshAsLegacyVTK ( const MeshTopologySaver & Topology, const VectorType& Geometry, const std::vector<VectorType>& Data, int dataType, std::string filename ) {

    std::ofstream out ( filename.c_str(), std::ios::out );
    int numV = Topology.getNumVertices();
    int numF = Topology.getNumFaces();
    
    out << "# vtk DataFile Version 2.0" << std::endl
        << "written by saveMeshAsLegacyVTK()" << std::endl
        << "ASCII" << std::endl
        << "DATASET POLYDATA" << std::endl
        << "POINTS " << numV << " float" << std::endl;

    // vertex coordinates
    for ( int i = 0; i < numV; i++ ) {
      for ( short j = 0; j < 3; ++j )
        out << ( j == 0 ? "" : " " ) << Geometry[j*numV+i];
      out << std::endl;
    }

    out << "POLYGONS " << numF << " " << 4 * numF << std::endl;
    // triangles' vertex indices
    for (  int i = 0; i < numF; i++ ) {
      out << "3";
      for ( short j = 0; j < 3; ++j )
        out << " "  << Topology.getNodeOfTriangle(i,j);
      out << std::endl;
    }

    // is there any data?
    if ( dataType > 0 ){
      if( dataType == 1 )
        out << "POINT_DATA " << numV << std::endl;
      if( dataType == 2 )
        out << "CELL_DATA " << numF << std::endl;
      // add scalar data 
      for (uint k = 0; k < Data.size(); k++ ) {
        int size = (dataType == 1) ? numV : numF;
        if( Data[k].size() != size )
          throw BasicException("saveMeshAsLegacyVTK(): wrong size of data container!");
        out << "SCALARS Data" << k+1 << " float" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < size; ++i)
          out << Data[k][i] << std::endl;
      }
    }

    out.close();
}


#endif