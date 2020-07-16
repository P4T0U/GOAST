// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//=================================================================================
//=================================================================================
class TreeCotreeDecomposition {
protected:
  const MeshTopologySaver &_Topology;

  const int _numVertices;
  const int _numEdges;
  const int _numFaces;

  std::vector<int> leftoverEdges;
  std::vector<std::vector<int>> dualBasis;

public:
  explicit TreeCotreeDecomposition( const MeshTopologySaver &Topology ) : _Topology( Topology ),
                                                                          _numVertices( Topology.getNumVertices()),
                                                                          _numEdges( Topology.getNumEdges()),
                                                                          _numFaces( Topology.getNumFaces()) {
    /**
     * Outline:
     * 1. Construct spanning tree T of the primal graph
     * 2. Construct dual spanning tree T* not intersecting T
     * 3. Determine left-over dual edges and connect to root of T*
     */
    int genus = -(_numVertices - _numEdges + _numFaces) / 2 + 1;

    std::cout << "TCD - Genus: " << genus << std::endl;

    //*** Step 0: Create Vertex -> edges list ***//
    std::vector<std::vector<int>> vertexEdges( _numVertices );


    for ( int edgeIdx = 0; edgeIdx < _numEdges; edgeIdx++ ) {
      int v1 = _Topology.getAdjacentNodeOfEdge( edgeIdx, 0 );
      int v2 = _Topology.getAdjacentNodeOfEdge( edgeIdx, 1 );
      vertexEdges[v1].push_back( edgeIdx );
      vertexEdges[v2].push_back( edgeIdx );
    }

    //*** Step 1: Dual spanning tree ***//
    std::vector<bool> faceVisited( _numFaces, false );
    std::vector<bool> edgeUsed( _numEdges, false );
    std::vector<int> previousFace( _numFaces, -1 );

    int rootFace = -1;

    // Queue of the next faces to visit
    std::queue<int> edgeQueue;


    // determine initial face and edge
    for ( int edgeIdx = 0; edgeIdx < _numEdges; edgeIdx++ ) {
      if ( !edgeUsed[edgeIdx] ) {
        if ( _Topology.getAdjacentTriangleOfEdge( edgeIdx, 0 ) == -1 ||
             _Topology.getAdjacentTriangleOfEdge( edgeIdx, 1 ) == -1 )
          continue;

        rootFace = _Topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );

        faceVisited[rootFace] = true;
        edgeQueue.emplace( edgeIdx );

        break;
      }
    }

    std::cout << "TCD - Root Face: " << rootFace << std::endl;

    while ( !edgeQueue.empty()) {
      const int currentEdgeIdx = edgeQueue.front();

      // Determine its edges
      std::array<int, 2> faces = { _Topology.getAdjacentTriangleOfEdge( currentEdgeIdx, 0 ),
                                   _Topology.getAdjacentTriangleOfEdge( currentEdgeIdx, 1 ) };


      for ( int i = 0; i < 2; i++ ) {
        const int faceIdx = faces[i];

        if ( faceIdx == -1 || faceVisited[faceIdx] )
          continue;

        faceVisited[faceIdx] = true;
        previousFace[faceIdx] = faces[(i + 1) % 2];
        edgeUsed[currentEdgeIdx] = true;

        std::array<int, 3> edges = { _Topology.getEdgeOfTriangle( faceIdx, 0 ),
                                     _Topology.getEdgeOfTriangle( faceIdx, 1 ),
                                     _Topology.getEdgeOfTriangle( faceIdx, 2 ) };

        for ( const int edgeIdx : edges )
          if ( !edgeUsed[edgeIdx] )
            edgeQueue.emplace( edgeIdx );

      }
      // Remove the currently visited edge from the queue
      edgeQueue.pop();
    }

    for ( int faceIdx = 0; faceIdx < _numFaces; faceIdx++) {
      if (!faceVisited[faceIdx])
        std::cout << "TCD - not visited face " << faceIdx << std::endl;
    }

    //*** Step 2: Primal spanning tree ***//
    // Keep track of already visited faces and vertices
    std::vector<bool> vertexVisited( _numVertices, false );

    for ( int edgeIdx = 0; edgeIdx < _numEdges; edgeIdx++ ) {
      if ( !edgeUsed[edgeIdx] ) {
        edgeQueue.emplace( edgeIdx );
        vertexVisited[_Topology.getAdjacentNodeOfEdge( edgeIdx, 0 )] = true;
        break;
      }
    }


    while ( !edgeQueue.empty()) {
      // Get the next dual vertex = primal face
      const int currentEdgeIdx = edgeQueue.front();

      // Determine its edges
      std::array<int, 2> vertices = { _Topology.getAdjacentNodeOfEdge( currentEdgeIdx, 0 ),
                                      _Topology.getAdjacentNodeOfEdge( currentEdgeIdx, 1 ) };


      for ( const auto &vertexIdx : vertices ) {
        if ( vertexVisited[vertexIdx] )
          continue;

        vertexVisited[vertexIdx] = true;
        edgeUsed[currentEdgeIdx] = true;

        for ( const int edgeIdx : vertexEdges[vertexIdx] )
          if ( !edgeUsed[edgeIdx] )
            edgeQueue.emplace( edgeIdx );

      }
      // Remove the currently visited edge from the queue
      edgeQueue.pop();
    }

    //*** Step 3: Left over edges ***//
    for ( int edgeIdx = 0; edgeIdx < _numEdges; edgeIdx++ )
      if ( !edgeUsed[edgeIdx] )
        leftoverEdges.emplace_back(edgeIdx);


    std::cout << "TCD - Leftover Edges: ";
    for (const int edgeIdx : leftoverEdges)
      std::cout << edgeIdx << "("
                << _Topology.getAdjacentTriangleOfEdge( edgeIdx, 0 ) << ", "
                << _Topology.getAdjacentTriangleOfEdge( edgeIdx, 1 ) << ") ";
    std::cout << std::endl;

    for ( const int edgeIdx : leftoverEdges ) {
      if ( _Topology.getAdjacentTriangleOfEdge( edgeIdx, 0 ) == -1 ||
           _Topology.getAdjacentTriangleOfEdge( edgeIdx, 1 ) == -1 )
        continue;

      int f1 = _Topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      int f2 = _Topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

      std::vector<int> loop;
      int currentFace = f1;
      while (currentFace != rootFace) {
        loop.push_back(currentFace);
        currentFace = previousFace[currentFace];
      }


      currentFace = f2;
      while (currentFace != rootFace) {
        loop.insert(loop.begin(), currentFace);
        currentFace = previousFace[currentFace];
      }

      loop.insert(loop.begin(), rootFace);


      std::cout << "TCD - Loop " << edgeIdx << ": ";
      for (const int loopIdx : loop)
        std::cout << loopIdx << " ";
      std::cout << std::endl;

      dualBasis.push_back(loop);

    }


  }
};
 
