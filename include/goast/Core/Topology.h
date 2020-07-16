// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef TOPOLOGY_HH
#define TOPOLOGY_HH


//== INCLUDES =================================================================
#include "OpenMeshIncludes.h"
#include "Auxiliary.h"

#include <queue>

//typedef std::vector<int>      Vec3i;
//typedef std::vector<int>      Vec2i;

//=============================================================================

/**
 * \brief Class to store all topological information in a triangle mesh
 * \author Heeren
 *
 * Class which only stores all important topology information of a given mesh,
 * i.e. indices, neighbors and boundary information
 */
class MeshTopologySaver {

protected:  
  typedef TriMesh TriMeshType;
  static const int IndexNotSet = -1;  
  const TriMeshType& _mesh;

  const int _numOfElems;
  const int _numOfNodes;
  const int _numOfEdges;
  
  // saves global indices of nodes belonging to each triangle
  std::vector< Vec3i > _nodesOfTriangles;
  // saves global indices of neighbourin triangles for each triangle (these are -1 at he boundary )
  std::vector< Vec3i > _neighboursOfTriangles;
  // saves global indices of node in neighboring triangle which is not contained in this
  std::vector< Vec3i > _oppositeNodesOfNeighboringTriangles;
  // saves global indices of edges belonging to each triangle
  std::vector< Vec3i > _edgesOfTriangles;

  // saves global indices of neighbouring triangles for each edge (these are -1 at the boundary )
  std::vector< Vec2i > _neighboursOfEdges;  
  // saves global indices of nodes belonging to each edge
  std::vector< Vec2i > _nodesOfEdges;
  // saves global indices of nodes opposite of each edge (one of them is -1 for boundary edges)
  std::vector< Vec2i > _oppositeNodesOfEdges;
   
  // stores unvalid edges and faces
  BitVector _isEdgeValid, _isFaceValid;
  // fixed nodes (e.g. for Dirichlet boundary conditions) 
  BitVector _isFixedNode; // not related to topological boundary nodes!
  int _numOfFixedNodes;



public:
  MeshTopologySaver( const TriMeshType& mesh )
      : _mesh( mesh ),
        _numOfElems( mesh.n_faces() ),
        _numOfNodes( mesh.n_vertices() ),
        _numOfEdges( mesh.n_edges() ),
        _nodesOfTriangles( _numOfElems ),
        _neighboursOfTriangles( _numOfElems ),
        _oppositeNodesOfNeighboringTriangles( _numOfElems ),
        _edgesOfTriangles( _numOfElems ),
        _neighboursOfEdges( _numOfEdges ),
        _nodesOfEdges( _numOfEdges ),
        _oppositeNodesOfEdges( _numOfEdges ),
        _isEdgeValid( _numOfEdges, true ), 
        _isFaceValid( _numOfElems, true ),
        _isFixedNode( _numOfNodes ),
        _numOfFixedNodes( 0 )
  {    
    // run over all faces
    for ( TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it ){
      // global node indices of current face 
      int locIdx = 0;
      for ( TriMesh::ConstFaceVertexIter cfv_it = mesh.cfv_iter(*f_it); cfv_it.is_valid(); ++cfv_it)
        _nodesOfTriangles[f_it->idx()][locIdx++] = cfv_it->idx();
      
      // indices of neighbours of current face
      // ordering convention: i'th neighboring triangle is opposite of node with local index i
      locIdx = 1;
      for ( TriMesh::ConstFaceHalfedgeIter cfh_it = mesh.cfh_iter(*f_it); cfh_it.is_valid(); ++cfh_it){
        _neighboursOfTriangles[f_it->idx()][locIdx%3] = mesh.is_boundary(*cfh_it) ? IndexNotSet : mesh.face_handle( mesh.opposite_halfedge_handle(*cfh_it) ).idx();
        locIdx++;
      }
      
      // indices of vertex in neighboring face that does not belong to current face
      // ordering convention: i'th entry is node in i'th neighboring triangle that does not belong to faceIdx
      locIdx = 0;
      for (TriMesh::ConstFaceHalfedgeIter cfh_it = mesh.cfh_iter(*f_it); cfh_it.is_valid(); ++cfh_it)      
        _oppositeNodesOfNeighboringTriangles[f_it->idx()][locIdx++] = mesh.opposite_he_opposite_vh( *cfh_it ).idx();
      
      //
      locIdx = 1;
      for (TriMesh::ConstFaceEdgeIter cfe_it = mesh.cfe_iter(*f_it); cfe_it.is_valid(); ++cfe_it) {
        _edgesOfTriangles[f_it->idx()][locIdx] = cfe_it->idx();
        locIdx = (locIdx + 1) % 3;
      }
    }

    // run over all edges
    for (  TriMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it ){
       // define two corresponding half edges
       HEH heh0( mesh.halfedge_handle( *e_it, 0) );
       HEH heh1( mesh.halfedge_handle( *e_it, 1) );
       
       // indices of vertices at current edge
       _nodesOfEdges[e_it->idx()][0] = mesh.from_vertex_handle( heh0 ).idx();
       _nodesOfEdges[e_it->idx()][1] = mesh.to_vertex_handle( heh0 ).idx();
        
       // indices of adjacent faces
       _neighboursOfEdges[e_it->idx()][0] = mesh.face_handle( heh0 ).idx();
       _neighboursOfEdges[e_it->idx()][1] = mesh.face_handle( heh1 ).idx();       
       
       // indices of opposite vertices in adjacent faces
       _oppositeNodesOfEdges[e_it->idx()][0] = mesh.opposite_vh( heh0 ).idx();
       _oppositeNodesOfEdges[e_it->idx()][1] = mesh.opposite_vh( heh1 ).idx();    
     }

  }
  
  const TriMeshType& getGrid() const {
    return _mesh;
  }

  int getNodeOfTriangle( int globalFaceIndex, int localNodeIndex ) const {
    return _nodesOfTriangles[globalFaceIndex][localNodeIndex];
  }
  
  int getEdgeOfTriangle( int globalFaceIndex, int localEdgeIndex ) const {
    return _edgesOfTriangles[globalFaceIndex][localEdgeIndex];
  }

  int getNeighbourOfTriangle( int globalFaceIndex, int localNodeIndex ) const {
    return _neighboursOfTriangles[globalFaceIndex][localNodeIndex];
  }

  int getOppositeNodeOfNeighbouringTriangle( int globalFaceIndex, int localNodeIndex ) const {
    return _oppositeNodesOfNeighboringTriangles[globalFaceIndex][localNodeIndex];
  }

  int getAdjacentNodeOfEdge( int globalEdgeIndex, int localIndex ) const{
    return _nodesOfEdges[globalEdgeIndex][localIndex];
  }

  std::vector<Vec2i> getAdjacentNodesOfEdges( ) const{
    return _nodesOfEdges;
  }
  
  int getAdjacentTriangleOfEdge( int globalEdgeIndex, int localIndex ) const{
    return _neighboursOfEdges[globalEdgeIndex][localIndex];
  }

  int getOppositeNodeOfEdge( int globalEdgeIndex, int localIndex ) const{
    return _oppositeNodesOfEdges[globalEdgeIndex][localIndex];
  }

  std::vector<Vec2i> getOppositeNodesOfEdges( ) const{
    return _oppositeNodesOfEdges;
  }

  int getNumEdges() const {
    return _numOfEdges;
  }

  int getNumFaces() const {
    return _numOfElems;
  }

  int getNumVertices() const {
    return _numOfNodes;
  }

  bool isFixedNode( int i ) const {
    return _isFixedNode[i];
  }
  
  void setFixedNode( int i ) {
    if( !_isFixedNode[i] ){
        _isFixedNode[i] = true;
        _numOfFixedNodes++;
    }
  }
    
  void resetFixedNodes() {
    _isFixedNode.setAll( false );
    _numOfFixedNodes = 0;
  }

  void setFixedNodes( const BitVector& mask ) {
    if( mask.size() != _numOfNodes )
      throw BasicException ( "MeshTopologySaver::setFixedNodes(): mask has wrong size!" );  
    _isFixedNode = mask;
    // determine number of fixed nodes
    _numOfFixedNodes = 0;
    for( int i = 0; i < _numOfNodes; i++ )
        if( _isFixedNode[i] )
            _numOfFixedNodes++;
  } 

  const BitVector& getFixedNodes() const {
    return _isFixedNode;
  }
  
  int getNumFixedNodes() const {
    return _numOfFixedNodes;    
  }
  
  bool isEdgeValid( int i ) const {
    return _isEdgeValid[i];
  }
  
  template<typename IntVectorType>
  void addUnvalidEdges( const IntVectorType& unvalidEdges ) {
    for( int i = 0; i < unvalidEdges.size(); i++ )
      _isEdgeValid.set( unvalidEdges[i], false );
  }
  
  void resetUnvalidEdges() {
    _isEdgeValid.setAll( true );
  }
  
  bool isFaceValid( int i ) const {
    return _isFaceValid[i];
  }
  
  template<typename IntVectorType>
  void addUnvalidFaces( const IntVectorType& unvalidFaces ) {
    for( int i = 0; i < unvalidFaces.size(); i++ )
      _isFaceValid.set( unvalidFaces[i], false );
  }
  
  void resetUnvalidFaces() {
    _isFaceValid.setAll( true );
  }
  
  // only considers one connected component of boundary (containing vertexIndex)!
  void fillPartialBoundaryMask( std::vector<int>& mask, int vertexIndex = -1 ) const {
    mask.clear();
    // get vertex handle at boundary
    TriMesh::VertexHandle vh = TriMesh::InvalidVertexHandle;
    if( vertexIndex < 0 ){      
      for (TriMesh::ConstVertexIter vIter = _mesh.vertices_begin(); vIter != _mesh.vertices_end(); ++vIter)
          if( _mesh.is_boundary ( *vIter ) ){
              vh = *vIter;
              break;
          }      
    }
    else{
      vh = _mesh.vertex_handle(vertexIndex);  
    }
    
    if( vh == TriMesh::InvalidVertexHandle )
        return;
      
    // get hafledge handle at boundary
    TriMesh::ConstVertexIHalfedgeIter hIter ( _mesh, vh );
    for (; !_mesh.is_boundary ( *hIter ); ++hIter ) ;
    
    // run along boundary
    HEH startHEH = *hIter;
    mask.push_back(  _mesh.to_vertex_handle ( startHEH ).idx() );
    HEH bdryHEH = _mesh.next_halfedge_handle ( startHEH );
    while( startHEH != bdryHEH ){
      mask.push_back(  _mesh.to_vertex_handle( bdryHEH ).idx() );
      bdryHEH = _mesh.next_halfedge_handle ( bdryHEH );
    }

  }
  
  void fillFullBoundaryMask( std::vector<int>& mask ) const {
    mask.clear();
    // get vertex handle at boundary
    for (TriMesh::ConstVertexIter cv_it = _mesh.vertices_begin(); cv_it != _mesh.vertices_end(); ++cv_it )
        if( _mesh.is_boundary ( *cv_it ) )
            mask.push_back( cv_it->idx() );

  }
  
  void getVertex1RingVertices( const int idx, std::vector<int>& indices) const {
	indices.resize(0);
	for (TriMesh::ConstVertexVertexIter vv_it = _mesh.cvv_iter(_mesh.vertex_handle(idx)); vv_it.is_valid(); ++vv_it)
		indices.push_back(static_cast<int>(vv_it->idx()));
  }

  void getVertex1RingFaces( const int idx, std::vector<int>& indices) const {
    indices.resize(0);
    for (TriMesh::ConstVertexFaceIter vf_it = _mesh.cvf_iter(_mesh.vertex_handle(idx)); vf_it.is_valid(); ++vf_it)
      indices.push_back(static_cast<int>(vf_it->idx()));
  }

  void getVertexNRingVertices( const int N, const int idx, std::vector<int>& indices) const {
	std::vector< std::set<int> > hierachicalSet(N + 1);
	hierachicalSet[0].insert(idx);

	for (int n = 1; n <= N; n++){
		for (std::set<int>::iterator it = hierachicalSet[n - 1].begin(); it != hierachicalSet[n - 1].end(); ++it){
			std::vector<int> oneRingIndices;
			getVertex1RingVertices( *it, oneRingIndices);
			for (uint j = 0; j < oneRingIndices.size(); j++)
				hierachicalSet[n].insert(oneRingIndices[j]);
		}

		// erase old vertices
		for (std::set<int>::iterator it = hierachicalSet[n - 1].begin(); it != hierachicalSet[n - 1].end(); ++it)
			hierachicalSet[n].erase(*it);
		if (n > 1)
			for (std::set<int>::iterator it = hierachicalSet[n - 2].begin(); it != hierachicalSet[n - 2].end(); ++it)
				hierachicalSet[n].erase(*it);
	}

	// merge all hierachical sets
	indices.resize(0);
	for (int n = 0; n <= N; n++)
		for (std::set<int>::iterator it = hierachicalSet[n].begin(); it != hierachicalSet[n].end(); ++it)
			indices.push_back(*it);
  }
  
  // compute union of all vertex-1-rings of all vertices contained in the index list
  void getVertexList1RingVertices( std::vector<int>& indices ) const {
	std::set<int> setOfVertices;
	for (int n = 0; n < indices.size(); n++){
            std::vector<int> oneRingIndices;
            getVertex1RingVertices( indices[n], oneRingIndices);
            for (uint j = 0; j < oneRingIndices.size(); j++)
                setOfVertices.insert(oneRingIndices[j]);
	}

	// merge all hierachical sets
	indices.resize(0);
	for (std::set<int>::iterator it = setOfVertices.begin(); it != setOfVertices.end(); ++it)
	  indices.push_back(*it);
  }

  Eigen::MatrixXi getConnectivityMatrix() const {
    Eigen::MatrixXi meshF(getNumFaces(), 3);
    for (int faceIdx = 0; faceIdx < getNumFaces(); faceIdx++) {
      meshF(faceIdx, 0) = getNodeOfTriangle(faceIdx, 0);
      meshF(faceIdx, 1) = getNodeOfTriangle(faceIdx, 1);
      meshF(faceIdx, 2) = getNodeOfTriangle(faceIdx, 2);
    }
    return meshF;
  }
};



//=================================================================================
//=================================================================================

// geom = (x,y,z) = (x_1, ..., x_n, y_1, ..., y_n, z_1, ..., z_n)
template<typename VectorType>
void getGeometry( const TriMesh& Mesh, VectorType& geom ){
    
    geom.resize( 3 * Mesh.n_vertices() );
    
    // get pointer to first argument in vertex coordinates list
    typedef TriMesh::Point Point;
    const Point* iter = Mesh.points();
    
    for( uint i = 0; i < Mesh.n_vertices(); i++, iter++ )
        for( int j = 0; j < 3; j++ )
            geom[j*Mesh.n_vertices()+i] = (*iter)[j];            
}

// geom = (x,y,z) = (x_1, ..., x_n, y_1, ..., y_n, z_1, ..., z_n)
template<typename VectorType>
void setGeometry( TriMesh& Mesh, const VectorType& geom ){
    
    if( 3 * Mesh.n_vertices() != (uint)geom.size() )
        throw BasicException( "setGeometry: sizes don't match!" );

    for( TriMesh::VertexIter v_it = Mesh.vertices_begin(); v_it!=Mesh.vertices_end(); ++v_it) {
      TriMesh::Point p;
      for(int j = 0; j < 3; ++j) 
          p[j] = geom[ j*Mesh.n_vertices() + v_it->idx() ];
      Mesh.set_point( *v_it, p );
    }        
}


#endif // TOPOLOGY_HH defined
