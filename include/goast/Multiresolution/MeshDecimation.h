// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __MESHDECIMATION_H
#define __MESHDECIMATION_H

#include <goast/Core/OpenMeshIncludes.h>

//==========================================================================================================
// HELPER CLASSES FOR DECIMATION
//==========================================================================================================
/**
 * This class describes the collapse of a single edge e = ( vt, vs ) to one vertex vs_new.
 * The collapse is performed in all of the N meshes that are decimated simultaneously.
 * The index of the vertex vs in the original mesh(es) is stored in vs_orig.
 * The displacements delta_vs = vs_new - vs (and delta_vt = vs_new - vt) are stored as well for every mesh to perform the prolongation afterwards.
 *
 * \author Heeren
 */
template <typename ConfiguratorType>
class Collapse{

typedef typename ConfiguratorType::RealType  RealType;
typedef typename ConfiguratorType::VecType   VecType;
    
public:
Collapse( int N ) : delta_vs( N ), delta_vt( N ){}

Collapse( const Collapse& other ) : vsOrig( other.vsOrig ), vtOrig ( other.vtOrig ) {
  delta_vs.clear();
  delta_vt.clear();
  for( int i = 0; i < other.delta_vs.size(); i++ ){
    delta_vs.push_back( other.delta_vs[i] );
    delta_vt.push_back( other.delta_vt[i] );
  }
}

// index of vs in original meshes (ie. the meshes that are commited in SimulProgMesh<>::calculateDecimation () )
int vsOrig;
// index of vt in original meshes
int vtOrig;
// displacements (needed by prolongation)
std::vector<VecType> delta_vs;
std::vector<VecType> delta_vt;

};


/**
 * The class DecimationInfo stores the mesh decimation of an initial mesh M_0 to a coarser mesh M_k.
 * The decimation has been performed by k edge collapses e=(vt,vs) to vs
 * The indices of the collapsed edges are stores in colIndices, detailed describtion of the collapse is stored in collapses.
 * vertexMap stores the indices of the deleted vertices vt as well as the indices of all vertices in M_k w.r.t. M_0 (see above)
 *
 * \author Heeren
 */
template <typename ConfiguratorType>
class DecimationInfo{

    
public:
DecimationInfo(){}

DecimationInfo( int n ) : vertexMap(n), vertexProjection(n), numCollapses(0) {
  for( int i = 0; i < n; i++ ){
    vertexMap[i] = i;
    vertexProjection[i] = i;
  }
}

void resize( int n ){
  vertexMap.resize( n );
  vertexProjection.resize( n );
  colIndices.clear();
  collapses.clear();
  numCollapses = 0;
  for( int i = 0; i < n; i++ ){
    vertexMap[i] = i;
    vertexProjection[i] = i;
  }
}

void pushBack( const Collapse<ConfiguratorType>& col ){
    collapses.push_back( col );
    numCollapses++;
}

void buildProjection() {
    int numFineNodes = vertexProjection.size();
    int numCoarseNodes = numFineNodes - numCollapses;
    for( int i = 0; i < numFineNodes; i++ )
        vertexProjection[i] = -1;
    for( uint i = 0; i < numCoarseNodes; i++ )
        vertexProjection[vertexMap[i]] = i;
}

std::vector< Collapse<ConfiguratorType> > collapses;
// After k < n decimations vertexMap[i] is  a) the index of vertex i in the original mesh for i < n-k
//                                          b) the index of the vertex vt, that was removed in the jth decimation, for i = n-k+j, j=0,...,k-1
std::vector<int> vertexMap;
// vertexProjection[i] is the vertex index in the coarse mesh, whose position is (likely to be the) closest to vertex i in the original mesh
std::vector<int> vertexProjection;
// stores all indices of the removed edges
std::vector<int> colIndices;

int numCollapses;
};


/**
 * \brief Heap interface
 * \author
 *
 * Copied and edited from class DecimaterT<>::HeapInterface in OpenMesh/Tools/Decimater/DecimaterT.hh
 * Replaced mesh properties _prio and _pos by mesh independent vectors.
 *
 */
template<typename ConfiguratorType>
class HeapInterface {
  
  typedef typename ConfiguratorType::VectorType  VectorType;  
  typedef typename TriMesh::VertexHandle    VertexHandle;
  typedef typename TriMesh::HalfedgeHandle  HalfedgeHandle;  
   
   VectorType& _prio;
   std::vector<int>& _pos;    
    
  public:      
    HeapInterface( VectorType& prio, std::vector<int>& position ) : _prio(prio), _pos(position)  { }

    inline bool less( VertexHandle vh0, VertexHandle vh1 )
    { return _prio[ vh0.idx() ] < _prio[ vh1.idx() ]; }

    inline bool greater( VertexHandle vh0, VertexHandle vh1 )
    { return _prio[ vh0.idx() ] > _prio[ vh1.idx() ]; }

    inline int get_heap_position(VertexHandle vh)
    { return _pos[vh.idx()]; }

    inline void set_heap_position(VertexHandle vh, int pos)
    { _pos[vh.idx()] = pos; }
   
};


//==========================================================================================================
// ERROR METRICS
//==========================================================================================================

//! \brief Abstract base class for simultaneous decimation
//! \author Heeren
template<typename ConfiguratorType>
class SimultaneousErrorMetricBase {
  
protected:  
  typedef typename ConfiguratorType::RealType  RealType;  
    
  typedef typename OpenMesh::Decimater::ModBaseT<TriMesh>::CollapseInfo CollapseInfo;
  typedef typename OpenMesh::Geometry::QuadricT<RealType> QuadricType;
  typedef typename OpenMesh::Vec3d    OpenMeshVecType;
  typedef typename TriMesh::FaceHandle    FaceHandle;
  typedef typename TriMesh::VertexHandle VertexHandle;
  typedef typename TriMesh::HalfedgeHandle HalfedgeHandle;
  typedef typename TriMesh::Point Point;
  
  std::vector<TriMesh>& _opMeshes;
  
public:
  SimultaneousErrorMetricBase( std::vector<TriMesh>& opMeshes ) : _opMeshes( opMeshes ) { }
  
  virtual ~SimultaneousErrorMetricBase() {}
  
  virtual bool hasOptimalPositions() const { return false; }
  
  virtual bool isInitialized() const { return false; }
  
  virtual void initialize() = 0;
  
  virtual float collapsePriority( const CollapseInfo& /*Collapse*/ ) const = 0;

  virtual void preprocessCollapse( const CollapseInfo& /*Collapse*/, std::vector<VertexHandle>& /*support*/ ) = 0;

  virtual void postprocessCollapse( const CollapseInfo& /*Collapse*/, std::vector<VertexHandle>& /*support*/ ) = 0;

};

  
//! \brief Random collapses - does not make sense usually
//! \author Heeren
template<typename ConfiguratorType>
class SimultaneousRandomMetric : public SimultaneousErrorMetricBase<ConfiguratorType> {
  
  typedef typename ConfiguratorType::RealType  RealType;  
  typedef typename OpenMesh::Decimater::ModBaseT<TriMesh>::CollapseInfo CollapseInfo;
   typedef typename TriMesh::VertexHandle VertexHandle;

  int _numOfMeshes; 
  
public:
  SimultaneousRandomMetric( std::vector<TriMesh>& opMeshes ) : SimultaneousErrorMetricBase<ConfiguratorType>( opMeshes ), _numOfMeshes( opMeshes.size() ) { std::srand(std::time(0)); }    
  
  float collapsePriority( const CollapseInfo& Collapse ) const {    
    for( int i = 0; i < _numOfMeshes; i++ )      
       if( !(Collapse.vr.is_valid() && Collapse.vl.is_valid()) )
	 return DBL_MAX;  
    return (RealType)rand() / RAND_MAX;
  }  
  
  bool isInitialized() const { 
    return true;     
  }  
  
  void initialize() {}
  void preprocessCollapse( const CollapseInfo& /*Collapse*/, std::vector<VertexHandle>& /*support*/ ) {}
  void postprocessCollapse( const CollapseInfo& /*Collapse*/, std::vector<VertexHandle>& /*support*/ ) {}
};

/**
 * \brief Quadric error metric for simultaneous mesh decimation of n meshes.
 * \author Heeren
 *
 * Defined as sum of ModQuadricT<> as in OpenMesh/Tools/Decimater/ModQuadricT.hh,
 * which is based on the quadric error metric proposed in Garland and Heckbert, 1997, "Surface Simplification Using Quadric Error Metrics"
 */
template<typename ConfiguratorType>
class SimultaneousQuadricErrorMetric : public SimultaneousErrorMetricBase<ConfiguratorType> {
  
  typedef typename ConfiguratorType::RealType  RealType;  
  typedef typename ConfiguratorType::VecType  VecType;  
  typedef typename ConfiguratorType::MatType  MatType;  
  
  typedef typename OpenMesh::Decimater::ModBaseT<TriMesh>::CollapseInfo CollapseInfo;
  typedef typename OpenMesh::Geometry::QuadricT<RealType> QuadricType;
  typedef typename OpenMesh::Vec3d    OpenMeshVecType;
  typedef typename TriMesh::FaceHandle    FaceHandle;
  typedef typename TriMesh::VertexHandle VertexHandle;
  typedef typename TriMesh::HalfedgeHandle HalfedgeHandle;
  typedef typename TriMesh::Point Point; 
  
  int _numOfMeshes;  
  bool _initialized, _computeOptimalPositions, _extendedPostProceesing;
  RealType _edgeLengthFactor;

  std::vector< OpenMesh::VPropHandleT<QuadricType> > _quadrics;
  std::vector< OpenMesh::EPropHandleT<Point> > _optPos;  

public:

  SimultaneousQuadricErrorMetric( std::vector<TriMesh>& opMeshes ) 
  : SimultaneousErrorMetricBase<ConfiguratorType>( opMeshes ), 
    _numOfMeshes( opMeshes.size() ), 
    _initialized(false), 
    _computeOptimalPositions(false), 
    _extendedPostProceesing(false), 
    _edgeLengthFactor(0.),
    _quadrics( _numOfMeshes ), 
    _optPos( _numOfMeshes ){ }
  
  SimultaneousQuadricErrorMetric( std::vector<TriMesh>& opMeshes, RealType edgeLengthParam, bool computeOptimalPositions, bool extendedPostProceesing ) 
  : SimultaneousErrorMetricBase<ConfiguratorType>( opMeshes ), 
    _numOfMeshes( opMeshes.size() ), 
    _initialized(false), 
    _computeOptimalPositions(computeOptimalPositions), 
    _extendedPostProceesing(extendedPostProceesing), 
    _edgeLengthFactor(edgeLengthParam),
    _quadrics( _numOfMeshes ), 
    _optPos( _numOfMeshes ){ }
  
  void initialize() {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( int i = 0; i < _numOfMeshes; i++ )
      initialize( this->_opMeshes[i], _quadrics[i] );
    
    if( _computeOptimalPositions )
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( int i = 0; i < _numOfMeshes; i++ )
      computeOptimalPositions( this->_opMeshes[i], _quadrics[i], _optPos[i] );
    
    _initialized = true;
  }
  
  bool hasOptimalPositions() const { return _computeOptimalPositions; }
  
  bool isInitialized() const { return _initialized; }
  
  // compute sum of errors
  //NOTE Collapse is some reference CollapseInfo, as only indices of v0 and v1 are called in subfunctions
  float collapsePriority( const CollapseInfo& Collapse ) const {
     
     RealType err = 0.;

     for( int i = 0; i < _numOfMeshes; i++ ){
       QuadricType q = this->_opMeshes[i].property(_quadrics[i], Collapse.v0);
       q += this->_opMeshes[i].property( _quadrics[i], Collapse.v1);
       err += q( this->_opMeshes[i].point(Collapse.v1) );
     }
     
     // add a little bit of edge length
     if( _edgeLengthFactor > 1e-12 ){
         for( int i = 0; i < _numOfMeshes; i++ ){
             Point Pi ( this->_opMeshes[i].point( Collapse.v0 ) ), Pj( this->_opMeshes[i].point( Collapse.v1 ) ); 
             err += _edgeLengthFactor * std::sqrt( lengthSqr( Pi, Pj) );
         }
     }

     return float( (err < DBL_MAX) ? err : -1. );
  }

  // no action required
  void preprocessCollapse( const CollapseInfo& /*Collapse*/, std::vector<VertexHandle>& /*support*/ ) {}
  
  // add quadric of deleted vertex to remaining vertex
  //NOTE Collapse is some reference CollapseInfo, as only indices of v0 and v1 are called in subfunctions
  void postprocessCollapse( const CollapseInfo& Collapse, std::vector<VertexHandle>& support )  {
    if( _extendedPostProceesing )
      return postprocessCollapseExtended( Collapse, support );
      
    for( int i = 0; i < _numOfMeshes; i++ )
      this->_opMeshes[i].property( _quadrics[i], Collapse.v1 ) +=  this->_opMeshes[i].property( _quadrics[i], Collapse.v0 );
  }

protected:

RealType lengthSqr( const Point& p, const Point& q ) const {
  return sqrFnc(p[0]-q[0]) + sqrFnc(p[1]-q[1]) + sqrFnc(p[2]-q[2]);
}

inline RealType sqrFnc( RealType a ) const {
    return a*a;
}
 
  //! recompute all face and vertex  matrices in neighbourhood of collapse
  //NOTE Collapse is some reference CollapseInfo, as only indices of v0 and v1 are called in subfunctions
  void postprocessCollapseExtended( const CollapseInfo& Collapse, std::vector<VertexHandle>& support )  {
     
    typedef typename std::vector<VertexHandle>::iterator SupportIterator;
    
    for( int i = 0; i < _numOfMeshes; i++ ){
       // clear quadrics and recompute
       this->_opMeshes[i].property( _quadrics[i], Collapse.v1 ).clear();
       for ( TriMesh::ConstVertexFaceIter vf_it = this->_opMeshes[i].cvf_iter( Collapse.v1 ); vf_it.is_valid(); ++vf_it )
         updateFaceMatrix( this->_opMeshes[i], _quadrics[i], *vf_it, Collapse.v1 );

       // for 1-ring
       for (SupportIterator s_it = support.begin(), s_end = support.end(); s_it != s_end; ++s_it){
	 this->_opMeshes[i].property( _quadrics[i], *s_it).clear();
         for ( TriMesh::ConstVertexFaceIter vf_it = this->_opMeshes[i].cvf_iter( *s_it ); vf_it.is_valid(); ++vf_it )
           updateFaceMatrix( this->_opMeshes[i], _quadrics[i], *vf_it, *s_it );
      }     

    }
  }

  //TODO documentation how error quadric is built
  void initialize( TriMesh& opMesh, OpenMesh::VPropHandleT<QuadricType>& quadric ) {

    // alloc quadrics
    if (!quadric.is_valid())
      opMesh.add_property( quadric );

    // clear quadrics
    typename TriMesh::VertexIter  v_it  = opMesh.vertices_begin(), v_end = opMesh.vertices_end();
    for (; v_it != v_end; ++v_it)
      opMesh.property( quadric, *v_it).clear();

    // calc (normal weighted) quadric
    VertexHandle vh;
    vh.invalidate();    
    typename TriMesh::FaceIter  f_it  = opMesh.faces_begin(), f_end = opMesh.faces_end();    
    // run over all faces, compute Kp for the face and add to all adjacent node quadrics
    for ( ; f_it != f_end; ++f_it)
      updateFaceMatrix( opMesh, quadric, *f_it, vh );
  }
  
  // Optimal position v_bar is the last row of A^{-1}:
  //         | Q0 Q1 Q2 Q3 | 0 |      | Q0 Q1 Q2 0 | -Q3 |
  // [A|b] = | Q1 Q4 Q5 Q6 | 0 | <->  | Q1 Q4 Q5 0 | -Q6 |
  //         | Q2 Q5 Q7 Q8 | 0 |      | Q2 Q5 Q7 0 | -Q8 |
  //         | 0  0  0  1  | 1 |      | 0  0  0  1 |  1  |
  void computeOptimalPositions( TriMesh& opMesh, 
				const OpenMesh::VPropHandleT<QuadricType>& quadric,
				OpenMesh::EPropHandleT<Point>& optPoints ) {
    
    // alloc properties
    opMesh.add_property( optPoints );
    
    typename TriMesh::EdgeIter v_end(opMesh.edges_end());  
    for ( typename TriMesh::EdgeIter e_it = opMesh.edges_begin(); e_it != v_end; ++e_it ){
      
      HalfedgeHandle heh( 2 * e_it->idx() );      
      VertexHandle vh0 = opMesh.from_vertex_handle( heh );
      VertexHandle vh1 = opMesh.to_vertex_handle( heh );
      
      // define sum of two quadrics
      const QuadricType& q0 = opMesh.property(quadric, vh0);
      const QuadricType& q1 = opMesh.property(quadric, vh1);
      
      RealType xx( q0.xx() + q1.xx() ), xy( q0.xy() + q1.xy() ), xz( q0.xz() + q1.xz() ), xw( q0.xw() + q1.xw() );
      RealType yy( q0.yy() + q1.yy() ), yz( q0.yz() + q1.yz() ), yw( q0.yy() + q1.yy() );
      RealType zz( q0.zz() + q1.zz() ), zw( q0.zw() + q1.zw() );
  
      // compute optimal position
      RealType detA = xx*yy*zz + 2.*xy*xz*yz - xz*xz*yy - yz*yz*xx - xy*xy*zz;
      VecType v_bar;
      if( (std::abs( detA ) > std::numeric_limits<RealType>::epsilon()) && !opMesh.is_boundary(vh0) && !opMesh.is_boundary(vh1) ){       
        MatType A( xx, xy, xz, xy, yy, yz, xz, yz, zz );
        MatType Ainv(A.inverse());
        VecType b( -xw, -yw, -zw );
        A.mult( b, v_bar );
      } 
      else{
        Point p = opMesh.point( vh0 );
        Point q = opMesh.point( vh1 );
        RealType lambda = opMesh.is_boundary( vh0 ) ? 1.0 : (opMesh.is_boundary( vh1 ) ? 0.0 : 0.5 );
        for( int i = 0; i < 3; i++ )
          v_bar[i] = lambda * p[i] + (1. - lambda) * q[i];
      }
      
      for( int j = 0; j < 3; j++ )
        opMesh.property( optPoints, *e_it )[j] = v_bar[j];      

    }
 
  }
  
  //
  void updateFaceMatrix( TriMesh& opMesh, OpenMesh::VPropHandleT<QuadricType>& quadric, const FaceHandle& fh, const VertexHandle& vh ) {
    VertexHandle vh0, vh1, vh2;
    typename TriMesh::FaceVertexIter fv_it = opMesh.fv_iter( fh );
    vh0 = *fv_it;  ++fv_it;
    vh1 = *fv_it;  ++fv_it;
    vh2 = *fv_it;

    // get nodal positions of triangle nodes
    OpenMeshVecType v0, v1, v2;
    v0 = OpenMesh::vector_cast<OpenMeshVecType>(opMesh.point(vh0));
    v1 = OpenMesh::vector_cast<OpenMeshVecType>(opMesh.point(vh1));
    v2 = OpenMesh::vector_cast<OpenMeshVecType>(opMesh.point(vh2));

    // compute normal
    OpenMeshVecType n = (v1-v0) % (v2-v0);
    double area = n.norm();
    if (area > FLT_MIN) {
      n /= area;
      area *= 0.5;
    }

    // Kp ist determined by normal and distance: a|b := a*b
    double dist = -(OpenMesh::vector_cast<OpenMeshVecType>(opMesh.point(vh0))|n);
    OpenMesh::Geometry::Quadricd q( n[0], n[1], n[2], dist );
    q *= area;

    // add to all?
    if( vh.is_valid() ){
      opMesh.property(quadric, vh) += q;
    }
    else{
      // add Kp to all adjacent node quadrics
      opMesh.property(quadric, vh0) += q;
      opMesh.property(quadric, vh1) += q;
      opMesh.property(quadric, vh2) += q;
    }
  }


};

//! \brief Error metric for simultaneous mesh decimation of n meshes based on dihedral angles
//! \author Heeren
//! \todo need to add non-trivial postprocessCollapse function!
template<typename ConfiguratorType>
class SimultaneousCurvatureErrorMetric : public SimultaneousErrorMetricBase<ConfiguratorType> {
  
  typedef typename ConfiguratorType::RealType  RealType;    
  typedef typename OpenMesh::Decimater::ModBaseT<TriMesh>::CollapseInfo CollapseInfo;
  typedef typename TriMesh::VertexHandle VertexHandle;
  typedef typename TriMesh::Point Point; 
 
  int _numOfMeshes;  

public:
  SimultaneousCurvatureErrorMetric( std::vector<TriMesh>& opMeshes ) 
  : SimultaneousErrorMetricBase<ConfiguratorType>( opMeshes ), _numOfMeshes( opMeshes.size() ){ }
  
  bool hasOptimalPositions() const { return false; }
  
  bool isInitialized() const { return true; }
  
  // compute sum of errors
  //NOTE Collapse is some reference CollapseInfo, as only indices of v0 and v1 are called in subfunctions
  float collapsePriority( const CollapseInfo& Collapse ) const {
     
     RealType err = 0.;

     for( int i = 0; i < _numOfMeshes; i++ ){
       
       if( !(Collapse.vr.is_valid() && Collapse.vl.is_valid()) )
	 return DBL_MAX;
       
       Point Pi ( this->_opMeshes[i].point( Collapse.v0 ) ), Pj( this->_opMeshes[i].point( Collapse.v1 ) );
       Point Pk ( this->_opMeshes[i].point( Collapse.vr) ), Pl( this->_opMeshes[i].point( Collapse.vl ) );
       
       RealType vol = getArea( Pi, Pj, Pk );
       vol += getArea( Pi, Pj, Pl );
       
       err += sqrFnc( getDihedralAngle( Pi, Pj, Pk, Pl ) ) * lengthSqr( Pi, Pj) / vol;
     }

     return float( (err < DBL_MAX) ? err : -1. );
  }

  void initialize() {}
  
  // no pre/postprocess here
  void preprocessCollapse( const CollapseInfo& /*Collapse*/, std::vector<VertexHandle>& /*support*/ )  {  }
  //TODO should be filled actually!
  void postprocessCollapse( const CollapseInfo& /*Collapse*/, std::vector<VertexHandle>& /*support*/ )  {  }
 
protected:

RealType dotProduct( const Point& p, const Point& q ) const {
  return p[0]*q[0] + p[1]*q[1] + p[2]*q[2];
}

RealType lengthSqr( const Point& p, const Point& q ) const {
  return sqrFnc(p[0]-q[0]) + sqrFnc(p[1]-q[1]) + sqrFnc(p[2]-q[2]);
}

RealType normSqr( const Point& p ) const {
  return dotProduct( p, p);
}

void makeCrossProduct( const Point& a, const Point& b, Point& axb, bool normalize ) const {
  axb[0] = a[1] * b[2] - a[2] * b[1]; 
  axb[1] = a[2] * b[0] - a[0] * b[2];
  axb[2] = a[0] * b[1] - a[1] * b[0];
  if( normalize ){
    RealType norm = std::sqrt( normSqr(axb) );
    for( int i = 0; i < 3; i++ )
      axb[i] /= norm;
  }
}
  
RealType getArea( const Point& Pi, const Point& Pj, const Point& Pk ) const {
  Point temp;   
  makeCrossProduct( Pk-Pj, Pi-Pk, temp, false );
  return 0.5 * std::sqrt( normSqr(temp) );
}

// returns dihedral angle between n1 and n2 at an edge e
RealType getDihedralAngle( const Point& nk, const Point& nl, const Point& e ) const {
    Point crossprod;
    makeCrossProduct( nk, nl, crossprod, false );
    //return std::asin( temp*e / e.norm() );
    return dotProduct(crossprod, e) < 0. ? -std::acos( min( dotProduct(nk,nl), 1. ) ) : std::acos( min( dotProduct(nk,nl), 1. ) );
}

//
RealType getDihedralAngle( const Point& Pi, const Point& Pj, const Point& Pk, const Point& Pl ) const {
    Point nk, nl;
    makeCrossProduct( Pk-Pj, Pi-Pk, nk, true );
    makeCrossProduct( Pl-Pi, Pj-Pl, nl, true );
    return getDihedralAngle( nk, nl, Pj-Pi );
}

inline RealType sqrFnc( RealType a ) const {
    return a*a;
}

};


//! \brief Error metric for simultaneous mesh decimation of n meshes  based on edge lengths
//! \author Heeren
template<typename ConfiguratorType>
class SimultaneousEdgeLengthErrorMetric : public SimultaneousErrorMetricBase<ConfiguratorType> {
  
  typedef typename ConfiguratorType::RealType  RealType;    
  typedef typename OpenMesh::Decimater::ModBaseT<TriMesh>::CollapseInfo CollapseInfo;
  typedef typename TriMesh::VertexHandle VertexHandle;
  typedef typename TriMesh::Point Point; 
 
  int _numOfMeshes;  

public:
  SimultaneousEdgeLengthErrorMetric( std::vector<TriMesh>& opMeshes ) 
  : SimultaneousErrorMetricBase<ConfiguratorType>( opMeshes ), _numOfMeshes( opMeshes.size() ){ }
  
  bool hasOptimalPositions() const { return false; }
  
  bool isInitialized() const { return true; }
  
  // compute sum of errors
  //NOTE Collapse is some reference CollapseInfo, as only indices of v0 and v1 are called in subfunctions
  float collapsePriority( const CollapseInfo& Collapse ) const {
     
     RealType err = 0.;

     for( int i = 0; i < _numOfMeshes; i++ ){
       
       if( !(Collapse.vr.is_valid() && Collapse.vl.is_valid()) )
	 return DBL_MAX;
       
       Point Pi ( this->_opMeshes[i].point( Collapse.v0 ) ), Pj( this->_opMeshes[i].point( Collapse.v1 ) );       
       err += std::sqrt( lengthSqr( Pi, Pj) );
     }

     return float( (err < DBL_MAX) ? err : -1. );
  }

  void initialize() {}
  
  // no preprocess here
  void preprocessCollapse( const CollapseInfo& /*Collapse*/, std::vector<VertexHandle>& /*support*/ )  {  }
  
  //
  void postprocessCollapse( const CollapseInfo& Collapse, std::vector<VertexHandle>& support )  {  
//     for( int i = 0; i < _numOfMeshes; i++ ){
//         Point Pi ( this->_opMeshes[i].point( Collapse.v1 ) ), Pj( this->_opMeshes[i].point( Collapse.v0 )  );
//         double length = std::sqrt( lengthSqr( Pi, Pj) );
//         for (uint j = 0; j < support.size(); j++ ){
//             if( (support[j] == Collapse.vr) || (support[j] == Collapse.vl) )
//                 continue;
//             this->_opMeshes[i].property( _quadrics[i], support[j] ) +=  length;
//         }
//     }
  }
 
protected:

RealType lengthSqr( const Point& p, const Point& q ) const {
  return sqrFnc(p[0]-q[0]) + sqrFnc(p[1]-q[1]) + sqrFnc(p[2]-q[2]);
}

inline RealType sqrFnc( RealType a ) const {
    return a*a;
}

};

//==========================================================================================================
// SIMULTANEOUS DECIMATION OPERATOR
//==========================================================================================================

/**
 * \brief Simultaneous decimation of several TriMeshes
 * \author Heeren
 *
 * Meshes are supposed to be in dense correspondence which is kept throughout the decimation process.
 *
 * Decimation is based on subsequently collapsing edges e = (vs, vt) to vs (and update position of vs afterwards if required)
 *
 * In order to determine edge to be collapsed next, a vertex heap is constructed sorting all remaining vertices according to some ErrorMetricType.
 *
 * ErrorMetricType is specified as template argument and has to be derived from SimultaneousErrorMetricBase, see documentation above.
 *
 * Each vertex vs in the vertex heap specifies an edge e = (vs, vt) which is locally cheapest to remove.
 *
 * Usage: call calcDecimation() or calcDecimationTp(), see documentation below of these functions.
 *
 * Afterwards, decimated meshes can be accessed by getCoarseMesh() or getCoarseMeshes(), respectively.
 */
template < typename ConfiguratorType, typename QuadricType = SimultaneousQuadricErrorMetric<ConfiguratorType> >
class SimultaneousDecimater {
  
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  typedef typename TriMesh::VertexHandle    VertexHandle;
  typedef typename TriMesh::FaceHandle    FaceHandle;
  typedef typename TriMesh::HalfedgeHandle  HalfedgeHandle;  
  typedef typename TriMesh::Point Point;
  typedef typename OpenMesh::Decimater::CollapseInfoT<TriMesh>          CollapseInfo;
  typedef typename OpenMesh::Utils::HeapT<VertexHandle, HeapInterface<ConfiguratorType> >  DecimationHeap;

  int _numOfMeshes;
  int _numOfOriginalVertices;
  std::vector<TriMesh> _opMeshes;
  QuadricType _quadric; 
  
  std::unique_ptr<DecimationHeap> _heap;

  VectorType _prio;
  std::vector<int> _pos;
  std::vector<HalfedgeHandle> _collapseTarget;
  std::vector<bool> _fixedVertices;
  
  bool _quiet, _updatePositions;
  
public: 
  /// Constructor
  SimultaneousDecimater( const std::vector<TriMesh>& meshes ) 
  :  _numOfMeshes( meshes.size() ), _numOfOriginalVertices( meshes[0].n_vertices() ), _opMeshes( meshes ), _quadric( _opMeshes ), //_heap(NULL),
  _prio( _numOfOriginalVertices ), _pos( _numOfOriginalVertices ), _collapseTarget( _numOfOriginalVertices ), _fixedVertices(_numOfOriginalVertices), _quiet( false ), _updatePositions( false ) {   
    _quadric.initialize();
  }  
  
  SimultaneousDecimater( const std::vector<TriMesh>& meshes, const std::vector<int>& FixedVertices ) 
  :  _numOfMeshes( meshes.size() ), _numOfOriginalVertices( meshes[0].n_vertices() ), _opMeshes( meshes ), _quadric( _opMeshes ), //_heap(NULL),
  _prio( _numOfOriginalVertices ), _pos( _numOfOriginalVertices ), _collapseTarget( _numOfOriginalVertices ), _fixedVertices(_numOfOriginalVertices), _quiet( false ), _updatePositions( false ) {   
    _quadric.initialize();
    for( uint i = 0; i < FixedVertices.size(); i++ )
        _fixedVertices[FixedVertices[i]] = true;
  } 
  
  // after edge collapse (vs, vt) -> vs, we set the position p[vs] of vs by p[vs] = 0.5 * (p[vs] + p[vt]) 
  void updatePositions() {
    _updatePositions = true;
  }  

  // Decimate until #vertices = theta * n, where n is the number of vertices of the original mesh
  // Decimation information (needed e.g. for prolongation) is tored in decInfo.
  void calcDecimation( RealType theta, DecimationInfo<ConfiguratorType>& decInfo ) {
    if( !(theta < 1.) )
      throw BasicException ( "SimultaneousDecimater::calcDecimation(): theta should be in (0,1)" );
    calcDecimationTo( std::floor(theta * _numOfOriginalVertices), decInfo );
  }
  
  // decimate until #vertices = theta * n, where n is the number of vertices of the original mesh
  void calcDecimation( RealType theta ){
    DecimationInfo<ConfiguratorType> decInfo( _numOfOriginalVertices );
    calcDecimationTo( std::floor(theta * _numOfOriginalVertices), decInfo );
  }
  
  // decimate until #vertices = targetNumOfVertices, where n is the number of vertices of the original mesh
  // Decimation information (needed e.g. for prolongation) is tored in decInfo.
  void calcDecimationTo( int targetNumOfVertices, DecimationInfo<ConfiguratorType>& decInfo ) {
    decimate( _numOfOriginalVertices - targetNumOfVertices, decInfo );
  }
  
  // decimate until #vertices = targetNumOfVertices, where n is the number of vertices of the original mesh
  void calcDecimationTo( int targetNumOfVertices ){
    DecimationInfo<ConfiguratorType> decInfo( _numOfOriginalVertices );
    calcDecimationTo( targetNumOfVertices, decInfo );
  }  
    
  // construct MeshType from TriMesh
  void getCoarseMeshes( std::vector<TriMesh>& coarseMeshes ) const {
    coarseMeshes.clear();
    for( int j = 0; j < _numOfMeshes; j++ )
      coarseMeshes.push_back( getCoarseMesh(j) );
  }
  
  // construct TriMesh from TriMesh
  const TriMesh& getCoarseMesh( int idx ) const {
    if( !( idx < _opMeshes.size() ) )
      throw BasicException( "SimultaneousDecimater::getCoarseMesh(): index out of bounds!" );
    return _opMeshes[idx];   
  }
  
  // check for bad angles and perform edge flip
  int flipEdgesWithBadAngles( RealType angleThreshold = 3.0 ) {
    int numOfFlippedEdges = 0;
    // run over all meshes
    for( int idx = 0; idx < _numOfMeshes; idx++ ){      
      // run over all edges
      for(typename TriMesh::EdgeIter  e_it = _opMeshes[idx].edges_begin(); e_it != _opMeshes[idx].edges_end(); ++e_it ) {
        
	// no edge flipping possible at boundary
	if( _opMeshes[idx].is_boundary( *e_it ) )
	  continue;
	
        bool flipEdge = false;
        
	// get first half edge
        HalfedgeHandle heh1( 2 * e_it->idx() );
	RealType angle1 = _opMeshes[idx].calc_sector_angle( _opMeshes[idx].next_halfedge_handle(heh1) );	
        if( angle1 > angleThreshold )
            flipEdge = true;
        
        // get second half edge
        HalfedgeHandle heh2( 2 * e_it->idx() + 1 );
	RealType angle2 = _opMeshes[idx].calc_sector_angle( _opMeshes[idx].next_halfedge_handle(heh2) );
        if( angle2 > angleThreshold )
            flipEdge = true;
        
        // do nothing if angle is smaller then threshold
        if( !flipEdge )           
	  continue;
	
        // check whether flip is ok
        bool flipIsOk = true;
        for( int j = 0; j < _numOfMeshes; j++ )
	  if( !_opMeshes[j].is_flip_ok( *e_it ) )
	    flipIsOk = false;
	 
        // perform edge flip
        if( flipIsOk ){      
          // Flip edge in all meshes
	  for( int j = 0; j < _numOfMeshes; j++ )
	    _opMeshes[j].flip( *e_it );
	  numOfFlippedEdges++;
        }
        else{
	  std::cerr << "WARNING: flipping failed!" << std::endl;
        }
      } // end edges
    } // end meshes
    return numOfFlippedEdges;
  }
  
  // compute coarse boundary mask from fine boundary mask
  static void computeDecimatedBoundaryMask( const DecimationInfo<ConfiguratorType>& decInfo, const BitVector& fineMask, BitVector& coarseMask ) {
    int numOfCoarseVertices = fineMask.size() - decInfo.collapses.size();
    coarseMask.resize( numOfCoarseVertices );
    coarseMask.setAll( false );
    for( int i = 0; i < numOfCoarseVertices; i++ )
      if( fineMask[decInfo.vertexMap[i]] )
	coarseMask.set( i, true );
  }

protected:  
  // Decimate (perform maxNumOfCollapses collapses).
  // Copied and edited from DecimaterT<Mesh>::decimate() in OpenMesh/Tools/Decimater/DecimaterT.hh,
  void decimate( int maxNumOfCollapses, DecimationInfo<ConfiguratorType>& decInfo ) {
    
    if( !_quadric.isInitialized() )
      throw BasicException( "SimultaneousDecimater::decimate(): quadric not initialized" );
     
    VertexHandle vp;
    HalfedgeHandle v0v1;
  
    decInfo.resize( _numOfOriginalVertices );
    int numCollapses(0);
    std::vector<int> idxOfCollapsedNode( maxNumOfCollapses );

    typedef std::vector<VertexHandle> Support;
    typedef typename Support::iterator SupportIterator;

    // support of the vertex that is deleted is given by its vertex 1-ring
    Support support(15);
    SupportIterator s_it, s_end;    
    
    // Note: although we define a non-const reference on the first mesh here, this is not changed!
    TriMesh& refMesh = _opMeshes[0];
    
    // initialize heap
    HeapInterface<ConfiguratorType> HI( _prio, _pos );
    _heap = std::unique_ptr<DecimationHeap>( new DecimationHeap(HI) );
    _heap->reserve( refMesh.n_vertices() );

    // fill heap with vertices
    // Note: Although we perform edge collapse the heap is filled with vertices.
    //       Then for each vertex the cheapest edge is determined.  
    typename TriMesh::VertexIter v_end(refMesh.vertices_end());  
    for (typename TriMesh::VertexIter v_it = refMesh.vertices_begin(); v_it != v_end; ++v_it) {
      _heap->reset_heap_position( *v_it );
      if ( !refMesh.status(*v_it).deleted() )
        insertVertexIntoHeap( refMesh, *v_it );
    }
    if( _heap->empty() )
      throw BasicException ( "SimultaneousDecimater::decimate(): heap is empty!" );    

    // process heap
    if( !_quiet ) std::cerr << "Start decimating " << maxNumOfCollapses << " collapses..." << std::endl;   
    while ( (!_heap->empty()) && (numCollapses < maxNumOfCollapses) ) {
      
      // get 1st heap entry, which is a vertex
      vp = _heap->front();
      v0v1 = _collapseTarget[ vp.idx() ];
      _heap->pop_front();

      // setup collapse info
      CollapseInfo ci( refMesh, v0v1);

      // check topological correctness AGAIN !
      if( !isCollapseLegal( refMesh, ci) )
        continue;

      // store support (= one ring of *vp)      
      support.clear();
      for (typename TriMesh::VertexVertexIter vv_it = refMesh.vv_iter(ci.v0); vv_it.is_valid(); ++vv_it)
        support.push_back( *vv_it );
      
      // pre-process collapse (if necessary!)
      _quadric.preprocessCollapse( ci, support );
      
      // store prolongation info
      Collapse<ConfiguratorType> col( _numOfMeshes );
      col.vsOrig = ci.v1.idx();
      col.vtOrig = ci.v0.idx();
      idxOfCollapsedNode[numCollapses] = ci.v0.idx();
      
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for( int j = 0; j < _numOfMeshes; j++ ){
        typename TriMesh::Point vs ( _opMeshes[j].point(ci.v1) ), vt ( _opMeshes[j].point(ci.v0) );
        // delta_vt = vs_new - vt, delta_vs = vs_new - vs
        // => vs = vs_new + delta_vs, vt = vs_new + delta_vs

	if ( _updatePositions ){
	  // vs_new = (vs + vt) / 2.
	  typename TriMesh::Point vsNew;
	  for( int k = 0; k < 3; k++ ){
	    col.delta_vt[j][k] = 0.5 * (vs[k] - vt[k]);
	    col.delta_vs[j][k] = 0.5 * (vt[k] - vs[k]);
	    vsNew[k] = 0.5 * (vt[k] + vs[k]);
	  }
	  _opMeshes[j].set_point( ci.v1, vsNew );
	}
	else{
	  // vs_new = vs
	  col.delta_vs[j].setZero();
	  for( int k = 0; k < 3; k++ )
	    col.delta_vt[j][k] = vs[k] - vt[k];
	}
      }
      decInfo.pushBack( col );
      
      // perform collapse
      for( int i = 0; i < _numOfMeshes; i++ )
        _opMeshes[i].collapse( v0v1 );                  
      ++numCollapses;      
      
      // post-process collapse
      _quadric.postprocessCollapse( ci, support );

      // update heap (former one ring of decimated vertex)
      for (s_it = support.begin(), s_end = support.end(); s_it != s_end; ++s_it) {
        assert( !refMesh.status(*s_it).deleted() );
        insertVertexIntoHeap( refMesh, *s_it);
      }
    }

    // build vertex map:
    // After k < n decimations decInfo.vertexMap[i] is  a) the index of vertex i in the original mesh for i < n-k
    //                                                  b) the index of the vertex vt, that was removed in the jth decimation, for i = n-k+j, j=0,...,k-1
    int i0( 0 ), i1( _numOfOriginalVertices - 1 );    
    while (1){
      // find 1st deleted and last un-deleted
      while (!refMesh.status(VertexHandle(i0)).deleted() && i0 < i1) ++i0;     
      while ( refMesh.status(VertexHandle(i1)).deleted() && i0 < i1) --i1;      
      if (i0 >= i1) break;      
      std::swap( decInfo.vertexMap[i0++],  decInfo.vertexMap[i1--] );
    }    
    
    for( int i = 0; i < numCollapses; i++ )
      decInfo.vertexMap[_numOfOriginalVertices - numCollapses + i ] = idxOfCollapsedNode[ numCollapses - 1 - i ];
  
    // delete heap
    _heap.reset();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( int i = 0; i < _numOfMeshes; i++ )
      _opMeshes[i].garbage_collection();
 
  }


  /// Insert vertex in heap
  //NOTE Copied and edited from DecimaterT<Mesh>::heap_vertex() in OpenMesh/Tools/Decimater/DecimaterT.hh
  void insertVertexIntoHeap( TriMesh& refMesh, VertexHandle _vh ) {   

    float prio, best_prio(FLT_MAX);
    typename TriMesh::HalfedgeHandle heh, collapse_target;

    // if vertex is supposed to survive the decimation
    if( _fixedVertices[_vh.idx()] ){
      _collapseTarget[_vh.idx()] = collapse_target;
      _prio[_vh.idx()] = best_prio;
      return;
    }
    
    // find best target in one ring
    typename TriMesh::VertexOHalfedgeIter voh_it(refMesh, _vh);
    for (; voh_it.is_valid(); ++voh_it) {
      heh = *voh_it;
      CollapseInfo ci( refMesh, heh );

      if ( isCollapseLegal( refMesh, ci) ) {
        prio = _quadric.collapsePriority( ci );
        if (prio >= 0.0 && prio < best_prio) {
          best_prio = prio;
          collapse_target = heh;
        }
      }
    }

    // target found -> put vertex on heap
    if (collapse_target.is_valid()) {

      _collapseTarget[_vh.idx()] = collapse_target;
      _prio[_vh.idx()] = best_prio;

      if (_heap->is_stored(_vh))
        _heap->update(_vh);
      else
        _heap->insert(_vh);
    }
    // not valid -> remove from heap
    else {
      if (_heap->is_stored(_vh))
        _heap->remove(_vh);

      _collapseTarget[_vh.idx()] = collapse_target;
      _prio[_vh.idx()] = -1.;
    }
  }
  
  /**
  *   Local configuration:
  *
  *       v1
  *        *
  *       / \
  *      /   \
  *     /     \
  * v0 *-------* prev
  *     \     /
  *      \   / edge 
  *       \ /
  *        *
  *        curr
  *
  **/
  bool isSameOrientation( const TriMesh &mesh, HalfedgeHandle edge, HalfedgeHandle curr_to_v0, VertexHandle curr, VertexHandle v1 ) const {
      
    Point edge_vec = mesh.calc_edge_vector(edge);
    Point curr_to_v0_vec = mesh.calc_edge_vector(curr_to_v0);  
    
    // is current triangle degenerated?
    Point cross_curr = cross(edge_vec, curr_to_v0_vec);
    if ( cross_curr.norm() < 1e-7 )
      return false;
    
    // is edge parallel to (curr - v1)?
    Point curr_to_v1_vec = mesh.point(v1) - mesh.point(curr);
    Point cross_v1 = cross(edge_vec, curr_to_v1_vec);    
    if ( cross_v1.norm() < 1e-7 )
      return false;

    cross_curr.normalize(); 
    cross_v1.normalize();
    if ( dot(cross_curr, cross_v1) < 1e-3)
      return false;
    
    return true;
  }


  //
  bool isCollapseLegalGeometrically( const TriMesh &mesh, const CollapseInfo& _ci ) const {
    VertexHandle  v0 = _ci.v0, v1 = _ci.v1, prev = v0;
    HalfedgeHandle ring_edge;
    for ( typename TriMesh::ConstVertexVertexIter vv_it = mesh.cvv_iter(v0); vv_it.is_valid(); ++vv_it) {
      if (prev == v0){
          prev = *vv_it;
          continue;
      }
      if (*vv_it != v1 && prev != v1){
				ring_edge = mesh.find_halfedge(prev, *vv_it);
				HalfedgeHandle curr_to_v0 = mesh.find_halfedge(*vv_it, v0);
				if (!isSameOrientation(mesh, ring_edge, curr_to_v0, *vv_it, v1)) 
                                    return false;
      }
      prev = *vv_it;
    }
    
    if (*mesh.cvv_iter(v0) != v1 && prev != v1){
			ring_edge = mesh.find_halfedge(prev, *mesh.cvv_iter(v0));
			HalfedgeHandle curr_to_v0 = mesh.find_halfedge(*mesh.cvv_iter(v0), v0);
			if (!isSameOrientation(mesh, ring_edge, curr_to_v0, *mesh.cvv_iter(v0), v1)) 
                            return false;
    }
    
    //
    return true;
  }
 
  //NOTE Copied and edited from BaseDecimaterT<>::is_collapse_legal() in OpenMesh/Tools/Decimater/BaseDecimaterT.hh
  /**
  *   Local configuration:
  *
  *       vl
  *        *
  *       / \
  *      /   \
  *     / fl  \
  * v0 *------>* v1
  *     \ fr  /
  *      \   /
  *       \ /
  *        *
  *        vr
  *
  **/
  bool isCollapseLegal( TriMesh& mesh, const CollapseInfo& _ci) const {

    // locked ?
    if ( mesh.status(_ci.v0).locked() )
      return false;

    // this test checks:
    // is v0v1 deleted?
    // is v0 deleted?
    // is v1 deleted?
    // are both vlv0 and v1vl boundary edges?
    // are both v0vr and vrv1 boundary edges?
    // are vl and vr equal or both invalid?
    // one ring intersection test
    // edge between two boundary vertices should be a boundary edge
    if ( !mesh.is_collapse_ok(_ci.v0v1) )
      return false;

    if (_ci.vl.is_valid() && _ci.vr.is_valid()
        && mesh.find_halfedge(_ci.vl, _ci.vr).is_valid()
        && mesh.valence(_ci.vl) == 3 && mesh.valence(_ci.vr) == 3) {
      return false;
    }
  
    //--- feature test ---
    if ( mesh.status(_ci.v0).feature() && !mesh.status(mesh.edge_handle(_ci.v0v1)).feature() )
      return false;

    //--- test boundary cases ---
    if (mesh.is_boundary(_ci.v0)) {

      // don't collapse a boundary vertex to an inner one
      if (!mesh.is_boundary(_ci.v1))
        return false;

      // only one one ring intersection
      if (_ci.vl.is_valid() && _ci.vr.is_valid())
        return false;
    }

    // there have to be at least 2 incident faces at v0
    if (mesh.cw_rotated_halfedge_handle(
        mesh.cw_rotated_halfedge_handle(_ci.v0v1)) == _ci.v0v1)
      return false;

    // geometric test for convexity
    if ( !isCollapseLegalGeometrically(mesh, _ci) ) 
        return false;
    
    // collapse passed all tests -> ok
    return true;
  }

};

#endif