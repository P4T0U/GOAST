// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Elastic mean computation.
 * \author Heeren, Sassen
 * 
 */
#ifndef ELASTICMEAN_H
#define ELASTICMEAN_H

#include <vector>
#include <goast/Core.h>
#include <goast/Optimization.h>
#include <goast/DiscreteShells/AuxiliaryFunctions.h>
#include <goast/NRIC/LeastSquaresReconstruction.h>

/**
 * \brief Computes approximation of elastic mean by reconstruction from \f$ L\Theta \f$ mean
 * \author Heeren
 *
 * For a sequence of given meshes \f$ s_1, \ldots, s_N \f$, compute corresponding representations \f$ z_1, \ldots, z_N \f$ in \f$ L \Theta \f$ - space,
 * i.e. by computing the vectors \$f z_k = (l^k, \theta^k)\f$ of edge lengths \f$ (l^k_e)_e \in R^m \f$ and dihedral angles \f$ (\theta^k_e)_e \in R^m \f$,
 * where \f$ m \f$ is the number of edges.
 * Then, perform linear averaging in \f$ L \Theta \f$ - space, i.e. compute \f$ \bar z = \frac1{N} \sum_{k=1}^N z_k \f$,
 * and reconstruct the (approximative) shape mean \f$ \bar s \f$ from \f$ \bar z \f$ by least squares optimization.
 *
 * Method as presented in the paper \cite FrBo11
 */
template<typename ConfiguratorType>
void computeLThetaMean( const typename ConfiguratorType::VectorType& Shapes,
                        const MeshTopologySaver& Topology,
                        const double bendingWeight,
                        const int maxIterations,
                        const std::vector<int> * bdryMask,
                        typename ConfiguratorType::VectorType& initialization,
                        bool quiet = true ) {

    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;

    int numV = Topology.getNumVertices();
    int numDofs = 3*numV;
    bool hasInitialization( initialization.size() == numDofs );
    if( !hasInitialization )
        initialization.resize( numDofs );

    int numOfShapes = Shapes.size() / numDofs;
    VectorType LThetaMean( 2 * Topology.getNumEdges() );
    LThetaMean.setZero();
    for (int k = 0; k < numOfShapes; k++) {
        Eigen::Ref<const VectorType> currShapeRef = Shapes.segment(k * numDofs, numDofs);
        VectorType currLengthsAngles;
        computeLengthsAngles<ConfiguratorType>( Topology, currShapeRef, currLengthsAngles );
        LThetaMean += currLengthsAngles;
    }
    LThetaMean *= 1./numOfShapes;

    // integration weights in reconstruction functionals
    VectorType Weights;
    VectorType refGeom( Shapes.segment(0, numDofs) );
    getReconstructionWeights<ConfiguratorType>( Topology, refGeom, Weights );

    // define reconstruction functionals
    RealType edgeLengthWeight = 1.;
    LinearReconstructionFunctional<ConfiguratorType> L( Topology, LThetaMean, Weights, edgeLengthWeight, bendingWeight );
    LinearReconstructionDerivative<ConfiguratorType> DL( Topology, LThetaMean, Weights, edgeLengthWeight, bendingWeight );

    // define initialization and apply Gauss-Newton
    VectorType init( refGeom );
    GaussNewtonAlgorithm<ConfiguratorType> GNOp( 2 * Topology.getNumEdges(), L, DL, maxIterations );
    if( quiet )
      GNOp.setQuiet();
    if( bdryMask ){
            GNOp.setBoundaryMask( *bdryMask );
            // use original content as initialization (at least at boundary)
            if( hasInitialization )
                init = initialization;
    }
    GNOp.solve( init, initialization );

}

/**
 * \brief Elastic average functional, i.e. weighted sum of deformation energies
 * \author Heeren, Sassen
 *
 * Given \f$ n \f$ example shapes \f$ s_1,\dots,s_n \f$, weights  \f$ \alpha_1,\ldots,\alpha_n \f$,
 * this class implements the elastic average energy \f$ F[s] = \sum_{i=1}^n \alpha_i W[s_i, s] \f$.
 */
template<typename ConfiguratorType>
class ElasticMeanFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes, _alpha;
    const unsigned int _numOfShapes;
    std::vector<Eigen::Ref<const VectorType> > _argRefs;

public:
    ElasticMeanFunctional(const DeformationBase<ConfiguratorType> &W, const VectorType &shapes,
                          const VectorType alpha, const unsigned int numOfShapes) : _W(W), _shapes(shapes),
                                                                                    _alpha(alpha),
                                                                                    _numOfShapes(numOfShapes) {
        if (shapes.size() % numOfShapes != 0)
            throw BasicException("ElasticMeanFunctional: wrong number of dof!");
        if (alpha.size() % numOfShapes != 0)
            throw BasicException("ElasticMeanFunctional: wrong number of alphas!");

        // Create references for this different shapes bec. of convenience
        _argRefs.reserve(numOfShapes);

        const int numLocalDofs = _shapes.size() / numOfShapes;
        for (int k = 0; k < numOfShapes; k++)
            _argRefs.push_back(_shapes.segment(k * numLocalDofs, numLocalDofs));
    }

    void apply(const VectorType &Arg, RealType &Dest) const override {
        if (Arg.size() != _argRefs[0].size())
            throw BasicException("ElasticMeanFunctional::apply: wrong size of dofs!");

        Dest = 0.;
        RealType v;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int k = 0; k < _numOfShapes; k++) {
          RealType Energy;
            _W.applyEnergy(_argRefs[k], Arg, Energy);
#pragma omp critical
            {
                Dest += _alpha[k] * Energy;
            }
        }
    }
};

//! \brief Gradient of elastic average functional
//! \author Heeren, Sassen
template<typename ConfiguratorType>
class ElasticMeanFunctionalGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes, _alpha;
    const unsigned int _numOfShapes;
    std::vector<Eigen::Ref<const VectorType> > _argRefs;

public:
    ElasticMeanFunctionalGradient(const DeformationBase<ConfiguratorType> &W, const VectorType &shapes,
                                  const VectorType &alpha, const unsigned int numOfShapes) : _W(W), _shapes(shapes),
                                                                                             _alpha(alpha),
                                                                                             _numOfShapes(numOfShapes) {
        if (shapes.size() % numOfShapes != 0)
            throw BasicException("ElasticMeanFunctional: wrong number of dof!");
        if (alpha.size() % numOfShapes != 0)
            throw BasicException("ElasticMeanFunctional: wrong number of alphas!");

        // Create references for this different shapes bec. of convenience
        _argRefs.reserve(numOfShapes);

        const int numLocalDofs = _shapes.size() / _numOfShapes;
        for (int k = 0; k < numOfShapes; k++)
            _argRefs.push_back(_shapes.segment(k * numLocalDofs, numLocalDofs));
    }

    void apply(const VectorType &Arg, VectorType &Dest) const override {
        if (Arg.size() != _argRefs[0].size())
            throw BasicException("ElasticMeanFunctional::apply: wrong size of dofs!");
        if (Dest.size() != Arg.size())
            Dest.resize(Arg.size());

        Dest.setZero();
        VectorType Gradient;

#ifdef _OPENMP
#pragma omp parallel for private(Gradient)
#endif
        for (int k = 0; k < _numOfShapes; k++) {
            _W.applyDefGradient(_argRefs[k], Arg, Gradient);
#ifdef _OPENMP
#pragma omp critical
#endif
            Dest += _alpha[k] * Gradient;
        }
    }
};

//! \brief Hessian of elastic average functional
//! \author Heeren, Sassen
template<typename ConfiguratorType>
class ElasticMeanFunctionalHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef std::vector<TripletType> TripletListType;

    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes, _alpha;
    const unsigned int _numOfShapes;
    std::vector<Eigen::Ref<const VectorType> > _argRefs;


public:
    ElasticMeanFunctionalHessian(const DeformationBase<ConfiguratorType> &W, const VectorType &shapes,
                                 const VectorType alpha, const unsigned int numOfShapes) : _W(W), _shapes(shapes),
                                                                                           _alpha(alpha),
                                                                                           _numOfShapes(numOfShapes) {
        if (shapes.size() % numOfShapes != 0)
            throw BasicException("ElasticMeanFunctional: wrong number of dof!");
        if (alpha.size() % numOfShapes != 0)
            throw BasicException("ElasticMeanFunctional: wrong number of alphas!");

        // Create references for this different shapes bec. of convenience
        _argRefs.reserve(numOfShapes);

        const int numLocalDofs = _shapes.size() / _numOfShapes;
        for (int k = 0; k < numOfShapes; k++)
            _argRefs.push_back(_shapes.segment(k * numLocalDofs, numLocalDofs));
    }

    void apply(const VectorType &Arg, MatrixType &Dest) const override {
        if (Arg.size() != _argRefs[0].size())
            throw BasicException("ElasticMeanFunctional::apply: wrong size of dofs!");

        if( (Dest.cols() != Arg.size()) || (Dest.rows() != Arg.size()) )
            Dest.resize( Arg.size(), Arg.size() );

        Dest.setZero();

        // fill triplet lists
        TripletListType tripletList;
        // we have 3*numOfShapes deformations
        tripletList.reserve(3 * _numOfShapes * _W.numOfNonZeroHessianEntries());

        pushTriplets(Arg, tripletList);

        Dest.setFromTriplets(tripletList.cbegin(), tripletList.cend());
    }

    void pushTriplets( const VectorType& Arg, TripletListType& tripletList ) const {
        if (Arg.size() != _argRefs[0].size())
            throw BasicException("ElasticMeanFunctional::pushTriplets: wrong size of dofs!");

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int k = 0; k < _numOfShapes; k++) {
            TripletListType localTripletList;
            localTripletList.reserve(_W.numOfNonZeroHessianEntries());
            _W.pushTripletsDefHessian(_argRefs[k], Arg, localTripletList, 0, 0, _alpha[k]);
#ifdef _OPENMP
#pragma omp critical
#endif
            tripletList.insert(tripletList.end(), localTripletList.begin(), localTripletList.end());

        }
    }

};

//! \brief Computation of (weighted) elastic average by optimization of ElasticMeanFunctional
//! \author Heeren, Sassen
template<typename ConfiguratorType>
class ElasticMean {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const MeshTopologySaver &_topology;
    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes, _alpha;
    const unsigned int _numOfShapes;

    mutable VectorType _solution;

    const OptimizationParameters<ConfiguratorType>& _optPars;
    std::vector<int> const* _mask;
    bool _quiet;

public:
    ElasticMean(const MeshTopologySaver &Topology,
                const DeformationBase<ConfiguratorType> &W,
                const VectorType &shapes,
                const VectorType &alpha,
                const unsigned int numOfShapes,
                const OptimizationParameters<ConfiguratorType>& optPars,
                bool quiet = true )
            : _W(W), _shapes(shapes), _alpha(alpha), _numOfShapes(numOfShapes),  _optPars(optPars),
              _topology(Topology), _mask(NULL), _quiet(quiet)
    {

    }

    void setBoundaryMask(const std::vector<int> &localMask) {
        _mask = &localMask;
    }

    void execute(VectorType &mean) const {
        ElasticMeanFunctional<ConfiguratorType> E(_W, _shapes, _alpha, _numOfShapes);
        ElasticMeanFunctionalGradient<ConfiguratorType> DE(_W, _shapes, _alpha, _numOfShapes);
        ElasticMeanFunctionalHessian<ConfiguratorType> D2E(_W, _shapes, _alpha, _numOfShapes);

	if( mean.size() != 3 * _topology.getNumVertices() ){
          if(!_quiet) std::cerr << "ElasticMean::execute() initialization has wrong size -> resize!" << std::endl;
	  mean = _shapes.segment( 0, 3 *_topology.getNumVertices() );
	}

        VectorType initialization(mean);

      RealType energy;
        VectorType grad;
	if(!_quiet){
          E.apply(initialization, energy);
          std::cerr << "Initial functional  = " << energy << std::endl;
          DE.apply(initialization, grad);
          std::cerr << "Initial gradient norm = " << grad.norm() << std::endl << std::endl;
	}

        // Optimization with gradient descent
        if (_optPars.getGradientIterations() > 0) {
            if(!_quiet) std::cerr << "Start gradient descent... " << std::endl;
            GradientDescent<ConfiguratorType> GD(E, DE, _optPars);
            if (_mask) {
                if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                GD.setBoundaryMask(*_mask);
            }
            GD.solve(initialization, mean);
            initialization = mean;
        }

        // optimization with BFGS
        if (_optPars.getBFGSIterations() > 0) {
            if(!_quiet) std::cerr << "Start Quasi-Newton... " << std::endl;
            QuasiNewtonBFGS<ConfiguratorType> QNM(E, DE, _optPars);
            if (_mask) {
                if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                QNM.setBoundaryMask(*_mask);
            }
            QNM.solve(initialization, mean);
            initialization = mean;
        }

        if (_optPars.getNewtonIterations() > 0) {
            if (_mask) {
                if(!_quiet) std::cerr << "Start Newton with boundary mask... " << std::endl;
//                NewtonMethod<ConfiguratorType> NM(DE, D2E, _optPars);
                LineSearchNewton<ConfiguratorType> NM (E, DE, D2E, _optPars );
                if (_mask) {
                    if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                    NM.setBoundaryMask(*_mask);
                }
                NM.solve(initialization, mean);
            } else {
                if(!_quiet) std::cerr << "Start Newton with Lagrange... " << std::endl;
                int numConstraints = 6;
                int totalNumDofs = 3 * _topology.getNumVertices();

                // Somehow using _argRefs here leads to a segfault
                VectorType startGeom = _shapes.block(0,0,totalNumDofs,1);

                RigidBodyMotionsConstraintHandler<ConfiguratorType> constHandler(_topology, startGeom, 1);


                LagrangeFunctional<ConfiguratorType,
                ElasticMeanFunctional<ConfiguratorType>,
                RigidBodyMotionsConstraintHandler<ConfiguratorType> > L(E, constHandler, totalNumDofs);

                LagrangeGradient<ConfiguratorType, ElasticMeanFunctionalGradient<ConfiguratorType>,
                RigidBodyMotionsConstraintHandler<ConfiguratorType> > DL(DE, constHandler, totalNumDofs);

                LagrangeHessian<ConfiguratorType, ElasticMeanFunctionalHessian<ConfiguratorType>,
                RigidBodyMotionsConstraintHandler<ConfiguratorType> > D2L(D2E, constHandler, totalNumDofs);

                initialization.conservativeResize(totalNumDofs + numConstraints);
                mean.conservativeResize(totalNumDofs + numConstraints);

                NewtonMethod<ConfiguratorType> NM(DL, D2L, _optPars);
                NM.solve(initialization, mean);
                mean.conservativeResize(totalNumDofs);
            }
        }

        if(!_quiet){
          E.apply(mean, energy);
          std::cerr << "Final energy = " << energy << std::endl;
          DE.apply(mean, grad);
          std::cerr << "Final gradient norm = " << grad.norm() << std::endl << std::endl;
	}
    }

};

//! \brief Wrapper function to compute elastic average
//! \author Heere, Sassen
template<typename ConfiguratorType>
void computeElasticAverage(const MeshTopologySaver &Topology,
                           const DeformationBase<ConfiguratorType> &W,
                           const typename ConfiguratorType::VectorType &shapes,
                           const typename ConfiguratorType::VectorType &alpha,
                           const unsigned int numOfShapes,
                           const OptimizationParameters<ConfiguratorType> &optPars,
                           std::vector<int> const *localMask,
                           typename ConfiguratorType::VectorType &average) {
    ElasticMean<ConfiguratorType> elastMeanOp(Topology, W, shapes, alpha, numOfShapes, optPars);
    if (localMask)
        elastMeanOp.setBoundaryMask(*localMask);
    elastMeanOp.execute(average);
}


/**
 * \brief Fast version of elastic average functional
 * \author Heeren
 *
 * Given \f$ n \f$ example shapes \f$ s_1,\dots,s_n \f$, weights  \f$ \alpha_1,\ldots,\alpha_n \f$
 * this class implements the elastic average energy \f$ F[s] = \sum_{i=1}^n \alpha_i W[s_i, s] \f$.
 *
 * Different from class ElasticMeanFunctional above, here we pre-compute and store all geometric quantities on the fixed (!) input shapes
 * that are needed in every energy evaluation
 */
template<typename ConfiguratorType>
class FastElasticMeanFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const MeshTopologySaver& _topology;
    const VectorType &_shapes;
    mutable VectorType _alpha;
    RealType _bendWeight, _mu, _lambdaQuarter, _muHalfPlusLambdaQuarter;
    const unsigned int _numOfShapes;
    std::vector<VectorType> _refDihedralAngles, _refFaceAreas, _refSqrEdgeLengths;

public:
    FastElasticMeanFunctional( const MeshTopologySaver& Topology,
                               const VectorType& shapes,
                               const VectorType& alpha,
                               RealType BendingWeight,
                               const unsigned int numOfShapes) : _topology(Topology), _shapes(shapes),
                                                                                    _alpha(alpha),
                                                                                    _bendWeight(BendingWeight),
                                                                                    _mu(1.),
                                                                                    _lambdaQuarter(0.25),
                                                                                    _muHalfPlusLambdaQuarter(0.5*_mu + _lambdaQuarter),
                                                                                    _numOfShapes(numOfShapes),
                                                                                    _refDihedralAngles(numOfShapes), _refFaceAreas(numOfShapes), _refSqrEdgeLengths(numOfShapes){
        if (shapes.size() % numOfShapes != 0)
            throw BasicException("FastElasticMeanFunctional: wrong number of dof!");
        if (alpha.size() != numOfShapes )
            throw BasicException("FastElasticMeanFunctional: wrong number of alphas!");

        // Create references for this different shapes bec. of convenience
        std::vector<Eigen::Ref<const VectorType> > argRefs;
        argRefs.reserve(numOfShapes);

        const int numLocalDofs = _shapes.size() / numOfShapes;
        for (int k = 0; k < numOfShapes; k++)
            argRefs.push_back(_shapes.segment(k * numLocalDofs, numLocalDofs));

        // precompute reference quantities
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for( int k = 0; k < numOfShapes; k++ ){
          VectorType dummy;
          computeReferenceQuantitiesBending<ConfiguratorType>( _topology, argRefs[k], _refDihedralAngles[k], dummy, _refSqrEdgeLengths[k]);
          getFaceAreas<ConfiguratorType>( _topology, argRefs[k], _refFaceAreas[k] );
        }
    }

    //
    void setAlpha( const VectorType& Alpha ) const {
        if( Alpha.size() != _alpha.size() )
             throw BasicException("FastElasticMeanFunctional::setAlpha: wrong size!");
        _alpha = Alpha;
    }

    //
    void apply(const VectorType &Argument, RealType &Dest) const override {

      if (Argument.size() != 3 * _topology.getNumVertices() )
          throw BasicException("FastElasticMeanFunctional::apply: wrong size of dofs!");

      // precompute quantities of argument
      VectorType dihedralAngles, sqrFaceAreas;
      getDihedralAngles<ConfiguratorType>( _topology, Argument, dihedralAngles );
      std::vector<VecType> dotProducts;
      getDotProductsAndSquaredFaceAreas<ConfiguratorType>( _topology, Argument, dotProducts, sqrFaceAreas );

      // set up single energies and parallelize
      VectorType singleEnergies( _numOfShapes );
      singleEnergies.setZero();

#ifdef _OPENMP
#pragma omp parallel for
#endif
     for (int k = 0; k < _numOfShapes; k++){
        // bending contribution
        for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

          if( !(_topology.isEdgeValid(edgeIdx)) )
	    continue;

          int pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ), pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

          // no bending at boundary edges
          if( std::min( pl, pk) < 0 )
            continue;

          int fl(_topology.getAdjacentTriangleOfEdge(edgeIdx,0) ), fr(_topology.getAdjacentTriangleOfEdge(edgeIdx,1) );
          RealType delTheta = dihedralAngles[edgeIdx] - _refDihedralAngles[k][edgeIdx];
          singleEnergies[k] += _bendWeight * _alpha[k] * delTheta * delTheta * _refSqrEdgeLengths[k][edgeIdx] / (_refFaceAreas[k][fl] + _refFaceAreas[k][fr]);
        }

        // membrane contribution
        int numFaces = _topology.getNumFaces();
        for ( int faceIdx = 0; faceIdx < numFaces; ++faceIdx ){
          RealType volRef( _refFaceAreas[k][faceIdx] );
          // trace term = -1. *  \sum_{i =0,1,2} <e_{i+1}, e_{i+2}> |\bar e_i|^2, where \bar e_i are reference edges
          // note the signs! This is since we have dotProducts[faceIdx] = { <e_1, -e_2>, <-e_2, e_0>, <e_0, e_1> }
          RealType traceTerm( dotProducts[faceIdx][0] * _refSqrEdgeLengths[k][ _topology.getEdgeOfTriangle(faceIdx,0)] + dotProducts[faceIdx][1] * _refSqrEdgeLengths[k][_topology.getEdgeOfTriangle(faceIdx,1)] - dotProducts[faceIdx][2] * _refSqrEdgeLengths[k][_topology.getEdgeOfTriangle(faceIdx,2)] );
          singleEnergies[k] += _alpha[k] * ( (_mu/8. *  traceTerm + _lambdaQuarter * sqrFaceAreas[faceIdx]) / volRef -  ( _muHalfPlusLambdaQuarter * std::log( sqrFaceAreas[faceIdx] / (volRef*volRef) ) + _mu + _lambdaQuarter) * volRef);
        }
     }

     //
     Dest = singleEnergies.sum();

    }

};

//! \brief Fast gradient of FastElasticMeanFunctional
//! \author Heeren
template<typename ConfiguratorType>
class FastElasticMeanGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const MeshTopologySaver& _topology;
    const VectorType &_shapes;
    mutable VectorType _alpha;
    RealType _bendWeight, _mu, _lambdaQuarter, _muHalfPlusLambdaQuarter;
    const unsigned int _numOfShapes;
    std::vector<VectorType> _refDihedralAngles, _refFaceAreas, _refSqrEdgeLengths, _refFactors;

public:
    FastElasticMeanGradient( const MeshTopologySaver& Topology,
                               const VectorType& shapes,
                               const VectorType& alpha,
                               RealType BendingWeight,
                               const unsigned int numOfShapes) : _topology(Topology), _shapes(shapes),
                                                                                    _alpha(alpha),
                                                                                    _bendWeight(BendingWeight),
                                                                                    _mu(1.),
                                                                                    _lambdaQuarter(0.25),
                                                                                    _muHalfPlusLambdaQuarter(0.5*_mu + _lambdaQuarter),
                                                                                    _numOfShapes(numOfShapes),
                                                                                    _refDihedralAngles(numOfShapes), _refFaceAreas(numOfShapes), _refSqrEdgeLengths(numOfShapes), _refFactors(numOfShapes){
        if (shapes.size() % numOfShapes != 0)
            throw BasicException("FastElasticMeanGradient: wrong number of dof!");
        if (alpha.size() != numOfShapes )
            throw BasicException("FastElasticMeanGradient: wrong number of alphas!");

        // Create references for this different shapes bec. of convenience
        std::vector<Eigen::Ref<const VectorType> > argRefs;
        argRefs.reserve(numOfShapes);

        const int numLocalDofs = _shapes.size() / numOfShapes;
        for (int k = 0; k < numOfShapes; k++)
            argRefs.push_back(_shapes.segment(k * numLocalDofs, numLocalDofs));

        // precompute reference quantities
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for( int k = 0; k < numOfShapes; k++ ){
          VectorType dummy;
          computeReferenceQuantitiesBending<ConfiguratorType>( _topology, argRefs[k], _refDihedralAngles[k], dummy, _refSqrEdgeLengths[k]);
          getFaceAreasAndFactors<ConfiguratorType>( _topology, argRefs[k], _refFaceAreas[k], _refFactors[k] );
        }
    }

    //
    void setAlpha( const VectorType& Alpha ) const {
        if( Alpha.size() != _alpha.size() )
             throw BasicException("FastElasticMeanGradient::setAlpha: wrong size!");
        _alpha = Alpha;
    }

    //
    void apply(const VectorType &Argument, VectorType &Dest) const override {

      if (Argument.size() != 3 * _topology.getNumVertices() )
          throw BasicException("FastElasticMeanGradient::apply: wrong size of dofs!");
      if (Dest.size() != 3 * _topology.getNumVertices() )
          Dest.resize(  3 * _topology.getNumVertices()  );
      Dest.setZero();

      // precompute quantities of argument
      std::vector<VecType> gradPi, gradPj, gradPk, gradPl;
      VectorType dihedralAngles;
      getDihedralAnglesGradients<ConfiguratorType>( _topology, Argument, dihedralAngles, gradPi, gradPj, gradPk, gradPl );

      // single gradients for parallelization
      std::vector<VectorType> gradients( _numOfShapes );
      for (int k = 0; k < _numOfShapes; k++){
        gradients[k].resize( 3 * _topology.getNumVertices()  );
        gradients[k].setZero();
      }

      // bending contribution
#ifdef _OPENMP
#pragma omp parallel for
#endif
     for (int k = 0; k < _numOfShapes; k++){
        for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

          if( !(_topology.isEdgeValid(edgeIdx)) )
	    continue;

	  int pi(_topology.getAdjacentNodeOfEdge(edgeIdx, 0)),
              pj(_topology.getAdjacentNodeOfEdge(edgeIdx, 1)),
	      pk(_topology.getOppositeNodeOfEdge(edgeIdx, 0)),
	      pl(_topology.getOppositeNodeOfEdge(edgeIdx, 1));

          // no bending at boundary edges
          if( std::min( pl, pk) < 0 )
            continue;

          // factor
          int fl(_topology.getAdjacentTriangleOfEdge(edgeIdx,0) ), fr(_topology.getAdjacentTriangleOfEdge(edgeIdx,1) );
          RealType delTheta( -2. *  _refSqrEdgeLengths[k][edgeIdx] * (_refDihedralAngles[k][edgeIdx] - dihedralAngles[edgeIdx] ) / (_refFaceAreas[k][fl] + _refFaceAreas[k][fr]) );

          // assemble in global vector
	  for (int i = 0; i < 3; i++){
            gradients[k][i*_topology.getNumVertices() + pi] += _bendWeight * _alpha[k] * delTheta * gradPi[edgeIdx][i];
            gradients[k][i*_topology.getNumVertices() + pj] += _bendWeight * _alpha[k] * delTheta * gradPj[edgeIdx][i];
            gradients[k][i*_topology.getNumVertices() + pk] += _bendWeight * _alpha[k] * delTheta * gradPk[edgeIdx][i];
            gradients[k][i*_topology.getNumVertices() + pl] += _bendWeight * _alpha[k] * delTheta * gradPl[edgeIdx][i];
	  }
        }
     }

     // get face areas and area gradients
     VectorType faceAreas;
     getAreaGradients<ConfiguratorType>( _topology, Argument, faceAreas, gradPi, gradPj, gradPk, gradPl );

     // membrane contribution I
#ifdef _OPENMP
#pragma omp parallel for
#endif
     for (int k = 0; k < _numOfShapes; k++){
        for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx ){
            std::vector<int> nodesIdx(3);
            for (int j = 0; j < 3; j++)
              nodesIdx[j] = _topology.getNodeOfTriangle(faceIdx, j);

            RealType factor = 2. * (_lambdaQuarter * faceAreas[faceIdx] / _refFaceAreas[k][faceIdx] - _muHalfPlusLambdaQuarter * _refFaceAreas[k][faceIdx] / faceAreas[faceIdx]);
            for (int j = 0; j < 3; j++){
              gradients[k][j*_topology.getNumVertices() + nodesIdx[0]] += _alpha[k] * factor * gradPi[faceIdx][j];
              gradients[k][j*_topology.getNumVertices() + nodesIdx[1]] += _alpha[k] * factor * gradPj[faceIdx][j];
              gradients[k][j*_topology.getNumVertices() + nodesIdx[2]] += _alpha[k] * factor * gradPk[faceIdx][j];
            }
        }
     }

     // membrane contribution II
     //! TODO
#ifdef _OPENMP
#pragma omp parallel for
#endif
     for (int k = 0; k < _numOfShapes; k++){
        RealType scalFactor = -0.25 * _mu * _alpha[k];
        for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx ){
            for (int j = 0; j < 3; j++)
                 for (int i = 0; i < 3; i++){
                     int nodeIdx = _topology.getNodeOfTriangle(faceIdx, i);
                     int nextIdx = (i + 1) % 3;
                     int prevIdx = (i + 2) % 3;
                     gradients[k][j*_topology.getNumVertices() + nodeIdx] += scalFactor * _refFactors[k][3*faceIdx + nextIdx] * gradPl[3*faceIdx + nextIdx][j];
                     gradients[k][j*_topology.getNumVertices() + nodeIdx] -= scalFactor * _refFactors[k][3*faceIdx + prevIdx] * gradPl[3*faceIdx + prevIdx][j];
                 }
        }
     }

     // sum all single gradients
     for (int k = 0; k < _numOfShapes; k++)
         Dest += gradients[k];

    }

};

/**
 * \brief Computation of fast (weighted) elastic average
 * \author Heeren
 * \note Since there is no fast version of the Hessian, we make use of the standard Hessian implemented in ElasticMeanFunctionalHessian here.
 */
template<typename ConfiguratorType>
class FastElasticMeanOperator {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const MeshTopologySaver &_topology;
    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes;
    mutable VectorType _alpha;
    const unsigned int _numOfShapes;

    FastElasticMeanFunctional<ConfiguratorType> _fastE;
    FastElasticMeanGradient<ConfiguratorType> _fastDE;

    mutable VectorType _solution;

    const OptimizationParameters<ConfiguratorType>& _optPars;
    std::vector<int> const* _mask;
    bool _quiet;

public:
    FastElasticMeanOperator( const MeshTopologySaver &Topology,
                             const DeformationBase<ConfiguratorType> &W,
                             const VectorType &shapes,
                             const unsigned int numOfShapes,
                             const OptimizationParameters<ConfiguratorType>& optPars,
                             bool quiet = true )
            : _W(W),
            _shapes(shapes),
            _alpha(VectorType::Constant(numOfShapes, 1.)),
            _numOfShapes(numOfShapes),
            _fastE( Topology, shapes, _alpha, W.getBendingWeight(), numOfShapes ),
            _fastDE( Topology, shapes, _alpha, W.getBendingWeight(), numOfShapes ),
            _optPars(optPars),
            _topology(Topology),
            _mask(NULL),
            _quiet(quiet){ }

    //
    void setAlpha( const VectorType& Alpha ) const {
        if( Alpha.size() != _alpha.size() )
             throw BasicException("FastElasticMeanOperator::setAlpha: wrong size!");
        _alpha = Alpha;
        _fastE.setAlpha( _alpha );
        _fastDE.setAlpha( _alpha );
    }

    void setBoundaryMask(const std::vector<int> &localMask) {
        _mask = &localMask;
    }

    void execute(VectorType &mean) const {
        ElasticMeanFunctionalHessian<ConfiguratorType> D2E(_W, _shapes, _alpha, _numOfShapes);

	if( mean.size() != 3 * _topology.getNumVertices() ){
          if(!_quiet) std::cerr << "FastElasticMeanOperator::execute() initialization has wrong size -> resize!" << std::endl;
	  mean = _shapes.segment( 0, 3 *_topology.getNumVertices() );
	}

        VectorType initialization(mean);

      RealType energy;
        VectorType grad;
	if(!_quiet){
          _fastE.apply(initialization, energy);
          std::cerr << "Initial functional  = " << energy << std::endl;
          _fastDE.apply(initialization, grad);
          std::cerr << "Initial gradient norm = " << grad.norm() << std::endl << std::endl;
	}

        // Optimization with gradient descent
        if (_optPars.getGradientIterations() > 0) {
            if(!_quiet) std::cerr << "Start gradient descent... " << std::endl;
            GradientDescent<ConfiguratorType> GD(_fastE, _fastDE, _optPars);
            if (_mask) {
                if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                GD.setBoundaryMask(*_mask);
            }
            GD.solve(initialization, mean);
            initialization = mean;
        }

        // optimization with BFGS
        if (_optPars.getBFGSIterations() > 0) {
            if(!_quiet) std::cerr << "Start Quasi-Newton... " << std::endl;
            QuasiNewtonBFGS<ConfiguratorType> QNM(_fastE, _fastDE, _optPars);
            if (_mask) {
                if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                QNM.setBoundaryMask(*_mask);
            }
            QNM.solve(initialization, mean);
            initialization = mean;
        }

        if (_optPars.getNewtonIterations() > 0) {
            if (_mask) {
                if(!_quiet) std::cerr << "Start Newton with boundary mask... " << std::endl;
                NewtonMethod<ConfiguratorType> NM(_fastDE, D2E, _optPars);
                if (_mask) {
                    if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                    NM.setBoundaryMask(*_mask);
                }
                NM.solve(initialization, mean);
            } else {
                if(!_quiet) std::cerr << "Start Newton with Lagrange... " << std::endl;
                int numConstraints = 6;
                int totalNumDofs = 3 * _topology.getNumVertices();

                // Somehow using _argRefs here leads to a segfault
                VectorType startGeom = _shapes.block(0,0,totalNumDofs,1);

                RigidBodyMotionsConstraintHandler<ConfiguratorType> constHandler(_topology, startGeom, 1);


                LagrangeFunctional<ConfiguratorType, FastElasticMeanFunctional<ConfiguratorType>,
                RigidBodyMotionsConstraintHandler<ConfiguratorType> > L(_fastE, constHandler, totalNumDofs);

                LagrangeGradient<ConfiguratorType, FastElasticMeanGradient<ConfiguratorType>,
                RigidBodyMotionsConstraintHandler<ConfiguratorType> > DL(_fastDE, constHandler, totalNumDofs);

                LagrangeHessian<ConfiguratorType, ElasticMeanFunctionalHessian<ConfiguratorType>,
                RigidBodyMotionsConstraintHandler<ConfiguratorType> > D2L(D2E, constHandler, totalNumDofs);

                initialization.conservativeResize(totalNumDofs + numConstraints);
                mean.conservativeResize(totalNumDofs + numConstraints);

                NewtonMethod<ConfiguratorType> NM(DL, D2L, _optPars);
                NM.solve(initialization, mean);
                mean.conservativeResize(totalNumDofs);
            }
        }

        if(!_quiet){
          _fastE.apply(mean, energy);
          std::cerr << "Final energy = " << energy << std::endl;
          _fastDE.apply(mean, grad);
          std::cerr << "Final gradient norm = " << grad.norm() << std::endl << std::endl;
	}
    }

};
#endif //ELASTICMEAN_H
