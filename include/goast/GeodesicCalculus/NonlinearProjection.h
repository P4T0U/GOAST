// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Nonlinear projection on submanifold spanned by principal geodesic variations.
 * \author Heeren
 * Defined in Heeren et al. \cite HeZhRu18
 * 
 */
#ifndef NONLINEARPROJECTION_H
#define NONLINEARPROJECTION_H

#include <goast/Core.h>
#include <goast/Optimization.h>

#include "DiscreteGeodesicCalculus.h"
#include "ElasticMean.h"

//===============================================================================================================================

/**
 * \brief Helper funciton to check whether submanifold coordinates are valid
 * \author Heeren
 *
 * An arbitrary shape can be represented in the submanifold spanned by \f$ J \f$ principal variations
 * by mean of \f$ 2J+1 \f$ coefficients.
 *
 * Here, we check whether a general vector \f$ \alpha \in \R^{2J} \f$ can be extended to a valid coefficient vector,
 * which is true if its entries are non-negative and sum to a value less than 1.
 *
 * See Heeren et al. 2018, \cite HeZhRu18
 */
template<typename ConfiguratorType>
bool isAlphaExtensionValid(const typename ConfiguratorType::VectorType &Alpha,
                           typename ConfiguratorType::VectorType &extendedAlpha) {
    int twoJ = Alpha.size();
    extendedAlpha.resize(twoJ + 1);
    bool isValid = true;
    
    extendedAlpha[0] = 1.;
    for (int j = 0; j < twoJ; j++) {
        if (Alpha[j] < 0.)
           isValid = false;
        extendedAlpha[0] -= Alpha[j];
        extendedAlpha[j + 1] = Alpha[j];
    }
    return !(extendedAlpha[0] < 0.) && isValid;
}

//! \brief Penalty function
template<typename RealType>
RealType penaltyFunction( const RealType& t, const RealType& eps ) {
    return (t > 0) ? 0. : t*t / eps;
}

//! \brief Derivative of penalty function
template<typename RealType>
RealType penaltyDerivative( const RealType& t, const RealType& eps ) {
    return (t > 0) ? 0. : 2.* t / eps;
}

//===============================================================================================================================

/**
 * \brief Nonlinear local projection energy
 * \author Heeren
 *
 * Given a geodesic mean \f$ p_0 = \bar s \f$, as well as \f$ J \f$ principal variations \f$ \{p_j\}_{j=1,...,J}\f$
 * and corresponding geodesic reflections \f$ \{p_{-j}\}_{j = 1,...,J}\f$  about \f$ \bar s \f$.
 *
 * For some unseen shape \f$ s \f$ (which is supposed to be close to the mean),
 * this functional represents the nonlinear local projection energy \f$ \mathcal{J}: \alpha \mapsto W[s, q[\alpha]]\f$
 * where \f$ q[\alpha] = \arg\min_q \sum_j \alpha_j W[p_j, q]\f$  is a local elastic average.
 *
 * \note Compare to eq. (18) in \cite HeZhRu18
 */
template<typename ConfiguratorType>
class NonlinearProjectionFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const MeshTopologySaver &_topology;
    const DeformationBase<ConfiguratorType> &_W;
    const OptimizationParameters<ConfiguratorType> &_optPars;
    const VectorType &_shape;
    const VectorType &_princVariations;
    const int _numDOFs;    
    FastElasticMeanOperator<ConfiguratorType> _fastElasticMeanOp;     
    mutable VectorType _q;
    //!TODO take care of mask!!!!
    std::vector<int> const *_mask;
    RealType _penaltyWeight;
    bool _quiet;
    mutable int _counter;

public:
    //! NOTE ordering convention in principalVariations: \f$ (p_j, p_0, p_{j+J}) \f$ is geodesic for \f$ j = 1,...,J \f$
    //! In particular, p[0] is the geodesic mean.
    NonlinearProjectionFunctional(const MeshTopologySaver &Topology,
                                  const DeformationBase<ConfiguratorType> &W,
                                  const VectorType &principalVariations,
                                  const VectorType &shape,
                                  const OptimizationParameters<ConfiguratorType> &optPars,
                                  RealType PenaltyWeight = 0 )
            : _topology(Topology),
              _W(W),
              _shape(shape),
              _optPars(optPars),
              _princVariations(principalVariations),
              _numDOFs(_princVariations.size() / (3 * _topology.getNumVertices()) - 1),
              _fastElasticMeanOp( _topology, _W, _princVariations, _numDOFs + 1, _optPars ),
              _mask(NULL),
              _penaltyWeight(PenaltyWeight), 
              _quiet(true),
              _counter(0){ getGeometry(_topology.getGrid(), _q); }

    
    ~NonlinearProjectionFunctional(){
        std::cerr << "Energy evaluations: " << _counter << std::endl;
    }
    
    void setBoundaryMask(const std::vector<int> &localMask) {
        _mask = &localMask;
        _fastElasticMeanOp.setBoundaryMask( *_mask );
    }
    
    const VectorType& getShape() const {
        return _q;
    }

    //
    void apply(const VectorType &Alpha, RealType &Dest) const override {
        _counter++;
        if (!_quiet) std::cerr << "Eval local proj. functional." << std::endl;
        auto t_start = std::chrono::high_resolution_clock::now();
        
        if (Alpha.size() != _numDOFs)
            throw BasicException("NonlinearProjectionFunctional::apply: wrong size of dofs!");

        VectorType extendedAlpha;
        bool validExtension = isAlphaExtensionValid<ConfiguratorType>(Alpha, extendedAlpha);
        
        // compute the local elastic average $q[\alpha] = \argmin_q \sum_j \alpha_j W[p_j, q]$
        if (!_quiet) std::cerr << "Solve elastic average for q." << std::endl;
        _fastElasticMeanOp.setAlpha( extendedAlpha );
        _fastElasticMeanOp.execute( _q );
 
        // compute W[s,q]
        _W.applyEnergy(_shape, _q, Dest);
       
        // penalty  
        if( _penaltyWeight > 0.){
          for( int i = 0; i < _numDOFs; i++ )
            Dest += penaltyFunction( Alpha[i], 1./_penaltyWeight );
          Dest += penaltyFunction( extendedAlpha[0], 1./_penaltyWeight );
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        if(!_quiet) std::cout << std::fixed << "Energy computed in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << "seconds." << std::endl << std::endl;
    }
};

/**
 * \brief Gradient of NonlinearProjectionFunctional
 * \author Heeren
 *
 * Let \f$ A[\alpha, q] = \sum_j \alpha_j W[p_j, q] \f$, where \f$ {p_j}_j \f$ are the \f$ 2J+1 \f$ principal variations.
 *
 * Let \f$ G[\alpha, q] = d_2 A[\alpha, q] \f$ and \f$ L[q,\alpha,\mu] = W[s,q] + G[\alpha,q] \cdot \mu \f$.
 *
 * Then the necessary condition \f$ 0 = DL \f$ can be written in three sets of equations:
 *
 * (1) \f$ 0 = D_q L = d_2 W[s,q] + d_2^2 A[\alpha,q] \cdot \mu \f$
 *
 * (2) \f$ 0 = D_{\alpha_j} L = (d_2 W[p_j, q] - d_2 W[p_0, q]) \cdot \mu \f$, for \f$ j = 1, ..., 2J \f$
 *
 * (3) \f$ 0 = D_\mu L = G[\alpha, q] \f$
 *
 * The evaluation of the gradient of J (ie the nonlinear local projection energy) for some coefficient vector \f$ \alpha \f$ can then be computed as follows:
 *
 * a) compute \f$ q \f$ by solving the nonlinear system (3) (using the argument \f$ \alpha \f$)
 *
 * b) compute \f$ \mu \f$ by solving the linear system (1)  (using the argument \f$ \alpha \f$ and  \f$ q \f$ from a))
 *
 * c) \f$ D_{\alpha_j} J = D_{\alpha_j} L \f$ given by (2)  (using the argument \f$ \alpha, q \f$ from a) and \f$ \mu \f$ from b))
 */
template<typename ConfiguratorType>
class NonlinearProjectionGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef std::vector<TripletType> TripletListType;

    const MeshTopologySaver &_topology;
    const DeformationBase<ConfiguratorType> &_W;
    const OptimizationParameters<ConfiguratorType> &_optPars;
    const VectorType &_shape;
    const VectorType &_princVariations;
    const int _numDOFs;
    FastElasticMeanOperator<ConfiguratorType> _fastElasticMeanOp;  
    mutable VectorType _q;
    std::vector<int> const *_mask;
    RealType _penaltyWeight;
    mutable TripletListType _constraintTriplets;
    bool _quiet;
    mutable int _counter;

public:
    //! NOTE ordering convention in principalVariations: \f$ (p_j, p_0, p_{j+J}) \f$ is geodesic for  \f$ j = 1,...,J \f$
    //! In particular, p[0] is the geodesic mean.
    NonlinearProjectionGradient(const MeshTopologySaver &Topology,
                                const DeformationBase<ConfiguratorType> &W,
                                const VectorType &principalVariations,
                                const VectorType &shape,
                                const OptimizationParameters<ConfiguratorType> &optPars,
                                RealType PenaltyWeight = 0.,
                                bool quiet = true )
            : _topology(Topology),
              _W(W),
              _optPars(optPars),
              _shape(shape),
              _princVariations(principalVariations),
              _numDOFs(_princVariations.size() / (3 * _topology.getNumVertices()) - 1),
              _fastElasticMeanOp( _topology, _W, _princVariations, _numDOFs + 1, _optPars ),
              _mask(NULL),
              _penaltyWeight(PenaltyWeight), 
              _quiet(quiet),
              _counter(0){
        // initialize q variable
        getGeometry(_topology.getGrid(), _q);
    }
    
    ~NonlinearProjectionGradient(){
        std::cerr << "Grad evaluations: " << _counter << std::endl;
    }

    void setBoundaryMask(const std::vector<int> &localMask) {
        _mask = &localMask;
    }

    //
    void apply(const VectorType &Alpha, VectorType &Dest) const override {
        _counter++;
        if (!_quiet) std::cerr << "Eval local proj. gradient." << std::endl;
        auto t_start = std::chrono::high_resolution_clock::now();
        
        if (Alpha.size() != _numDOFs)
            throw BasicException("NonlinearProjectionGradient::apply: wrong size of dofs!");
        if (Dest.size() != Alpha.size())
            Dest.resize(Alpha.size());

        VectorType extendedAlpha;
        bool validExtension = isAlphaExtensionValid<ConfiguratorType>(Alpha, extendedAlpha);
        
        if (extendedAlpha.size() != _numDOFs + 1)
            throw BasicException("NonlinearProjectionGradient::apply: extendedAlpha has wrong size, should not happen!");

        // bring into more convenient form
        std::vector<Eigen::Ref<const VectorType> > princVarsRefs;
        princVarsRefs.reserve(_numDOFs + 1);
        int numLocalDofs = 3 * _topology.getNumVertices();
        int dimMatrix = _mask ? numLocalDofs : numLocalDofs + 6;
        for (int k = 0; k <= _numDOFs; k++)
            princVarsRefs.push_back(_princVariations.segment(k * numLocalDofs, numLocalDofs));

        // compute the local elastic average $q[\alpha] = \argmin_q \sum_j \alpha_j W[p_j, q]$
        if (!_quiet) std::cerr << "Solve elastic average for q." << std::endl;
        _fastElasticMeanOp.setAlpha( extendedAlpha );
        _fastElasticMeanOp.execute( _q );

        // solve for dual variable
        VectorType mu, negativeRhs;
        _W.applyDefGradient(_shape, _q, negativeRhs);

        // fill triplet lists and assembla matrix
        if (!_quiet) std::cerr << "Assemble matrix to solve for dual variable mu." << std::endl;
        TripletListType tripletList;
        tripletList.reserve(extendedAlpha.size() * _W.numOfNonZeroHessianEntries());
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (uint j = 0; j < extendedAlpha.size(); j++) {
            TripletListType localTripletList;
            localTripletList.reserve(_W.numOfNonZeroHessianEntries());
            _W.pushTripletsDefHessian(princVarsRefs[j], _q, localTripletList, 0, 0, -1. * extendedAlpha[j]);
#ifdef _OPENMP
#pragma omp critical
#endif
            tripletList.insert(tripletList.end(), localTripletList.begin(), localTripletList.end());
        }

        // if there is no boundary mask, we make use of Lagrange setting and push constraint triplets
        if (!_mask) {
            if (!_quiet) std::cerr << "Account for Lagrange setup." << std::endl;
            // push triplets to fix zeroth and first momentum
            if (_constraintTriplets.size() == 0) {
                if (!_quiet) std::cerr << "Push constraint triplets." << std::endl;
                RigidBodyMotionsConstraintHandler<ConfiguratorType>(_topology, _shape, 1).addConstraintHessian(
                        _constraintTriplets);
            }
            tripletList.insert(tripletList.end(), _constraintTriplets.begin(), _constraintTriplets.end());

            // extend rhs
            negativeRhs.conservativeResize(dimMatrix);
            for (int i = 0; i < 6; i++)
                negativeRhs[numLocalDofs + i] = 0.;
        }

        // assemble matrix
        if (!_quiet) std::cerr << "Assemble matrix." << std::endl;
        MatrixType negativeHessian(dimMatrix, dimMatrix);
        negativeHessian.setFromTriplets(tripletList.cbegin(), tripletList.cend());

        // mask matrix and rhs
        if (_mask) {
            if (!_quiet) std::cerr << "Mask system matrix and rhs." << std::endl;
            applyMaskToSymmetricMatrixAndVector(*_mask, negativeHessian, negativeRhs);
        }

        // solve linear system
        if (!_quiet) std::cerr << "Solve linear system for mu." << std::endl;
        LinearSolver<ConfiguratorType>().solve(negativeHessian, negativeRhs, mu);

        // resize solution
        if (!_mask)
            mu.conservativeResize(numLocalDofs);

        // compute gradient
        if (!_quiet) std::cerr << "Compute final gradient." << std::endl;
        VectorType constGradPart;
        _W.applyDefGradient(princVarsRefs[0], _q, constGradPart);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j < _numDOFs; j++) {
            VectorType tempGradPart;
            _W.applyDefGradient(princVarsRefs[j + 1], _q, tempGradPart);
            tempGradPart -= constGradPart;
            Dest[j] = mu.dot(tempGradPart);
        }
     
        // penalty gradient
        if( _penaltyWeight > 0. ){
          RealType constPenaltyDeriv = penaltyDerivative( extendedAlpha[0], 1./_penaltyWeight );
          for( int i = 0; i < _numDOFs; i++ )
            Dest[i] += penaltyDerivative( Alpha[i], 1./_penaltyWeight ) - constPenaltyDeriv;
        }

        //std::cerr << "\n Grad NL proh = " << Dest << std::endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        if(!_quiet) std::cout << std::fixed << "Gradient computed in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << "seconds." << std::endl << std::endl;
    }
};

//! \brief Perform nonlinear local projection.
//! \author Heeren
//! \note Compare to eq. (18) in \cite HeZhRu18
template<typename ConfiguratorType>
class NonlinearLocalProjection {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const MeshTopologySaver &_topology;
    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shape;
    const VectorType &_principalVariations;


    const OptimizationParameters<ConfiguratorType>& _optParsOuter, _optParsInner;
    std::vector<int> const *_mask;
    bool _quiet;
    RealType _penaltyWeight;
    const int _numDOFs;

public:
    NonlinearLocalProjection(const MeshTopologySaver &topology,
                        const DeformationBase<ConfiguratorType> &W,
                        const VectorType &principalVariations,
                        const VectorType &shape,
                        const OptimizationParameters<ConfiguratorType> &optParsOuter,
                        const OptimizationParameters<ConfiguratorType> &optParsInner,
                        bool quiet = true)
            : _topology(topology), 
              _W(W),  
              _optParsOuter(optParsOuter), 
              _optParsInner(optParsInner), 
              _shape(shape),
              _principalVariations(principalVariations), 
              _quiet(quiet), 
              _penaltyWeight(0.),
              _numDOFs(_principalVariations.size() / (3 * _topology.getNumVertices()) - 1),
              _mask(NULL) { }

    void setBoundaryMask(const std::vector<int> &localMask) {
        _mask = &localMask;
    }

    void execute( VectorType &Alpha, VectorType& optShape ) const {
        NonlinearProjectionFunctional<ConfiguratorType> F(_topology, _W, _principalVariations, _shape, _optParsInner, _penaltyWeight);
        NonlinearProjectionGradient<ConfiguratorType> DF(_topology, _W, _principalVariations, _shape, _optParsInner, _penaltyWeight);

        if( Alpha.size() != _numDOFs){
            Alpha.resize( _numDOFs );
            Alpha.setConstant(1. / _numDOFs);
        }

        VectorType initialization(Alpha);

        RealType energy;
        VectorType grad;
        if(!_quiet){
            F.apply(initialization, energy);
            std::cerr << "Initial nonlinear projection functional  = " << energy << std::endl;
            DF.apply(initialization, grad);
            std::cerr << "Initial nonlinear projection gradient norm = " << grad.norm() << std::endl << std::endl;
        }
        
        // gradient test
        if( false ){
          std::cerr << "Start gradient test!" << std::endl;
          ScalarValuedDerivativeTester<ConfiguratorType> ( F, DF, 1e-3 ).plotAllDirections ( initialization, "initial_gradTest" );
        }

        // Optimization with gradient descent
        if (_optParsOuter.getGradientIterations() > 0) {
            if(!_quiet) std::cerr << "Start gradient descent... " << std::endl;
            GradientDescent<ConfiguratorType> GD(F, DF, _optParsOuter);
            if (_mask) {
                if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                GD.setBoundaryMask(*_mask);
            }
            GD.solve(initialization, Alpha);
            initialization = Alpha;
        }

        // optimization with BFGS
        if (_optParsOuter.getBFGSIterations() > 0) {
            if(!_quiet) std::cerr << "Start Quasi-Newton... " << std::endl;
            QuasiNewtonBFGS<ConfiguratorType> QNM(F, DF, _optParsOuter);
            if (_mask) {
                if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                QNM.setBoundaryMask(*_mask);
            }
            QNM.solve(initialization, Alpha);
            initialization = Alpha;
        }

        if(!_quiet){
            F.apply(Alpha, energy);
            std::cerr << "Final energy = " << energy << std::endl;
            DF.apply(Alpha, grad);
            std::cerr << "Final gradient norm = " << grad.norm() << std::endl << std::endl;
        }
        
        // gradient test
        if( false ){
          std::cerr << "Start gradient test!" << std::endl;
          ScalarValuedDerivativeTester<ConfiguratorType> ( F, DF, 1e-4 ).plotAllDirections ( initialization, "final_gradTest" );
        }
        
        // 
        optShape = F.getShape();

    }

};


/*
//! \brief Perform nonlinear global projection.
//! \author Heeren
//! \note Compare to Fig. 7, right hand side, in \cite HeZhRu18.
template<typename ConfiguratorType, typename ShellDeformationType>
class NonlinearGlobalProjection {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;
    
    const MeshTopologySaver &_topology;
    const ParameterParser& _pparser;
    const DiscretePGAModel<ConfiguratorType, ShellDeformationType>& _pgaModel;
    ShellDeformationType _W;

    std::vector<int> const *_mask;
    bool _quiet;
    RealType _penaltyWeight;
    const int _numLocalDofs, _numPrincVars;

public:
    NonlinearGlobalProjection( const MeshTopologySaver &topology,
                               const ParameterParser& pparser,
                               const DiscretePGAModel<ConfiguratorType, ShellDeformationType>& pgaModel,                              
                               bool quiet = true )
            : _topology(topology), 
              _pparser(pparser),
              _pgaModel(pgaModel),
              _W(_topology, pparser.getDouble("bendWeight") ),  
              _quiet(quiet), 
              _penaltyWeight(0.),
              _numLocalDofs( 3 * _topology.getNumVertices() ),
              _numPrincVars( _pgaModel.getNumPrincipalVariations() ),
              _mask(NULL) { }
              
              
    void execute( const VectorType &unseenGeom, VectorType &Alpha, VectorType& globalProjection ) const {
        
        if( unseenGeom.size() != _numLocalDofs )
            throw BasicException("NonlinearGlobalProjection::execute(): argument has wrong size!");
        
        // #### Principal Variations ####
        if(!_quiet) std::cerr << "\n\n====================================================================================================" << std::endl;   
        if(!_quiet) std::cerr << "START NONLINEAR PROJECTION." << std::endl; 
        auto t_start_proj = std::chrono::high_resolution_clock::now();        
   
        // local projection optim. params
        OptimizationParameters<ConfiguratorType> optParsProjInner;
        optParsProjInner.setGradientIterations( _pparser.getInt("localProjInnerGradDescIters"));
        optParsProjInner.setNewtonIterations( _pparser.getInt("localProjInnerNewtonIters") );
        optParsProjInner.setQuiet();
        //optParsProj.setVerbose();
   
        OptimizationParameters<ConfiguratorType> optParsProjOuter;
        optParsProjOuter.setGradientIterations( _pparser.getInt("localProjOuterGradDescIters"));
        optParsProjOuter.setBFGSIterations( _pparser.getInt("localProjOuterBFGSIters") );
        optParsProjOuter.setVerbose();
 
        // compute scaling factor
        RealType rho = computeScalingFactor( unseenGeom );
        int K = _pgaModel.getK();
        
        // set boundary mask
        std::vector<int> mask;
        if( _mask )
            mask = *_mask;
   
        // optim. params for scaling
        OptimizationParameters<ConfiguratorType> optParsScale;
        optParsScale.setGradientIterations( _pparser.getInt("geodScalingGradDescIters") );
        optParsScale.setBFGSIterations( _pparser.getInt("geodScalingBFGSIters") );
        optParsScale.setNewtonIterations( _pparser.getInt("geodScalingNewtonIters") );
        optParsScale.setInitializationScheme( _pparser.getInt("geodScalingInitializationScheme") );
        //optParsScale.setVerbose();
        optParsScale.setQuiet();   
   
        // scale geometry  
        VectorType geodesicPathToUnseenShape;
        if(!_quiet) std::cerr << "1) SCALING." << std::endl;
        VectorType scaledGeom( _pgaModel.getMean() );
        geodesicScaling<ConfiguratorType>( _topology, _pgaModel.getMean(), unseenGeom, mask, _W, optParsScale, K + 1, rho, geodesicPathToUnseenShape, scaledGeom, false );


        auto t_start_loc_proj = std::chrono::high_resolution_clock::now(); 
        if(!_quiet) std::cerr << "\n2) LOCAL PROJECTION." << std::endl; 
        
        VectorType localProjection( _pgaModel.getMean() );
        NonlinearLocalProjection<ConfiguratorType> localProjOp( _topology, _W, _pgaModel.getPrincipalVariations(), scaledGeom, optParsProjOuter, optParsProjInner, _penaltyWeight, _quiet );   
        Alpha.resize( 2 * _numPrincVars );
        for( int i = 0; i < 2 * _numPrincVars; i++ )
          Alpha[i] = 1. / (2.*_numPrincVars + 1.); 
        localProjOp.execute( Alpha, localProjection );
        if(!_quiet) std::cerr << "Optimal weights = " << Alpha << std::endl;
        if(!_quiet) std::cerr << "sum = " << Alpha.sum() << std::endl << std::endl;
        
        RealType tempE;
        _W.applyEnergy( scaledGeom, localProjection, tempE);
        if(!_quiet) std::cerr << "Local projection error = " << tempE[0] << std::endl;
        
        auto t_end_loc_proj = std::chrono::high_resolution_clock::now();
        if(!_quiet) std::cout << std::fixed << "Local projection done in " << std::chrono::duration<double, std::ratio<1> >(t_end_loc_proj - t_start_loc_proj).count() << "seconds." << std::endl;
   
        // rescale 
        if(!_quiet) std::cerr << "3) RESCALING." << std::endl;
        globalProjection = localProjection;
        geodesicScalingInverse<ConfiguratorType>( _topology, _pgaModel.getMean(), localProjection, mask, _W, optParsScale, K + 1, rho, globalProjection );
        _W.applyEnergy( globalProjection, unseenGeom, tempE);
        if(!_quiet) std::cerr << "Global projection error = " << tempE[0] << std::endl;
   
        auto t_end_proj = std::chrono::high_resolution_clock::now();
        if(!_quiet) std::cout << std::fixed << "Nonlinear projection done in " << std::chrono::duration<double, std::ratio<1> >(t_end_proj - t_start_proj).count() << "seconds." << std::endl;
    }          
                        
protected:
    //
    RealType computeScalingFactor( const VectorType& unseenGeom ) const {
        if(!_quiet) std::cerr << "Compute scaling factor" << std::endl;
        RealType kappa = _pparser.getDoubleOrDefault("kappa", 0.5);
        RealType rho = FLT_MAX;
   
        // compare against all active pcs
        RealType tempE;
        for( int i = 0; i < _numPrincVars; i++ ){       
          _W.applyEnergy ( _pgaModel.getMean(), _pgaModel.getPrincipalVariations().segment( (i + 1) * _numLocalDofs, _numLocalDofs), tempE);
          if(!_quiet) std::cerr << "dist[ mean, pc_" << i << "] approx " << std::sqrt( tempE[0] ) << std::endl;
          if( rho > std::sqrt( tempE[0] ) )
            rho = std::sqrt( tempE[0] );
        }
        _W.applyEnergy ( _pgaModel.getMean(), unseenGeom, tempE);
        if(!_quiet) std::cerr << "dist[ mean, s]  approx " << std::sqrt( tempE[0] ) << std::endl;
        return rho * kappa / std::sqrt( tempE[0] );
    }
                        
};
*/

#endif