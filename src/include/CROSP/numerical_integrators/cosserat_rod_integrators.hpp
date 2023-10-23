#ifndef COSSERAT_ROD_INTEGRATORS_HPP
#define COSSERAT_ROD_INTEGRATORS_HPP





#include <Eigen/Dense>
#include <memory>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"

#include <ATORS/ATORS.hpp>

#include "full_state_reconstructor.hpp"


namespace CROSP::numerical_integrators {



template<class NumericalIntegrator>
concept CosseratIntegrator = requires(NumericalIntegrator integrator){

        NumericalIntegrator(polynomial_representation::PolynomialRepresentation(),
                            int(),
                            std::shared_ptr<rod_properties::RodProperties>());


        integrator.addInternalActuation(::ATORS::distributed_actuation::DistributedActuationUptr());


        integrator.printIntegratorProperties();


        integrator.updateParameterisation(Eigen::VectorXd(),
                                          Eigen::VectorXd(),
                                          Eigen::VectorXd());


        integrator.forwardKinematics();


        integrator.forwardKinematics(Eigen::Vector4d(),
                                       Eigen::Vector3d(),
                                       Eigen::Vector3d(),
                                       Eigen::Vector3d(),
                                       Eigen::Vector3d(),
                                       Eigen::Vector3d());


       integrator.updateDeltaParameterisation(Eigen::VectorXd(),
                                                Eigen::VectorXd(),
                                                Eigen::VectorXd());

        integrator.forwardTangentKinematics();


        integrator.forwardTangentKinematics(Eigen::Vector3d(),
                                              Eigen::Vector3d(),
                                              Eigen::Vector3d(),
                                              Eigen::Vector3d(),
                                              Eigen::Vector3d(),
                                              Eigen::Vector3d());


        ::LieAlgebra::Kinematics() = integrator.getKinematicsAtTip();

        ::LieAlgebra::TangentKinematics() = integrator.getTangentKinematicsAtTip();



        integrator.backwardDynamics(::LieAlgebra::Vector6d());


        integrator.updateInternalActuation(double());

        Eigen::VectorXd() = integrator.getQad();



        integrator.backwardTangentDynamics(::LieAlgebra::Vector6d());




        ::LieAlgebra::Vector6d() = integrator.getLambdaAtBase();



        LieAlgebra::Vector6d() = integrator.getDeltaLambdaAtBase();


        ::LieAlgebra::Vector6d() = integrator.getQaAtBase();



        LieAlgebra::Vector6d() = integrator.getDeltaQaAtBase();


        integrator.updateIntegrationDomain(double());


        Eigen::MatrixXd() = integrator.getRodPositions();


        FullODEStatesObservations() = integrator.getIDMStatesObservations();

        FullODEStatesObservations() = integrator.getTIDMStatesObservations();

        Eigen::MatrixXd() = integrator.integratePhiTPhi();

};




}   //  namespace CROSP::numerical_integrators





#endif // COSSERAT_ROD_INTEGRATORS_HPP
