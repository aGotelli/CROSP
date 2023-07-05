#ifndef COSSERAT_ROD_INTEGRATORS_HPP
#define COSSERAT_ROD_INTEGRATORS_HPP





#include <Eigen/Dense>
#include <memory>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"




namespace CROSP::numerical_integrators {


struct FullODEStatesObservations {


    Eigen::MatrixXd orientation_stack;
    Eigen::MatrixXd r_stack;

    Eigen::MatrixXd Omega_stack;
    Eigen::MatrixXd V_stack;

    Eigen::MatrixXd dot_Omega_stack;
    Eigen::MatrixXd dot_V_stack;

    Eigen::MatrixXd C_stack;
    Eigen::MatrixXd N_stack;

    Eigen::MatrixXd Qa_stack;
    Eigen::MatrixXd Qad_stack;

};





template<class NumericalIntegrator>
concept CosseratIntegrator = requires(NumericalIntegrator integrator){

        NumericalIntegrator(polynomial_representation::PolynomialRepresentation(),
                            int(),
                            std::shared_ptr<rod_properties::RodProperties>());

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


        FullODEStatesObservations() = integrator.getFullODEStatesObservations();

};




}   //  namespace CROSP::numerical_integrators





#endif // COSSERAT_ROD_INTEGRATORS_HPP
