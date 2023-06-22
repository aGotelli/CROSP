#ifndef COSSERAT_ROD_INTEGRATORS_HPP
#define COSSERAT_ROD_INTEGRATORS_HPP





#include <Eigen/Dense>
#include <memory>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"


namespace CROSP::numerical_integrators {



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



        integrator.backwardTangentDynamics(::LieAlgebra::Vector6d());




        ::LieAlgebra::Vector6d() = integrator.getLambdaAtBase();



        LieAlgebra::Vector6d() = integrator.getDeltaLambdaAtBase();


        ::LieAlgebra::Vector6d() = integrator.getQaAtBase();



        LieAlgebra::Vector6d() = integrator.getDeltaQaAtBase();


        integrator.updateIntegrationDomain(double());


        Eigen::MatrixXd() = integrator.getRodPositions();

};


struct CosseratRodIntegrators {


    virtual ~CosseratRodIntegrators()=default;


    virtual std::string printIntegratorProperties()const=0;


    virtual void updateParameterisation(const Eigen::VectorXd &t_qe,
                                        const Eigen::VectorXd &t_dot_qe,
                                        const Eigen::VectorXd &t_ddot_qe)=0;


    virtual void forwardKinematics()=0;


    virtual void forwardKinematics(const Eigen::Vector4d &t_initial_quaternion,
                                   const Eigen::Vector3d &t_initial_position,
                                   const Eigen::Vector3d &t_initial_angular_velocity,
                                   const Eigen::Vector3d &t_initial_linear_velocity,
                                   const Eigen::Vector3d &t_initial_angular_acceleration,
                                   const Eigen::Vector3d &t_initial_linear_acceleration)=0;


   virtual void updateDeltaParameterisation(const Eigen::VectorXd &t_Delta_qe,
                                            const Eigen::VectorXd &t_Delta_dot_qe,
                                            const Eigen::VectorXd &t_Delta_ddot_qe)=0;

    virtual void forwardTangentKinematics()=0;


    virtual void forwardTangentKinematics(const Eigen::Vector3d &t_initial_Delta_orientation,
                                          const Eigen::Vector3d &t_initial_Delta_position,
                                          const Eigen::Vector3d &t_initial_Delta_angular_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_linear_velocity,
                                          const Eigen::Vector3d &t_initial_Delta_angular_acceleration,
                                          const Eigen::Vector3d &t_initial_Delta_linear_acceleration)=0;


    virtual ::LieAlgebra::Kinematics getKinematicsAtTip()const=0;

    virtual ::LieAlgebra::TangentKinematics getTangentKinematicsAtTip()const=0;



    virtual void backwardDynamics(const ::LieAlgebra::Vector6d &t_Lambda_X1)=0;



    virtual void backwardTangentDynamics(const ::LieAlgebra::Vector6d &t_Delta_Lambda_X1)=0;




    virtual ::LieAlgebra::Vector6d getLambdaAtBase()const=0;



    virtual LieAlgebra::Vector6d getDeltaLambdaAtBase()const=0;


    virtual ::LieAlgebra::Vector6d getQaAtBase()const=0;



    virtual LieAlgebra::Vector6d getDeltaQaAtBase()const=0;


    virtual void updateIntegrationDomain(const double &t_rod_lenght)=0;


    virtual Eigen::MatrixXd getRodPositions()const=0;



};

typedef std::unique_ptr<CosseratRodIntegrators> CosseratIntegratorUPtr;

typedef std::shared_ptr<CosseratRodIntegrators> CosseratIntegratorSPtr;



}   //  namespace CROSP::numerical_integrators





#endif // COSSERAT_ROD_INTEGRATORS_HPP
