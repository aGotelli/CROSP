#ifndef STRAIN_PARAMETERISATION_STACK_HPP
#define STRAIN_PARAMETERISATION_STACK_HPP

#include <memory>
#include "CROSP/polynomial_representation/polynomial_representation.hpp"

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

namespace CROSP::strain_parameterisation_stack {


typedef std::function<::LieAlgebra::Vector6d(const double&)> StrainFunction;


static StrainFunction default_constrained_strain = [](const double&){return ::LieAlgebra::Vector6d::Unit(3);};

static StrainFunction zero_constrained_strain = [](const double&){return ::LieAlgebra::Vector6d::Zero();};


struct StrainParameterisationStack {

    StrainParameterisationStack(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                const unsigned int t_number_of_Chebyshev_point,
                                const StrainFunction &t_Xi_c=default_constrained_strain);

    StrainParameterisationStack(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                const std::vector<double> &t_observation_points,
                                const StrainFunction &t_Xi_c=default_constrained_strain);


    void updateStrainParameterisation(const Eigen::VectorXd &t_q,
                                      const Eigen::VectorXd &t_dot_q,
                                      const Eigen::VectorXd &t_ddot_q);


    unsigned int m_number_of_points;

    StrainFunction m_Xi_c;


    std::vector<Eigen::MatrixXd> m_Phi_stack;

    Eigen::MatrixXd m_B_Phi_stack;

    Eigen::MatrixXd m_Xi_stack;
    Eigen::MatrixXd m_dot_Xi_stack;
    Eigen::MatrixXd m_ddot_Xi_stack;
    Eigen::MatrixXd m_Xi_c_stack;


    std::unique_ptr<std::vector<Eigen::Ref<Eigen::Vector3d>>> m_k_stack_;
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>()
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>()
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_K_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>()
    };


    std::shared_ptr<std::vector<Eigen::Vector3d>> m_Gamma_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>()
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_dot_Gamma_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>()
    };
    std::shared_ptr<std::vector<Eigen::Vector3d>> m_ddot_Gamma_stack {
        std::make_shared<std::vector<Eigen::Vector3d>>()
    };


    std::shared_ptr<std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        std::make_shared<std::vector<Eigen::Matrix3d>>()
    };
    std::shared_ptr<std::vector<Eigen::Matrix3d>> m_hat_dot_K_stack {
        std::make_shared<std::vector<Eigen::Matrix3d>>()
    };
    std::shared_ptr<std::vector<Eigen::Matrix3d>> m_hat_ddot_K_stack {
        std::make_shared<std::vector<Eigen::Matrix3d>>()
    };


    std::shared_ptr<std::vector<Eigen::Matrix3d>> m_hat_Gamma_stack {
        std::make_shared<std::vector<Eigen::Matrix3d>>()
    };
    std::shared_ptr<std::vector<Eigen::Matrix3d>> m_hat_dot_Gamma_stack {
        std::make_shared<std::vector<Eigen::Matrix3d>>()
    };
    std::shared_ptr<std::vector<Eigen::Matrix3d>> m_hat_ddot_Gamma_stack {
        std::make_shared<std::vector<Eigen::Matrix3d>>()
    };


};

}   //  namespace CROSP::strain_parameterisation_stack

#endif // STRAIN_PARAMETERISATION_STACK_HPP
