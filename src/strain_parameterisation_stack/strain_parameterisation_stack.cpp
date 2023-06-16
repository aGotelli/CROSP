#include "CROSP/strain_parameterisation_stack/strain_parameterisation_stack.hpp"

namespace CROSP::strain_parameterisation_stack {





StrainParameterisationStack::StrainParameterisationStack(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                                         const unsigned int t_number_of_Chebyshev_point,
                                                         const StrainFunction &t_Xi_c)
    : m_number_of_points(t_number_of_Chebyshev_point),
      m_Xi_c(t_Xi_c)
{

    const Eigen::MatrixXd B = t_polynomial_representation.m_B;
    Eigen::MatrixXd Phi = t_polynomial_representation.getPhi( 0 );
    Eigen::MatrixXd BPhi = B*Phi;


    m_B_Phi_stack =
            Eigen::MatrixXd(t_number_of_Chebyshev_point*BPhi.rows(), BPhi.cols());

    m_Xi_stack =
            Eigen::MatrixXd(t_number_of_Chebyshev_point*6, 1);

    m_dot_Xi_stack =
            Eigen::MatrixXd(t_number_of_Chebyshev_point*6, 1);

    m_ddot_Xi_stack =
            Eigen::MatrixXd(t_number_of_Chebyshev_point*6, 1);

    m_Xi_c_stack =
            Eigen::MatrixXd(t_number_of_Chebyshev_point*6, 1);


    ::LieAlgebra::Vector6d Xi_c = ::LieAlgebra::Vector6d::Zero();
    Xi_c(3) = 1;

    unsigned int row;
    const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(t_number_of_Chebyshev_point);

    for(unsigned int step=0; const auto point : Chebyshev_points){
            Phi = t_polynomial_representation.getPhi( point );

            row = step * BPhi.rows();

            m_B_Phi_stack.block(row, 0, BPhi.rows(), BPhi.cols()) = B*Phi;

            m_Xi_c_stack.block(row, 0, BPhi.rows(), 1) = m_Xi_c(point);


            m_K_stack->push_back( Eigen::Vector3d::Zero() );
            m_dot_K_stack->push_back( Eigen::Vector3d::Zero() );
            m_ddot_K_stack->push_back( Eigen::Vector3d::Zero() );

            m_Gamma_stack->push_back( Eigen::Vector3d::Zero() );
            m_dot_Gamma_stack->push_back( Eigen::Vector3d::Zero() );
            m_ddot_Gamma_stack->push_back( Eigen::Vector3d::Zero() );




            m_hat_K_stack->push_back( Eigen::Matrix3d::Zero() );
            m_hat_dot_K_stack->push_back( Eigen::Matrix3d::Zero() );
            m_hat_ddot_K_stack->push_back( Eigen::Matrix3d::Zero() );

            m_hat_Gamma_stack->push_back( Eigen::Matrix3d::Zero() );
            m_hat_dot_Gamma_stack->push_back( Eigen::Matrix3d::Zero() );
            m_hat_ddot_Gamma_stack->push_back( Eigen::Matrix3d::Zero() );

            step++;
    }


}

StrainParameterisationStack::StrainParameterisationStack(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                                         const std::vector<double> &t_observation_points,
                                                         const StrainFunction &t_Xi_c)
    : m_number_of_points(t_observation_points.size()),
      m_Xi_c(t_Xi_c)
{

    const Eigen::MatrixXd B = t_polynomial_representation.m_B;
    Eigen::MatrixXd Phi = t_polynomial_representation.getPhi( 0 );
    Eigen::MatrixXd BPhi = B*Phi;


    m_B_Phi_stack =
            Eigen::MatrixXd(m_number_of_points*BPhi.rows(), BPhi.cols());

    m_Xi_stack =
            Eigen::MatrixXd(m_number_of_points*6, 1);

    m_dot_Xi_stack =
            Eigen::MatrixXd(m_number_of_points*6, 1);

    m_ddot_Xi_stack =
            Eigen::MatrixXd(m_number_of_points*6, 1);

    m_Xi_c_stack =
            Eigen::MatrixXd(m_number_of_points*6, 1);


    ::LieAlgebra::Vector6d Xi_c = ::LieAlgebra::Vector6d::Zero();
    Xi_c(3) = 1;

    unsigned int row;
    for(unsigned int step=0; const auto point : t_observation_points){
            Phi = t_polynomial_representation.getPhi( point );

            row = step * BPhi.rows();

            m_B_Phi_stack.block(row, 0, BPhi.rows(), BPhi.cols()) = B*Phi;

            m_Xi_c_stack.block(row, 0, BPhi.rows(), 1) = m_Xi_c(point);


            m_K_stack->push_back( Eigen::Vector3d::Zero() );
            m_dot_K_stack->push_back( Eigen::Vector3d::Zero() );
            m_ddot_K_stack->push_back( Eigen::Vector3d::Zero() );

            m_Gamma_stack->push_back( Eigen::Vector3d::Zero() );
            m_dot_Gamma_stack->push_back( Eigen::Vector3d::Zero() );
            m_ddot_Gamma_stack->push_back( Eigen::Vector3d::Zero() );




            m_hat_K_stack->push_back( Eigen::Matrix3d::Zero() );
            m_hat_dot_K_stack->push_back( Eigen::Matrix3d::Zero() );
            m_hat_ddot_K_stack->push_back( Eigen::Matrix3d::Zero() );

            m_hat_Gamma_stack->push_back( Eigen::Matrix3d::Zero() );
            m_hat_dot_Gamma_stack->push_back( Eigen::Matrix3d::Zero() );
            m_hat_ddot_Gamma_stack->push_back( Eigen::Matrix3d::Zero() );

            step++;
    }


}


void StrainParameterisationStack::updateStrainParameterisation(const Eigen::VectorXd &t_q,
                                                               const Eigen::VectorXd &t_dot_q,
                                                               const Eigen::VectorXd &t_ddot_q)
{
    m_Xi_stack = m_B_Phi_stack * t_q + m_Xi_c_stack;
    m_dot_Xi_stack = m_B_Phi_stack * t_dot_q;
    m_ddot_Xi_stack = m_B_Phi_stack * t_ddot_q;

    Eigen::Vector3d K;
    Eigen::Vector3d dot_K;
    Eigen::Vector3d ddot_K;

    Eigen::Vector3d Gamma;
    Eigen::Vector3d dot_Gamma;
    Eigen::Vector3d ddot_Gamma;

    for(unsigned int step=0; step<m_number_of_points; step++){

        K      =      m_Xi_stack.block<3, 1>(step*6, 0);
        dot_K  =  m_dot_Xi_stack.block<3, 1>(step*6, 0);
        ddot_K = m_ddot_Xi_stack.block<3, 1>(step*6, 0);

        Gamma      =      m_Xi_stack.block<3, 1>(step*6 + 3, 0);
        dot_Gamma  =  m_dot_Xi_stack.block<3, 1>(step*6 + 3, 0);
        ddot_Gamma = m_ddot_Xi_stack.block<3, 1>(step*6 + 3, 0);


        m_K_stack->at(step)      = K;
        m_dot_K_stack->at(step)  = dot_K;
        m_ddot_K_stack->at(step) = ddot_K;

        m_Gamma_stack->at(step)      = Gamma;
        m_dot_Gamma_stack->at(step)  = dot_Gamma;
        m_ddot_Gamma_stack->at(step) = ddot_Gamma;



        m_hat_K_stack->at(step)      = ::LieAlgebra::skew( K );
        m_hat_dot_K_stack->at(step)  = ::LieAlgebra::skew( dot_K );
        m_hat_ddot_K_stack->at(step) = ::LieAlgebra::skew( ddot_K );

        m_hat_Gamma_stack->at(step)      = ::LieAlgebra::skew( Gamma );
        m_hat_dot_Gamma_stack->at(step)  = ::LieAlgebra::skew( dot_Gamma );
        m_hat_ddot_Gamma_stack->at(step) = ::LieAlgebra::skew( ddot_Gamma );

    }
}





}   //  namespace CROSP::strain_parameterisation_stack
