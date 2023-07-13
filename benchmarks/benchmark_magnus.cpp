 #include "OMNI/OMNI.hpp"
#include "memory"

#include "CROSP/CROSP/cosserat_rod.hpp"

#include "benchmark/benchmark.h"


struct IntegratorRotationMatrix : public ::OMNI::ODEA<3, 3>{

    IntegratorRotationMatrix(const unsigned int t_number_of_Chebyshev_points,
                             std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack)
        : ::OMNI::ODEA<3, 3>((t_number_of_Chebyshev_points-1)*3),
          m_strain_parameterisation_stack(t_strain_parameterisation_stack),
          m_interpolator(::Chebyshev::ChebyshevInterpolator(t_number_of_Chebyshev_points, this->m_quadrature_points))
    {}


    virtual Eigen::Matrix3d computeCofficientsMatrixAtQuadraturePoint(const unsigned int t_Chebyshev_point,
                                                                      const unsigned int t_quadrature_point)override
    {
        const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;

        return m_strain_parameterisation_stack->m_hat_K_stack->at(stack_index);
//        return Eigen::Matrix<double, 3, 3>::Zero();

    }




    virtual Eigen::Matrix3d recursiveExponentialMapping(const Eigen::Matrix<double, 3, 3> &t_previous_state,
                                                                    const MatrixDd &t_exp_Omega) override
    {
        return t_previous_state * t_exp_Omega;
    }


    virtual Eigen::Matrix3d exponentialOmega(const MatrixDd &t_Omega)override
    {
        return ::LieAlgebra::expRodigues( t_Omega );
    }

//    virtual void postProcessInterpolation()override
//    {

//        Eigen::Quaterniond q;
//        for(unsigned int point=0; point<this->m_number_of_Chebyshev_points; point++){
//            q = Eigen::Quaterniond(m_states_stack[point]);
//            m_quaternion_stack.col(point) << q.w(), q.x(), q.y(), q.z();
//        }

//        m_quaternions_at_quadrature_points = m_interpolator(m_quaternion_stack);
//    }

    virtual Eigen::MatrixXd getStateAtQuadraturePoint(const unsigned int t_quadrature_point)const override
    {
        return m_quaternions_at_quadrature_points.col(t_quadrature_point);
    }




    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;


    Eigen::MatrixXd m_quaternion_stack { Eigen::MatrixXd::Zero(4, this->m_number_of_Chebyshev_points) };

    ::Chebyshev::ChebyshevInterpolator m_interpolator;

    Eigen::MatrixXd m_quaternions_at_quadrature_points { Eigen::MatrixXd::Zero(4, m_interpolator.getNumberOfInterpolationPoints()) };

};










struct IntegratorSO3 : public ::OMNI::MagnusIntegratorSO3{

    IntegratorSO3(const unsigned int t_number_of_Chebyshev_points,
                std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack)
        : ::OMNI::MagnusIntegratorSO3((t_number_of_Chebyshev_points-1)*3),
          m_strain_parameterisation_stack(t_strain_parameterisation_stack),
          m_interpolator(::Chebyshev::ChebyshevInterpolator(t_number_of_Chebyshev_points, this->m_quadrature_points))
    {

        for(unsigned int point=0; point<m_number_of_Chebyshev_points; point++){

//            const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;

            this->m_A1_stack.push_back( m_strain_parameterisation_stack->m_K_stack->at(point*4 + 1) );
            this->m_A2_stack.push_back( m_strain_parameterisation_stack->m_K_stack->at(point*4 + 2) );
            this->m_A3_stack.push_back( m_strain_parameterisation_stack->m_K_stack->at(point*4 + 3) );
        }

    }





    virtual Eigen::MatrixXd getStateAtQuadraturePoint(const unsigned int t_quadrature_point)const override
    {
        return m_quaternions_at_quadrature_points.col(t_quadrature_point);
    }




    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;


    Eigen::MatrixXd m_quaternion_stack { Eigen::MatrixXd::Zero(4, this->m_number_of_Chebyshev_points) };

    ::Chebyshev::ChebyshevInterpolator m_interpolator;

    Eigen::MatrixXd m_quaternions_at_quadrature_points { Eigen::MatrixXd::Zero(4, m_interpolator.getNumberOfInterpolationPoints()) };

};













struct PositionIntegrator : public OSNI::ODEb {

    PositionIntegrator(std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                       const unsigned int t_number_of_Chebyshev_points,
                       std::shared_ptr<const IntegratorRotationMatrix> t_rotation_matrix_integrator,
                       const double &t_upper_integration_limit=1.0f,
                       const Eigen::Vector3d &t_initial_condition=Eigen::Vector3d::Zero())
        : OSNI::ODEb(3,
                     ::Chebyshev::INTEGRATION_DIRECTION::FORWARD,
                     t_number_of_Chebyshev_points,
                     t_upper_integration_limit),
          m_strain_parameterisation_stack(t_strain_parameterisation_stack),
          m_rotation_matrix_integrator(t_rotation_matrix_integrator)
    {
        this->setInitialConditions(t_initial_condition);
    }


    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) final
    {
        const unsigned int stack_index = t_point*4;

        const auto Nc = m_rotation_matrix_integrator->getChebyshevPointsNumber();

        const auto R = m_rotation_matrix_integrator->getStateAtPoint(Nc - 1 - t_point);

//        std::cout << "R at point : " << t_point << "\n" << R << "\n\n";

        return R * m_Gamma_stack->at(stack_index);
//        return Eigen::Matrix<double, 3, 1>::Zero();

    }


    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_Gamma_stack {
        m_strain_parameterisation_stack->m_Gamma_stack
    };

    std::shared_ptr<const IntegratorRotationMatrix> m_rotation_matrix_integrator;


};


struct IntegratorOmega : public ::OMNI::ODEAb<3>{

    IntegratorOmega(const unsigned int t_number_of_Chebyshev_points,
                    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack)
        : ::OMNI::ODEAb<3>(t_number_of_Chebyshev_points),
          m_strain_parameterisation_stack(t_strain_parameterisation_stack),
          m_interpolator(::Chebyshev::ChebyshevInterpolator(t_number_of_Chebyshev_points, this->m_quadrature_points))
    {}



    virtual Eigen::Matrix<double, 3, 3> computeCofficientsMatrixAtQuadraturePoint(const unsigned int t_Chebyshev_point,
                                                                                  const unsigned int t_quadrature_point)override
    {
        const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;
        return -m_strain_parameterisation_stack->m_hat_K_stack->at(stack_index);
//        return Eigen::Matrix<double, 3, 3>::Zero();

    }


    virtual Eigen::Matrix<double, 3, 1> computeParametersVectorAtQuadraturePoint(const unsigned int t_Chebyshev_point,
                                                                                  const unsigned int t_quadrature_point)override
    {
        const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;

        return m_strain_parameterisation_stack->m_dot_K_stack->at(stack_index);
//        return Eigen::Matrix<double, 3, 1>::Zero();

    }


    virtual void postProcessInterpolation()override
    {
        m_states_at_quadrature_points = m_interpolator( this->getStackAsMatrix() );
    }


    virtual Eigen::MatrixXd getStackAsMatrix()const override
    {
        for(unsigned int Chebyshev_point=0; Chebyshev_point<this->m_number_of_Chebyshev_points; Chebyshev_point++){
            this->m_states_stack_as_Eigen_matrix.block<3, 1>(0, 1*Chebyshev_point) =
                    this->m_states_stack[Chebyshev_point].block<3, 1>(0, 0);
        }

        return this->m_states_stack_as_Eigen_matrix;
    }

    virtual Eigen::MatrixXd getStateAtQuadraturePoint(const unsigned int t_quadrature_point)const override
    {
        return m_states_at_quadrature_points.col(t_quadrature_point);
    }



    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;


    ::Chebyshev::ChebyshevInterpolator m_interpolator;

    Eigen::MatrixXd m_states_at_quadrature_points { Eigen::MatrixXd::Zero(3, m_interpolator.getNumberOfInterpolationPoints()) };
};


struct IntegratorV : public ::OMNI::ODEAb<3>{

    IntegratorV(const unsigned int t_number_of_Chebyshev_points,
                std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                std::shared_ptr<const IntegratorOmega> t_angular_velocity)
        : ::OMNI::ODEAb<3>(t_number_of_Chebyshev_points),
          m_strain_parameterisation_stack(t_strain_parameterisation_stack),
          m_angular_velocity(t_angular_velocity),
          m_interpolator(::Chebyshev::ChebyshevInterpolator(t_number_of_Chebyshev_points, this->m_quadrature_points))
    {}



    virtual Eigen::Matrix<double, 3, 3> computeCofficientsMatrixAtQuadraturePoint(const unsigned int t_Chebyshev_point,
                                                                                  const unsigned int t_quadrature_point)override
    {
        const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;
        return -m_strain_parameterisation_stack->m_hat_K_stack->at(stack_index);
//        return Eigen::Matrix<double, 3, 3>::Zero();

    }


    virtual Eigen::Matrix<double, 3, 1> computeParametersVectorAtQuadraturePoint(const unsigned int t_Chebyshev_point,
                                                                                  const unsigned int t_quadrature_point)override
    {
        const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;
        const auto Nc = m_angular_velocity->getChebyshevPointsNumber();
        const auto offset = (Nc-1)*3 - 1;
        const auto quadrature_index = offset - (t_Chebyshev_point*m_number_of_quadrature_points + t_quadrature_point);

        Eigen::Vector3d dot_Gamma = m_dot_Gamma_stack->at(stack_index);
        Eigen::Vector3d Omega = m_angular_velocity->getStateAtQuadraturePoint(quadrature_index);

//        std::cout << "Omega at Chebyshev point : " << t_Chebyshev_point << ", quadrature point : " << t_quadrature_point << "\n" << Omega << "\n\n";


        return dot_Gamma - m_hat_Gamma_stack->at(stack_index)*Omega;
//        return Eigen::Matrix<double, 3, 1>::Zero();


    }


    virtual void postProcessInterpolation()override
    {
        m_states_at_quadrature_points = m_interpolator( this->getStackAsMatrix() );
    }


    virtual Eigen::MatrixXd getStackAsMatrix()const override
    {
        for(unsigned int Chebyshev_point=0; Chebyshev_point<this->m_number_of_Chebyshev_points; Chebyshev_point++){
            this->m_states_stack_as_Eigen_matrix.block<3, 1>(0, 1*Chebyshev_point) =
                    this->m_states_stack[Chebyshev_point].block<3, 1>(0, 0);
        }

        return this->m_states_stack_as_Eigen_matrix;
    }

    virtual Eigen::MatrixXd getStateAtQuadraturePoint(const unsigned int t_quadrature_point)const override
    {
        return m_states_at_quadrature_points.col(t_quadrature_point);
    }



    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        m_strain_parameterisation_stack->m_hat_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_Gamma_stack {
        m_strain_parameterisation_stack->m_hat_Gamma_stack
    };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_dot_Gamma_stack {
        m_strain_parameterisation_stack->m_dot_Gamma_stack
    };

    std::shared_ptr<const IntegratorOmega> m_angular_velocity;


    ::Chebyshev::ChebyshevInterpolator m_interpolator;

    Eigen::MatrixXd m_states_at_quadrature_points { Eigen::MatrixXd::Zero(4, m_interpolator.getNumberOfInterpolationPoints()) };
};


struct IntegratordotOmega : public ::OMNI::ODEAb<3>{

    IntegratordotOmega(const unsigned int t_number_of_Chebyshev_points,
                    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                       std::shared_ptr<const IntegratorOmega> t_angular_velocity)
        : ::OMNI::ODEAb<3>(t_number_of_Chebyshev_points),
          m_strain_parameterisation_stack(t_strain_parameterisation_stack),
          m_angular_velocity(t_angular_velocity),
          m_interpolator(::Chebyshev::ChebyshevInterpolator(t_number_of_Chebyshev_points, this->m_quadrature_points))
    {}



    virtual Eigen::Matrix<double, 3, 3> computeCofficientsMatrixAtQuadraturePoint(const unsigned int t_Chebyshev_point,
                                                                                  const unsigned int t_quadrature_point)override
    {
        const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;
        return -m_strain_parameterisation_stack->m_hat_K_stack->at(stack_index);
//        return Eigen::Matrix<double, 3, 3>::Zero();

    }


    virtual Eigen::Matrix<double, 3, 1> computeParametersVectorAtQuadraturePoint(const unsigned int t_Chebyshev_point,
                                                                                  const unsigned int t_quadrature_point)override
    {
        const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;

        const auto Nc = m_angular_velocity->getChebyshevPointsNumber();
        const auto offset = (Nc-1)*3 - 1;
        const auto quadrature_index = offset - (t_Chebyshev_point*m_number_of_quadrature_points + t_quadrature_point);

        return m_ddot_K_stack->at(stack_index)
                - m_hat_dot_K_stack->at(stack_index) * m_angular_velocity->getStateAtQuadraturePoint(quadrature_index);
//        return Eigen::Matrix<double, 3, 1>::Zero();

    }


    virtual void postProcessInterpolation()override
    {
        m_states_at_quadrature_points = m_interpolator( this->getStackAsMatrix() );
    }


    virtual Eigen::MatrixXd getStackAsMatrix()const override
    {
        for(unsigned int Chebyshev_point=0; Chebyshev_point<this->m_number_of_Chebyshev_points; Chebyshev_point++){
            this->m_states_stack_as_Eigen_matrix.block<3, 1>(0, 1*Chebyshev_point) =
                    this->m_states_stack[Chebyshev_point].block<3, 1>(0, 0);
        }

        return this->m_states_stack_as_Eigen_matrix;
    }

    virtual Eigen::MatrixXd getStateAtQuadraturePoint(const unsigned int t_quadrature_point)const override
    {
        return m_states_at_quadrature_points.col(t_quadrature_point);
    }



    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;


    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        m_strain_parameterisation_stack->m_hat_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_dot_K_stack {
        m_strain_parameterisation_stack->m_hat_dot_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_K_stack {
        m_strain_parameterisation_stack->m_ddot_K_stack
    };

    std::shared_ptr<const IntegratorOmega> m_angular_velocity;


    ::Chebyshev::ChebyshevInterpolator m_interpolator;

    Eigen::MatrixXd m_states_at_quadrature_points { Eigen::MatrixXd::Zero(3, m_interpolator.getNumberOfInterpolationPoints()) };
};


struct IntegratordotV : public ::OMNI::ODEAb<3>{

    IntegratordotV(const unsigned int t_number_of_Chebyshev_points,
                std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> t_strain_parameterisation_stack,
                std::shared_ptr<const IntegratorOmega> t_angular_velocity,
                   std::shared_ptr<const IntegratorV> t_linear_velocity,
                   std::shared_ptr<const IntegratordotOmega> t_angular_acceleration)
        : ::OMNI::ODEAb<3>(t_number_of_Chebyshev_points),
          m_strain_parameterisation_stack(t_strain_parameterisation_stack),
          m_angular_velocity(t_angular_velocity),
          m_linear_velocity(t_linear_velocity),
          m_angular_acceleration(t_angular_acceleration),
          m_interpolator(::Chebyshev::ChebyshevInterpolator(t_number_of_Chebyshev_points, this->m_quadrature_points))
    {}



    virtual Eigen::Matrix<double, 3, 3> computeCofficientsMatrixAtQuadraturePoint(const unsigned int t_Chebyshev_point,
                                                                                  const unsigned int t_quadrature_point)override
    {
        const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;
        return -m_strain_parameterisation_stack->m_hat_K_stack->at(stack_index);

//        return Eigen::Matrix<double, 3, 3>::Zero();

    }


    virtual Eigen::Matrix<double, 3, 1> computeParametersVectorAtQuadraturePoint(const unsigned int t_Chebyshev_point,
                                                                                  const unsigned int t_quadrature_point)override
    {
        const unsigned int stack_index = t_Chebyshev_point*4 + t_quadrature_point+1;

        const auto Nc = m_angular_velocity->getChebyshevPointsNumber();
        const auto offset = (Nc-1)*3 - 1;
        const auto quadrature_index = offset - (t_Chebyshev_point*m_number_of_quadrature_points + t_quadrature_point);

        const Eigen::Vector3d ddot_Gamma = m_ddot_Gamma_stack->at(stack_index);
        const Eigen::Matrix3d hat_Gamma = m_hat_Gamma_stack->at(stack_index);
        const Eigen::Matrix3d hat_dot_Gamma = m_hat_dot_Gamma_stack->at(stack_index);
        const Eigen::Matrix3d hat_dot_K = m_hat_dot_K_stack->at(stack_index);
        const Eigen::Vector3d Omega = m_angular_velocity->getStateAtQuadraturePoint(quadrature_index);
        const Eigen::Vector3d V = m_linear_velocity->getStateAtQuadraturePoint(quadrature_index);
        const Eigen::Vector3d dot_Omega = m_angular_acceleration->getStateAtQuadraturePoint(quadrature_index);

        return ddot_Gamma
                - hat_Gamma * dot_Omega
                - hat_dot_Gamma * Omega
                - hat_dot_K * V;

//        return Eigen::Matrix<double, 3, 1>::Zero();

    }


    virtual void postProcessInterpolation()override
    {
        m_states_at_quadrature_points = m_interpolator( this->getStackAsMatrix() );
    }


    virtual Eigen::MatrixXd getStackAsMatrix()const override
    {
        for(unsigned int Chebyshev_point=0; Chebyshev_point<this->m_number_of_Chebyshev_points; Chebyshev_point++){
            this->m_states_stack_as_Eigen_matrix.block<3, 1>(0, 1*Chebyshev_point) =
                    this->m_states_stack[Chebyshev_point].block<3, 1>(0, 0);
        }

        return this->m_states_stack_as_Eigen_matrix;
    }

    virtual Eigen::MatrixXd getStateAtQuadraturePoint(const unsigned int t_quadrature_point)const override
    {
        return m_states_at_quadrature_points.col(t_quadrature_point);
    }



    std::shared_ptr<const ::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack;

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_K_stack {
        m_strain_parameterisation_stack->m_hat_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_dot_K_stack {
        m_strain_parameterisation_stack->m_hat_dot_K_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_Gamma_stack {
        m_strain_parameterisation_stack->m_hat_Gamma_stack
    };

    std::shared_ptr<const std::vector<Eigen::Matrix3d>> m_hat_dot_Gamma_stack{
        m_strain_parameterisation_stack->m_hat_dot_Gamma_stack
    };

    std::shared_ptr<const std::vector<Eigen::Vector3d>> m_ddot_Gamma_stack{
        m_strain_parameterisation_stack->m_ddot_Gamma_stack
    };

    std::shared_ptr<const IntegratorOmega> m_angular_velocity;

    std::shared_ptr<const IntegratorV> m_linear_velocity;

    std::shared_ptr<const IntegratordotOmega> m_angular_acceleration;


    ::Chebyshev::ChebyshevInterpolator m_interpolator;

    Eigen::MatrixXd m_states_at_quadrature_points { Eigen::MatrixXd::Zero(4, m_interpolator.getNumberOfInterpolationPoints()) };
};


int main(int argc, char *argv[])
{

    const unsigned int ne = 4;

    static constexpr std::array<bool, 6> admitted_deformations = {
        true,
        true,
        true,

        false,
        false,
        false
    };

    static constexpr unsigned int na = std::count(admitted_deformations.begin(),
                                                  admitted_deformations.end(),
                                                  true);

    constexpr unsigned int number_of_Chebyshev_points = 31;


    const unsigned int coordinated_dimension = na * ne;




    const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(number_of_Chebyshev_points);
    auto quadrature_points = ::OMNI::defineQuadrarturePoints(Chebyshev_points);
    auto double_quadrature_points = ::OMNI::defineQuadrarturePoints(quadrature_points);

    std::vector<double> obervation_points = Chebyshev_points;

    obervation_points.insert( obervation_points.end(), quadrature_points.begin(), quadrature_points.end() );
    obervation_points.insert( obervation_points.end(), double_quadrature_points.begin(), double_quadrature_points.end() );

    std::sort(obervation_points.begin(), obervation_points.end());




    auto polynomial_representation =
            std::make_shared<::CROSP::polynomial_representation::PolynomialRepresentation>(admitted_deformations, ne);


    std::shared_ptr<::CROSP::strain_parameterisation_stack::StrainParameterisationStack> m_strain_parameterisation_stack =
            std::make_shared<::CROSP::strain_parameterisation_stack::StrainParameterisationStack>(polynomial_representation, obervation_points);





    std::shared_ptr<IntegratorRotationMatrix> integrate_rotation_matrix =
            std::make_shared<IntegratorRotationMatrix>(number_of_Chebyshev_points, m_strain_parameterisation_stack);

    std::shared_ptr<PositionIntegrator> integrate_position =
            std::make_shared<PositionIntegrator>(m_strain_parameterisation_stack, number_of_Chebyshev_points, integrate_rotation_matrix);

    std::shared_ptr<IntegratorOmega> Omega_integrator =
            std::make_shared<IntegratorOmega>(number_of_Chebyshev_points, m_strain_parameterisation_stack);

    std::shared_ptr<IntegratorV> V_integrator =
            std::make_shared<IntegratorV>(number_of_Chebyshev_points, m_strain_parameterisation_stack, Omega_integrator);

    std::shared_ptr<IntegratordotOmega> dot_Omega_integrator =
            std::make_shared<IntegratordotOmega>(number_of_Chebyshev_points, m_strain_parameterisation_stack, Omega_integrator);

    std::shared_ptr<IntegratordotV> dot_V_integrator =
            std::make_shared<IntegratordotV>(number_of_Chebyshev_points, m_strain_parameterisation_stack, Omega_integrator, V_integrator, dot_Omega_integrator);


    Eigen::VectorXd q      = Eigen::VectorXd::Random(coordinated_dimension);
    const Eigen::VectorXd dot_q  = Eigen::VectorXd::Random(coordinated_dimension);
    const Eigen::VectorXd ddot_q = Eigen::VectorXd::Random(coordinated_dimension);


    q << 0, 0, 0, 0,
          0, 0, 0, 0,
          0.75, -2.5, 0.25, 0.3333;


    const Eigen::Matrix3d R_X0 = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d r_X0 = Eigen::Vector3d::Zero();
    const Eigen::Vector3d Omega = Eigen::Vector3d::Zero();
    const Eigen::Vector3d V_X0 = Eigen::Vector3d::Zero();
    const Eigen::Vector3d dot_Omega = Eigen::Vector3d::Zero();
    const Eigen::Vector3d dot_V_X0 = Eigen::Vector3d::Zero();




//    ::benchmark::RegisterBenchmark("Serial", [&](::benchmark::State &t_state){


//    });

//    ::benchmark::RegisterBenchmark("Serial", [&](::benchmark::State &t_state){

//        unsigned int coord_dim = ne*na;

//        Eigen::MatrixXd BPhi = Eigen::MatrixXd::Random(number_of_Chebyshev_points*na, coord_dim);

//        Eigen::VectorXd q = Eigen::VectorXd::Random(coord_dim);
//        Eigen::VectorXd dot_q = Eigen::VectorXd::Random(coord_dim);
//        Eigen::VectorXd ddot_q = Eigen::VectorXd::Random(coord_dim);

//        Eigen::VectorXd xi = BPhi * q;
//        Eigen::VectorXd dot_xi = BPhi * dot_q;
//        Eigen::VectorXd ddot_xi = BPhi * ddot_q;

//        while(t_state.KeepRunning()){

//            xi = BPhi * q;
//            dot_xi = BPhi * dot_q;
//            ddot_xi = BPhi * ddot_q;

//        }
//    })->Unit(::benchmark::kMicrosecond)->Repetitions(5);



//    ::benchmark::RegisterBenchmark("Parallel", [&](::benchmark::State &t_state){

//        unsigned int coord_dim = ne*na;

//        Eigen::MatrixXd BPhi = Eigen::MatrixXd::Random(number_of_Chebyshev_points*na, coord_dim);

//        Eigen::MatrixXd q = Eigen::VectorXd::Random(coord_dim, 3);


//        Eigen::VectorXd xi_dot_x_ddot_xi = BPhi * q;


//        while(t_state.KeepRunning()){

//            xi_dot_x_ddot_xi = BPhi * q;

//        }
//    })->Unit(::benchmark::kMicrosecond)->Repetitions(5);





//    ::benchmark::RegisterBenchmark("Update", [&](::benchmark::State &t_state){

//        t_state.counters = {
//            {"ne", ne},
//            {"na", na},
//            {"Nc", number_of_Chebyshev_points}
//        };

//        while(t_state.KeepRunning()){


//            m_strain_parameterisation_stack->updateStrainParameterisation(q, dot_q, ddot_q);

//        }
//    })->Unit(::benchmark::kMicrosecond)->Repetitions(5);





    ::benchmark::RegisterBenchmark("Integration R", [&](::benchmark::State &t_state){

        t_state.counters = {
            {"ne", ne},
            {"na", na},
            {"Nc", number_of_Chebyshev_points}
        };






        m_strain_parameterisation_stack->updateStrainParameterisation(q, dot_q, ddot_q);



        while(t_state.KeepRunning()){


            integrate_rotation_matrix->computesCoefficientsMatricesAtQuadraturePoints();
            integrate_rotation_matrix->integrate(R_X0);

        }
    })->Unit(::benchmark::kMicrosecond)->Repetitions(5);


    std::shared_ptr<IntegratorSO3> integrate_SO3 =
            std::make_shared<IntegratorSO3>(number_of_Chebyshev_points, m_strain_parameterisation_stack);

    ::benchmark::RegisterBenchmark("Integration SO3", [&](::benchmark::State &t_state){

        t_state.counters = {
            {"ne", ne},
            {"na", na},
            {"Nc", number_of_Chebyshev_points}
        };







        m_strain_parameterisation_stack->updateStrainParameterisation(q, dot_q, ddot_q);



        while(t_state.KeepRunning()){


            integrate_SO3->computesCoefficientsMatricesAtQuadraturePoints();
            integrate_SO3->integrate(R_X0);

        }
    })->Unit(::benchmark::kMicrosecond)->Repetitions(5);


//    ::benchmark::RegisterBenchmark("Integrations FK", [&](::benchmark::State &t_state){

//        t_state.counters = {
//            {"ne", ne},
//            {"na", na},
//            {"Nc", number_of_Chebyshev_points}
//        };






//        m_strain_parameterisation_stack->updateStrainParameterisation(q, dot_q, ddot_q);

//        integrate_rotation_matrix->computesCoefficientsMatricesAtQuadraturePoints();
//        integrate_rotation_matrix->integrate(R_X0);

//        while(t_state.KeepRunning()){


//            integrate_position->integrate(r_X0);

//            Omega_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//            Omega_integrator->integrate(Omega);

//            V_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//            V_integrator->integrate(V_X0);

//            dot_Omega_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//            dot_Omega_integrator->integrate(dot_Omega);

//            dot_V_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//            dot_V_integrator->integrate(dot_V_X0);
//        }
//    })->Unit(::benchmark::kMicrosecond)->Repetitions(5);



//    ::benchmark::RegisterBenchmark("Fake BD", [&](::benchmark::State &t_state){

//        t_state.counters = {
//            {"ne", ne},
//            {"na", na},
//            {"Nc", number_of_Chebyshev_points}
//        };



//        m_strain_parameterisation_stack->updateStrainParameterisation(q, dot_q, ddot_q);

//        integrate_rotation_matrix->computesCoefficientsMatricesAtQuadraturePoints();
//        integrate_rotation_matrix->integrate(R_X0);


//        while(t_state.KeepRunning()){



//            integrate_position->integrate(r_X0);

//            Omega_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//            Omega_integrator->integrate(Omega);

//            V_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//            V_integrator->integrate(V_X0);

//        }
//    })->Unit(::benchmark::kMicrosecond)->Repetitions(5);



//    m_strain_parameterisation_stack->updateStrainParameterisation(q, dot_q, ddot_q);



//    integrate_rotation_matrix->computesCoefficientsMatricesAtQuadraturePoints();
//    integrate_rotation_matrix->integrate(R_X0);


//    for(unsigned int i=0; i<number_of_Chebyshev_points; i++){
//        std::cout << "c" << i << ", point X=" << Chebyshev_points[i] << "\n";
//        std::cout << "R : \n" << integrate_rotation_matrix->getStateAtPoint(i) << "\n\n";

//    }

//    ::benchmark::Initialize(&argc, argv);

//    ::benchmark::RunSpecifiedBenchmarks();

//    return -64;

//    integrate_position->integrate(r_X0);

//    Omega_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//    Omega_integrator->integrate(Omega);

//    V_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//    V_integrator->integrate(V_X0);

//    dot_Omega_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//    dot_Omega_integrator->integrate(dot_Omega);

//    dot_V_integrator->computesCoefficientsMatricesAtQuadraturePoints();
//    dot_V_integrator->integrate(dot_V_X0);



//    auto strain_parameterisation =
//            std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(polynomial_representation, number_of_Chebyshev_points);

//    ::CROSP::rod_properties::CircularCrossSection cs;
//    const double length = 1.0;
//    ::CROSP::rod_properties::RodDimensions rod_dimensions(&cs, length);

//    //  Now use it in rod properties
//    auto rod_properties =
//            std::make_shared<::CROSP::rod_properties::RodProperties>(strain_parameterisation,
//                                                                     rod_dimensions,
//                                                                     ::CROSP::rod_properties::MaterialProperties());




//    //  Variables related the perturbation of the strain parameterisation
//    std::shared_ptr<::CROSP::strain_parameterisation::StrainParameterisation> strain_parameterisation_Delta {
//        std::make_shared<::CROSP::strain_parameterisation::StrainParameterisation>(strain_parameterisation,
//                                                                          ::LieAlgebra::Vector6d::Zero())
//    };

//    ::CROSP::numerical_integrators::spectral_method::IntegratorsSPtr m_cosserat_rod_integrators {
//      std::make_shared<::CROSP::numerical_integrators::spectral_method::SpectralIntegrators>(strain_parameterisation,
//                                                                                             strain_parameterisation_Delta,
//                                                                                             rod_properties)
//    };

//    m_cosserat_rod_integrators->updateParameterisation(q, dot_q, ddot_q);

//    m_cosserat_rod_integrators->forwardKinematics();


//    const Eigen::MatrixXd Q_spectral = m_cosserat_rod_integrators->m_idm_integrators->m_quaternion->getStackAsMatrix();
//    const Eigen::MatrixXd R_magnus   = integrate_rotation_matrix->getStackAsMatrix();
//    const Eigen::MatrixXd R_SO3      = integrate_SO3->getStackAsMatrix();
//    Eigen::MatrixXd R_sidebyside(R_magnus.rows(), 2*R_magnus.cols());
//    R_sidebyside << R_magnus.transpose(), R_SO3.transpose();

//    std::cout << "error orientation : \n" << R_sidebyside << "\n\n";

    m_strain_parameterisation_stack->updateStrainParameterisation(q, dot_q, ddot_q);


    integrate_rotation_matrix->computesCoefficientsMatricesAtQuadraturePoints();
    integrate_rotation_matrix->integrate(R_X0);

    integrate_SO3->computesCoefficientsMatricesAtQuadraturePoints();
    integrate_SO3->integrate(R_X0);

    for(unsigned int i=0; i<number_of_Chebyshev_points; i++){
        std::cout << "error orientation : \n" << integrate_rotation_matrix->getStateAtPoint(i) << "\n" << integrate_SO3->getStateAtPoint(i) << "\n\n";
    }


//    Eigen::MatrixXd Q_magnus = Eigen::MatrixXd::Zero(4, R_magnus.cols()/3);
//    for(unsigned int i=0; i<Q_magnus.cols(); i++){
//        Eigen::Quaterniond q(R_magnus.block<3, 3>(0,i*3));
//        Q_magnus.col(i) << q.w(), q.x(), q.y(), q.z();
//    }
//    std::cout << "Q_spectral:\n" << Q_spectral << "\n\n";
//    std::cout << "Q_magnus:\n" << Q_magnus << "\n\n";
//    std::cout << "error orientation : \n" << (Q_spectral - Q_magnus).norm() << "\n\n";


//    const Eigen::MatrixXd r_spectral = m_cosserat_rod_integrators->m_idm_integrators->m_position->getStackAsMatrix();
//    const Eigen::MatrixXd r_magnus   = integrate_position->getStackAsMatrix();
//    std::cout << "r_spectral:\n" << r_spectral << "\n\n";
//    std::cout << "r_magnus:\n" << r_magnus << "\n\n";
//    std::cout << "error position : \n" << (r_spectral - r_magnus).norm() << "\n\n";
//    std::cout.flush();


//    const Eigen::MatrixXd Omega_spectral = m_cosserat_rod_integrators->m_idm_integrators->m_angular_velocity->getStackAsMatrix();
//    const Eigen::MatrixXd Omega_magnus   = Omega_integrator->getStackAsMatrix();
//    std::cout << "error angular velocity : \n" << (Omega_spectral - Omega_magnus).norm() << "\n\n";
//    std::cout.flush();


//    const Eigen::MatrixXd V_spectral = m_cosserat_rod_integrators->m_idm_integrators->m_linear_velocity->getStackAsMatrix();
//    const Eigen::MatrixXd V_magnus   = V_integrator->getStackAsMatrix();
//    std::cout << "error linear velocity : \n" << (V_spectral - V_magnus).norm() << "\n\n";
//    std::cout.flush();


//    const Eigen::MatrixXd dot_Omega_spectral = m_cosserat_rod_integrators->m_idm_integrators->m_angular_acceleration->getStackAsMatrix();
//    const Eigen::MatrixXd dot_Omega_magnus   = dot_Omega_integrator->getStackAsMatrix();
//    std::cout << "error angular acceleration : \n" << (dot_Omega_spectral - dot_Omega_magnus).norm() << "\n\n";
//    std::cout.flush();


//    const Eigen::MatrixXd dot_V_spectral = m_cosserat_rod_integrators->m_idm_integrators->m_linear_acceleration->getStackAsMatrix();
//    const Eigen::MatrixXd dot_V_magnus   = dot_V_integrator->getStackAsMatrix();
//    std::cout << "error linear acceleration : \n" << (dot_V_spectral - dot_V_magnus).norm() << "\n\n";
//    std::cout.flush();

//    ::benchmark::Initialize(&argc, argv);

//    ::benchmark::RunSpecifiedBenchmarks();


    return 0;
}
