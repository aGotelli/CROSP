#include <iostream>

#include "CROSP/CROSP/cosserat_rod.hpp"

#include <benchmark/benchmark.h>


#include <math_tools/Chebyshev/chebyshev_differentiation.hpp>

#include <utilities/Eigen/eigen_io.hpp>



void toMatrixForward(const Eigen::VectorXd &t_stack, Eigen::MatrixXd &t_matrix, const unsigned int t_state_dim, const unsigned int t_Nc)
{
    for(unsigned int it=0; it<t_state_dim; it++){
        unsigned int start = it*(t_Nc-1);

        for(unsigned int c=1; c<t_Nc; c++){

            unsigned int point = start + c - 1;
            t_matrix(it, t_Nc - c) = t_stack(point);
        }
    }
}


void toMatrixBackward(const Eigen::VectorXd &t_stack, Eigen::MatrixXd &t_matrix, const unsigned int t_state_dim, const unsigned int t_Nc)
{

    // Eigen::MatrixXd numbered_stack(t_stack.size(), 2);
    // for(unsigned int it=0; it<t_stack.size(); it++){
    //     numbered_stack(it, 0) = it;
    //     numbered_stack(it, 1) = t_stack(it);
    // }

    // std::cout << "numbered_stack: \n" << numbered_stack << std::endl;

    for(unsigned int it=0; it<t_state_dim; it++){
        unsigned int start = it*(t_Nc-1);
        // std::cout << "start: " << start << std::endl;
        for(unsigned int c=0; c<t_Nc-1; c++){

            unsigned int point = start + c;
            // std::cout << "\tpoint: " << point << std::endl;
            t_matrix(it, t_Nc - 2 - c) = t_stack(point);
        }
    }
}





Eigen::Vector3d skew_bracket(const Eigen::Vector3d &t_u,const Eigen::Vector3d &t_v)
{
    return Eigen::Vector3d(t_u.y()*t_v.z() - t_v.y()*t_u.z(),
                           t_v.x()*t_u.z() - t_u.x()*t_v.z(),
                           t_u.x()*t_v.y() - t_v.x()*t_u.y());

}



/*!
 * \brief The INTEGRATION_DOMAIN enum define the domain for the numerical integration.
 *
 * In this enumeration, we refer with BEGIN the lower limit of the domain, typically X=0, and
 * with END the upper limit of the domain, typically X=1.
 *
 * This enum is defined for being used to access directly and more clearly the begin and the end of the
 * integration domain for reading the corresponding value of the state.
 */
enum class INTEGRATION_DOMAIN {
    BEGIN,
    END
};

/*!
 * \brief The MagnusIntegrator class implements a generic templated object for the numerical integration via magnus Expansion
 * \tparam state_rows is the compiled time known number of rows of the state vector
 * \tparam state_cols is the compiled time known number of columns of the state vector
 */
class MagnusIntegratorSO3 /*: public ODESolverInterface*/{

public:



    //virtual ~MagnusIntegratorSO3()=default;

    /*!
     * \brief MagnusIntegrator constructs the object allocating all the needed memory
     * \param t_number_of_Chebyshev_points the number of Chebyshev points for the numerical integration
     */
    MagnusIntegratorSO3(const unsigned int t_number_of_Chebyshev_points)
        : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
    {

        //  Define the Chenyshev points we use
        const auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(m_number_of_Chebyshev_points);

        //  Some loop variables
        double begin, end;
        double h, step;

        double x1, x2, x3;

        //  Define quadrature points between every Chebyshev point
        for(unsigned int i=0; i<m_number_of_Chebyshev_points-1; i++){

            //  Take the two points
            begin = Chebyshev_points[ (m_number_of_Chebyshev_points-1) - i];
            end = Chebyshev_points[ (m_number_of_Chebyshev_points-1) - (i + 1)];


            //  Compute and save the distance
            h = end - begin;
            m_h_stack.push_back(h);


            step = sqrt(15.0)/10.0;

            x1 = begin + (0.5 - step)*h;
            x2 = begin + 0.5*h;
            x3 = begin + (0.5 + step)*h;


            //  Save the quadrature points
            m_quadrature_points.push_back( x1 );
            m_quadrature_points.push_back( x2 );
            m_quadrature_points.push_back( x3 );

        }


    }


    virtual void computesCoefficientsMatricesAtQuadraturePoints(){};

    /*!
     * \brief solveSystem solves the system of linear ODEs
     *
     * This function implements the integration via magnus expansion as Explained in https://doi.org/10.1016/j.physrep.2008.11.001
     */
    void solveSystem()
    // {
    //     double h;
    //     for(unsigned int i=0; i<m_number_of_Chebyshev_points-1; i++){

    //         h = m_h_stack[i];

    //         //  Compute contants
    //         alpha1 = h*m_A2_stack[i];
    //         alpha2 = (m_A3_stack[i] - m_A1_stack[i])*sqrt(15.0)*h/3.0;
    //         alpha3 = (m_A3_stack[i] - 2*m_A2_stack[i] + m_A1_stack[i])*10.0*h/3.0;


    //         C1 = skew_bracket(alpha1, alpha2);
    //         C2 = 2*alpha3 + C1;
    //         C3 = C2.cross(alpha1)/60.0;
    //         C4 = -20*alpha1 - alpha3 + C1;
    //         C5 = alpha2 + C3;
    //         C6 = skew_bracket(C4, C5);

    //         Omega = alpha1
    //                 + alpha3/12.0
    //                 + C6/240.0;

    //         //            std::cout << "Omega : \n" << ::LieAlgebra::skew(Omega) << "\n\n";
    //         //            std::cout.flush();


    //         m_exp_Omega_stack[i] = exponentialOmega(Omega);

    //     }

    //     for(unsigned int i=0; i<m_number_of_Chebyshev_points-1; i++)
    //         m_states_stack[i+1] = recursiveExponentialMapping( m_states_stack[i], m_exp_Omega_stack[i]);

    //     // postProcessInterpolation();
    // }
    {
        double h;
        for(unsigned int i=0; i<m_number_of_Chebyshev_points-1; i++){

            h = m_h_stack[i];

            //  Compute contants
            alpha1 = h*m_A2_stack[i];
            alpha2 = (m_A3_stack[i] - m_A1_stack[i])*sqrt(15.0)*h/3.0;
            alpha3 = (m_A3_stack[i] - 2*m_A2_stack[i] + m_A1_stack[i])*10.0*h/3.0;


            C1 = alpha1*alpha2 - alpha2*alpha1;
            C2 = 2*alpha3 + C1;
            C3 = -( alpha1*C2 - C2*alpha1 )/60.0;
            C4 = -20*alpha1 - alpha3 + C1;
            C5 = alpha2 + C3;


            Omega = alpha1
                    + alpha3/12.0
                    + ( C4*C5  - C5*C4 )/240.0;

            //            std::cout << "old Omega : \n" << Omega << "\n\n";
            //            std::cout.flush();


            m_exp_Omega_stack[i] = exponentialOmega(::LieAlgebra::antiSkew(Omega) );

        }

        for(unsigned int i=0; i<m_number_of_Chebyshev_points-1; i++)
            m_states_stack[i+1] = recursiveExponentialMapping( m_states_stack[i], m_exp_Omega_stack[i]);

        // postProcessInterpolation();
    }


    /*!
     * \brief recursiveExponentialMapping this function defines the shape of the Homogeneous ODE we want to solve
     * \param t_previous_state the so far last known state
     * \param t_exp_Omega the exponential mapping to reach the next state
     * \return the next state
     *
     * This function is used to shape the homogenous ODE. It could be in two ways
     *
     * y' = Ay
     * or
     * y' = yA
     *
     *
     * By default is y' = Ay
     */
    Eigen::Matrix3d recursiveExponentialMapping(const Eigen::Matrix3d &t_previous_state,
                                                        const Eigen::Matrix3d &t_exp_Omega)
    {
        //        return t_exp_Omega * t_previous_state;
        return t_previous_state * t_exp_Omega;
    }


    /*!
     * \brief exponentialOmega computes the exponential of the matrix Omega
     * \param t_Omega the matrix Omega
     * \return the exponential of the matrix Omega
     *
     * This function defines how we want to compute the exponential of the matrix.
     *
     * By defalt  we use the Cayley expression of the exponential as it is faster.
     * This expression requires to compute the inverse of a small matrix.
     *
     * Depending on the size of the matrix, we use different Eigen routines :
     * -> If the size is at the most 4, the .inverse() routine is called which is faster on such small matrices
     * -> For bigger matrices, the LLT decomposition is used (Check that your matrix can be used in this case).
     *
     * More details of these routines are given here https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
     *
     * The differentiation happensa at compile time, with a constexpr if.
     *
     *
     */
    Eigen::Matrix3d exponentialOmega(const Eigen::Vector3d &t_Omega)
    {

        return ::LieAlgebra::expRodigues(t_Omega);

    }



    void updateAMatrix(std::shared_ptr<std::vector<Eigen::Vector3d>> m_K_stack)
    // {
    //     for(unsigned int it=0; it<m_number_of_Chebyshev_points-1; it++){
    //         m_A1_stack.push_back( m_K_stack->at(3*it + 0) );
    //         m_A2_stack.push_back( m_K_stack->at(3*it + 1) );
    //         m_A3_stack.push_back( m_K_stack->at(3*it + 2) );
    //     }
    // }
    {
        for(unsigned int it=0; it<m_number_of_Chebyshev_points-1; it++){
            m_A1_stack.push_back( ::LieAlgebra::skew( m_K_stack->at(3*it + 0) ) );
            m_A2_stack.push_back( ::LieAlgebra::skew( m_K_stack->at(3*it + 1) ) );
            m_A3_stack.push_back( ::LieAlgebra::skew( m_K_stack->at(3*it + 2) ) );
        }
    }


    // void postProcessInterpolation(){}

    Eigen::MatrixXd getStateAtQuadraturePoint(const unsigned int t_quadrature_point)const{}



    void integrate(const Eigen::MatrixXd &t_initial_state)
    {

        setInitialConditions( t_initial_state );

        solveSystem();
    }

    void setInitialConditions(const Eigen::MatrixXd &t_initial_state)
    {
        m_states_stack[0] = t_initial_state;
    }



    unsigned int getChebyshevPointsNumber()const
    {
        return m_number_of_Chebyshev_points;
    }



    Eigen::MatrixXd getStateAtPoint(const unsigned int t_point)const
    {
        return m_states_stack[t_point];
    }


    Eigen::MatrixXd getStateAtPoint(INTEGRATION_DOMAIN t_position)const
    {
        //  Check if its base
        if(t_position == INTEGRATION_DOMAIN::BEGIN)
            return getStateAtPoint(m_number_of_Chebyshev_points - 1);

        // If not its at tip
        return getStateAtPoint(0);
    }


    Eigen::MatrixXd getStack()const
    {
        return getStackAsMatrix();
    }


    Eigen::MatrixXd getStackAsMatrix()const
    {
        for(unsigned int Chebyshev_point=0; Chebyshev_point<m_number_of_Chebyshev_points; Chebyshev_point++){
            m_states_stack_as_Eigen_matrix.block<3, 3>(0, 3*Chebyshev_point) =
                m_states_stack[Chebyshev_point];
        }

        return m_states_stack_as_Eigen_matrix;
    }


    // void setUpperIntegrationDomain(const double &t_upper_integration_limit){};



protected:

    unsigned int m_number_of_Chebyshev_points;

    unsigned int m_number_of_quadrature_points { 3 };

    // std::vector<double> m_quadrature_points;


    // std::vector<Eigen::Vector3d> m_A1_stack;
    // std::vector<Eigen::Vector3d> m_A2_stack;
    // std::vector<Eigen::Vector3d> m_A3_stack;
    std::vector<Eigen::Matrix3d> m_A1_stack;
    std::vector<Eigen::Matrix3d> m_A2_stack;
    std::vector<Eigen::Matrix3d> m_A3_stack;


    //private:
public:


    std::vector<double> m_quadrature_points;


    // Eigen::Vector3d alpha1 { Eigen::Vector3d::Zero() };
    // Eigen::Vector3d alpha2 { Eigen::Vector3d::Zero() };
    // Eigen::Vector3d alpha3 { Eigen::Vector3d::Zero() };

    // Eigen::Vector3d C1 { Eigen::Vector3d::Zero() };
    // Eigen::Vector3d C2 { Eigen::Vector3d::Zero() };
    // Eigen::Vector3d C3 { Eigen::Vector3d::Zero() };
    // Eigen::Vector3d C4 { Eigen::Vector3d::Zero() };
    // Eigen::Vector3d C5 { Eigen::Vector3d::Zero() };
    // Eigen::Vector3d C6 { Eigen::Vector3d::Zero() };

    // Eigen::Vector3d Omega  { Eigen::Vector3d::Zero() };
    Eigen::Matrix3d alpha1 { Eigen::Matrix3d::Zero() };
    Eigen::Matrix3d alpha2 { Eigen::Matrix3d::Zero() };
    Eigen::Matrix3d alpha3 { Eigen::Matrix3d::Zero() };

    Eigen::Matrix3d C1 { Eigen::Matrix3d::Zero() };
    Eigen::Matrix3d C2 { Eigen::Matrix3d::Zero() };
    Eigen::Matrix3d C3 { Eigen::Matrix3d::Zero() };
    Eigen::Matrix3d C4 { Eigen::Matrix3d::Zero() };
    Eigen::Matrix3d C5 { Eigen::Matrix3d::Zero() };
    Eigen::Matrix3d C6 { Eigen::Matrix3d::Zero() };

    Eigen::Matrix3d Omega  { Eigen::Matrix3d::Zero() };


    std::vector<Eigen::Matrix3d> m_exp_Omega_stack { std::vector<Eigen::Matrix3d>(m_number_of_Chebyshev_points-1) };

    std::vector<Eigen::Matrix3d> m_states_stack { std::vector<Eigen::Matrix3d>(m_number_of_Chebyshev_points) };

    mutable Eigen::MatrixXd m_states_stack_as_Eigen_matrix { Eigen::MatrixXd::Zero(3, (m_number_of_Chebyshev_points)*3) };

    std::vector<double> m_h_stack;


};



void compareSolutionsQandR(const Eigen::MatrixXd &t_Q_stack, const Eigen::MatrixXd &t_R_stack, const unsigned int Nc)
{
    for(unsigned int it=0; it<Nc-1; ++it){
        Eigen::Vector4d Q = t_Q_stack.col(it);
        Eigen::Matrix3d R_Q = Eigen::Quaterniond(Q(0), Q(1), Q(2), Q(3)).toRotationMatrix();
        Eigen::Matrix3d R = t_R_stack.block<3, 3>(0, 3*it);

        std::cout << "R from quaternion : \n" << R_Q << "\nR : \n" << R <<"\n\n";
    }
}



int main(int argc, char *argv[])
{

    //  Rod
    ::CROSP::rod_properties::RodProperties rod_properties;

    //  Integration parameters
    const unsigned int Nc = 7;
    const unsigned int state_dim = 3;


    //  Integration points
    const auto forward_points  = ::Chebyshev::defineIntegrationPoints(Nc, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD);
    const auto backward_points = ::Chebyshev::defineIntegrationPoints(Nc, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD);

    //  VSA
    const std::array<bool, 6> admitted_deformations = {
        true,
        true,
        true,
        0,
        0,
        0
    };

    const unsigned int na = std::count(admitted_deformations.begin(), admitted_deformations.end(), true);
    const unsigned int modes = 5;
    const unsigned int ne = modes*na;

    ::CROSP::polynomial_representation::PolynomialRepresentation polynomial_representation(admitted_deformations, modes);

    ::CROSP::strain_parameterisation_stack::StrainParameterisationStack strain_parameterisation_stack(
        polynomial_representation,
        Nc
    );




    ::CROSP::CosseratRod<> rod(polynomial_representation, Nc, rod_properties.m_rod_dimensions, rod_properties.m_material_properties);



    //  Generalized elastic coordinates

    srand( 7 );
    Eigen::VectorXd q      = Eigen::VectorXd::Random(ne);
    Eigen::VectorXd dot_q  = Eigen::VectorXd::Random(ne)*10;
    Eigen::VectorXd ddot_q = Eigen::VectorXd::Random(ne)*100;

    //  Print values
    std::cout << "q : \n" << q << "\n\n";
    std::cout << "dot_q : \n" << dot_q << "\n\n";
    std::cout << "ddot_q : \n" << ddot_q << "\n\n";






    auto computeMatrixAtChebyshevPoint = [&](const unsigned int t_point) -> Eigen::Matrix3d
    {
        return - strain_parameterisation_stack.m_hat_K_stack->at(t_point);
    };


    auto computeMatrixAtChebyshevPointQuaternion = [&](const unsigned int t_point) -> Eigen::Matrix4d
    {
        Eigen::Matrix4d A;
        const auto m_K_stack = strain_parameterisation_stack.m_K_stack;

        A  <<              0             ,   -m_K_stack->at(t_point)(0),   -m_K_stack->at(t_point)(1),   -m_K_stack->at(t_point)(2),
            m_K_stack->at(t_point)(0) ,                0            ,    m_K_stack->at(t_point)(2),   -m_K_stack->at(t_point)(1),
            m_K_stack->at(t_point)(1) ,   -m_K_stack->at(t_point)(2),                0            ,    m_K_stack->at(t_point)(0),
            m_K_stack->at(t_point)(2) ,    m_K_stack->at(t_point)(1),   -m_K_stack->at(t_point)(0),                0            ;


        return 0.5*A;
    };





    MagnusIntegratorSO3 integratorR(Nc);
    ::CROSP::strain_parameterisation_stack::StrainParameterisationStack strain_stack_magnus(
        polynomial_representation,
        integratorR.m_quadrature_points
    );








    //  Differentiation matrices
    const Eigen::MatrixXd Df = ::Chebyshev::getD_NN(Nc, state_dim, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD);
    const Eigen::MatrixXd Df_quat = ::Chebyshev::getD_NN(Nc, 4, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD);
    const Eigen::MatrixXd Db = ::Chebyshev::getD_NN(Nc, state_dim, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD);

    //  Initial conditions
    const Eigen::Matrix3d R_X0 = Eigen::Matrix3d::Identity();
    const Eigen::Vector4d Q_X0 = Eigen::Vector4d(1, 0, 0, 0);
    const Eigen::Vector3d r_X0 = Eigen::Vector3d::Zero();

    const Eigen::Vector3d Omega_X0 = Eigen::Vector3d::Zero();//Eigen::Vector3d(.5, .25, -.33);
    const Eigen::Vector3d V_X0 = Eigen::Vector3d::Zero();

    const Eigen::Vector3d dot_Omega_X0 = Eigen::Vector3d::Zero();
    const Eigen::Vector3d dot_V_X0 = Eigen::Vector3d::Zero();

    const Eigen::Vector3d N_X0 = Eigen::Vector3d::Zero();
    const Eigen::Vector3d C_X0 = Eigen::Vector3d::Zero();


    //  Influence of the initial conditions
    Eigen::MatrixXd DIf = ::Chebyshev::getD_IN(Nc, state_dim, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD);
    Eigen::MatrixXd DIf_quat = ::Chebyshev::getD_IN(Nc, 4, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD);
    Eigen::MatrixXd DIb = ::Chebyshev::getD_IN(Nc, state_dim, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD);



    //  Solution IVP
    const auto ivp_Quaternion = DIf_quat*Q_X0;
    const auto ivp_r = DIf*r_X0;

    const auto ivp_Omega = DIf*Omega_X0;
    const auto ivp_V = DIf*V_X0;

    const auto ivp_dot_Omega = DIf*dot_Omega_X0;
    const auto ivp_dot_V = DIf*dot_V_X0;

    const auto ivp_N = DIf*N_X0;
    const auto ivp_C = DIf*C_X0;


    //  Kinematics stacks
    Eigen::VectorXd Quaternion_stack = Eigen::VectorXd::Zero((Nc-1)*4);

    Eigen::VectorXd Omega_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);
    Eigen::VectorXd V_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);

    Eigen::VectorXd dot_Omega_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);
    Eigen::VectorXd dot_V_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);

    Eigen::VectorXd N_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);
    Eigen::VectorXd C_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);




    //  Coefficient matrices
    Eigen::MatrixXd Cf = Df;
    Eigen::MatrixXd Cf_quat = Df_quat;
    Eigen::MatrixXd Cb = Db;


    //  Inversions
    Eigen::PartialPivLU<Eigen::MatrixXd> Cf_inv(Cf);
    Eigen::PartialPivLU<Eigen::MatrixXd> Cb_inv(Cb);
    // Eigen::HouseholderQR<Eigen::MatrixXd> Cf_inv(Cf);
    // Eigen::HouseholderQR<Eigen::MatrixXd> Cb_inv(Cb);


    Eigen::MatrixXd Quaternion_matrix(4, Nc);
    Quaternion_matrix.col(0) << Q_X0;

    Eigen::MatrixXd Omega_optimized(3, Nc);
    Omega_optimized.col(0) << Omega_X0;
    Eigen::MatrixXd V_optimized(3, Nc);
    V_optimized.col(0) << V_X0;

    Eigen::MatrixXd dot_Omega_optimized(3, Nc);
    dot_Omega_optimized.col(0) << dot_Omega_X0;
    Eigen::MatrixXd dot_V_optimized(3, Nc);
    dot_V_optimized.col(0) << dot_V_X0;

    Eigen::MatrixXd N_optimized(3, Nc);
    N_optimized.col(Nc-1) << N_X0;
    Eigen::MatrixXd C_optimized(3, Nc);
    C_optimized.col(Nc-1) << C_X0;





    auto solve_IDM = [&](const bool optimized=true){



        strain_parameterisation_stack.updateStrainParameterisation(q, dot_q, ddot_q);

        strain_stack_magnus.updateStrainParameterisation(q, dot_q, ddot_q);
        integratorR.updateAMatrix(strain_stack_magnus.m_K_stack);



        //  Compute coefficient matrix forward integration
        for(int index = 0; const auto it : forward_points) {

            Eigen::Matrix3d A = computeMatrixAtChebyshevPoint(it);




            for (unsigned int row = 0; row < state_dim; ++row) {
                for (unsigned int col = 0; col < state_dim; ++col) {
                    int row_index = row*(Nc-1) + index;
                    int col_index = col*(Nc-1) + index;
                    Cf(row_index, col_index) = Df(row_index, col_index) - A(row, col);
                }
            }

            index++;

        }

        //  Compute coefficient matrix backward integration
        for(int index = 0; const auto it : backward_points) {

            Eigen::Matrix3d A = computeMatrixAtChebyshevPoint(it);


            for (unsigned int row = 0; row < state_dim; ++row) {
                for (unsigned int col = 0; col < state_dim; ++col) {
                    int row_index = row*(Nc-1) + index;
                    int col_index = col*(Nc-1) + index;
                    Cb(row_index, col_index) = Db(row_index, col_index) - A(row, col);
                }
            }

            index++;

        }





        if(optimized){
            integratorR.integrate(R_X0);
        } else {

            //  Compute coefficient matrix forward integration
            for(int index = 0; const auto it : forward_points) {

                Eigen::Matrix4d A = computeMatrixAtChebyshevPointQuaternion(it);


                for (unsigned int row = 0; row < 4; ++row) {
                    for (unsigned int col = 0; col < 4; ++col) {
                        int row_index = row*(Nc-1) + index;
                        int col_index = col*(Nc-1) + index;
                        Cf_quat(row_index, col_index) = Df_quat(row_index, col_index) - A(row, col);
                    }
                }

                index++;

            }


            Quaternion_stack = Cf_quat.partialPivLu().solve(-ivp_Quaternion);

        }





        if(optimized){
            //  Inversions
            Cf_inv.compute(Cf);
            Cb_inv.compute(Cb);
        }




        //                      Omega Integration
        //  Coefficients
        Eigen::Vector3d b;
        Eigen::VectorXd b_Omega = 0*ivp_Omega;
        for(int index = 0; const auto it : forward_points) {

            b = strain_parameterisation_stack.m_dot_K_stack->at(it);

            for (unsigned int row = 0; row < state_dim; ++row)
                b_Omega(index + row * (Nc-1)) = b(row);
            index++;
        }
        // Solve
        if(optimized)
            Omega_stack = Cf_inv.solve(b_Omega - ivp_Omega);
        else
            Omega_stack = Cf.partialPivLu().solve(b_Omega - ivp_Omega);

        toMatrixForward(Omega_stack, Omega_optimized, state_dim, Nc);






        //                      V Integration
        //  Coefficients
        Eigen::VectorXd b_V = 0*ivp_V;
        Eigen::Vector3d dot_Gamma;
        Eigen::Vector3d Omega;
        for(int index = 0; const auto it : forward_points) {

            dot_Gamma = strain_parameterisation_stack.m_dot_Gamma_stack->at(it);

            for(int i=0; i<state_dim; i++)
                Omega(i) = Omega_stack(it + i*(Nc-1));

            b = dot_Gamma - strain_parameterisation_stack.m_hat_Gamma_stack->at(it)*Omega;

            for (unsigned int row = 0; row < state_dim; ++row)
                b_V(index + row * (Nc-1)) = b(row);
            index++;
        }
        //  Solve
        if(optimized)
            V_stack = Cf_inv.solve(b_V - ivp_V);
        else
            V_stack = Cf.partialPivLu().solve(b_V - ivp_V);

        toMatrixForward(V_stack, V_optimized, state_dim, Nc);


        //                      dot_Omega Integration
        //  Coefficients
        Eigen::VectorXd b_dot_Omega = 0*ivp_dot_Omega;
        Eigen::Vector3d dot_Gamma_dot;
        for(int index = 0; const auto it : forward_points) {

            for(int i=0; i<state_dim; i++)
                Omega(i) = Omega_stack(it + i*(Nc-1));

            b = strain_parameterisation_stack.m_ddot_K_stack->at(it)
                - strain_parameterisation_stack.m_hat_dot_K_stack->at(it)*Omega;

            for (unsigned int row = 0; row < state_dim; ++row)
                b_dot_Omega(index + row * (Nc-1)) = b(row);
            index++;
        }
        //  Solve
        if(optimized)
            dot_Omega_stack = Cf_inv.solve(b_dot_Omega - ivp_dot_Omega);
        else
            dot_Omega_stack = Cf.partialPivLu().solve(b_dot_Omega - ivp_dot_Omega);

        toMatrixForward(dot_Omega_stack, dot_Omega_optimized, state_dim, Nc);


        //                      dot_V Integration
        //  Coefficients
        Eigen::VectorXd b_dot_V = 0*ivp_dot_V;
        Eigen::Vector3d dot_Omega;
        Eigen::Vector3d V;
        for(int index = 0; const auto it : forward_points) {

            for(int i=0; i<state_dim; i++)
                Omega(i) = Omega_stack(it + i*(Nc-1));

            for(int i=0; i<state_dim; i++)
                dot_Omega(i) = dot_Omega_stack(it + i*(Nc-1));

            for(int i=0; i<state_dim; i++)
                V(i) = V_stack(it + i*(Nc-1));



            b = strain_parameterisation_stack.m_ddot_Gamma_stack->at(it)
                - strain_parameterisation_stack.m_hat_Gamma_stack->at(it) * dot_Omega
                - strain_parameterisation_stack.m_hat_dot_Gamma_stack->at(it) * Omega
                - strain_parameterisation_stack.m_hat_dot_K_stack->at(it) * V;

            for (unsigned int row = 0; row < state_dim; ++row)
                b_dot_V(index + row * (Nc-1)) = b(row);
            index++;
        }
        //  Solve
        if(optimized)
            dot_V_stack = Cf_inv.solve(b_dot_V - ivp_dot_V);
        else
            dot_V_stack = Cf.partialPivLu().solve(b_dot_V - ivp_dot_V);

        toMatrixForward(dot_V_stack, dot_V_optimized, state_dim, Nc);

        // std::cout << "dot_V_optimized = [\n " << dot_V_optimized << "\n]" << std::endl;

        //                      N Integration
        //  Coefficients
        Eigen::VectorXd b_N = 0*ivp_N;
        const auto M_linear = rod_properties.getMLinear();
        Eigen::Vector3d dot_V;
        for(int index = 0; const auto it : backward_points) {


            // for(int i=0; i<state_dim; i++)
            //     Omega(i) = Omega_stack(it-1 + i*(Nc-1));

            // for(int i=0; i<state_dim; i++)
            //     dot_Omega(i) = dot_Omega_stack(it-1 + i*(Nc-1));

            // for(int i=0; i<state_dim; i++)
            //     V(i) = V_stack(it-1 + i*(Nc-1));

            // for(int i=0; i<state_dim; i++)
            //     dot_V(i) = dot_V_stack(it + i*(Nc-1));
            Omega     = Omega_optimized.col(Nc - 1 - it);
            dot_Omega = dot_Omega_optimized.col(Nc - 1 - it);
            V         = V_optimized.col(Nc - 1 - it);
            dot_V     = dot_V_optimized.col(Nc - 1 - it);


            // std::cout << "At point : " << it << "\n" <<
            //     "\t dot_V (sec) : \n" << dot_V << "\n" << std::endl;
            //     // "\t inertial_velocities : \n" << ::LieAlgebra::skew( Omega ).transpose() * M_linear * V << std::endl;


            b = M_linear*dot_V - ::LieAlgebra::skew( Omega ).transpose() * M_linear * V;

            for (unsigned int row = 0; row < state_dim; ++row)
                b_N(index + row * (Nc-1)) = b(row);
            index++;
        }

        // std::cout << "b_N = [\n " << b_N << "\n]" << std::endl;
        // std::cout << "ivp_N = [\n " << ivp_N << "\n]" << std::endl;
        // std::cout << "Cb = [\n " << Cb << "\n]" << std::endl;
        //  Solve
        if(optimized)
            N_stack = Cb_inv.solve(b_N - ivp_N);
        else
            N_stack = Cb.partialPivLu().solve(b_N - ivp_N);

        toMatrixBackward(N_stack, N_optimized, state_dim, Nc);

        //                      C Integration
        //  Coefficients
        Eigen::VectorXd b_C = 0*ivp_C;
        const auto M_angular = rod_properties.getMAngular();
        Eigen::Vector3d N;
        Eigen::Vector3d Gamma;
        for(int index = 0; const auto it : backward_points) {

            Gamma = strain_parameterisation_stack.m_Gamma_stack->at(it);
            dot_Gamma = strain_parameterisation_stack.m_dot_Gamma_stack->at(it);

            // for(int i=0; i<state_dim; i++)
            //     Omega(i) = Omega_stack(it-1 + i*(Nc-1));

            // for(int i=0; i<state_dim; i++)
            //     dot_Omega(i) = dot_Omega_stack(it-1 + i*(Nc-1));

            // for(int i=0; i<state_dim; i++)
            //     V(i) = V_stack(it-1 + i*(Nc-1));

            // for(int i=0; i<state_dim; i++)
            //     dot_V(i) = dot_V_stack(it-1 + i*(Nc-1));

            for(int i=0; i<state_dim; i++)
                N(i) = N_stack(it-1 + i*(Nc-1));


            Omega     = Omega_optimized.col(Nc - 1 - it);
            dot_Omega = dot_Omega_optimized.col(Nc - 1 - it);
            V         = V_optimized.col(Nc - 1 - it);
            dot_V     = dot_V_optimized.col(Nc - 1 - it);


            b = ::LieAlgebra::skew( Gamma ).transpose()*N
                + M_angular*dot_Omega
                - ::LieAlgebra::skew( Omega).transpose() * M_angular * Omega
                - ::LieAlgebra::skew( V ).transpose()* M_linear * V;

            for (unsigned int row = 0; row < state_dim; ++row)
                b_C(index + row * (Nc-1)) = b(row);

            index++;
        }
        //  Solve
        if(optimized)
            C_stack = Cb_inv.solve(b_C - ivp_C);
        else
            C_stack = Cb.partialPivLu().solve(b_C - ivp_C);

        toMatrixBackward(C_stack, C_optimized, state_dim, Nc);
    };


    solve_IDM();


    solve_IDM(false);
    toMatrixForward(Quaternion_stack, Quaternion_matrix, 4, Nc);


    // //  Compare rotation to quaternion
    // compareSolutionsQandR(Quaternion_matrix, integratorR.getStack(), Nc);


    rod.m_rod_properties->m_gravity.setZero();
    rod.updateParameterisation(q, dot_q, ddot_q);
    rod.forwardKinematics(Eigen::Vector4d::UnitX(), r_X0,
                            Omega_X0, V_X0,
                          dot_Omega_X0, dot_V_X0);
    // rod.backwardDynamics(::LieAlgebra::Vector6d::Zero());

    const auto Omega_reference = rod.m_cosserat_rod_integrators->m_idm_integrators->m_angular_velocity->getStackAsMatrix();
    // const auto Q_reference = rod.m_cosserat_rod_integrators->m_idm_integrators->m_quaternion->getStackAsMatrix();


    std::cout << "Omega_reference : \n" << Omega_reference << "\n\n\n";


    return -1;

    //  Save data
    const std::string path = "benchmarks/data/optimized_IDM";


    //  Save data as .csv file
    ::writeToFile("q", q, path);
    ::writeToFile("dot_q", dot_q, path);
    ::writeToFile("ddot_q", ddot_q, path);


    ::writeToFile("Quaternion", Quaternion_matrix, path);
    ::writeToFile("Omega", Omega_optimized, path);
    ::writeToFile("V", V_optimized, path);
    ::writeToFile("dot_Omega", dot_Omega_optimized, path);
    ::writeToFile("dot_V", dot_V_optimized, path);
    ::writeToFile("N", N_optimized, path);
    ::writeToFile("C", C_optimized, path);


    const unsigned int repetitions = 10;


    ::benchmark::RegisterBenchmark("Solution IDM", [&](::benchmark::State &t_state){
        const bool optimized = false;
        while(t_state.KeepRunning())
            solve_IDM(optimized);
    })->Unit(::benchmark::kMicrosecond)->Repetitions(repetitions);


    ::benchmark::RegisterBenchmark("Solution IDM optimized", [&](::benchmark::State &t_state){
        while(t_state.KeepRunning())
            solve_IDM();
    })->Unit(::benchmark::kMicrosecond)->Repetitions(repetitions);



    ::CROSP::strain_parameterisation_stack::StrainParameterisationStack Delta_strain_parameterisation_stack(
        polynomial_representation,
        Nc,
        ::CROSP::strain_parameterisation_stack::zero_constrained_strain
        );






    Eigen::VectorXd Delta_Pi_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);
    Eigen::VectorXd Delta_D_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);

    Eigen::VectorXd Delta_Omega_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);
    Eigen::VectorXd Delta_V_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);

    Eigen::VectorXd Delta_dot_Omega_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);
    Eigen::VectorXd Delta_dot_V_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);

    Eigen::MatrixXd Delta_Pi = Eigen::MatrixXd::Zero(state_dim, Nc);
    Eigen::MatrixXd Delta_D = Eigen::MatrixXd::Zero(state_dim, Nc);

    Eigen::MatrixXd Delta_Omega = Eigen::MatrixXd::Zero(state_dim, Nc);
    Eigen::MatrixXd Delta_V = Eigen::MatrixXd::Zero(state_dim, Nc);

    Eigen::MatrixXd Delta_dot_Omega = Eigen::MatrixXd::Zero(state_dim, Nc);
    Eigen::MatrixXd Delta_dot_V = Eigen::MatrixXd::Zero(state_dim, Nc);


    Eigen::VectorXd Delta_N_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);
    Eigen::VectorXd Delta_C_stack = Eigen::VectorXd::Zero((Nc-1)*state_dim);

    Eigen::MatrixXd Delta_N = Eigen::MatrixXd::Zero(state_dim, Nc);
    Eigen::MatrixXd Delta_C = Eigen::MatrixXd::Zero(state_dim, Nc);








    Eigen::MatrixXd Delta_Pi_optimized = Eigen::MatrixXd::Zero((Nc-1)*state_dim, ne);
    Eigen::MatrixXd Delta_D_optimized  = Eigen::MatrixXd::Zero((Nc-1)*state_dim, ne);

    Eigen::MatrixXd Delta_Omega_optimized = Eigen::MatrixXd::Zero((Nc-1)*state_dim, ne);
    Eigen::MatrixXd Delta_V_optimized     = Eigen::MatrixXd::Zero((Nc-1)*state_dim, ne);

    Eigen::MatrixXd Delta_dot_Omega_optimized = Eigen::MatrixXd::Zero((Nc-1)*state_dim, ne);
    Eigen::MatrixXd Delta_dot_V_optimized     = Eigen::MatrixXd::Zero((Nc-1)*state_dim, ne);


    Eigen::MatrixXd Delta_N_optimized = Eigen::MatrixXd::Zero((Nc-1)*state_dim, ne);
    Eigen::MatrixXd Delta_C_optimized     = Eigen::MatrixXd::Zero((Nc-1)*state_dim, ne);


    Eigen::VectorXd Delta_q   = 0*q;
    Eigen::VectorXd Delta_dq  = 0*dot_q;
    Eigen::VectorXd Delta_ddq = 0*ddot_q;

    Eigen::VectorXd Delta = Eigen::VectorXd::Zero((Nc-1)*state_dim);
    Eigen::MatrixXd Delta_Matrix = Eigen::MatrixXd::Zero((Nc-1)*state_dim, ne);



    Eigen::MatrixXd b_Pi = Eigen::MatrixXd(3*(Nc-1), ne);
    for(unsigned int dne=0; dne<ne; ++dne){

        Delta_q.setZero();
        Delta_q(dne) = 1;
        Delta_dq     = 10*Delta_q;
        Delta_ddq    = 100*Delta_q;

        Delta_strain_parameterisation_stack.updateStrainParameterisation(Delta_q, Delta_dq, Delta_ddq);

        for(int index = 0; const auto it : forward_points) {

            const auto Delta_K = Delta_strain_parameterisation_stack.m_K_stack->at(it);
            const auto Delta_dot_K = Delta_strain_parameterisation_stack.m_K_stack->at(it);
            const auto Delta_ddot_K = Delta_strain_parameterisation_stack.m_K_stack->at(it);

            for (unsigned int row = 0; row < state_dim; ++row)
                b_Pi(index + row * (Nc-1)) = Delta_K(row);
            index++;
        }
    }

    auto solve_TIDM = [&](const bool optimized=true){


        const auto M_linear = rod_properties.getMLinear();
        const auto M_angular = rod_properties.getMAngular();
        if(optimized){


            // for(unsigned int it=0; it<ne; it++){

            //     Delta_q.setZero();
            //     Delta_q(it) = 1;
            //     Delta_dq    = 10*Delta_q;
            //     Delta_ddq   = 100*Delta_q;


            //     Delta_strain_parameterisation_stack.updateStrainParameterisation(Delta_q, Delta_dq, Delta_ddq);

            // }    No need to do this every time, we can just make it once at the beginning of the simulation

            Delta_Pi_optimized = Cf_inv.solve(Delta_Matrix);
            Delta_D_optimized = Cf_inv.solve(Delta_Matrix);

            Delta_Omega_optimized = Cf_inv.solve(Delta_Matrix);
            Delta_V_optimized = Cf_inv.solve(Delta_Matrix);


            Delta_dot_Omega_optimized = Cf_inv.solve(Delta_Matrix);
            Delta_dot_V_optimized = Cf_inv.solve(Delta_Matrix);


            Delta_N_optimized = Cb_inv.solve(Delta_Matrix);
            Delta_C_optimized = Cb_inv.solve(Delta_Matrix);



        } else {

            for(unsigned int it=0; it<ne; it++){


                // Delta_q.setZero();
                // Delta_q(it) = 1;
                // Delta_dq    = 10*Delta_q;
                // Delta_ddq   = 100*Delta_q;


                // Delta_strain_parameterisation_stack.updateStrainParameterisation(Delta_q, Delta_dq, Delta_ddq);
                //  No need to do this every time, we can just make it once at the beginning of the simulation




                //                      Delta Pi Integration
                //  Coefficients
                Eigen::Vector3d b;
                Eigen::VectorXd b_stack = 0*ivp_Omega;
                for(int index = 0; const auto it : forward_points) {

                    b = Delta_strain_parameterisation_stack.m_K_stack->at(it);

                    for (unsigned int row = 0; row < state_dim; ++row)
                        b_stack(index + row * (Nc-1)) = b(row);
                    index++;
                }
                // Solve
                Delta_Pi_stack = Cf.partialPivLu().solve(b_stack);
                toMatrixForward(Delta_Pi_stack, Delta_Pi, state_dim, Nc);



                //                      Delta D Integration
                //  Coefficients
                for(int index = 0; const auto it : forward_points) {


                    const auto Delta_Gamma = Delta_strain_parameterisation_stack.m_Gamma_stack->at(it);
                    const auto Gamma = strain_parameterisation_stack.m_Gamma_stack->at(it);
                    const auto Delta_Pi_ = Delta_Pi.col(it);

                    b = - ::LieAlgebra::skew( Gamma )*Delta_Pi_ + Delta_Gamma;

                    for (unsigned int row = 0; row < state_dim; ++row)
                        b_stack(index + row * (Nc-1)) = b(row);
                    index++;
                }
                // Solve
                Delta_D_stack = Cf.partialPivLu().solve(b_stack);
                toMatrixForward(Delta_D_stack, Delta_D, state_dim, Nc);



                //                      Delta Omega Integration
                //  Coefficients
                for(int index = 0; const auto it : forward_points) {


                    const auto Delta_K = Delta_strain_parameterisation_stack.m_K_stack->at(it);
                    const auto Delta_dot_K = strain_parameterisation_stack.m_dot_K_stack->at(it);
                    const auto Omega = Omega_optimized.col(it);

                    b = - ::LieAlgebra::skew( Delta_K )*Omega + Delta_dot_K;

                    for (unsigned int row = 0; row < state_dim; ++row)
                        b_stack(index + row * (Nc-1)) = b(row);
                    index++;
                }
                // Solve
                Delta_Omega_stack = Cf.partialPivLu().solve(b_stack);
                toMatrixForward(Delta_Omega_stack, Delta_Omega, state_dim, Nc);




                //                      Delta V Integration
                //  Coefficients
                for(int index = 0; const auto it : forward_points) {


                    const auto Delta_K          = Delta_strain_parameterisation_stack.m_K_stack->at(it);
                    const auto Delta_dot_K      = strain_parameterisation_stack.m_dot_K_stack->at(it);
                    const auto Gamma            = strain_parameterisation_stack.m_Gamma_stack->at(it);
                    const auto Delta_Gamma      = Delta_strain_parameterisation_stack.m_Gamma_stack->at(it);
                    const auto Delta_dot_Gamma  = Delta_strain_parameterisation_stack.m_dot_Gamma_stack->at(it);

                    const auto Omega = Omega_optimized.col(it);
                    const auto V = V_optimized.col(it);

                    const auto Delta_Omega_ = Delta_Omega.col(it);

                    b = - ::LieAlgebra::skew( Gamma )*Delta_Omega_ - ::LieAlgebra::skew( Delta_Gamma )*Omega
                        - ::LieAlgebra::skew( Delta_K )*V + Delta_dot_Gamma;

                    for (unsigned int row = 0; row < state_dim; ++row)
                        b_stack(index + row * (Nc-1)) = b(row);
                    index++;
                }
                // Solve
                Delta_V_stack = Cf.partialPivLu().solve(b_stack);
                toMatrixForward(Delta_V_stack, Delta_V, state_dim, Nc);





                //                      Delta dot Omega Integration
                //  Coefficients
                for(int index = 0; const auto it : forward_points) {

                    const auto dot_K  = strain_parameterisation_stack.m_dot_K_stack->at(it);

                    const auto Delta_K      = Delta_strain_parameterisation_stack.m_K_stack->at(it);
                    const auto Delta_dot_K  = Delta_strain_parameterisation_stack.m_dot_K_stack->at(it);
                    const auto Delta_ddot_K = Delta_strain_parameterisation_stack.m_ddot_K_stack->at(it);

                    const auto Omega = Omega_optimized.col(it);
                    const auto dot_Omega = dot_Omega_optimized.col(it);
                    const auto Delta_Omega_ = Delta_Omega.col(it);

                    b = - ::LieAlgebra::skew( Delta_K )*dot_Omega
                        - ::LieAlgebra::skew( Delta_dot_K )*Omega
                        - ::LieAlgebra::skew( dot_K )*Delta_Omega_
                        - ::LieAlgebra::skew( Delta_dot_K )*Omega
                        + Delta_ddot_K;

                    for (unsigned int row = 0; row < state_dim; ++row)
                        b_stack(index + row * (Nc-1)) = b(row);
                    index++;
                }
                // Solve
                Delta_dot_Omega_stack = Cf.partialPivLu().solve(b_stack);
                toMatrixForward(Delta_dot_Omega_stack, Delta_dot_Omega, state_dim, Nc);







                //                      Delta dot V Integration
                //  Coefficients
                for(int index = 0; const auto it : forward_points) {



                    const auto K     = strain_parameterisation_stack.m_K_stack->at(it);
                    const auto dot_K = strain_parameterisation_stack.m_dot_K_stack->at(it);


                    const auto Delta_K          = Delta_strain_parameterisation_stack.m_K_stack->at(it);
                    const auto Delta_dot_K      = Delta_strain_parameterisation_stack.m_dot_K_stack->at(it);

                    const auto Gamma            = strain_parameterisation_stack.m_Gamma_stack->at(it);
                    const auto dot_Gamma        = strain_parameterisation_stack.m_dot_Gamma_stack->at(it);
                    const auto Delta_Gamma      = Delta_strain_parameterisation_stack.m_Gamma_stack->at(it);
                    const auto Delta_dot_Gamma  = Delta_strain_parameterisation_stack.m_dot_Gamma_stack->at(it);
                    const auto Delta_ddot_Gamma = Delta_strain_parameterisation_stack.m_ddot_Gamma_stack->at(it);

                    const auto Omega = Omega_optimized.col(it);
                    const auto V = V_optimized.col(it);

                    const auto dot_Omega = dot_Omega_optimized.col(it);
                    const auto dot_V     = dot_V_optimized.col(it);


                    const auto Delta_Omega_ = Delta_Omega.col(it);
                    const auto Delta_V_     = Delta_V.col(it);

                    const auto Delta_dot_Omega_ = Delta_dot_Omega.col(it);



                    b = - ::LieAlgebra::skew( Delta_K )*dot_V
                        - ::LieAlgebra::skew( Delta_Gamma )*dot_Omega
                        - ::LieAlgebra::skew( Gamma )*Delta_dot_Omega_
                        - ::LieAlgebra::skew( Delta_dot_Gamma )*Omega
                        - ::LieAlgebra::skew( dot_Gamma )*Delta_Omega_
                        - ::LieAlgebra::skew( Delta_dot_K )*V
                        - ::LieAlgebra::skew( dot_K )*Delta_V_
                        + Delta_ddot_Gamma;

                    for (unsigned int row = 0; row < state_dim; ++row)
                        b_stack(index + row * (Nc-1)) = b(row);
                    index++;
                }
                // Solve
                Delta_dot_V_stack = Cf.partialPivLu().solve(b_stack);
                toMatrixForward(Delta_dot_V_stack, Delta_dot_V, state_dim, Nc);










                //                      Delta N Integration
                //  Coefficients
                for(int index = 0; const auto it : forward_points) {



                    // const auto K     = strain_parameterisation_stack.m_K_stack->at(it);
                    // const auto dot_K = strain_parameterisation_stack.m_dot_K_stack->at(it);


                    const auto Delta_K          = Delta_strain_parameterisation_stack.m_K_stack->at(it);
                    // const auto Delta_dot_K      = Delta_strain_parameterisation_stack.m_dot_K_stack->at(it);

                    // const auto Gamma            = strain_parameterisation_stack.m_Gamma_stack->at(it);
                    // const auto dot_Gamma        = strain_parameterisation_stack.m_dot_Gamma_stack->at(it);
                    // const auto Delta_Gamma      = Delta_strain_parameterisation_stack.m_Gamma_stack->at(it);
                    const auto Delta_dot_Gamma  = Delta_strain_parameterisation_stack.m_dot_Gamma_stack->at(it);
                    // const auto Delta_ddot_Gamma = Delta_strain_parameterisation_stack.m_ddot_Gamma_stack->at(it);

                    const auto Omega = Omega_optimized.col(it);
                    const auto V = V_optimized.col(it);

                    // const auto dot_Omega = dot_Omega_optimized.col(it);
                    // const auto dot_V     = dot_V_optimized.col(it);


                    const auto Delta_Omega_ = Delta_Omega.col(it);
                    const auto Delta_V_     = Delta_V.col(it);
                    const auto Delta_dot_V_ = Delta_dot_V.col(it);

                    const auto Delta_dot_Omega_ = Delta_dot_Omega.col(it);


                    // const auto C = C_optimized.col(it);
                    const auto N = N_optimized.col(it);



                    b = - ::LieAlgebra::skew( Delta_K )*N
                        - M_linear*Delta_dot_V_
                        - ::LieAlgebra::skew( Delta_Omega_ )*M_linear*V
                        - ::LieAlgebra::skew( Omega )*M_linear*Delta_V_
                        - ::LieAlgebra::skew( Delta_dot_Gamma )*Omega;

                    for (unsigned int row = 0; row < state_dim; ++row)
                        b_stack(index + row * (Nc-1)) = b(row);
                    index++;
                }
                // Solve
                Delta_N_stack = Cf.partialPivLu().solve(b_stack);
                toMatrixForward(Delta_N_stack, Delta_N, state_dim, Nc);








                //                      Delta C Integration
                //  Coefficients
                for(int index = 0; const auto it : forward_points) {



                    // const auto K     = strain_parameterisation_stack.m_K_stack->at(it);
                    // const auto dot_K = strain_parameterisation_stack.m_dot_K_stack->at(it);


                    const auto Delta_K          = Delta_strain_parameterisation_stack.m_K_stack->at(it);
                    // const auto Delta_dot_K      = Delta_strain_parameterisation_stack.m_dot_K_stack->at(it);

                    const auto Gamma            = strain_parameterisation_stack.m_Gamma_stack->at(it);
                    // const auto dot_Gamma        = strain_parameterisation_stack.m_dot_Gamma_stack->at(it);
                    const auto Delta_Gamma      = Delta_strain_parameterisation_stack.m_Gamma_stack->at(it);
                    const auto Delta_dot_Gamma  = Delta_strain_parameterisation_stack.m_dot_Gamma_stack->at(it);
                    // const auto Delta_ddot_Gamma = Delta_strain_parameterisation_stack.m_ddot_Gamma_stack->at(it);

                    const auto Omega = Omega_optimized.col(it);
                    const auto V = V_optimized.col(it);

                    // const auto dot_Omega = dot_Omega_optimized.col(it);
                    // const auto dot_V     = dot_V_optimized.col(it);


                    const auto Delta_Omega_ = Delta_Omega.col(it);
                    const auto Delta_V_     = Delta_V.col(it);
                    const auto Delta_dot_V_ = Delta_dot_V.col(it);

                    const auto Delta_dot_Omega_ = Delta_dot_Omega.col(it);


                    const auto C = C_optimized.col(it);
                    const auto N = N_optimized.col(it);

                    const auto Delta_N_ = Delta_N.col(it);



                    b = - ::LieAlgebra::skew( Delta_K )*C
                        - ::LieAlgebra::skew( Delta_Gamma )*N
                        - ::LieAlgebra::skew( Gamma )*Delta_N_
                        - M_angular*Delta_dot_Omega_
                        - ::LieAlgebra::skew( Delta_Omega_ )*M_angular*Omega
                        - ::LieAlgebra::skew( Omega )*M_angular*Delta_Omega_
                        - ::LieAlgebra::skew( Delta_V_ )*M_linear*V
                        - ::LieAlgebra::skew( V )*M_linear*Delta_V_
                        ;

                    for (unsigned int row = 0; row < state_dim; ++row)
                        b_stack(index + row * (Nc-1)) = b(row);
                    index++;
                }
                // Solve
                Delta_C_stack = Cf.partialPivLu().solve(b_stack);
                toMatrixForward(Delta_C_stack, Delta_C, state_dim, Nc);


            }

        }

    };



    solve_TIDM();

    // ::benchmark::RegisterBenchmark("Solution TIDM", [&](::benchmark::State &t_state){
    //     const bool optimized = false;
    //     while(t_state.KeepRunning())
    //         solve_TIDM(optimized);
    // })->Unit(::benchmark::kMicrosecond)->Repetitions(repetitions);


    // ::benchmark::RegisterBenchmark("Solution TIDM optimized", [&](::benchmark::State &t_state){
    //     while(t_state.KeepRunning())
    //         solve_TIDM();
    // })->Unit(::benchmark::kMicrosecond)->Repetitions(repetitions);





    // ::benchmark::Initialize(&argc, argv);


    // ::benchmark::RunSpecifiedBenchmarks();



    return 0;
}
