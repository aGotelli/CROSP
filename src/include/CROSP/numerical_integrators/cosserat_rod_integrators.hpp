#ifndef COSSERAT_ROD_INTEGRATORS_HPP
#define COSSERAT_ROD_INTEGRATORS_HPP





#include <Eigen/Dense>
#include <memory>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include "CROSP/polynomial_representation/polynomial_representation.hpp"
#include "CROSP/rod_properties/rod_properties.hpp"




namespace CROSP::numerical_integrators {


struct FullODEStatesObservations {

    FullODEStatesObservations()=default;

    FullODEStatesObservations(const unsigned int t_number_of_cols,
                              const unsigned int t_generalised_coordinates_dimension,
                              const std::pair<unsigned int, unsigned int> t_orientation_dimensions={4, 1})
        : m_orientation_dimensions(t_orientation_dimensions),
          m_generalised_coordinates_dimension(t_generalised_coordinates_dimension)
    {
        orientation_stack.resize(m_orientation_dimensions.first, m_orientation_dimensions.second*t_number_of_cols);
        r_stack.resize(3, t_number_of_cols);

        Omega_stack.resize(3, t_number_of_cols);
        V_stack.resize(3, t_number_of_cols);

        dot_Omega_stack.resize(3, t_number_of_cols);
        dot_V_stack.resize(3, t_number_of_cols);

        C_stack.resize(3, t_number_of_cols);
        N_stack.resize(3, t_number_of_cols);

        Qa_stack.resize(m_generalised_coordinates_dimension, t_number_of_cols);
        Qad_stack.resize(m_generalised_coordinates_dimension, t_number_of_cols);
    }



    void writeState(const Eigen::VectorXd t_state, const unsigned int t_col)
    {

        unsigned int start_row = 0;
        const unsigned int start_col = 0;
        unsigned int block_rows = m_orientation_dimensions.first;
        const unsigned int block_cols = 1;

        orientation_stack.block(0, t_col, block_rows, 1) =
                t_state.block(start_row, start_col, block_rows, block_cols);

        //  Update row initial
        start_row += block_rows;

        block_rows = 3;

        r_stack.block(0, t_col, block_rows, 1) =
                t_state.block(start_row, start_col, block_rows, block_cols);
        //  Update row initial
        start_row += block_rows;


        Omega_stack.block(0, t_col, block_rows, 1) =
                t_state.block(start_row, start_col, block_rows, block_cols);
        //  Update row initial
        start_row += block_rows;


        V_stack.block(0, t_col, block_rows, 1) =
                t_state.block(start_row, start_col, block_rows, block_cols);
        //  Update row initial
        start_row += block_rows;


        dot_Omega_stack.block(0, t_col, block_rows, 1) =
                t_state.block(start_row, start_col, block_rows, block_cols);
        //  Update row initial
        start_row += block_rows;


        dot_V_stack.block(0, t_col, block_rows, 1) =
                t_state.block(start_row, start_col, block_rows, block_cols);
        //  Update row initial
        start_row += block_rows;


        C_stack.block(0, t_col, block_rows, 1) =
                t_state.block(start_row, start_col, block_rows, block_cols);
        //  Update row initial
        start_row += block_rows;


        N_stack.block(0, t_col, block_rows, 1) =
                t_state.block(start_row, start_col, block_rows, block_cols);
        //  Update row initial
        start_row += block_rows;


        block_rows = m_generalised_coordinates_dimension;

        Qa_stack.block(0, t_col, block_rows, 1) =
                t_state.block(start_row, start_col, block_rows, block_cols);

    }

    std::pair<unsigned int, unsigned int> m_orientation_dimensions;

    unsigned int m_generalised_coordinates_dimension;

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


        FullODEStatesObservations() = integrator.getIDMStatesObservations();

        FullODEStatesObservations() = integrator.getTIDMStatesObservations();

        Eigen::MatrixXd() = integrator.integratePhiTPhi();

};




}   //  namespace CROSP::numerical_integrators





#endif // COSSERAT_ROD_INTEGRATORS_HPP
