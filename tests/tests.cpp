
#include "cosserat_rod.hpp"

#include <fstream>


/*!
 * \brief writeToFile writes a Eigen matrix into file
 * \param t_name    name of the file
 * \param t_matrix  the Eigen matrix to write into the file
 * \param t_relative_path_from_build the relative path from the build folder to the file location. Default is none so the file is written in the build directory)
 * \param t_format  the specification for writing. (Default in column major allignment, with comma column separator and 8 digits precision)
 */
void writeToFile(std::string t_name,
                 const Eigen::MatrixXd &t_matrix,
                 std::string t_relative_path_from_build="",
                 const Eigen::IOFormat &t_format=Eigen::IOFormat(16, 0, ","))
{
    if(not t_relative_path_from_build.empty()
            ){
        //  Ensure relative path ends with a backslash only if a path is given
        if(not t_relative_path_from_build.ends_with('/'))
            t_relative_path_from_build.append("/");
    }


    //  Ensure it ends with .csv
    if(t_name.find(".csv") == std::string::npos)
        t_name.append(".csv");

    //  The file will be created in the location given by the realtive path and with the given name
    const auto file_name_and_location = t_relative_path_from_build + t_name;

    //  Create file in given location with given name
    std::ofstream file(file_name_and_location.c_str());

    //  Put matrix in this file
    file << t_matrix.format(t_format);

    //  Close the file
    file.close();
 }


void LoadEigenMatrixFromFile(Eigen::MatrixXd &t_matrix,
                             std::string t_file_name,
                             std::string t_relative_file_path_from_build="",
                             const Eigen::IOFormat &t_format=Eigen::IOFormat())
{
    // the matrix entries are stored in this variable row-wise. For example if we have the matrix:
    // M=[a b c
    //    d e f]
    // the entries are stored as matrixEntries=[a,b,c,d,e,f], that is the variable "matrixEntries" is a row vector
    // later on, this vector is mapped into the Eigen matrix format
    std::vector<double> matrixEntries;

    std::string file_to_open = t_relative_file_path_from_build + t_file_name;

    // in this object we store the data from the matrix
    std::ifstream matrixDataFile;
    matrixDataFile.open(file_to_open);

    if(not matrixDataFile.is_open()){
        std::cout << "Could not open the file" << std::endl;
        return;
    }


    // this variable is used to store the row of the matrix that contains commas
    std::string matrixRowString;

    // this variable is used to store the matrix entry;
    std::string matrixEntry;

    // this variable is used to track the number of rows
    int matrixRowNumber = 0;

    // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
    while (std::getline(matrixDataFile, matrixRowString))
    {
        //convert matrixRowString that is a string to a stream variable.
        std::stringstream matrixRowStringStream(matrixRowString);

        // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
        while (std::getline(matrixRowStringStream, matrixEntry, ','))
        {
            //here we convert the string to double and fill in the row vector storing all the matrix entries
            matrixEntries.push_back(stod(matrixEntry));
        }
        //update the column numbers
        matrixRowNumber++;
    }

    // here we convet the vector variable into the matrix and return the resulting object,
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
    t_matrix = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}



void testRelativeComputations()
{
    const unsigned int number_of_chebyshev_points = 17;

    const auto chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(number_of_chebyshev_points);

    for(unsigned int point=0; point<chebyshev_points.size(); point++)
        std::cout << "Point : " << point << " is at " << chebyshev_points[point] << std::endl;

    CROSP::CosseratRod cosserat_rod(number_of_chebyshev_points);

    constexpr unsigned int ne = 3;
    constexpr unsigned int na = 3;
    Eigen::VectorXd qe, dot_qe, ddot_qe;

    Eigen::MatrixXd generalised_coordinates;

    LoadEigenMatrixFromFile(generalised_coordinates, "qe-dot_qe-ddotqe.csv", "../../tests/data/");

    qe = generalised_coordinates(Eigen::all, 0);
    dot_qe = generalised_coordinates(Eigen::all, 1);
    ddot_qe = generalised_coordinates(Eigen::all, 2);

    std::cout << "qe : \n" << qe << "\n\n\n";
    std::cout << "dot_qe : \n" << dot_qe << "\n\n\n";
    std::cout << "ddot_qe : \n" << ddot_qe << "\n\n\n";

    cosserat_rod.updateParameterisation(qe, dot_qe, ddot_qe);
    std::cout << "\n\n\n";

    Eigen::MatrixXd init_quaternion;
    Eigen::MatrixXd init_position;
    Eigen::MatrixXd init_eta;
    Eigen::MatrixXd init_angular_velocity;
    Eigen::MatrixXd init_linear_velocity;
    Eigen::MatrixXd init_dot_eta;
    Eigen::MatrixXd init_angular_acceleration;
    Eigen::MatrixXd init_linear_acceleration;

    LoadEigenMatrixFromFile(init_quaternion, "Q0.csv", "../../tests/data/");
    LoadEigenMatrixFromFile(init_position, "r0.csv", "../../tests/data/");
    LoadEigenMatrixFromFile(init_eta, "eta0.csv", "../../tests/data/");
    LoadEigenMatrixFromFile(init_dot_eta, "dot_eta0.csv", "../../tests/data/");

    init_angular_velocity = init_eta.block<3,1>(0, 0);
    init_linear_velocity = init_eta.block<3,1>(3, 0);

    init_angular_acceleration = init_dot_eta.block<3,1>(0, 0);
    init_linear_acceleration = init_dot_eta.block<3,1>(3, 0);

    cosserat_rod.forwardKinematics(init_quaternion,
                                   init_position,
                                   init_angular_velocity,
                                   init_linear_velocity,
                                   init_angular_acceleration,
                                   init_linear_acceleration);


    Eigen::MatrixXd quaternion_stack_from_ode45;
    Eigen::MatrixXd position_stack_from_ode45;
    Eigen::MatrixXd angular_velocity_stack_from_ode45;
    Eigen::MatrixXd linear_velocity_stack_from_ode45;
    Eigen::MatrixXd angular_acceleration_stack_from_ode45;
    Eigen::MatrixXd linear_acceleration_stack_from_ode45;

    LoadEigenMatrixFromFile(quaternion_stack_from_ode45, "quaternion_stack.csv", "../../tests/data/");
    LoadEigenMatrixFromFile(position_stack_from_ode45, "position_stack.csv", "../../tests/data/");
    LoadEigenMatrixFromFile(angular_velocity_stack_from_ode45, "angular_velocity_stack.csv", "../../tests/data/");
    LoadEigenMatrixFromFile(linear_velocity_stack_from_ode45, "linear_velocity_stack.csv", "../../tests/data/");
    LoadEigenMatrixFromFile(angular_acceleration_stack_from_ode45, "angular_acceleration_stack.csv", "../../tests/data/");
    LoadEigenMatrixFromFile(linear_acceleration_stack_from_ode45, "linear_acceleration_stack.csv", "../../tests/data/");


    const unsigned int number_of_Chebyshev_points = position_stack_from_ode45.rows();

    const Eigen::MatrixXd quaternions_relative_error = quaternion_stack_from_ode45 - cosserat_rod.m_quaternion_integrator->getStack();
    const Eigen::MatrixXd positions_relative_error = position_stack_from_ode45 - cosserat_rod.m_position_integrator->getStack();
    const Eigen::MatrixXd angular_velocity_relative_error = angular_velocity_stack_from_ode45 - cosserat_rod.m_angular_velocity_integrator->getStack();
    const Eigen::MatrixXd linear_velocity_relative_error = linear_velocity_stack_from_ode45 - cosserat_rod.m_linear_velocity_integrator->getStack();
    const Eigen::MatrixXd angular_acceleration_relative_error = angular_acceleration_stack_from_ode45 - cosserat_rod.m_angular_acceleration_integrator->getStack();
    const Eigen::MatrixXd linear_acceleration_relative_error = linear_acceleration_stack_from_ode45 - cosserat_rod.m_linear_acceleration_integrator->getStack();

//    writeToFile("quaternions_relative_error", quaternions_relative_error, "../../tests/data/");
//    writeToFile("positions_relative_error", positions_relative_error, "../../tests/data/");
//    writeToFile("angular_velocity_relative_error", angular_velocity_relative_error, "../../tests/data/");
//    writeToFile("linear_velocity_relative_error", linear_velocity_relative_error, "../../tests/data/");
//    writeToFile("angular_acceleration_relative_error", angular_acceleration_relative_error, "../../tests/data/");
//    writeToFile("linear_acceleration_relative_error", linear_acceleration_relative_error, "../../tests/data/");

//    std::cout << "Relative errors in quaternions : \n" << quaternions_relative_error << "\n\n" << std::endl;
//    std::cout << "Relative errors in positions : \n" << positions_relative_error << "\n\n" << std::endl;
//    std::cout << "Relative errors in angular velocities : \n" << angular_velocity_relative_error << "\n\n" << std::endl;
//    std::cout << "Relative errors in linear velocities : \n" << linear_velocity_relative_error << "\n\n" << std::endl;
//    std::cout << "Relative errors in angular accelerations : \n" << angular_acceleration_relative_error << "\n\n" << std::endl;
//    std::cout << "Relative errors in linear accelerations : \n" << linear_acceleration_relative_error << "\n\n" << std::endl;



    Eigen::MatrixXd init_F1;
    Eigen::MatrixXd init_force;
    Eigen::MatrixXd init_couple;

    LoadEigenMatrixFromFile(init_F1, "F1.csv", "../../tests/data/");
    init_couple = init_F1.block<3,1>(0, 0);
    init_force = init_F1.block<3,1>(3, 0);

    cosserat_rod.backwardDynamics(init_force, init_couple);

    Eigen::MatrixXd force_stack_from_ode45;
    Eigen::MatrixXd couple_stack_from_ode45;

    LoadEigenMatrixFromFile(force_stack_from_ode45, "force_stack.csv", "../../tests/data/");
    LoadEigenMatrixFromFile(couple_stack_from_ode45, "couple_stack.csv", "../../tests/data/");


    Eigen::MatrixXd forces(number_of_Chebyshev_points*3, 2);
    forces << force_stack_from_ode45, cosserat_rod.m_internal_forces_integrator->getStack();
    std::cout << "Forces stakcs : \n" << forces << std::endl << std::endl << std::endl;

    const Eigen::MatrixXd force_relative_error = force_stack_from_ode45 - cosserat_rod.m_internal_forces_integrator->getStack();
    const Eigen::MatrixXd couple_relative_error = couple_stack_from_ode45 - cosserat_rod.m_internal_couples_integrator->getStack();



//    for(unsigned int point = 1; point<=number_of_Chebyshev_points; point++)
//        std::cout << "At point : " << point << " that is " << cosserat_rod.m_position_integrator->getStateAtPoint(point)(0) << " the couples : " << cosserat_rod.m_internal_couples_integrator->getStateAtPoint(point).transpose() << std::endl;


//    writeToFile("force_relative_error", force_relative_error, "../../tests/data/");
//    writeToFile("couple_relative_error", couple_relative_error, "../../tests/data/");

    std::cout << "Relative errors in forces : \n" << force_relative_error << "\n\n" << std::endl;
    std::cout << "Relative errors in couples : \n" << couple_relative_error << "\n\n" << std::endl;

    Eigen::MatrixXd couples(number_of_Chebyshev_points*3, 2);
    couples << couple_stack_from_ode45, cosserat_rod.m_internal_couples_integrator->getStack();
    std::cout << "Couples stakcs : \n" << couples << std::endl << std::endl << std::endl;
}





void testRelativePrecision(const unsigned int t_number_of_Chebyshev_points)
{


    CROSP::CosseratRod cosserat_rod(t_number_of_Chebyshev_points);

    constexpr unsigned int ne = 3;
    constexpr unsigned int na = 3;
    Eigen::VectorXd qe, dot_qe, ddot_qe;

    Eigen::MatrixXd generalised_coordinates;

    LoadEigenMatrixFromFile(generalised_coordinates, "qe-dot_qe-ddotqe_"+std::to_string(t_number_of_Chebyshev_points)+".csv", "../../tests/data/");

    qe = generalised_coordinates(Eigen::all, 0);
    dot_qe = generalised_coordinates(Eigen::all, 1);
    ddot_qe = generalised_coordinates(Eigen::all, 2);


    cosserat_rod.updateParameterisation(qe, dot_qe, ddot_qe);

    Eigen::MatrixXd init_quaternion;
    Eigen::MatrixXd init_position;
    Eigen::MatrixXd init_eta;
    Eigen::MatrixXd init_angular_velocity;
    Eigen::MatrixXd init_linear_velocity;
    Eigen::MatrixXd init_dot_eta;
    Eigen::MatrixXd init_angular_acceleration;
    Eigen::MatrixXd init_linear_acceleration;

    LoadEigenMatrixFromFile(init_quaternion, "Q0_"+std::to_string(t_number_of_Chebyshev_points)+".csv", "../../tests/data/");
    LoadEigenMatrixFromFile(init_position, "r0_"+std::to_string(t_number_of_Chebyshev_points)+".csv", "../../tests/data/");
    LoadEigenMatrixFromFile(init_eta, "eta0_"+std::to_string(t_number_of_Chebyshev_points)+".csv", "../../tests/data/");
    LoadEigenMatrixFromFile(init_dot_eta, "dot_eta0_"+std::to_string(t_number_of_Chebyshev_points)+".csv", "../../tests/data/");


    init_angular_velocity = init_eta.block<3,1>(0, 0);
    init_linear_velocity = init_eta.block<3,1>(3, 0);

    init_angular_acceleration = init_dot_eta.block<3,1>(0, 0);
    init_linear_acceleration = init_dot_eta.block<3,1>(3, 0);

    cosserat_rod.forwardKinematics(init_quaternion,
                                   init_position,
                                   init_angular_velocity,
                                   init_linear_velocity,
                                   init_angular_acceleration,
                                   init_linear_acceleration);




    Eigen::MatrixXd final_state_from_ode45;


    LoadEigenMatrixFromFile(final_state_from_ode45, "final_state_"+std::to_string(t_number_of_Chebyshev_points)+".csv", "../../tests/data/");


    Eigen::MatrixXd final_state = final_state_from_ode45;


    final_state.block<4,1>(0, 0) = cosserat_rod.m_quaternion_integrator->getStateAtPoint(0);
    final_state.block<3,1>(4, 0) = cosserat_rod.m_position_integrator->getStateAtPoint(0);
    final_state.block<3,1>(7, 0) = cosserat_rod.m_angular_velocity_integrator->getStateAtPoint(0);
    final_state.block<3,1>(10, 0) = cosserat_rod.m_linear_velocity_integrator->getStateAtPoint(0);
    final_state.block<3,1>(13, 0) = cosserat_rod.m_angular_acceleration_integrator->getStateAtPoint(0);
    final_state.block<3,1>(16, 0) = cosserat_rod.m_linear_acceleration_integrator->getStateAtPoint(0);


//    std::cout << "final_state :\n" << final_state << std::endl;

    const Eigen::MatrixXd final_state_relative_error = final_state_from_ode45 - final_state;

    writeToFile("final_state_relative_error_"+std::to_string(t_number_of_Chebyshev_points), final_state_relative_error, "../../tests/data/");

    std::cout << "For : " << t_number_of_Chebyshev_points << ", the errors are : \n" << final_state_relative_error << "\n\n" << std::endl;
}


void processData(const std::vector<unsigned int> &t_tested_points)
{

    Eigen::MatrixXd precision_benchmack(t_tested_points.size(), 2);

    for(unsigned int i=0; i<t_tested_points.size(); i++){
        Eigen::MatrixXd final_state_relative_error;
        LoadEigenMatrixFromFile(final_state_relative_error, "final_state_relative_error_"+std::to_string(t_tested_points[i])+".csv", "../../tests/data/");

        precision_benchmack(i, 0) = t_tested_points[i];
        precision_benchmack(i, 1) = final_state_relative_error.norm();
    }

    writeToFile("precision_benchmack", precision_benchmack, "../../tests/data/");
}

#include <iomanip>
void processDataAdavanced(const std::vector<unsigned int> &t_tested_points)
{

    Eigen::MatrixXd precision_benchmack_advanced(t_tested_points.size(), 7);
//    std::ofstream myfile;
//    myfile.open ("precision_benchmack_advanced.csv");
//    myfile << "points, quaternions_norm, positions_norm, angular_velocities_norm, linear_velocites_norm, angular_acceleration_norm, linear_acceleration_norm, time, standard_deviation \n";

//    myfile << std::fixed << std::setprecision(16) << std::endl;
    for(unsigned int i=0; i<t_tested_points.size(); i++){
        Eigen::MatrixXd final_state_relative_error;
        LoadEigenMatrixFromFile(final_state_relative_error, "final_state_relative_error_"+std::to_string(t_tested_points[i])+".csv", "../../tests/data/");

        precision_benchmack_advanced.row(i) <<  t_tested_points[i],
                                                final_state_relative_error.block<4,1>(0, 0).norm(),
                                                final_state_relative_error.block<3,1>(4, 0).norm(),
                                                final_state_relative_error.block<3,1>(7, 0).norm(),
                                                final_state_relative_error.block<3,1>(10, 0).norm(),
                                                final_state_relative_error.block<3,1>(13, 0).norm(),
                                                final_state_relative_error.block<3,1>(16, 0).norm();
    }

    writeToFile("precision_benchmack_advanced", precision_benchmack_advanced, "../../tests/data/");
}


int main(int argc, char *argv[])
{


    testRelativeComputations();

//    std::vector<unsigned int> points_to_test {
//        10,
//        15,
//        20,
//        25,
//        30,
//        35
//    };

//    for(const auto points : points_to_test)
//        testRelativePrecision( points );

//    processDataAdavanced( points_to_test );

    return 0;
}
