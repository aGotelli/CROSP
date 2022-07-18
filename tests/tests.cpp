
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


int main(int argc, char *argv[])
{

    CROSP::CosseratRod cosserat_rod;

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

    writeToFile("quaternions_relative_error", quaternions_relative_error, "../../tests/data/");
    writeToFile("positions_relative_error", positions_relative_error, "../../tests/data/");
    writeToFile("angular_velocity_relative_error", angular_velocity_relative_error, "../../tests/data/");
    writeToFile("linear_velocity_relative_error", linear_velocity_relative_error, "../../tests/data/");
    writeToFile("angular_acceleration_relative_error", angular_acceleration_relative_error, "../../tests/data/");
    writeToFile("linear_acceleration_relative_error", linear_acceleration_relative_error, "../../tests/data/");

    std::cout << "Relative errors in quaternions : \n" << quaternions_relative_error << "\n\n" << std::endl;
    std::cout << "Relative errors in positions : \n" << positions_relative_error << "\n\n" << std::endl;
    std::cout << "Relative errors in angular velocities : \n" << angular_velocity_relative_error << "\n\n" << std::endl;
    std::cout << "Relative errors in linear velocities : \n" << linear_velocity_relative_error << "\n\n" << std::endl;
    std::cout << "Relative errors in angular accelerations : \n" << angular_acceleration_relative_error << "\n\n" << std::endl;
    std::cout << "Relative errors in linear accelerations : \n" << linear_acceleration_relative_error << "\n\n" << std::endl;

    return 0;
}
