#ifndef CHEBYSHEV_RECONSTRUCTOR_HPP
#define CHEBYSHEV_RECONSTRUCTOR_HPP


#include <Eigen/Dense>

#include <vector>
#include <boost/math/special_functions/chebyshev.hpp>

#include <iostream>


class ChebyshevReconstructor {
public:


    ChebyshevReconstructor(const unsigned int t_number_of_Chebyshev_points,
                           const std::vector<double> &t_reconstruction_points,
                           const unsigned int t_state_dimension)
        : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
          m_reconstruction_points(t_reconstruction_points),
          m_state_dimension(t_state_dimension)
    {}


    Eigen::MatrixXd reconstructVector(const Eigen::MatrixXd &t_vector_on_Chebyshev_grid)const
    {
        //  Perform the cosine transform
        m_CN = m_DDCT * t_vector_on_Chebyshev_grid.transpose();

        std::cout << "CN : \n" << m_CN << "\n\n";

        double T0 = 1.0/sqrt(2.0);  //  First Chebyshev polynomial (Normalized)
        double Th;                  //  Other Chebyshev polynomial (Normalized)
        double s;                   //  Coordinate in the [-1, 1] domain


        Eigen::VectorXd reconstructed_point(m_state_dimension);
        for(unsigned int i=0; const auto &x : m_reconstruction_points) {

            //  Compute for the first point
            reconstructed_point = m_CN.row(0) * T0;

            //  Obtain the others Chebyshev points
            for(unsigned int h=1; h<=m_number_of_Chebyshev_points-1;h++) {

                //  Define the corresponding point on the unit circle
                s = 2*x - 1;
                Th = boost::math::chebyshev_t(h, s);

                reconstructed_point += m_CN.row(h) * Th;

            }

            m_reconstructed_points.row(i++) = reconstructed_point;

        }

        return m_reconstructed_points;
    }



    inline unsigned int getNumberOfReconstructionPoints()const{return m_reconstruction_points.size();}


    const unsigned int m_state_dimension { 3 };

    //  Number of Chebyshev points
    const unsigned int m_number_of_Chebyshev_points;

    const std::vector<double> m_reconstruction_points;

    //  Number of points describing the rod centerline
    const unsigned int m_number_of_reconstruction_points { static_cast<unsigned int>(m_reconstruction_points.size()) };



    //  Direct Discrete Cosine Transform
    const Eigen::MatrixXd m_DDCT { [&](){
            //  Initialize the matrix
            Eigen::MatrixXd DDCT(m_number_of_Chebyshev_points, m_number_of_Chebyshev_points);

            const auto k = [](const unsigned int h){ return h==0 ? 1.0/sqrt(2.0) : pow((-1), h);};

            for(unsigned int h=0; h<=m_number_of_Chebyshev_points-1;h++) {
                for(unsigned int j=0; j<=m_number_of_Chebyshev_points-1;j++){
//                    DDCT(h, j) = (2.0/m_number_of_Chebyshev_points)
//                                * k(h)*cos( (2.0*j + 1.0)*h*M_PI/(2.0*m_number_of_Chebyshev_points) );

                    const auto pm = [this](const unsigned int m){
                        if(m == 0 || m == m_number_of_Chebyshev_points)
                            return 1.0;

                        return 2.0;
                    };
                    DDCT(h, j) = (2.0/m_number_of_Chebyshev_points)
                                * k(h)*cos( (2.0*j + 1.0)*h*M_PI/(2.0*m_number_of_Chebyshev_points) );
                }
            }


            //  Return the matrix
            return DDCT;

          }()};




    //  Cosine transform of the fucntion (marked mutable, here we want just to allocate the memory)
    mutable Eigen::MatrixXd m_CN { Eigen::MatrixXd(m_number_of_Chebyshev_points, m_state_dimension) };


    //  Matrix of reconstructed points (marked mutable, here we want just to allocate the memory)
    mutable Eigen::MatrixXd m_reconstructed_points { Eigen::MatrixXd(m_number_of_reconstruction_points, m_state_dimension) };


};


#endif // CHEBYSHEV_RECONSTRUCTOR_HPP
