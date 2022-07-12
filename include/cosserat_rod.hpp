/**
 * \file cosserat_rod.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the definition of strain based parameterised Cosserat rod
 * \date 12-07-2022
 *
 * \copyright Copyright (c) 2022 Andrea Gotelli
 *
 *
 */



#ifndef COSSERAT_ROD_HPP
#define COSSERAT_ROD_HPP


#include "OSNI/include/OSNI.hpp"
#include <Eigen/Dense>


//struct QuaternionIntegrator : public OSNI::ODEA {

//    QuaternionIntegrator() : OSNI::ODEA(4) {}

//    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(unsigned int t_point)
//    {
//        return Eigen::MatrixXd();
//    }
//};



class CosseratRod
{
public:
    CosseratRod();



};

#endif // COSSERAT_ROD_HPP
