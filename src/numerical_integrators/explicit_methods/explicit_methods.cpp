#include "CROSP/numerical_integrators/explicit_methods/explicit_methods.hpp"


#include <regex>

namespace CROSP::numerical_integrators::explicit_methods {




std::string extractBoostSignatureFromtypeid(std::string t_method_type)
{

    std::regex lowercase_regex("[a-z_]+");
    std::smatch match;

    std::stringstream name;
    while (std::regex_search(t_method_type, match, lowercase_regex)) {


        if (match[0].str().find("igen") != std::string::npos)
            break;
        else
            name << "::";


        name << match[0];
        t_method_type = match.suffix().str();
    }

    return name.str();

}




Eigen::Matrix4d getA(const Eigen::Vector3d t_k){
    Eigen::Matrix4d A;
    A   <<     0  , -t_k(0),  -t_k(1),  -t_k(2),
            t_k(0),     0  ,   t_k(2),  -t_k(1),
            t_k(1), -t_k(2),      0  ,   t_k(0),
            t_k(2),  t_k(1),  -t_k(0),      0  ;

    return A;
}


ExplicitIntegrationODEs::ExplicitIntegrationODEs(const polynomial_representation::PolynomialRepresentation t_polynomial_representation,
                                                     std::shared_ptr<const rod_properties::RodProperties> t_rod_properties)
    : m_rod_properties(t_rod_properties),
      m_polynomial_representation(t_polynomial_representation)
{}



void ExplicitIntegrationODEs::forwardStaticODEs(const Eigen::VectorXd &t_y,
                                                Eigen::VectorXd &t_dyds,
                                                const double t_X) const
{
    /*  Preprocessing    */
    Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

    //  Get the strains for the rod
    const Eigen::VectorXd Xi = BPhi*m_qe + m_constant_strain;

    t_dyds = forwardStaticStep(t_y, Xi);
}

ExplicitIntegrationODEs::PoseState ExplicitIntegrationODEs::forwardStaticStep(const PoseState &t_state,
                                                      const ::LieAlgebra::Vector6d &t_Xi) const
{
    //  Decompose the strain
    const Eigen::Vector3d k = t_Xi.block<3,1>(0,0);
    const Eigen::Vector3d gamma = t_Xi.block<3,1>(3,0);


    /*  The state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     */


    //  Unpack state vector
    const Eigen::Quaterniond Q(t_state[0], t_state[1],t_state[2], t_state[3]);
    [[maybe_unused]] const Eigen::Vector3d r = t_state.block<3,1>(4,0);


    const Eigen::Matrix3d R = Q.toRotationMatrix();


    //  Actual ODE
    const Eigen::Vector4d Q_prime = 0.5*getA(k)*Eigen::Vector4d(Q.w(),Q.x(),Q.y(),Q.z());
    const Eigen::Vector3d r_prime = R*gamma;


    //  Packing state vector derivative
    ExplicitIntegrationODEs::PoseState dyds;
    dyds << Q_prime,
            r_prime;

    dyds *= m_rod_length;

    return dyds;

}



void ExplicitIntegrationODEs::forwardODEs(const Eigen::VectorXd &t_y,
                                          Eigen::VectorXd &t_dyds,
                                          const double t_X) const
{
    /*  Preprocessing    */
    Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

    //  Get the strains for the rod
    const ::LieAlgebra::Vector6d Xi      = BPhi*m_qe + m_constant_strain;
    const ::LieAlgebra::Vector6d dot_Xi  = BPhi*m_dot_qe;
    const ::LieAlgebra::Vector6d ddot_Xi = BPhi*m_ddot_qe;

    t_dyds = forwardStep(t_y, Xi, dot_Xi, ddot_Xi);
}



ExplicitIntegrationODEs::ForwardKinematicState ExplicitIntegrationODEs::forwardStep(const ExplicitIntegrationODEs::ForwardKinematicState &t_state,
                                                                              const ::LieAlgebra::Vector6d &t_Xi,
                                                                              const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                                              const ::LieAlgebra::Vector6d &t_ddot_Xi) const
{
    /*  The state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     *  | η |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 7-12
     *  | η̇ |   Ω1 , Ω2 , Ω3 , V1 , V2 , V3 13-18
     */


    //  Unpack state vector
    const ExplicitIntegrationODEs::PoseState g = t_state.block<7,1>(0,0);
    const ::LieAlgebra::Vector6d eta = t_state.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d eta_dot = t_state.block<6,1>(13,0);


    const auto ad_Xi = ::LieAlgebra::ad(t_Xi);
    const auto ad_dot_Xi = ::LieAlgebra::ad(t_dot_Xi);

    //  Actual ODE
    const ExplicitIntegrationODEs::PoseState g_prime = forwardStaticStep(g, t_Xi);
    const ::LieAlgebra::Vector6d eta_prime = - ad_Xi*eta + t_dot_Xi;
    const ::LieAlgebra::Vector6d eta_dot_prime = - ad_Xi*eta_dot - ad_dot_Xi*eta + t_ddot_Xi;



    ExplicitIntegrationODEs::ForwardKinematicState dydx;

    //  Packing state vector derivative
    dydx << g_prime,
            eta_prime,
            eta_dot_prime;

    dydx *= m_rod_length;

    return dydx;
}






::LieAlgebra::Vector6d ExplicitIntegrationODEs::getLambdaPrime(const Eigen::Quaterniond &t_Q,
                                                            const ::LieAlgebra::Vector6d &t_Lambda,
                                                            const ::LieAlgebra::Vector6d &t_eta,
                                                            const ::LieAlgebra::Vector6d &t_dot_eta,
                                                            const ::LieAlgebra::Matrix6d &t_ad_Xi,
                                                            const double &t_X)const
{

    //  Some needed variables
    const Eigen::Matrix3d R = t_Q.toRotationMatrix();
    ::LieAlgebra::Vector6d F_bar = ::LieAlgebra::Vector6d::Zero();
    F_bar.block<3, 1>(3, 0) = R.transpose() * m_rod_properties->distributedGravitationalForce();


    const auto M = m_rod_properties->m_M;

    const Eigen::VectorXd Lambda_prime = t_ad_Xi.transpose()*t_Lambda
                                            + M*t_dot_eta
                                            - ::LieAlgebra::ad(t_eta).transpose()*M*t_eta
                                            - F_bar;

    return Lambda_prime;
}




void ExplicitIntegrationODEs::backwardODEs(const Eigen::VectorXd &t_y,
                                          Eigen::VectorXd &t_dyds,
                                          const double t_X) const
{
    /*  Preprocessing    */
    Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

    //  Get the strains for the rod
    const Eigen::VectorXd Xi      = BPhi*m_qe + m_constant_strain;
    const Eigen::VectorXd dot_Xi  = BPhi*m_dot_qe;
    const Eigen::VectorXd ddot_Xi = BPhi*m_ddot_qe;

    t_dyds = backwardStep(t_y, Xi, dot_Xi, ddot_Xi, BPhi, t_X);
}




Eigen::VectorXd ExplicitIntegrationODEs::backwardStep(const Eigen::VectorXd &t_y,
                                                   const ::LieAlgebra::Vector6d &t_Xi,
                                                   const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                   const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                                   const Eigen::MatrixXd &t_BPhi,
                                                   const double &t_X)const
{
    /*  The state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     *  | η |   Ω1, Ω2, Ω3, V1, V2, V3      7-12
     *  | η̇ |   Ω1, Ω2, Ω3, V1, V2, V3     13-18
     *  | Λ |   C1, C2, C3, N1, N2, N3     19-24
     *  | Qa|                              25-(25+ne*na)
     */

    //  Unpack state vector
    const ExplicitIntegrationODEs::ForwardKinematicState kinematic_state = t_y.block<19,1>(0,0);

    const ::LieAlgebra::Vector6d eta = t_y.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d dot_eta = t_y.block<6,1>(13,0);
    const ::LieAlgebra::Vector6d Lambda = t_y.block<6,1>(19,0);




    const auto ad_Xi = ::LieAlgebra::ad(t_Xi);



    //  Actual ODE
    const ExplicitIntegrationODEs::ForwardKinematicState kinematic_state_prime = forwardStep(kinematic_state, t_Xi, t_dot_Xi, t_ddot_Xi);


    const Eigen::VectorXd Lambda_prime = getLambdaPrime(Eigen::Quaterniond(t_y[0], t_y[1],t_y[2], t_y[3]),
                                                                           Lambda, eta, dot_eta, ad_Xi, t_X);
    const Eigen::VectorXd Qa_prime = - t_BPhi.transpose()*Lambda;

    //  Packing state vector derivative
    Eigen::VectorXd dyds(25 + m_generalised_coordinates_dimension);
    dyds <<  kinematic_state_prime,
             Lambda_prime,
             Qa_prime;

    dyds *= m_rod_length;

    return dyds;
}

void ExplicitIntegrationODEs::tangentKinematicsODEs(const Eigen::VectorXd &t_y,
                                          Eigen::VectorXd &t_dyds,
                                          const double t_X)const
{
    /*  Preprocessing    */
    Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

    //  Get the strains for the rod
    const ::LieAlgebra::Vector6d Xi      = BPhi*m_qe + m_constant_strain;
    const ::LieAlgebra::Vector6d dot_Xi  = BPhi*m_dot_qe;
    const ::LieAlgebra::Vector6d ddot_Xi = BPhi*m_ddot_qe;


    //  Get the strains for the rod
    const ::LieAlgebra::Vector6d Delta_Xi      = BPhi*m_Delta_qe;
    const ::LieAlgebra::Vector6d Delta_dot_Xi  = BPhi*m_Delta_dot_qe;
    const ::LieAlgebra::Vector6d Delta_ddot_Xi = BPhi*m_Delta_ddot_qe;


    t_dyds = tangentKinematicsStep(t_y, Xi, dot_Xi, ddot_Xi, Delta_Xi, Delta_dot_Xi, Delta_ddot_Xi);
}



ExplicitIntegrationODEs::TangentKinematicState ExplicitIntegrationODEs::tangentKinematicsStep(const ExplicitIntegrationODEs::TangentKinematicState &t_state,
                                                                                        const ::LieAlgebra::Vector6d &t_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_Delta_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_Delta_dot_Xi,
                                                                                        const ::LieAlgebra::Vector6d &t_Delta_ddot_Xi) const
{
    /*  The state has the form
     *  | Q |   w, x, y, z                      0-3
     *  | r |   x, y, z                         4-6
     *  | η |   Ω1, Ω2, Ω3, V1, V2, V3         7-12
     *  | η̇ |   Ω1, Ω2, Ω3, V1, V2, V3         13-18
     *  | ∆ζ |  ∆K1, ∆K2, ∆K3, ∆Γ1, ∆Γ2, ∆Γ3    19-24
     *  | ∆η |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   25-30
     *  | ∆η̇ |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   31-37
     */

    //  Unpack state vector
    const ExplicitIntegrationODEs::ForwardKinematicState kinematic_state = t_state.block<19,1>(0,0);
    const ::LieAlgebra::Vector6d eta = t_state.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d dot_eta = t_state.block<6,1>(13,0);

    const ::LieAlgebra::Vector6d Delta_zeta    = t_state.block<6,1>(19,0);
    const ::LieAlgebra::Vector6d Delta_eta     = t_state.block<6,1>(25,0);
    const ::LieAlgebra::Vector6d Delta_dot_eta = t_state.block<6,1>(31,0);



    const auto ad_Xi = ::LieAlgebra::ad(t_Xi);

    const auto ad_Delta_dot_Xi = ::LieAlgebra::ad(t_Delta_dot_Xi);

    const auto ad_eta = ::LieAlgebra::ad(eta);
    const auto ad_dot_eta = ::LieAlgebra::ad(dot_eta);


    const ExplicitIntegrationODEs::ForwardKinematicState kinematic_state_prime = forwardStep(kinematic_state, t_Xi, t_dot_Xi, t_ddot_Xi);
    const ::LieAlgebra::Vector6d Delta_zeta_prime =
            - ad_Xi*Delta_zeta + t_Delta_Xi;
    const ::LieAlgebra::Vector6d Delta_eta_prime =
            - ad_Xi*Delta_eta + ad_eta*t_Delta_Xi + t_Delta_dot_Xi;
    const ::LieAlgebra::Vector6d Delta_dot_eta_prime =
            - ad_Xi*Delta_dot_eta - ad_Delta_dot_Xi*Delta_eta + ad_eta*t_Delta_dot_Xi + ad_dot_eta*t_Delta_Xi + t_Delta_ddot_Xi;


    //  Packing state vector derivative
    ExplicitIntegrationODEs::TangentKinematicState dydx;
    dydx <<  kinematic_state_prime,
             Delta_zeta_prime,
             Delta_eta_prime,
             Delta_dot_eta_prime;

    dydx *= m_rod_length;

    return dydx;
}



void ExplicitIntegrationODEs::tangentDynamicsODEs(const Eigen::VectorXd &t_y,
                                          Eigen::VectorXd &t_dyds,
                                          const double t_X)const
{
    /*  The state has the form
     *  | Q  |   w, x, y, z                       0-3
     *  | r  |   x, y, z                          4-6
     *  | η  |   Ω1, Ω2, Ω3, V1, V2, V3          7-12
     *  | η̇  |   Ω1, Ω2, Ω3, V1, V2, V3         13-18
     *  | Λ  |   C1, C2, C3, N1, N2, N3          19-24
     *  | ∆ζ |  ∆K1, ∆K2, ∆K3, ∆Γ1, ∆Γ2, ∆Γ3    25-30
     *  | ∆η |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   31-36
     *  | ∆η̇ |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   37-42
     *  | ∆Λ |   C1, C2, C3, N1, N2, N3          43-48
     *  | ∆Qa|                                   59-59+ne
     */

    /*  Preprocessing    */
    Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

    //  Get the strains for the rod
    const ::LieAlgebra::Vector6d Xi      = BPhi*m_qe + m_constant_strain;
    const ::LieAlgebra::Vector6d dot_Xi  = BPhi*m_dot_qe;
    const ::LieAlgebra::Vector6d ddot_Xi = BPhi*m_ddot_qe;


    //  Get the strains for the rod
    const ::LieAlgebra::Vector6d Delta_Xi      = BPhi*m_Delta_qe;
    const ::LieAlgebra::Vector6d Delta_dot_Xi  = BPhi*m_Delta_dot_qe;
    const ::LieAlgebra::Vector6d Delta_ddot_Xi = BPhi*m_Delta_ddot_qe;

    t_dyds = tangentDynamicsStep(t_y,
                                       Xi, dot_Xi, ddot_Xi,
                                       Delta_Xi, Delta_dot_Xi, Delta_ddot_Xi,
                                       BPhi, t_X);
}



Eigen::VectorXd ExplicitIntegrationODEs::tangentDynamicsStep(const Eigen::VectorXd &t_state,
                                                          const ::LieAlgebra::Vector6d &t_Xi,
                                                          const ::LieAlgebra::Vector6d &t_dot_Xi,
                                                          const ::LieAlgebra::Vector6d &t_ddot_Xi,
                                                          const ::LieAlgebra::Vector6d &t_Delta_Xi,
                                                          const ::LieAlgebra::Vector6d &t_Delta_dot_Xi,
                                                          const ::LieAlgebra::Vector6d &t_Delta_ddot_Xi,
                                                          const Eigen::MatrixXd &t_BPhi,
                                                          const double &t_X)const
{

    /*  The state has the form
     *  | Q  |   w, x, y, z                       0-3
     *  | r  |   x, y, z                          4-6
     *  | η  |   Ω1, Ω2, Ω3, V1, V2, V3          7-12
     *  | η̇  |   Ω1, Ω2, Ω3, V1, V2, V3         13-18
     *  | Λ  |   C1, C2, C3, N1, N2, N3          19-24
     *  | ∆ζ |  ∆K1, ∆K2, ∆K3, ∆Γ1, ∆Γ2, ∆Γ3    25-30
     *  | ∆η |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   31-36
     *  | ∆η̇ |  ∆Ω1, ∆Ω2, ∆Ω3, ∆V1, ∆V2, ∆V3   37-42
     *  | ∆Λ |   C1, C2, C3, N1, N2, N3          43-48
     *  | ∆Qa|                                   49-49+ne
     */


    //  Unpack state vector
    ExplicitIntegrationODEs::TangentKinematicState tangent_kinematic_state;
    tangent_kinematic_state << t_state.block<19,1>(0, 0), t_state.block<18,1>(25, 0);


    const ::LieAlgebra::Vector6d Lambda = t_state.block<6, 1>(19,0);

    const ::LieAlgebra::Vector6d eta = t_state.block<6,1>(7,0);
    const ::LieAlgebra::Vector6d dot_eta = t_state.block<6,1>(13,0);

    const ::LieAlgebra::Vector6d Delta_zeta    = t_state.block<6,1>(25,0);
    const ::LieAlgebra::Vector6d Delta_eta     = t_state.block<6,1>(31,0);
    const ::LieAlgebra::Vector6d Delta_dot_eta = t_state.block<6,1>(37,0);


    const ::LieAlgebra::Vector6d Delta_Lambda = t_state.block<6, 1>(43,0);



    const auto ad_Xi = ::LieAlgebra::ad(t_Xi);
    const auto ad_Delta_Xi = ::LieAlgebra::ad(t_Delta_Xi);


    const auto ad_eta = ::LieAlgebra::ad(eta);
    const auto ad_Delta_eta = ::LieAlgebra::ad(Delta_eta);

    const auto M = m_rod_properties->m_M;








    const Eigen::Matrix3d R = Eigen::Quaterniond(t_state[0], t_state[1],t_state[2], t_state[3]).toRotationMatrix();

    Eigen::Vector3d Delta_rotation = Delta_zeta.block<3,1>(0, 0);

    Eigen::Vector3d Delta_N_bar = ::LieAlgebra::skew( Delta_rotation ).transpose()
                                    *R.transpose()
                                    *m_rod_properties->distributedGravitationalForce();

    ::LieAlgebra::Vector6d Delta_F_bar = ::LieAlgebra::Vector6d::Zero();
    Delta_F_bar.block<3,1>(3, 0) = Delta_N_bar;





    const ExplicitIntegrationODEs::TangentKinematicState tangent_kinematic_state_prime = tangentKinematicsStep(tangent_kinematic_state,
                                                                                      t_Xi, t_dot_Xi, t_ddot_Xi,
                                                                                      t_Delta_Xi, t_Delta_dot_Xi, t_Delta_ddot_Xi);
    const Eigen::VectorXd Lambda_prime = getLambdaPrime(Eigen::Quaterniond(t_state[0], t_state[1],t_state[2], t_state[3]),
                                                                           Lambda, eta, dot_eta, ad_Xi, t_X);
    const ::LieAlgebra::Vector6d Delta_Lambda_prime =
            M*Delta_dot_eta - ad_eta.transpose()*M*Delta_eta - ad_Delta_eta.transpose()*M*eta + ad_Xi.transpose()*Delta_Lambda + ad_Delta_Xi.transpose()*Lambda - Delta_F_bar;

    const Eigen::VectorXd Delta_Qa_prime = -t_BPhi.transpose()*Delta_Lambda;


    //  Packing state vector derivative
    Eigen::VectorXd dydx(49 + m_generalised_coordinates_dimension);
    dydx <<  tangent_kinematic_state_prime.block<19,1>(0,0),
             Lambda_prime,
             tangent_kinematic_state_prime.block<18,1>(19,0),
             Delta_Lambda_prime,
             Delta_Qa_prime;


    dydx *= m_rod_length;

    return dydx;
}



void ExplicitIntegrationODEs::distributedActuationODE(const Eigen::VectorXd &,
                                                     Eigen::VectorXd &t_dyds,
                                                     const double t_X)const
{
    /*  Preprocessing    */
    Eigen::MatrixXd BPhi = m_polynomial_representation.m_B*m_polynomial_representation.getPhi(t_X);

    //  Get the strains for the rod
    const Eigen::VectorXd Xi = BPhi*m_qe + m_constant_strain;


    //  Decompose the strain
    const Eigen::Vector3d K = Xi.block<3,1>(0,0);
    const Eigen::Vector3d Gamma = Xi.block<3,1>(3,0);


    t_dyds = m_distributed_actuation.getInternalActuatedStresses(K,
                                                        Gamma,
                                                        BPhi);
}



}   //  namespace CROSP::numerical_integrators::explicit_methods
