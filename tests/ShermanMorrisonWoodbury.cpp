#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "utilities/Eigen/eigen_io.hpp"

#include <benchmark/benchmark.h>

using namespace std;
using namespace Eigen;

constexpr int blockSize = 20;
constexpr int matrixSize = 60;
constexpr int numBlocks = 3;
constexpr int numUpdates = 6; // off-diagonal blocks
constexpr int updateRankPerBlock = blockSize; // 20 per block
constexpr int totalRank = numUpdates * updateRankPerBlock; // 120 columns in U and V


int main(int argc, char *argv[])
{
    // Load full matrix A from CSV
//    MatrixXd A = loadCSV("data/A.csv");

    MatrixXd A;
    LoadEigenMatrixFromFile(A, "A.csv", "tests/data/");

    // Extract the 20x20 diagonal block D from top-left
    MatrixXd D = A.block(0, 0, blockSize, blockSize);

    // Construct block-diagonal matrix A0 = blkdiag(D, D, D)
    MatrixXd A0 = MatrixXd::Zero(matrixSize, matrixSize);
    for (int i = 0; i < numBlocks; ++i) {
        A0.block(i * blockSize, i * blockSize, blockSize, blockSize) = D;
    }

    // Compute ΔA = A - A0
    MatrixXd DeltaA = A - A0;


    cout << "Rank DeltaA : " << ColPivHouseholderQR<MatrixXd>(DeltaA).rank() << "\n";

//    // Build U and V from ΔA using block-wise decomposition
//    MatrixXd U(matrixSize, totalRank);
//    MatrixXd V(matrixSize, totalRank);
//    int colIndex = 0;

//    for (int i = 0; i < numBlocks; ++i) {
//        for (int j = 0; j < numBlocks; ++j) {
//            if (i == j) continue; // skip diagonal

//            MatrixXd Bij = DeltaA.block(i * blockSize, j * blockSize, blockSize, blockSize);

//            // Decompose each 20x20 block into 20 rank-1 updates
//            for (int k = 0; k < blockSize; ++k) {
//                VectorXd u = VectorXd::Zero(matrixSize);
//                VectorXd v = VectorXd::Zero(matrixSize);

//                // u holds the k-th row of Bij placed in the correct global row block
//                u.segment(i * blockSize, blockSize) = Bij.row(k).transpose();
//                // v is a unit vector with 1 at the corresponding global column
//                v(j * blockSize + k) = 1.0;

//                U.col(colIndex) = u;
//                V.col(colIndex) = v;
//                ++colIndex;
//            }
//        }
//    }

//    constexpr int totalRank_v2 = 6; // One per off-diagonal block
//    MatrixXd U = MatrixXd::Zero(matrixSize, totalRank_v2);
//    MatrixXd V = MatrixXd::Zero(matrixSize, totalRank_v2);

//    int colIndex = 0;
//    for (int i = 0; i < numBlocks; ++i) {
//        for (int j = 0; j < numBlocks; ++j) {
//            if (i == j) continue;

//            MatrixXd Bij = DeltaA.block(i * blockSize, j * blockSize, blockSize, blockSize);

//            // Construct u_k from row-sums of block B_{ij}
//            VectorXd u = VectorXd::Zero(matrixSize);
//            u.segment(i * blockSize, blockSize) = Bij.rowwise().sum();

//            // Construct v_k as a flat weight vector (e.g., ones)
//            VectorXd v = VectorXd::Zero(matrixSize);
//            v.segment(j * blockSize, blockSize).setOnes();

//            U.col(colIndex) = u;
//            V.col(colIndex) = v;
//            ++colIndex;
//        }
//    }

    writeToFile("DeltaA", DeltaA, "tests/data");


    //  IMPLEMENTATION FROM GPT WITH SVD DECOMPOSITION
    JacobiSVD<MatrixXd> svd(DeltaA, ComputeThinU | ComputeThinV);

    int r = svd.rank();

    cout << "rank sVD : " << r << endl;

    VectorXd sigma = svd.singularValues().head(r);
    MatrixXd sqrtSigma = sigma.cwiseSqrt().asDiagonal();

    MatrixXd U = svd.matrixU().leftCols(r) * sqrtSigma;
    MatrixXd V = svd.matrixV().leftCols(r) * sqrtSigma;

    cout << "U size: " << U.rows() << " x " << U.cols() << endl;
    cout << "V size: " << V.rows() << " x " << V.cols() << endl;

    writeToFile("U", U, "tests/data");
    writeToFile("V", V, "tests/data");



    // Verify ΔA ≈ UVᵗ
    MatrixXd DeltaA_check = U * V.transpose();
    double reconstruction_error = (DeltaA - DeltaA_check).norm();
    cout << "Frobenius norm of reconstruction error: " << reconstruction_error << endl;

    // ---- Solve Ax = b using SMW ----
    VectorXd b = VectorXd::Random(matrixSize);

    // Compute A0 inverse blockwise
    MatrixXd A0inv = MatrixXd::Zero(matrixSize, matrixSize);
    MatrixXd Dinv = D.inverse();
    for (int i = 0; i < numBlocks; ++i) {
        A0inv.block(i * blockSize, i * blockSize, blockSize, blockSize) = Dinv;
    }

    // SMW formula
    VectorXd A0inv_b = A0inv * b;
    MatrixXd A0inv_U = A0inv * U;
    MatrixXd M = MatrixXd::Identity(r, r) + V.transpose() * A0inv_U;
    writeToFile("M", M, "tests/data");
//    MatrixXd M_inv = M.inverse();
    VectorXd x = A0inv_b - A0inv_U * M.householderQr().solve(V.transpose() * A0inv_b);

    // Direct solve for validation
    VectorXd x_direct = A.householderQr().solve(b);
    double rel_error = (x - x_direct).norm() / x_direct.norm();
    cout << "Relative error between SMW and direct solve: " << rel_error << endl;

    ::benchmark::RegisterBenchmark("householderQr", [&](::benchmark::State &t_state){
        while(t_state.KeepRunning())
            x_direct = A.householderQr().solve(b);
    });

//    ::benchmark::RegisterBenchmark("SMW", [&](::benchmark::State &t_state){
//        while(t_state.KeepRunning()){
//            A0inv_b = A0inv * b;
//            A0inv_U = A0inv * U;
//            M = MatrixXd::Identity(r, r) + V.transpose() * A0inv_U;
//            M_inv = M.inverse();
//            x = A0inv_b - A0inv_U * M_inv * (V.transpose() * A0inv_b);
//        }
//    });


    ::benchmark::RegisterBenchmark("SMW + QR", [&](::benchmark::State &t_state){
        while(t_state.KeepRunning()){
            A0inv_b = A0inv * b;
            A0inv_U = A0inv * U;
            M = MatrixXd::Identity(r, r) + V.transpose() * A0inv_U;
//            M_inv = M.inverse();
            x = A0inv_b - A0inv_U * M.householderQr().solve(V.transpose() * A0inv_b);
        }
    });


    ::benchmark::RegisterBenchmark("SVD + SMW + QR", [&](::benchmark::State &t_state){
        while(t_state.KeepRunning()){

            svd  = JacobiSVD<MatrixXd>(DeltaA, ComputeThinU | ComputeThinV);

            sigma = svd.singularValues().head(r);
            sqrtSigma = sigma.cwiseSqrt().asDiagonal();

            U = svd.matrixU().leftCols(r) * sqrtSigma;
            V = svd.matrixV().leftCols(r) * sqrtSigma;

            A0inv_b = A0inv * b;
            A0inv_U = A0inv * U;
            M = MatrixXd::Identity(r, r) + V.transpose() * A0inv_U;
//            M_inv = M.inverse();
            x = A0inv_b - A0inv_U * M.householderQr().solve(V.transpose() * A0inv_b);
        }
    });



//    ::benchmark::RegisterBenchmark("SMW + QR no A0inv_b", [&](::benchmark::State &t_state){
//        while(t_state.KeepRunning()){
////            A0inv_b = A0inv * b;
//            A0inv_U = A0inv * U;
//            M = MatrixXd::Identity(r, r) + V.transpose() * A0inv_U;
////            M_inv = M.inverse();
//            x = A0inv_b - A0inv_U * M.householderQr().solve(V.transpose() * A0inv_b);
//        }
//    });

//    ::benchmark::RegisterBenchmark("SMW + QR no A0inv_U", [&](::benchmark::State &t_state){
//        while(t_state.KeepRunning()){
//            A0inv_b = A0inv * b;
////            A0inv_U = A0inv * U;
//            M = MatrixXd::Identity(r, r) + V.transpose() * A0inv_U;
////            M_inv = M.inverse();
//            x = A0inv_b - A0inv_U * M.householderQr().solve(V.transpose() * A0inv_b);
//        }
//    });


//    ::benchmark::RegisterBenchmark("SMW + QR no M", [&](::benchmark::State &t_state){
//        while(t_state.KeepRunning()){
//            A0inv_b = A0inv * b;
//            A0inv_U = A0inv * U;
////            M = MatrixXd::Identity(r, r) + V.transpose() * A0inv_U;
////            M_inv = M.inverse();
//            x = A0inv_b - A0inv_U * M.householderQr().solve(V.transpose() * A0inv_b);
//        }
//    });


//    ::benchmark::RegisterBenchmark("SMW + QR no A0inv_b", [&](::benchmark::State &t_state){
//        while(t_state.KeepRunning()){
//            A0inv_b = A0inv * b;
//            A0inv_U = A0inv * U;
//            M = MatrixXd::Identity(r, r) + V.transpose() * A0inv_U;
////            M_inv = M.inverse();
//            x = A0inv_b - A0inv_U * M.householderQr().solve(V.transpose() * A0inv_b);
//        }
//    });





    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();

    return 0;
}


