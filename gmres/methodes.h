#ifndef METHODES_H

#include <vector>
#include <Eigen/Dense>

void initMat(Eigen::MatrixXd & A, double alpha);
void initb(Eigen::VectorXd &b);
void arnoldi(Eigen::MatrixXd &A, Eigen::VectorXd &b, Eigen::MatrixXd &V, Eigen::MatrixXd &H);
double normeTordue(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);
void givens_heissen(Eigen::MatrixXd &H, Eigen::MatrixXd &Q, Eigen::MatrixXd &R);
void argmin(Eigen::MatrixXd &Q, Eigen::MatrixXd &R, double beta, Eigen::VectorXd &y);
void test(Eigen::MatrixXd V,Eigen::MatrixXd H,Eigen::MatrixXd Q,Eigen::MatrixXd R); 
#define METHODES_H
#endif
