#ifndef METHODES_CPP

#include "methodes.h"
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <assert.h>

using namespace std;
using namespace Eigen;

void givens_heissen(MatrixXd &H, MatrixXd &Q, MatrixXd &R){

  int const m = H.rows();
  Q = MatrixXd::Identity(m,m);
  R = H;

  // Déclaration de la matrice de rotation à l'étape i
  MatrixXd rot_i(m,m);

  // Déclaration des 2 coeficient de la rotation
  double c(0.),s(0.);

  // Boucle pour enlever les termes sous diagonaux
  for(int i=0;i<m-1;i++){
    rot_i = MatrixXd::Identity(m,m);

    c=R(i,i)/sqrt(R(i,i)*R(i,i)+R(i+1,i)*R(i+1,i));
    s=R(i+1,i)/sqrt(R(i,i)*R(i,i)+R(i+1,i)*R(i+1,i));

    rot_i(i,i) = c;
    rot_i(i+1,i) = -s;
    rot_i(i,i+1) = s;
    rot_i(i+1,i+1) = c;


    // cout << endl << "---------- ieme Rotation" << endl;
    // cout << rot_i;
    // cout << endl << "----------" << endl;
    Q = Q*rot_i.transpose();
    R = rot_i*R;
    /*
    cout << "--------------------------------------\n" << endl;
    cout << "R = \n" << R << endl;
    cout << "--------------------------------------\n" << endl;

    cout << "--------------------------------------\n" << endl;
    cout << "Q = \n" << Q << endl;
    cout << "--------------------------------------\n" << endl;
    cout << "\n Matrice de rotation = \n" << rot_i << endl;
    */
  }

}

void argmin(MatrixXd &Q, MatrixXd &R, double beta, VectorXd &y)
{

  // Résolution de R*y = Q.transpose()*beta*e1
  // avec R triangulaire supérieure et e1 =(1,0,0,0...)
  int const m = R.cols();
  VectorXd e1 = VectorXd::Zero(m);
  e1(0) = 1.;

  // On met dans e1 Q.transpose()*beta*e1 pour avoir un systeme R*y= e1
  e1 = beta*Q.transpose()*e1;

  for(int i=m-1;i>=0;i--){
    y(i) = e1(i);
    for(int j=i+1;j<m;j++){
      y(i) -= y(j)*R(i,j);
    }
    y(i)=y(i)/R(i,i);
  }

}

void arnoldi(MatrixXd &A, VectorXd &b, MatrixXd &V, MatrixXd &H)
{
  // On vérifie les tailles des matrices d'entrée
  assert (A.cols() == A.rows());
  assert (A.rows() == b.size());
  assert (H.rows() == H.cols());
  assert ((V.rows() == A.cols()) and (V.cols()==H.cols()));

  // Déclartion des constantes
  int const n = V.rows();
  int const m = V.cols();

  // On crée les matrice H_barre de taille m,m+1
  // et  V_barre de taille m+1,m
  MatrixXd V_barre = MatrixXd::Zero(n,m+1);
  MatrixXd H_barre = MatrixXd::Zero(m+1,m);

  // Création d'un vecteur temporaire pour stoker A*V_barre.cols(i-1)
  VectorXd AV_barre = VectorXd::Zero(A.rows());

  // On met le premier vecteur de la base dans V_barre
  V_barre.col(0) = b/b.norm();

  for(int i=1;i<V_barre.cols();i++){

    // On commence a créer le vecteur i à partir de A*le vecteur i-1
    AV_barre = A*V_barre.col(i-1);
    V_barre.col(i) = AV_barre;

    // On othogonalise le vecteur avec tous les vecteurs précédents en stockant au passage
    // les différents produits scalaires dans H_barre
    for(int j = 0;j<i;j++){
      H_barre(j,i-1) = AV_barre.dot(V_barre.col(j));
      V_barre.col(i) = V_barre.col(i) -H_barre(j,i-1)*V_barre.col(j);
    }
    // Calcul du dernier produit scalaire qui correspond en fait à la norme (pas la norme au carré)
    // c'est le produit scalaire (vi|A*vi-1) avec vi de norme 1 ... magie
    H_barre(i,i-1) = V_barre.col(i).norm();
    V_barre.col(i) = V_barre.col(i)/H_barre(i,i-1);



  }

  V = V_barre.block(0,0,V.rows(),V.cols());
  H = H_barre.block(0,0,H.rows(),H.cols());
}

void test(MatrixXd V,MatrixXd H,MatrixXd Q,MatrixXd R){

  cout << endl << "--------V--------" << endl;
  cout << V << endl;
  cout << endl << " Produit scalaire de V(:,-2) avec V(:,-1) pour tester l'orthogonalité : " << V.col(V.cols()-2).dot(V.col(V.cols()-1)) << endl << endl;
  cout << endl << "--------H--------" << endl;
  cout << H << endl;
  cout << endl << "--------Q--------" << endl;
  cout << Q << endl;
  cout << endl << "--------R--------" << endl;
  cout << R << endl;
  cout << endl << "--------QR--------" << endl;
  cout << Q*R << endl;
  cout << endl << "------------------" << endl;
}

#define METHODES_CPP
#endif
