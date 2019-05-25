#include <vector>
#include "methodes.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;

int main(){
  const int n=6;
  const int m=3;
  const double alpha=8.;

  // Définition de la matrice A du systeme linéaire de base de taille n*n
  MatrixXd A = MatrixXd::Random(n,n);
  for(int i=0;i<n;i++){
    A(i,i) += alpha;
  }
  // Définition du vecteur b du systeme linéaire de base de taille n
  VectorXd b = VectorXd::Random(n);

	cout << "------------------------------------" << endl;
	cout << "Affichage de la matrice A" << endl;
	cout << "------------------------------------" << endl;
	cout << "A = \n" << A << endl;

	cout << "------------------------------------" << endl;
	cout << "Affichage du vecteur b" << endl;
	cout << "------------------------------------" << endl;
	cout << "b = \n" << b << endl;




  // Définition de la matrice V qui contient les vecteurs de la base de K orthonormalisé
  MatrixXd V = MatrixXd::Zero(n,m);

  // Définition de la matrice H tq, H(i,j) = (Avj,vi)
  MatrixXd H = MatrixXd::Zero(m,m);
  // Définition des matrices Q et R tq H=QR
  MatrixXd Q(m,m), R(m,m);

  // Vecteur y=argmin ...
  VectorXd y = VectorXd::Zero(m);

  // Beta = norme de r a chaque étape
  double beta(0.);

  // Déclaration et initialisation du vecteur x
  VectorXd x = VectorXd::Zero(n);

  // Déclaration et initialisation du résidu
  VectorXd r = VectorXd::Zero(n);
  r = b-A*x;
  beta = r.norm();

  VectorXd e1 = VectorXd::Zero(m);
  e1(0)=1.;

  ofstream flux("gmres.dat");

  int k=0;
  while((k<50)and(beta>0.00000000000001)){

    // Algorithme d'arnoldi sur V avec la matrice A et le vecteur b
    arnoldi(A,r,V,H);
    // Décomposition QR de H grace à la méthode de givens optimisé pour les matrice d'Heissenberg


    givens_heissen(H,Q,R);
    //test(V,H,Q,R);

		cout << "--------------------------------------\n" << endl;
		cout << "Q = \n" << Q << endl;
		cout << "--------------------------------------\n" << endl;

		cout << "--------------------------------------\n" << endl;
		cout << "H = \n" << H << endl;
		cout << "--------------------------------------\n" << endl;

		cout << "--------------------------------------\n" << endl;
		cout << "R = \n" << R << endl;
		cout << "--------------------------------------\n" << endl;

    cout << "--------------------------------------\n" << endl;
    cout << "Q * H - R = \n" << H-Q*R << endl;
    cout << "--------------------------------------\n" << endl;

    // Résolution de argmin de .. pour avoir y
    argmin(Q,R,beta,y);
    //cout << endl << "argmin " << (beta*e1-H*y).norm() << endl;

		cout << "--------------------------------------\n" << endl;
		cout << "y = \n" << y << endl;
		cout << "--------------------------------------\n" << endl;

    // On actualise x
    x = x + V*y;
    // On actualise r et beta
    r = b-A*x;
    beta = r.norm();
    flux << k << " " << beta << endl;

    k+=1;
  }

  cout << endl << endl << "----------" << endl;
  cout << "Après " << k << " itérations, r.norm() = "<<  r.norm() << endl << endl;
  cout << "On obtient x =  " << endl << x << endl << "--------------";

  return 0;
}
