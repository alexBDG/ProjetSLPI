#include "MethodesCreux.h"
#include <string>
#include <iostream>

using namespace std;
using namespace Eigen;


//constructeur par défaut
Methodes::Methodes(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax)
{
	_A=A;
	_b=b;
	_erreur=erreur;
	_kmax=kmax;
	_x0=x0;
}

//destructeur par défaut
Methodes::~Methodes()
{}

void Methodes::InitialisationFichier(const string nom_fichier)
{
	_fichier.open(nom_fichier);
}

void Methodes::Ecriture(int k, double erreur)
{
	_fichier << k << " " << erreur << endl;
}

////////////////////////////////////////////////////////////////////////////////
//Sous classe Jacobi
////////////////////////////////////////////////////////////////////////////////

//constructeur
Jacobi::Jacobi(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax) : Methodes(A, b, x0, erreur, kmax)
{}

//Algorithme générique
void Jacobi::AlgoGene()
{
	InitialisationFichier("Jacobi.txt");
	SparseVector<double> r = _b - _A*_x0;
	int k = 0;
	_x = _x0;
	SparseMatrix<double> D(_b.size(),_b.size()), Dinv(_b.size(),_b.size());
	D.setZero();
	Dinv.setZero();
	for (int i = 0 ; i < _b.size() ; i++)
	{
		D.coeffRef(i,i) = _A.coeffRef(i,i);
		if (_A.coeffRef(i,i) == 0)
		{
			cout << "------------------------------------" << endl;
			cout << "ERREUR : la matrice A est mauvaise" << endl;
			cout << "------------------------------------" << endl;
			break;
		}
		Dinv.coeffRef(i,i) = 1./_A.coeffRef(i,i);
	}
	Ecriture(0, r.norm());
	while ((r.norm() > _erreur) && (k <= _kmax))
	{
		k += 1;
		_x = Dinv*(D-_A)*_x + Dinv*_b;
		r = _b - _A*_x;
		Ecriture(k, r.norm());
	}
	cout << "La méthode de Jacobi à AX=b converge avec une erreur de : " << r.norm() << endl;
	if (k > _kmax)
	{
		cout << "       /!\ /!\ /!\\" << endl;
		cout << "Tolérance non atteinte : " << r.norm() << endl;
		cout << "       /!\ /!\ /!\\" << endl;
	}
}



////////////////////////////////////////////////////////////////////////////////
//Sous classe Gradient à pas Optimal
////////////////////////////////////////////////////////////////////////////////

//constructeur
GradientOptimal::GradientOptimal(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax) : Methodes(A, b, x0, erreur, kmax)
{}

//Algorithme générique
void GradientOptimal::AlgoGene()
{
	InitialisationFichier("GradientOptimal.txt");
	SparseVector<double> r = _b - _A*_x0;
	int k = 0;
	_x = _x0;
	SparseVector<double> z;
	double alpha;
	Ecriture(0, r.norm());
	while ((r.norm() > _erreur) && (k <= _kmax))
	{
		k += 1;
		z = _A*r;
		alpha = r.dot(r)/(z.dot(r));
		_x += alpha*r;
		r -= alpha*z;
		Ecriture(k, r.norm());
	}
	cout << "La méthode du gradient à pas optimal à AX=b converge avec une erreur de : " << r.norm() << endl;
	if (k > _kmax)
	{
		cout << "       /!\ /!\ /!\\" << endl;
		cout << "Tolérance non atteinte : " << r.norm() << endl;
		cout << "       /!\ /!\ /!\\" << endl;
	}
}


////////////////////////////////////////////////////////////////////////////////
//Sous classe Résidu minimum
////////////////////////////////////////////////////////////////////////////////

//constructeur
ResiduMinimum::ResiduMinimum(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax) : Methodes(A, b, x0, erreur, kmax)
{}

//Algorithme du résidu minimum
void ResiduMinimum::AlgoGene()
{
	InitialisationFichier("ResiduMinimum.txt");
	SparseVector<double> r = _b - _A*_x0;
	int k = 0;
	_x = _x0;
	SparseVector<double> z(_b.size());
	z.setZero();
	double alpha = 0.;
	Ecriture(0, r.norm());
	while ((r.norm() > _erreur) && (k <= _kmax))
	{
		k += 1;
		z = _A*r;
		alpha = r.dot(z)/z.dot(z);
		_x += alpha*r;
		r -= alpha*z;
		Ecriture(k, r.norm());
	}
	cout << "La méthode du résidu minimum à AX=b converge avec une erreur de : " << r.norm() << endl;
	if (k > _kmax)
	{
		cout << "       /!\ /!\ /!\\" << endl;
		cout << "Tolérance non atteinte : " << r.norm() << endl;
		cout << "       /!\ /!\ /!\\" << endl;
	}
}


////////////////////////////////////////////////////////////////////////////////
//Sous classe GMRes
////////////////////////////////////////////////////////////////////////////////

//constructeur
GMRes::GMRes(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax) : Methodes(A, b, x0, erreur, kmax)
{}

//Algorithme d'Arnoldi (Gram-Schmidt) - version naïve + méthode de G

void GMRes::arnoldi(SparseMatrix<double> & A, SparseVector<double> & b, SparseMatrix<double> & V, SparseMatrix<double> & H)
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
  SparseMatrix<double> V_barre(n,m+1);
  SparseMatrix<double> H_barre(m+1,m);
  V_barre.setZero();
  H_barre.setZero();

  // Création d'un vecteur temporaire pour stoker A*V_barre.cols(i-1)
  SparseVector<double> AV_barre(A.rows());
  AV_barre.setZero();

  // On met le premier vecteur de la base dans V_barre
  V_barre.col(0) = b/b.norm();

  for(int i=1;i<V_barre.cols();i++){

    // On commence a créer le vecteur i à partir de A*le vecteur i-1
    AV_barre = A*V_barre.col(i-1);
    V_barre.col(i) = AV_barre;

    // On othogonalise le vecteur avec tous les vecteurs précédents en stockant au passage
    // les différents produits scalaires dans H_barre
    for(int j = 0;j<i;j++){
      H_barre.coeffRef(j,i-1) = AV_barre.dot(V_barre.col(j));
      V_barre.col(i) = V_barre.col(i) -H_barre.coeffRef(j,i-1)*V_barre.col(j);
    }
    // Calcul du dernier produit scalaire qui correspond en fait à la norme (pas la norme au carré)
    // c'est le produit scalaire (vi|A*vi-1) avec vi de norme 1 ... magie
    H_barre.coeffRef(i,i-1) = V_barre.col(i).norm();
    V_barre.col(i) = V_barre.col(i)/H_barre.coeffRef(i,i-1);



  }

  V = V_barre.block(0,0,V.rows(),V.cols());
  H = H_barre.block(0,0,H.rows(),H.cols());
}

//Méthode de givens
void GMRes::givens_heissen(SparseMatrix<double> &H, SparseMatrix<double> &Q, SparseMatrix<double> &R){

  Q.setIdentity();
  //cout << " Q1 = \n" << Q<< endl;
  R = H;
  double tmp(0.);

  // Déclaration des 2 coeficient de la rotation
  double c(0.),s(0.);

  // Boucle pour enlever les termes sous diagonaux
  for(int i=0;i<m-1;i++){
    c=R.coeffRef(i,i)/sqrt(R.coeffRef(i,i)*R.coeffRef(i,i)+R.coeffRef(i+1,i)*R.coeffRef(i+1,i));
    s=R.coeffRef(i+1,i)/sqrt(R.coeffRef(i,i)*R.coeffRef(i,i)+R.coeffRef(i+1,i)*R.coeffRef(i+1,i));
    //cout << " Q = \n" << i << endl<< Q<< endl;
    for(int j=0;j<Q.rows();j++){
      tmp = Q.coeffRef(j,i);
      Q.coeffRef(j,i) = Q.coeffRef(j,i)*c+Q.coeffRef(j,i+1)*s;
      Q.coeffRef(j,i+1) = -s*tmp + Q.coeffRef(j,i+1)*c;
      tmp = R.coeffRef(i,j);
      R.coeffRef(i,j) = c*R.coeffRef(i,j) + s*R.coeffRef(i+1,j);
      R.coeffRef(i+1,j) = -s*tmp + c*R.coeffRef(i+1,j);
    }
  }

}

//Calcul de Argmin
void GMRes::ArgMin(double beta, SparseMatrix<double> & Q, SparseVector<double> & y, SparseMatrix<double> & R)
{

	for (int i = m-1 ; i >= 0  ; i--)
	{
		y.coeffRef(i) =beta*Q.coeffRef(0,i);
		for (int j = i+1 ; j < m ; j++)
		{
			y.coeffRef(i) -= y.coeffRef(j)*R.coeffRef(i,j);
		}
		y.coeffRef(i) = y.coeffRef(i)/R.coeffRef(i,i);
	}
}


//Algorithme de GMRes
void GMRes::AlgoGene()
{
	//tailles pour les Matrices
	cout << "------------------------------------" << endl;
	cout << "Valeur du paramètre de GMRes, m = " << endl;
	cout << "------------------------------------" << endl;
	cin >> m;

	//dimensionnement
	_Vm.resize(n,m+1);
	_H.resize(m+1,m);
	_y.resize(m);
	_R.resize(m+1,m);

	//variables d'Arnoldi
	_H.setZero();
	_Vm.setZero();
	_R.setZero();
	_y.setZero();
	SparseVector<double> _g(m+1);
	_g.setZero();
	_g.coeffRef(0) = 1.;

	//variables de la résolution
	SparseVector<double> _e(m+1);
	_e.setZero();
	_e.coeffRef(0) = 1.;

	//redimentionnement
	SparseMatrix<double> _R_(m,m), _H_(m,m), _Vm_(n,m), _Q_(m,m);
	_R_.setZero();
	_H_.setZero();
	_Vm_.setZero();
	_Q_.setZero();
	SparseVector<double> _g_(m);

	InitialisationFichier("GMRes.txt");
	SparseVector<double> r = _b - _A*_x0;
	double beta = r.norm();
	int k = 0;
	_x = _x0;
	Ecriture(0, r.norm());

	SparseVector<double> e1(m);
	e1.coeffRef(0) = 1.;

	while ((beta > _erreur) && (k <= _kmax))
	{
		cout.flush();
		cout << "Progression générale : " << (double)k/(double)_kmax*100 << "%" << " \r";
		k += 1;
		//Arnoldi(r,_g,_Vm,_H,_R);
		arnoldi(_A,r,_Vm_,_H_);

		givens_heissen(_H_,_Q_,_R_);
		// _R_ = _R.block(0,0,m,m);
		// _Vm_ = _Vm.block(0,0,n,m);
		e1.coeffRef(0) = 1.;

		ArgMin(beta,_Q_,_y,_R_);

	//	cout << " R = \n" << _R_ << endl;
	//	cout << "y = \n" << _y << endl;
	//	cout << "z = \n" << _g_ << endl;


		_x +=  _Vm_*_y;
		cout << "x = \n" << _x << endl;

		r = _b - _A*_x;
		Ecriture(k, r.norm());
		_g = _e;
		_y.setZero();
		beta = r.norm();

	}
	cout << "La méthode GMRes à AX=b converge avec une erreur de : " << beta << endl;
	if (k > _kmax)
	{
		cout << "       /!\ /!\ /!\\" << endl;
		cout << "Tolérance non atteinte : " << r.norm() << endl;
		cout << "       /!\ /!\ /!\\" << endl;
	}
}
