#include "PreconditionneurCreux.h"
#include <string>
#include <iostream>

using namespace std;
using namespace Eigen;

//constructeur par défaut
MethodesPreconditionnees::MethodesPreconditionnees(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax)
{
	_A=A;
	_b=b;
	_erreur=erreur;
	_kmax=kmax;
	_x0=x0;
}

//destructeur par défaut
MethodesPreconditionnees::~MethodesPreconditionnees()
{}

void MethodesPreconditionnees::InitialisationFichier(const string nom_fichier)
{
	_fichier.open(nom_fichier);
}

void MethodesPreconditionnees::Ecriture(int k, double erreur)
{
	_fichier << k << " " << erreur << endl;
}

////////////////////////////////////////////////////////////////////////////////
//Sous classe Résidu minimum préconditionné à gauche avec Jacobi
////////////////////////////////////////////////////////////////////////////////


//constructeur
RMPJ::RMPJ(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax) : MethodesPreconditionnees(A, b, x0, erreur, kmax)
{}

//Algorithme du résidu minimum
void RMPJ::AlgoGene()
{
	InitialisationFichier("ResiduMinimumPregaucheJacobi.txt");
	SparseVector<double> r = _b - _A*_x0;
	SparseVector<double> q(r.size());
	int k = 0;
	//Resoudre
	for (int i=0; i<q.rows(); i++)
	{
		q.coeffRef(i) = r.coeffRef(i)/_A.coeffRef(i,i);
	}
	_x = _x0;
	SparseVector<double> z(_b.size()), w(_b.size());
	z.setZero();
	double alpha = 0.;
	Ecriture(0, r.norm());
	while ((r.norm() > _erreur) && (k <= _kmax))
	{
		k += 1;
		w = _A*q;
		//Resoudre
		for (int i=0; i<z.rows(); i++)
		{
			z.coeffRef(i) = w.coeffRef(i)/_A.coeffRef(i,i);
		}
		alpha = q.dot(z)/z.dot(z);
		_x += alpha*q;
		r -= alpha*w;
		q -= alpha*z;
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
//Sous classe Résidu minimum préconditionné à gauche avec SSOR
////////////////////////////////////////////////////////////////////////////////


//constructeur
RMPSSOR::RMPSSOR(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax) : MethodesPreconditionnees(A, b, x0, erreur, kmax)
{}

//Algorithme du résidu minimum
void RMPSSOR::AlgoGene()
{
	InitialisationFichier("ResiduMinimumPregaucheSSOR.txt");
	SparseVector<double> r = _b - _A*_x0;
	SparseVector<double> q(r.size());
	int k = 0;
	//Resoudre
	resoudre(_A,q,r,1.5);

	_x = _x0;
	SparseVector<double> z(_b.size()), w(_b.size());
	z.setZero();
	double alpha = 0.;
	Ecriture(0, r.norm());
	while ((r.norm() > _erreur) && (k <= _kmax))
	{
		k += 1;
		w = _A*q;
		//Resoudre
		resoudre(_A,z,w,1.5);

		alpha = q.dot(z)/z.dot(z);
		_x += alpha*q;
		r -= alpha*w;
		q -= alpha*z;
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

void RMPSSOR::resoudre(SparseMatrix<double> & A, SparseVector<double> & x,SparseVector<double> & b, double w)
{
	SparseVector<double> y(x.size());
	for(int i=0;i<x.size();i++)
	{
		y.coeffRef(i)=b.coeffRef(i);
		for(int j=0;j<i;j++)
		{
			if(_A.coeffRef(i,j)!=0)
			{
				y.coeffRef(i) = y.coeffRef(i)-w*_A.coeffRef(i,j)*y.coeffRef(j)/_A.coeffRef(i,i);
			}
		}
	}

	for(int i=x.size()-1;i>=0;i=i-1)
	{
		x.coeffRef(i)= y.coeffRef(i);
		for(int j=i+1;j<x.size();j++)
		{
			if(_A.coeffRef(i,j)!=0)
			{
				x.coeffRef(i) = x.coeffRef(i)-w*_A.coeffRef(j,i)*x.coeffRef(j);
			}
		}
		x.coeffRef(i) = x.coeffRef(i)/_A.coeffRef(i,i);
	}
}


////////////////////////////////////////////////////////////////////////////////
//Sous classe Résidu minimum préconditionné à droite
////////////////////////////////////////////////////////////////////////////////


//constructeur
RMPD::RMPD(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax) : MethodesPreconditionnees(A, b, x0, erreur, kmax)
{}

//Algorithme du résidu minimum
void RMPD::AlgoGene()
{
	InitialisationFichier("ResiduMinimumPredroiteJacobi.txt");
	SparseVector<double> r = _b - _A*_x0;
	int k = 0;
	_x = _x0;
	SparseVector<double> z(_b.size()), w(_b.size());
	z.setZero();
	double alpha = 0.;
	Ecriture(0, r.norm());
	while ((r.norm() > _erreur) && (k <= _kmax))
	{
		k += 1;
		//Resoudre
		for (int i=0; i<z.rows(); i++)
		{
			z.coeffRef(i) = r.coeffRef(i)/_A.coeffRef(i,i);
		}
		//
		w = _A*z;
		alpha = r.dot(w)/w.dot(w);
		_x += alpha*z;
		r -= alpha*w;
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
//Sous classe Résidu minimum préconditionné à droite version flexible avec SSOR
////////////////////////////////////////////////////////////////////////////////


//constructeur
RMPDflexSSOR::RMPDflexSSOR(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax) : MethodesPreconditionnees(A, b, x0, erreur, kmax)
{}

//Algorithme du résidu minimum
void RMPDflexSSOR::AlgoGene()
{
	InitialisationFichier("ResiduMinimumPredroiteSSORflex.txt");
	SparseVector<double> r = _b - _A*_x0;
	int k = 0;
	_x = _x0;
	SparseVector<double> z(_b.size()), w(_b.size());
	z.setZero();
	double alpha = 0.;
	Ecriture(0, r.norm());
	while ((r.norm() > _erreur) && (k <= _kmax))
	{
		k += 1;
		//Resoudre
		resoudre(_A,z,r,0.5+k%2);
		//
		w = _A*z;
		alpha = r.dot(w)/w.dot(w);
		_x += alpha*z;
		r -= alpha*w;
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

void RMPDflexSSOR::resoudre(SparseMatrix<double> & A, SparseVector<double> & x,SparseVector<double> & b, double w)
{
	SparseVector<double> y(x.size());
	for(int i=0;i<x.size();i++)
	{
		y.coeffRef(i)=b.coeffRef(i);
		for(int j=0;j<i;j++)
		{
			if(_A.coeffRef(i,j)!=0)
			{
				y.coeffRef(i) = y.coeffRef(i)-w*_A.coeffRef(i,j)*y.coeffRef(j)/_A.coeffRef(i,i);
			}
		}
	}

	for(int i=x.size()-1;i>=0;i=i-1)
	{
		x.coeffRef(i)= y.coeffRef(i);
		for(int j=i+1;j<x.size();j++)
		{
			if(_A.coeffRef(i,j)!=0)
			{
				x.coeffRef(i) = x.coeffRef(i)-w*_A.coeffRef(j,i)*x.coeffRef(j);
			}
		}
		x.coeffRef(i) = x.coeffRef(i)/_A.coeffRef(i,i);
	}
}


////////////////////////////////////////////////////////////////////////////////
//Sous classe Résidu minimum préconditionné à droite version flexible avec le Résidu Minimum
////////////////////////////////////////////////////////////////////////////////


//constructeur
RMPDflexResMin::RMPDflexResMin(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, double erreur, int kmax) : MethodesPreconditionnees(A, b, x0, erreur, kmax)
{}

//Algorithme du résidu minimum
void RMPDflexResMin::AlgoGene()
{
	InitialisationFichier("ResiduMinimumPredroiteResMinflex.txt");
	SparseVector<double> r = _b - _A*_x0;
	int k = 0;
	_x = _x0;
	SparseVector<double> z(_b.size()), w(_b.size());
	z.setZero();
	double alpha = 0.;
	Ecriture(0, r.norm());
	while ((r.norm() > _erreur) && (k <= _kmax))
	{
		k += 1;
		//Resoudre
		ResiduMinimum(_A,r,z,10);
		z = GetSol();
		//
		w = _A*z;
		alpha = r.dot(w)/w.dot(w);
		cout << "alpha =\n" << alpha << endl;
		_x += alpha*z;
		r -= alpha*w;
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


//Algorithme du résidu minimum
void RMPDflexResMin::ResiduMinimum(SparseMatrix<double> & A, SparseVector<double> & b, SparseVector<double> & x0, int kmax)
{
	SparseVector<double> r = _b - _A*_x0;
	int k = 0;
	_x = _x0;
	SparseVector<double> z(_b.size());
	z.setZero();
	double alpha = 0.;
	while (k <= _kmax)
	{
		k += 1;
		z = _A*r;
		alpha = r.dot(z)/z.dot(z);
		_x += alpha*r;
		r -= alpha*z;
	}
}
