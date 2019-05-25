#include "Dense"
#include "Sparse"
#include <fstream>

class Methodes
{
	protected:
		//matrice d'entrée
		Eigen::SparseMatrix<double> _A;
		//vecteur d'entrée
		Eigen::SparseVector<double> _b;
		//vecteur solution
		Eigen::SparseVector<double> _x;
		//vecteur initial
		Eigen::SparseVector<double> _x0;
		//valeur d'erreur minimale
		double _erreur;
		//valeur d'itération maximale
		int _kmax;
		//Ecriture du fichier
		std::ofstream _fichier;
	public:
		//constructeur
		Methodes(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);
		//destructeur
		virtual ~Methodes();
		// Initialiser le nom du fichier solution
		void InitialisationFichier(const std::string nom_fichier);
		// Sauvegarde la solution
		void Ecriture(int k, double erreur);
		//algortithme générique de la méthode appelée
		virtual void AlgoGene() = 0;
};

class Jacobi : public Methodes
{
	public:
		//constructeur
		Jacobi(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);
		//algortithme générique d'une méthode itérative - version améliorée
		void AlgoGene();
};

class GradientOptimal : public Methodes
{
	public:
		//constructeur
		GradientOptimal(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);
		//algortithme générique du gradient à pas optimal
		void AlgoGene();
};

class ResiduMinimum : public Methodes
{
	public:
		//constructeur
		ResiduMinimum(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);
		//algorithme du résidu minimum
		void AlgoGene();
};

class GMRes : public Methodes
{
	private:
		int n = _b.size();
		int m=0;
		//Matrices et Vecteurs pour GMRes
		Eigen::SparseMatrix<double> _Vm;
		Eigen::SparseMatrix<double> _H;
		Eigen::SparseVector<double> _y;
		//Vecteur des cos et sin de rotation
		Eigen::SparseVector<double> _Cos;
		Eigen::SparseVector<double> _Sin;
		//Vecteur de la dernière colonne de H
		Eigen::SparseVector<double> _Hc;
		//Décomposition QR avec g=Q*beta*e1
		Eigen::SparseMatrix<double> _R;
		//pour GMRes
	//	Eigen::VectorXd _g;
	public:
		//constructeur
		GMRes(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);
		void givens_heissen(Eigen::SparseMatrix<double> &H, Eigen::SparseMatrix<double> &Q, Eigen::SparseMatrix<double> &R);
		//méthode d'Arnoldi
		void arnoldi(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseMatrix<double> & V, Eigen::SparseMatrix<double> & H);
		//fonction qui calcule _y minimum
		void ArgMin(double beta, Eigen::SparseMatrix<double> & Q, Eigen::SparseVector<double> & y, Eigen::SparseMatrix<double> & R);
		//algorithme du GRMes
		void AlgoGene();
	//	Eigen::VectorXd Getg(){return _g;};
};
