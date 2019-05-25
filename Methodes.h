#include "Dense"
#include <fstream>

class Methodes
{
	protected:
		//matrice d'entrée
		Eigen::MatrixXd _A;
		//vecteur d'entrée
		Eigen::VectorXd _b;
		//vecteur solution
		Eigen::VectorXd _x;
		//vecteur initial
		Eigen::VectorXd _x0;
		//valeur d'erreur minimale
		double _erreur;
		//valeur d'itération maximale
		int _kmax;
		//Ecriture du fichier
		std::ofstream _fichier;
	public:
		//constructeur
		Methodes(Eigen::MatrixXd & A, Eigen::VectorXd & b, Eigen::VectorXd & x0, double erreur, int kmax);
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
		Jacobi(Eigen::MatrixXd & A, Eigen::VectorXd & b, Eigen::VectorXd & x0, double erreur, int kmax);
		//algortithme générique d'une méthode itérative - version améliorée
		void AlgoGene();
};

class GradientOptimal : public Methodes
{
	public:
		//constructeur
		GradientOptimal(Eigen::MatrixXd & A, Eigen::VectorXd & b, Eigen::VectorXd & x0, double erreur, int kmax);
		//algortithme générique du gradient à pas optimal
		void AlgoGene();
};

class ResiduMinimum : public Methodes
{
	public:
		//constructeur
		ResiduMinimum(Eigen::MatrixXd & A, Eigen::VectorXd & b, Eigen::VectorXd & x0, double erreur, int kmax);
		//algorithme du résidu minimum préconditionné à gauche
		void AlgoGene();

};

class GMRes : public Methodes
{
	private:
		//Matrices et Vecteurs pour GMRes
		Eigen::MatrixXd _Vm;
		Eigen::MatrixXd _H;
		Eigen::VectorXd _y;
		//Vecteur des cos et sin de rotation
		Eigen::VectorXd _Cos;
		Eigen::VectorXd _Sin;
		//Vecteur de la dernière colonne de H
		Eigen::VectorXd _Hc;
		//Décomposition QR avec g=Q*beta*e1
		Eigen::MatrixXd _R;
		//pour GMRes
	//	Eigen::VectorXd _g;
	public:
		//méthode d'Arnoldi
		void Arnoldi(Eigen::VectorXd & v, Eigen::VectorXd & g, Eigen::MatrixXd & Vm, Eigen::MatrixXd & H, Eigen::MatrixXd & R);
		//fonction qui calcule le produit d'une matrice rotation par un vecteur
		void rotation(int i, Eigen::VectorXd & v);
		//fonction qui calcule _y minimum
		void ArgMin(double beta, Eigen::VectorXd & g, Eigen::VectorXd & y, Eigen::MatrixXd & R);
		//constructeur
		GMRes(Eigen::MatrixXd & A, Eigen::VectorXd & b, Eigen::VectorXd & x0, double erreur, int kmax);
		//Destructeur
		virtual ~GMRes();
		//algorithme du GRMes
		void AlgoGene();
	//	Eigen::VectorXd Getg(){return _g;};
};
