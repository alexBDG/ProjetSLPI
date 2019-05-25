#include "Dense"
#include "Sparse"
#include <fstream>

class MethodesPreconditionnees
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
		MethodesPreconditionnees(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);
		//destructeur
		virtual ~MethodesPreconditionnees();
		// Initialiser le nom du fichier solution
		void InitialisationFichier(const std::string nom_fichier);
		// Sauvegarde la solution
		void Ecriture(int k, double erreur);
		//algortithme générique de la méthode appelée
		virtual void AlgoGene() = 0;
};

class RMPJ: public MethodesPreconditionnees
{
	public:
		//constructeur
		RMPJ(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);

		//algortithme générique d'une méthode itérative
		void AlgoGene();
};

class RMPSSOR: public MethodesPreconditionnees
{
	public:
		//constructeur
		RMPSSOR(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);

		//algortithme générique d'une méthode itérative
		void AlgoGene();

		void resoudre(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & x,Eigen::SparseVector<double> & b, double w);
};

class RMPD: public MethodesPreconditionnees
{
	public:
		//constructeur
		RMPD(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);

		//algortithme générique d'une méthode itérative
		void AlgoGene();
};

class RMPDflexSSOR: public MethodesPreconditionnees
{
	public:
		//constructeur
		RMPDflexSSOR(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);

		//algortithme générique d'une méthode itérative
		void AlgoGene();

		void resoudre(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & x,Eigen::SparseVector<double> & b, double w);
};

class RMPDflexResMin: public MethodesPreconditionnees
{
	public:
		//constructeur
		RMPDflexResMin(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, double erreur, int kmax);

		//algortithme générique d'une méthode itérative
		void AlgoGene();

		void ResiduMinimum(Eigen::SparseMatrix<double> & A, Eigen::SparseVector<double> & b, Eigen::SparseVector<double> & x0, int kmax);
		Eigen::SparseVector<double> _x;
		Eigen::SparseVector<double> GetSol(){return _x;};
};
