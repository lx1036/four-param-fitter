#include<iostream>
#include<ctime>
using namespace std;
const double alpha = -1.0;     // reflection coefficient
const  double beta = 0.5;   // contraction coefficient
const  double gamma = 2.0;      // expansion coefficient
const double root2 = 1.414214; // square root of 2
const  double maxError = 1e-10; 

class CurveFitter {      

	const static int IterFactor = 500;
	const static int Len=5;

	const static int defaultRestarts = 2; 
	const static int getNumParams=4; // paremeter number
	double xData[Len], yData[Len];  // x,y data to fit
	int numPoints;          // number of data points
	int numParams;          // number of parametres
	int numVertices;        // numParams+1 (includes sumLocalResiduaalsSqrd)
	int worst;          // worst current parametre estimates
	int nextWorst;      // 2nd worst current parametre estimates
	int best;           // best current parametre estimates
	double simp[getNumParams+1][getNumParams+1];        // the simplex (the last element of the array at each vertice is the sum of the square of the residuals)
	double next[getNumParams+1];      // new vertex to be tested
	int numIter;        // number of iterations so far
	int maxIter;    // maximum number of iterations per restart
	int restarts;   // number of times to restart simplex after first soln.
	// default number of restarts
	int nRestarts;  // the number of restarts that occurred
	// maximum error tolerance
//	double initialParams[Len];  // user specified initial parameters
	time_t times;  //elapsed time in ms

public:
	CurveFitter (double x[], double y[]);
	~CurveFitter();
	void doFit() ;
//	void initialize();
	void restart(int n);
	double f( double p[], double x);
	double *getParams();
	double *getResiduals();
	double getSumResidualsSqr();
	double  getSD();
	double  getRSquared() ;
	void print();
	double sqr(double d);
	void  sumResiduals (double x[]) ;
	void newVertex();
	void order();
	int  getIterations() ;
	int getMaxIterations();
};