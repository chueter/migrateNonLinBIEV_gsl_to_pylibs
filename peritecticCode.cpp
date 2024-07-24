#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <list>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string.h>

#define numberGridPoints  90 
#define addGridPoints 10
//numberGridPoints is #of points when doing ONE grid from 0 to +inf.  
#define numberGridPoints_y 90
// numberGridPoints_y is number of points where x(y) solver; y=0-> x=0to save lambda, which is x(yMax) and which is not solved by generic equation but by saving at x(y=0) in x-solutionVector 

using namespace std;

// I keep the parameters strucs though they are redundant now

struct parameters{
double peclet1;
double peclet2;
double delta1;
double delta2;
double slope1;	// in entire code, suffix 1 means that the quantity is considered on the right branch.
double slope2;
};

struct integralParameters{
double sigma;
double xObserv;
double yObserv;
gsl_spline* splineG1loc;
gsl_interp_accel* accelG1loc;
gsl_spline* kappaRHP_splineG1;
gsl_interp_accel* kappaRHP_accelG1;
//gsl_spline* kappa_splineG12;
//gsl_interp_accel* kappa_accelG12;	
};

struct integralParameters_l1l2Interface{
double sigma;
double xObserv;
double yObserv;
gsl_spline* splineG12loc;
gsl_interp_accel* accelG12loc;
gsl_spline* kappa_splineG12;
gsl_interp_accel* kappa_accelG12;	
};
////////////////////////////////// testing variables
//bool shifted(0);
bool testedIntegrandAtxObs(0);
/////////////////////////////////// variables of standard types  ///////////////////////////
const  double PI=4.*atan(1.0);
//double Dc1S_Dc12(100.), Dc2S_Dc12(100.); // ratio of miscibility gaps
double DcLDelta_DcLGamma(2.);


double lambda(1e5),gammaD(1e5),pG1(1e5),pG2(1e5),deltaG1(1e5),deltaG2(1e5),slopeG1(1e5),slopeG2(1e5),slopeG12(1e5),initialDisc(1e5),A_N_L_end(1e5),B_N_L_end(1e5),A_N_R_end(1e5),B_N_R_end(1e5),C_N_L_end(1e5),C_N_R_end(1e5),tailFactor(1e5),angle_two_phi(100.),angle_delta(100.),expG2(1e5),expG1(1e5);
int BC_switch(10000);

//////////////////////////////////// gsl specific variables, arrays, lists, vectors, other STL type instances   
size_t iterations(0);
int errorCount(0);
const gsl_interp_type* splineT = gsl_interp_cspline;
double xGrid[numberGridPoints+2*addGridPoints],xGridG1[numberGridPoints+2*addGridPoints],yGridG12[numberGridPoints_y+2*addGridPoints];
double prolonguedYGridG1[numberGridPoints+2*addGridPoints],prolonguedXGridG12[numberGridPoints_y+2*addGridPoints];
vector <double> Discretization(numberGridPoints+2*addGridPoints,10.);
vector <double> A_N_R(numberGridPoints+addGridPoints,10.);
vector <double> B_N_R(numberGridPoints+addGridPoints,10.);
vector <double> C_N_R(numberGridPoints+addGridPoints,10.);
vector <double> locCurvature(numberGridPoints,10.);
vector <double> curvature(numberGridPoints+numberGridPoints_y,10.); 
vector <double> HMK_curvature_RHP((numberGridPoints+1)/2,10.);
vector <double> integrand(numberGridPoints+numberGridPoints_y,10.);
////////////////////////////////////////// functions which do not return value ////////////////////////////
void setNumericalParameters(double& initialDiscretization, double& error_tol, string& gridSwitch, double& p1_o_p2);
void initializeGridAndSolution(double& discretization, string& gridSwitch,string& loadValue, gsl_vector* x, double& p1_o_p2);
void printState (gsl_multiroot_fdfsolver* solver/*gsl_multiroot_fsolver* solver*/);
void printChangingParams(gsl_vector* x_solution);
void defineGrid(double& discretization);
void loadSolution(gsl_vector* x,double& discretization);
void initializeSolutionVector(gsl_vector* x,double& p1_o_p2);
void setFileNames(string& namePath, string& shapeFileName);
void displayInitialGuess(gsl_vector* x);
void testIvantsovRelationWithLocalSplines(gsl_vector* x);
void defineFittingParabola(const gsl_vector* x);
void displayFittingParabolaAndDisc();
void updateSplines(const gsl_vector* x, gsl_spline* localSplineG1,  gsl_interp_accel* localAccelG1, gsl_spline* kappaRHP_splineG1, gsl_interp_accel* kappaRHP_accelG1,  gsl_spline* localSplineG12, gsl_interp_accel* localAccelG12, gsl_spline* kappa_splineG12, gsl_interp_accel* kappa_accelG12);
void calc_curv_HMK_RHP(vector<double>& s,vector<double>& x,vector<double>&
y, vector<double>& nx, vector<double>& ny, vector<double>& curv);
void calc_sep_RHP(vector<double>& s,vector<double>& x,vector<double>&
y,vector<double>& xm,vector<double>& ym,vector<double>& sp);
void calc_curv_HMK_G12(vector<double>& s,vector<double>& x,vector<double>&
y, vector<double>& nx, vector<double>& ny, vector<double>& curv);
void calc_sep_LHP(vector<double>& s,vector<double>& x,vector<double>&
y,vector<double>& xm,vector<double>& ym,vector<double>& sp);
float sgn(float val);
//void saveResultsInFile(string& filename, gsl_vector* x);
void createDirectory();


////////////////////////////////////////// functions which do return value ////////////////////////////
int set_equationSystem_Meiron(const gsl_vector* x, void* params, gsl_vector* equationVector);
int set_equationSystem_Meiron_fdf(const gsl_vector* x, void* params, gsl_vector* equationVector, gsl_matrix* jacobian);
int set_jacobian_equationSystem_Meiron(const gsl_vector* x, void* params, gsl_matrix* jacobian);
double exactDiffusionIntegrand_conservationLaw(double x, void* params);
double exactDiffusionIntegrand_locEq(double x, void* params);
double exactDiffusionIntegrand_l1l2(double y, void* params);

string IntToString ( int number);
string DoubleToString ( double doubleNumber);
///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// MAIN ///////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	if ((argv[1]) == NULL)
	{	
		cout << " -------------------------------------------------------- " << endl;
		cout << " fix if load or noload of old result: " << endl;
		cout << " [program name] load " << endl;
		cout << " [program name] noload " << endl; 
		cout << " -------------------------------------------------------- " << endl;
		exit(1);
	}
	string loadSwitch(argv[1]);
	cout << " -------------------------------------------------------- " << endl;
	cout << " call: " << argv[0] << " " << argv[1] << endl;
	cout << " -------------------------------------------------------- " << endl;
///////////////////////////////////////////////////////////////////////////////////////////////
	double discretization(0.),error_tol(0.0001),p1_o_p2(0.);
	string gridSwitch;
	int status(0)/*,goodConvergeSwitch(0)*/;
		
	
	gsl_vector* x = gsl_vector_calloc(numberGridPoints+numberGridPoints_y-1);
	setNumericalParameters(discretization,error_tol,gridSwitch,p1_o_p2);
	struct parameters meiron_Parameters = {pG1,pG2,deltaG1,deltaG2,slopeG1,slopeG2};
	initializeGridAndSolution(discretization,gridSwitch,loadSwitch,x,p1_o_p2);
	
	defineFittingParabola(x);
	
	displayFittingParabolaAndDisc();
////////////////////////////////////// UP TO HERE EVERYTHING FINE. //////////////////////////////
//// also tested in setEq -delta-function*1/2 of integrand if proper -Delta/2 comes out for p=0.01,0.1,1.0 - works 
////////////////////////////////// asymm Meiron ///////////////////////////////////////////////////////////////
  	const gsl_multiroot_fdfsolver_type* fdfsolverT = gsl_multiroot_fdfsolver_gnewton;/*hybridsj*/
  	gsl_multiroot_fdfsolver* solver1= gsl_multiroot_fdfsolver_alloc(fdfsolverT,numberGridPoints+numberGridPoints_y-1) ;

 	gsl_multiroot_function_fdf model_Meiron = 	{&set_equationSystem_Meiron,&set_jacobian_equationSystem_Meiron,&set_equationSystem_Meiron_fdf, numberGridPoints+numberGridPoints_y-1, &meiron_Parameters};

	gsl_multiroot_fdfsolver_set(solver1,&model_Meiron,x);
	
	printState(solver1);

//////////////////////////////////////////////////////////////////////////////////////////////////
	do
	{
		iterations++;
		status = gsl_multiroot_fdfsolver_iterate(solver1);
	
		printState(solver1);

		printChangingParams((solver1->x));

		if (status) {break;}

		status = gsl_multiroot_test_residual(solver1->f, error_tol);
		
		//exit(1);
	}
	while (status == GSL_CONTINUE && iterations < 1000);
	
	cout << "status= " << gsl_strerror(status) << endl;

	for (int i=0; i<numberGridPoints+numberGridPoints_y-1;i++)
	{
		if (i < numberGridPoints)
		{
			cout << "x,root(x),y_Iv: " <<  xGridG1[addGridPoints+i] << " " << gsl_vector_get(solver1->x,i) << " " << -0.5*pow(xGridG1[addGridPoints+i],2.0)   << endl;
		}
		else 
		{		
			cout << " y,root(y): " << yGridG12[addGridPoints+i-numberGridPoints] << ", " << gsl_vector_get(solver1->x,i) << endl; 

		}
	}  
	for (int i=0; i<numberGridPoints+numberGridPoints_y-1;i++)
	{
		if (i < numberGridPoints)
		{
			cout <<  xGridG1[addGridPoints+i] << " " << gsl_vector_get(solver1->x,i)  << endl;
		}
		else 
		{		
			cout <<  gsl_vector_get(solver1->x,i) << " "  << yGridG12[addGridPoints+i-numberGridPoints] << endl; 
			
		}
	}  
	gammaD = gsl_vector_get(solver1->x,0);
	//angle_delta = gsl_vector_get(solver1->x,(numberGridPoints+1)/2);
	lambda = gsl_vector_get(solver1->x,numberGridPoints);

/////////////////////////START: ///////////////////////////// save results in file //////////////////////////////////
	string pathName, shapeFileName,shapeLoadFileName;
	setFileNames(pathName, shapeFileName);
	double xSave(1e5),ySave(1e7),d1SOR1Save(1e3) ,lambdaSave(1e5);	

	createDirectory();

	ofstream outShapeFile;   
	const string outShapeFilename = pathName+shapeFileName;
	cout << outShapeFilename << " = name of file with shape data " << endl;
	outShapeFile.open(outShapeFilename.c_str());
	outShapeFile << setprecision(15);
	
	ofstream outCheckFile;
	const string outCheckFilename = pathName+"checkFile_"+shapeFileName;
	outCheckFile.open(outCheckFilename.c_str());
	outCheckFile << setprecision(15);
	
	string loadPath("loadFiles/");
	ofstream outShapeLoadFile;   
	const string outShapeLoadFilename = loadPath+"shapePreviousResult.dat" ;
	cout << outShapeLoadFilename << " = location of LoadFile with shape data " << endl;
	outShapeLoadFile.open(outShapeLoadFilename.c_str());
	outShapeLoadFile << setprecision(15);

	ofstream autoGnuScriptFile;
	const string autoGnuScriptName = "autoGnuScript.gnu";
	autoGnuScriptFile.open(autoGnuScriptName.c_str());
	ofstream autoGnuScriptFileCloseUp;
	const string autoGnuScriptCloseUpName = "autoGnuScriptCloseUp.gnu";
	autoGnuScriptFileCloseUp.open(autoGnuScriptCloseUpName.c_str());
	ofstream autoGnuScriptCheckFile;
	const string autoGnuScriptCheckName = "autoGnuScriptCheck.gnu";
	autoGnuScriptCheckFile.open(autoGnuScriptCheckName.c_str());
	
///////////////////////////////////////////////////////////////
string gnuFileBasisShape,  gnuOutputShape, plot1,completeGnuFileShape,gnuClosing,gnuFileBasisCheck, gnuOutputCheck;
////////////////////////////////////////////////////////////////////////////// define gnu file BEGIN
gnuFileBasisShape = "set terminal postscript enhanced monochrome 'Helvetica' 25;set key top right;set style line 1 lt -1 lw 2 pt 1 ps 1;set style line 2 lt 2 lw 1 pt 1 ps 1;set style line 3 lt 0 lw 6 pt 1 ps 1;set style line 4 lt 5 lw 2 pt 1 ps 1;set style line 5 lt 4 lw 2 pt 1 ps 1;set xlabel '{x}';set ylabel '{/Symbol x}';set xrange [-20:20];set yrange [-20:20];";
plot1 = "plot '"+outShapeFilename/*+"' using (-$1):($2) ls 1 title '-{/Symbol x}(x)', '" + resultFilePosition*/ + "' using ($1):($2) ls 2 title '{/Symbol x}(x)';";
gnuOutputShape = "set out 'peritecticPaperPreparePics/epsOf_"+shapeFileName+"_.eps"+"';";
gnuClosing = "set terminal x11;set out";
completeGnuFileShape = gnuFileBasisShape+gnuOutputShape+plot1+gnuClosing;
autoGnuScriptFile << completeGnuFileShape << endl;
////////////////////////////////////////// close up gnufile ////////////////////////////
gnuFileBasisShape = "set terminal postscript enhanced monochrome 'Helvetica' 25;set key top right;set style line 1 lt -1 lw 2 pt 1 ps 1;set style line 2 lt 2 lw 1 pt 1 ps 1;set style line 3 lt 0 lw 6 pt 1 ps 1;set style line 4 lt 5 lw 2 pt 1 ps 1;set style line 5 lt 4 lw 2 pt 1 ps 1;set xlabel '{x}';set ylabel '{/Symbol x}';set xrange [-7.5:7.5];set yrange [-7.5:7.5];";
//plot1 = "plot 'firstConvergence/"+fileName_Result/*+"' using (-$1):($2) ls 1 title '-{/Symbol x}(x)', '" + resultFilePosition*/ + "' using ($1):($2) ls 2 title '{/Symbol x}(x)';";
gnuOutputShape = "set out 'peritecticPaperPreparePics/epsOf_closeUp_"+shapeFileName+"_.eps"+"';";
gnuClosing = "set terminal x11;set out";
completeGnuFileShape = gnuFileBasisShape+gnuOutputShape+plot1+gnuClosing;
autoGnuScriptFileCloseUp << completeGnuFileShape << endl;
//////////////////////////////////// checkFile gnufile //////////////////////////////////
gnuFileBasisCheck = "set terminal postscript enhanced color 'Helvetica' 25;set key top right;set style line 1 lt -1 lw 2 pt 1 ps 1;set style line 2 lt 2 lw 1 pt 1 ps 1;set style line 3 lt 0 lw 6 pt 1 ps 1;set style line 4 lt 5 lw 2 pt 1 ps 1;set style line 5 lt 4 lw 2 pt 1 ps 1;set xlabel '{x}';set ylabel 'Int,{/Symbol k}';";
plot1 = "plot '"+outCheckFilename/*+"' using (-$1):($2) ls 1 title '-{/Symbol x}(x)', '" + resultFilePosition*/ + "' using ($1):($3) ls 2 title '{/Symbol k}(x)', '"+outCheckFilename+"' using ($1):($4) ls 3 title 'Int';";
gnuOutputCheck = "set out 'peritecticPaperPreparePics/epsOf_Check_"+shapeFileName+"_.eps"+"';";
gnuClosing = "set terminal x11;set out";
completeGnuFileShape = gnuFileBasisCheck+gnuOutputCheck+plot1+gnuClosing;
autoGnuScriptCheckFile << completeGnuFileShape << endl;
///////////////////////////////////////////////////////////////////////// define gnu file END
//////////////////////////////////////////////////////// define system call : gnuplot xxx.gnu BEGIN
string commandPlotResult("gnuplot "+autoGnuScriptName);
string commandPlotCloseUpResult("gnuplot "+autoGnuScriptCloseUpName);
string commandPlotResultCheck("gnuplot "+autoGnuScriptCheckName);

const char* commandPlotShape = commandPlotResult.c_str();
const char* commandPlotShapeCloseUp = commandPlotCloseUpResult.c_str();
const char* commandPlotCheck = commandPlotResultCheck.c_str();
//////////////////////////////////////////////////////// define system call : gnuplot xxx.gnu END
//////////////////////////////////////////////////////////////// END automated eps-generation	
	
	
	
	for(int k = 0; k<(numberGridPoints+numberGridPoints_y-1) ; k++) 
	{
		if (k == 0)
		{
		outShapeFile << "0." << " " << "0." << endl;  
		outShapeLoadFile << "0." << " " << "0." << endl;  
		outCheckFile << "0." << " " << "0." << " " <<  "0." << " " << "0." << endl;
		}
		else if (k>0 && (k<numberGridPoints))
		{	
			xSave = xGridG1[addGridPoints+k]; 
			ySave = gsl_vector_get(solver1->x,k);
			outShapeFile << xSave << " " << ySave << endl; 
			outShapeLoadFile << xSave << " " << ySave << endl; 
			outCheckFile << xSave << " " << ySave << " " << curvature.at(k) << " " << integrand.at(k) << endl; 
		}
		else if (k == numberGridPoints)
		{
			outShapeFile << "0." << " " << "0." << endl;  
			outShapeLoadFile << "0." << " " << "0." << endl;  
			outCheckFile << "0." << " " << "0." << " " <<  "0." << " " << "0." << endl;
		}
		else if (k>(numberGridPoints) && (k<(numberGridPoints+numberGridPoints_y-1)))
		{	
			xSave = gsl_vector_get(solver1->x,k); 
			ySave = yGridG12[addGridPoints+k-(numberGridPoints)];
			outShapeFile << xSave << " " << ySave << endl; 
			outShapeLoadFile << xSave << " " << ySave << endl; 
			//outCheckFile << xSave << " " << ySave << " " << curvature.at(k) << " " << integrand.at(k) << endl; 
		}
		else 
		{
			cout << " index error in main " << endl;
			exit(1);
		}			
	
		
	}
	outShapeFile << " " << endl; 
	outShapeLoadFile << " " << endl; 
	outCheckFile << " " << endl;

	d1SOR1Save = gsl_vector_get(solver1->x,0);
	//deltaSave = gsl_vector_get(solver1->x,(numberGridPoints+1)/2);
	lambdaSave = gsl_vector_get(solver1->x,numberGridPoints);
	outShapeFile << " ### d1S/R1 | delta | lambda | dG1 | dG2 ### " << endl;
 	outShapeLoadFile << " ### d1S/R1 | delta | lambda ### " << endl;
	outCheckFile << " ### d1S/R1 | delta | lambda ### " << endl;
	outShapeFile << " ### " << d1SOR1Save  << "     " << lambdaSave << "   " << deltaG1 << "   " << deltaG2  <<" ### " << endl;
	outShapeLoadFile << " ### " << d1SOR1Save  << "     " << lambdaSave << " ### " << endl;
	outCheckFile << " ### " << d1SOR1Save  << "     " << lambdaSave << " ### " << endl;
	
////////////END: ////////////////////////////save results in file //////////////////////////////////
 //////// autoplot figures BEGIN
  system(commandPlotShape);
  system(commandPlotShapeCloseUp);
  system(commandPlotCheck);
  cout << " plot check: " << commandPlotCheck << endl;
 //////// autoplot figures END
 


	gsl_multiroot_fdfsolver_free(solver1);
	gsl_vector_free(x);

	cout << " hier samma !" << endl;
////////////////////////////////////////////////////////////////////////////////////////////////////

return 0;
}	

///////////////////////////////////////////////////////////////////////
void initializeGridAndSolution(double& discretization, string& gridSwitch, string& loadValue, gsl_vector* x, double& p1_o_p2)
{
	if (loadValue == "load")
	{
		loadSolution(x,discretization);
	}
	else if (loadValue == "noload")
	{
		defineGrid(discretization);
		initializeSolutionVector(x, p1_o_p2);
	}
	else 
	{
		cout << " -------------------------------------------------------- " << endl;
		cout << " please fix if you want to: " << endl;
		cout << " -------------------------------------------------------- " << endl;
		cout << " load old solutions: enter [program name] load " << endl;
		cout << " start from asymptotical Ivantsov solution: enter [program name] noload " << endl;
		cout << " -------------------------------------------------------- " << endl;
		cout << " now please call program again with described flag values" << endl;
		exit(1);
	}
}
/////////////////////////////////////////////////////////////////////////////////
void loadSolution(gsl_vector* x,double& discretization)
{
	string previousSolution("shapePreviousResult.dat");
	string loadFilePath("loadFiles/");
			
	//loadPath+"YofXPreviousResult.dat" is entire load file generated in main after result obtained,
	//where loadPath = "loadFiles/"
	
	double readXValue(0), readYValue(0);
	double initialGuess[numberGridPoints+numberGridPoints_y];
	double sigmaInitialGuess(0.5);
	vector<double> oldGridVector, solutionVector;
			
	fstream loadPreviousSolutionFile;
	const string loadPreviousSolutionFileName = loadFilePath+previousSolution;
	loadPreviousSolutionFile.open(loadPreviousSolutionFileName.c_str());
	loadPreviousSolutionFile<<setprecision(15);	
	//cout << " filename: " << loadPreviousSolutionFileName << endl;
	
	int counter(0);
	while (loadPreviousSolutionFile >> readXValue >> readYValue)
	{
			oldGridVector.push_back(readXValue);
			solutionVector.push_back(readYValue);
			//cout << " commence avec entrainer les valeurs " << endl;	
			cout << " old Data: " << oldGridVector.at(counter) << " " << solutionVector.at(counter) << " " << counter << endl;
			++counter;
	}
			
	if ( (solutionVector.size() != numberGridPoints+numberGridPoints_y-1) || (oldGridVector.size() != solutionVector.size()))	// the solution does not explicitly contain y= 0 as first element, but Grid_y does !
	{
		std::cout << " les nombres du points de la grille ne coincident pas, prends la grille sans solution anterieure " << endl;
		cout << " solutionVector.size(): " << solutionVector.size() << " num. Gridpoints: " <<  (numberGridPoints+numberGridPoints_y) << endl;

		cout << " incorrect number of points, define grid " << endl; 
		exit(1);
		defineGrid(discretization);	
				
	}
			
	else 
	{	
		counter = 0;
		for ( vector<double>::iterator iterY = solutionVector.begin(), iterX = oldGridVector.begin(); iterY != solutionVector.end(); ++iterX, ++iterY )
		{
			int xGridIndex(counter);
			if (xGridIndex < (numberGridPoints))
			{
				initialGuess[xGridIndex] = *iterY; 
				cout << " @1.1 in loadSolution: i= " << xGridIndex << " iniGuess= " << initialGuess[xGridIndex] << " *iter= " << *iterY << endl;
			}
		
			if (xGridIndex >= (numberGridPoints))
			{
				initialGuess[xGridIndex] = *iterX; 
				cout << " @1.2 in loadSolution: i= " << xGridIndex << " iniGuess= " << initialGuess[xGridIndex] << " *iter= " << *iterX << endl;
			}
			++counter;
		}	 
		
		

		for (int i = 0; i<numberGridPoints+numberGridPoints_y-2; i++)
		{
			gsl_vector_set(x,i,initialGuess[i]);
			cout << " @2 in loadSolution: i= " << i << " iniGuess= " << initialGuess[i] << " x(i)= " << gsl_vector_get(x,i) << endl;
		}
		gsl_vector_set(x,0,sigmaInitialGuess);
		gsl_vector_set(x,numberGridPoints,oldGridVector.at(numberGridPoints+numberGridPoints_y-2));

		cout << " correct number of points, ini with old solution " << endl; 

		defineGrid(discretization);
	}

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initializeSolutionVector(gsl_vector* x, double& p1_o_p2)
{
	double initialGuess[numberGridPoints+numberGridPoints_y-1], sigmaInitialGuess(0.6),lambdaInitialGuess(lambda);
	list<double> solutionG1_G2,shortGrid;
	//int counter(0);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////	
	for (int xGridIndex = addGridPoints; xGridIndex < numberGridPoints+addGridPoints; xGridIndex++)
	{
		
		shortGrid.push_back(xGridG1[xGridIndex]);
		
	}	 
	
	for (list<double>::iterator iter = shortGrid.begin(); iter != shortGrid.end(); ++iter)
	{
		
		solutionG1_G2.push_back(0.5*lambdaInitialGuess*lambdaInitialGuess-0.5*(*iter-lambdaInitialGuess)*(*iter-lambdaInitialGuess));
		
	}
	
	for (list<double>::iterator iter = solutionG1_G2.begin(); iter != solutionG1_G2.end(); ++iter)         
	{                                                                                                      
		static int cnt(0);                                                                                    
		//cout <<  " xGrid1, *iter in solutionG1_G2 " << xGrid[cnt+addGridPoints] << " " << *iter << endl;      
		cnt++;                                                                                                
	}                                                                                                      
	
	for ( list<double>::iterator iter = solutionG1_G2.begin(); iter != solutionG1_G2.end(); ++iter )
	{
		static int cnt(0);
		initialGuess[cnt] = *iter;
		cnt++;
		//cout << "ini guess: cnt, iniGuess[cnt], *iter solutionG1_G2 " <<  cnt << " " << initialGuess[cnt-1] << " " << *iter  << endl;   // works up to here !!
	}	
	
	//exit(1);
	initialGuess[0] = sigmaInitialGuess;
	initialGuess[numberGridPoints] = lambdaInitialGuess; // standard initial value
	
	for (int i=1; i< numberGridPoints_y-1;i++)
	{
		//initialGuess[numberGridPoints+i] = -(-exp(-yGridG12[addGridPoints+i])+1.)*lambdaInitialGuess+exp(-yGridG12[addGridPoints+i])*slopeG12*yGridG12[addGridPoints+i]; // standard initial guess
		initialGuess[numberGridPoints+i] = exp(-10.*yGridG12[addGridPoints+i])*slopeG12*yGridG12[addGridPoints+i] + lambdaInitialGuess*(1.+exp(-yGridG12[addGridPoints+i]) - 2.*exp(-2./5. * yGridG12[addGridPoints+i]));
		//initialGuess[numberGridPoints+i] = sqrt(2.*yGridG12[addGridPoints+i]);
		//initialGuess[numberGridPoints+i] = 0.5*yGridG12[addGridPoints+i]*yGridG12[addGridPoints+i];
		//-1.*sqrt(2.*yGridG12[addGridPoints+i]); // testIvantsovrelation 	
		//cout << " initialGuess[(numberGridPoints+1)+i]= " << initialGuess[(numberGridPoints+1)+i] << endl;
	}
	
	
	for (int i = 0; i<numberGridPoints+numberGridPoints_y-1; i++)
	{
		//cout << " in set ini guess loop: i, iniguess: " << i << " " << initialGuess[i] << endl;
		gsl_vector_set(x,i,initialGuess[i]);
	}
	// up to here the initialGuess is just as it should be !!!
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
void defineGrid(double& discretization)
{
	
	list<double> rightGrid,yGrid;
	double xValue(1e5),yValue(1e+5);
	int x_counter(0),y_counter(0);
	// Discretization is vector <double> with n+2*add elmnts

// 	if (gridSwitch=="tanh")
	{
		double x1(1e5),x2(1e5),h1(1e5),h2(1e5);
		
		for (int extGridIndex = 0; extGridIndex < addGridPoints; extGridIndex++)
		{
			rightGrid.push_back(-discretization*(addGridPoints-extGridIndex));	
		}
		
		rightGrid.push_back(0.);
		for (int xGridIndex=1; xGridIndex < numberGridPoints+addGridPoints; ++xGridIndex)
		{
			h1 = 5./((numberGridPoints+addGridPoints-1)/2);
			h2 = (((numberGridPoints+addGridPoints-1))*discretization)/(((numberGridPoints+addGridPoints-1))*(10.*tanh(0.)+11.));
			x1 = xGridIndex*h1;
			x2 = 10.*tanh(x1-5.)+11.;
			xValue = 2.*xGridIndex*h2*x2;
			xValue = 0.0075 + pow(xValue,expG1);
			rightGrid.push_back(xValue);
			
		}
	
		for ( list<double>::iterator iter = rightGrid.begin(); iter != rightGrid.end(); ++iter )
		{
			int xGridIndex(x_counter);
			if (xGridIndex < (numberGridPoints+2*addGridPoints))
			{
				xGrid[xGridIndex] = *iter; 
			}
			++x_counter;
		}	 
	
//////////////////////////////////////////////////////////////
		double y1(1e+4),y2(1e+6);

		yGrid.push_back(0.);
		for (int yGridIndex=1; yGridIndex < numberGridPoints_y+addGridPoints; ++yGridIndex)
		{
			h1 = 5./((numberGridPoints_y-1));
			h2 = ((numberGridPoints_y-1)*discretization)/((numberGridPoints_y-1)*(10.*tanh(0.)+11.));
			y1 = yGridIndex*h1;
			y2 = 10.*tanh(y1-5.)+11.;
			yValue = 2.*yGridIndex*h2*y2;
			yValue = 0.0025 + pow(yValue,1.);
			yGrid.push_back(4.5*pow(yValue,1.375));
		}
		for ( list<double>::iterator iter = yGrid.begin(); iter != yGrid.end(); ++iter )
		{
			int yGridIndex(y_counter);
			if (yGridIndex < (numberGridPoints_y+addGridPoints))
			{
				yGridG12[yGridIndex+addGridPoints] = *iter; 
				//cout << " i= " << counter << " yGridG12[i]= " << yGridG12[yGridIndex+addGridPoints] << endl;  
			}
			++y_counter;
		}	 



/////////////////////////////////////////////////////////////
	}
///////////////////////////////////////   Discretization   //////////////////////////////////////////////////////////////////////
	
	for (int gridIndex = 1; gridIndex < numberGridPoints+2*addGridPoints; gridIndex++)
	{
		Discretization.at(gridIndex) = xGrid[gridIndex]-xGrid[gridIndex-1];
	}
	Discretization.at(0) = Discretization.at(1);
	
/////////////////////////////////////////// G1 spline (usually right branch)///////////////////////////////////
	for (int i =0; i<numberGridPoints+2*addGridPoints;i++)
	{
		xGridG1[i] = xGrid[i]; 
	}
/////////////////////////////////////////// G1 spline (usually right branch)///////////////////////////////////
	for (int i =0; i<addGridPoints;i++)
	{
		yGridG12[i] = -yGridG12[2*addGridPoints-i]; 
	}
/////////////////////////////////////////// G1 spline (usually right branch)///////////////////////////////////

}
////////////////////////////////////////////////////////////////////////
void setNumericalParameters(double& initialDiscretization, double& error_tol, string& gridSwitch, double& p1_o_p2)
{

	initialDiscretization = 0.08125;
	expG1 = 1.4;
	//expG2 = expG1;
	gridSwitch = "tanh";
	error_tol = 5.*1e-3;
	initialDisc = initialDiscretization; // initialDisc is used as globally declared parameter for normalShift and naming 
	pG1 = 0.00825;
	//pG2 = 0.1;
	deltaG1 = sqrt(PI*pG1)*exp(pG1)*erfc(sqrt(pG1));
	//deltaG2 = sqrt(PI*pG2)*exp(pG2)*erfc(sqrt(pG2));
	p1_o_p2 = pG1/pG2;

	angle_two_phi = 2.0944; // opening angle exactly 120 degrees
	//angle_two_phi = 1.0472; // opening angle exactly 60 degrees
	//angle_two_phi = 1.5708; //opening angle exactly 90 degrees
	//angle_two_phi = 1.7918; //slopeG1 = -0.8 
	//angle_two_phi = 2.7468; //slopeG1 = -0.2
	
	//angle_delta = -0.05;
	lambda =.5;
	//sigma = 0.5;
	//Dc1S_Dc12=-0.5;
	//Dc2S_Dc12=.5;
	DcLDelta_DcLGamma = -2.0;

	//slopeG2 = -1./tan(angle_delta-0.5*angle_two_phi);
	slopeG1 = 1.;//-1./tan(angle_delta+0.5*angle_two_phi);
	slopeG12 = -1.;//-1.*tan(angle_delta);
	
	//exit(1); 
	
	//slopeG1 = -0.0; // right branch
	//slopeG2 =  0.0; // left branch
	//slopeG12 = -0.2; //

	cout << " slopeG1: " << slopeG1  << " slopeG2: " << slopeG2 << " slopeG12: " << slopeG12 << endl;

	tailFactor = 1.00;

}
/////////////////////////////////////////////////////////////////////////
int set_equationSystem_Meiron(const gsl_vector* x, void* params,gsl_vector* equationVector)
{
	static size_t nr_access(0);
	//cout << " # accesses to setEqSystem " << nr_access << endl;
	++nr_access;
	
// for meiron problem, since slopes are not changed, not necessary yet to set prolongued splines for calculated slope in each step.
	double equation[numberGridPoints+numberGridPoints_y-1];	
	for (int i=0;i<numberGridPoints+numberGridPoints_y-1;i++) {equation[i] = 100.;}
	double dydx(0.),dxdy(1e6),yObs(0.),y_0(1e5),y_p1(1e5),x_p1,kappa(0.),/*kappaHMK(0.),*/d2ydx2(0.),d2xdy2(1e5),sigma(1e5);
	double xObs(0.);
	double Integral(100.),error(1e5),rightIntegral_locEq(10.),rightIntegral_conservationLaw(10.),conservationLawIntegral(10.),locEqIntegral(10.),l1l2Integral(1e5);

	gsl_set_error_handler_off();
	
	gsl_spline* kappaRHP_splineG1 = gsl_spline_alloc(splineT,numberGridPoints);
	gsl_spline* kappa_splineG12 = gsl_spline_alloc(splineT,numberGridPoints_y);
	//gsl_spline* kappaLHP_splineG2 = gsl_spline_alloc(splineT,(numberGridPoints+1)/2);
	gsl_interp_accel* kappaRHP_accelG1 = gsl_interp_accel_alloc();
	gsl_interp_accel* kappa_accelG12 = gsl_interp_accel_alloc();
	//gsl_interp_accel* kappaLHP_accelG2 = gsl_interp_accel_alloc();

	gsl_spline* current_splineG1 = gsl_spline_alloc(splineT,numberGridPoints+2*addGridPoints);
	//gsl_spline* current_splineG2 = gsl_spline_alloc(splineT,(numberGridPoints+1)/2+2*addGridPoints);
	gsl_spline* current_splineG12 = gsl_spline_alloc(splineT,numberGridPoints_y+2*addGridPoints);
	gsl_interp_accel* current_accelG1 = gsl_interp_accel_alloc();
	//gsl_interp_accel* current_accelG2 = gsl_interp_accel_alloc();
	gsl_interp_accel* current_accelG12 = gsl_interp_accel_alloc();

	sigma = gsl_vector_get(x,0);
	//angle_delta = gsl_vector_get(x,(numberGridPoints+1)/2);
	lambda = gsl_vector_get(x,numberGridPoints);

	//slopeG2 = -1./tan(angle_delta-0.5*angle_two_phi);
	slopeG1 = 1.;//-1./tan(angle_delta+0.5*angle_two_phi);
	slopeG12 = -1.;//-1.*tan(angle_delta);
	
	updateSplines(x, current_splineG1,  current_accelG1, kappaRHP_splineG1,  kappaRHP_accelG1,  current_splineG12, current_accelG12, kappa_splineG12, kappa_accelG12);
	
	//for (int i = 0; i < numberGridPoints+2*addGridPoints; i++)
//	{
//		cout << " i, , x, prol YGridG1 " << i << " " << xGridG1[i] << " " << prolonguedYGridG1[i] << endl;
//	}
//	for (int i = 0; i < numberGridPoints_y+2*addGridPoints; i++)
//	{
//		cout << " i, , x, prol XGridG12 " << i << " " << yGridG12[i] << " " << prolonguedXGridG12[i] << endl;
//	}
	
	
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(2000);

	struct integralParameters currentParams = {sigma,xObs,yObs,current_splineG1,current_accelG1,kappaRHP_splineG1, kappaRHP_accelG1};
	
	struct integralParameters_l1l2Interface current_l1l2Params = 
	{sigma, xObs,yObs,current_splineG12, current_accelG12, kappa_splineG12, kappa_accelG12};

	//////////////////////////////////////////////////////////////////////////////////////////////
	y_0 =  gsl_spline_eval(current_splineG1,0.,current_accelG1);
	y_p1 = gsl_spline_eval(current_splineG1,xGridG1[addGridPoints+1],current_accelG1);
	x_p1 = gsl_spline_eval(current_splineG12,yGridG12[addGridPoints+1],current_accelG12);	

// ////////////////////////////////////////////////// equationVectorSetting ////////////////////////	
	equation[0] = y_p1-slopeG1*xGridG1[addGridPoints+1];
	//equation[0] = gsl_spline_eval(current_splineG1,xGridG1[addGridPoints+1],current_accelG1) - xGridG1[addGridPoints+1];

	//cout << " i=0: eq[0], y_p1 " << equation[0] << " " << y_p1 << endl; 
	
	for (int i = 1; i< numberGridPoints;i++)
	{
		size_t n_evaluations;
		xObs = xGridG1[i+addGridPoints];	

		yObs = gsl_spline_eval(current_splineG1,xObs,current_accelG1);
		
		dydx = gsl_spline_eval_deriv(current_splineG1,xObs,current_accelG1);
		d2ydx2 = gsl_spline_eval_deriv2(current_splineG1,xObs,current_accelG1);
		kappa = gsl_spline_eval(kappaRHP_splineG1,xObs,kappaRHP_accelG1);
/////////////////////////////// arbitray p integral /////////////////////////////
		gsl_function integrateExactDiffKernel;
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_conservationLaw;
		integrateExactDiffKernel.params =  &currentParams;
		currentParams.xObserv = xObs;
		currentParams.yObserv = yObs;
		
 		gsl_integration_qng(&integrateExactDiffKernel,0.00/*xGridG1[addGridPoints+1]*/, tailFactor*xGridG1[addGridPoints+numberGridPoints-2],0,1e-7,&Integral,&error,&n_evaluations);
		rightIntegral_conservationLaw = Integral;
		conservationLawIntegral = rightIntegral_conservationLaw; 


////////////////////////////////////arbitrary p integral /////////////////////////////////////////////////////////////////////
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_locEq;
		integrateExactDiffKernel.params =  &currentParams;
		currentParams.xObserv = xObs;
		currentParams.yObserv = yObs;

//////////////////////////////////////////////////////////////////////////////////////////
		

		gsl_integration_qng(&integrateExactDiffKernel,/*xGridG1[addGridPoints+1]*/0.00, tailFactor*xGridG1[addGridPoints+numberGridPoints-2],0,1e-7,&Integral,&error,&n_evaluations);

		rightIntegral_locEq = Integral;

		locEqIntegral = rightIntegral_locEq;
///////////////////////////////////////////////////////////////
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_l1l2;
		integrateExactDiffKernel.params = &current_l1l2Params;
		current_l1l2Params.yObserv = yObs;
		current_l1l2Params.xObserv = xObs;
		
		// test curvature spline along L Delta interface
		//cout << " test curv spline L Delta @ yObs  " << gsl_spline_eval( kappa_splineG12,yObs, kappa_accelG12) << " " << yObs << endl;
		
		gsl_integration_qng(&integrateExactDiffKernel, 0., yGridG12[addGridPoints+numberGridPoints_y-2],0,1e-7,&l1l2Integral,&error,&n_evaluations);
		
		Integral = 1.*conservationLawIntegral + 1.*locEqIntegral + l1l2Integral; 


		integrand.at(i) = pG1*(1./(2.*PI))*Integral;
		curvature.at(i) = kappa;
		
		equation[i] = 0.5*(deltaG1-kappa*sigma)-(pG1/(2.*PI))*Integral;
		
		//cout << " in L Gamma interface block: i, eq[i] = " << i << " " << equation[i] << endl;

}
	
	//equation[0] = B_N_L_end; 
	equation[numberGridPoints-1] = B_N_R_end;
	//cout << " BNR end eq : eq, BNR_end " <<  equation[numberGridPoints] << " " << B_N_R_end << endl;
//	equation[0] = C_N_L_end+0.5; 
	//equation[numberGridPoints] = C_N_R_end+0.5;
	
	equation[numberGridPoints] = x_p1-slopeG12*yGridG12[addGridPoints+1];
	//equation[numberGridPoints] = gsl_spline_eval(current_splineG12,yGridG12[addGridPoints+1],current_accelG12) + yGridG12[addGridPoints+1];

	//cout << " slopeG12 eq: xp1, slopeG12, yGridG12, eq " << x_p1 << " " << slopeG12 << " " << yGridG12[addGridPoints+1] << " " <<  equation[(numberGridPoints+1)] << endl;
	
for (int i=1; i<numberGridPoints_y-1;i++)
{
		size_t n_evaluations;
		yObs = yGridG12[i+addGridPoints];	
		xObs = gsl_spline_eval(current_splineG12,yObs,current_accelG12);
		dxdy = gsl_spline_eval_deriv(current_splineG12,yObs,current_accelG12);
		d2xdy2 = gsl_spline_eval_deriv2(current_splineG12,yObs,current_accelG12);
		//kappa = -d2xdy2/pow((1+dxdy*dxdy),1.5);
		kappa = gsl_spline_eval(kappa_splineG12,yObs,kappa_accelG12);
	
/////////////////////////////////////////// arbitrary p integral ////////////////////////////////////////////////
		gsl_function integrateExactDiffKernel;
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_conservationLaw;
		integrateExactDiffKernel.params =  &currentParams;
		currentParams.yObserv = yObs;
		currentParams.xObserv = xObs;
//////////////////////////////////////////kappa from parabola at outer points ////////////////////////////
		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
 		gsl_integration_qng(&integrateExactDiffKernel,/*xGridG1[addGridPoints+1]*/0.00, tailFactor*xGridG1[addGridPoints+numberGridPoints-2],0,1e-7,&Integral,&error,&n_evaluations);
		rightIntegral_conservationLaw = Integral;
		conservationLawIntegral  =  rightIntegral_conservationLaw;
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////arbitrary p integral ////////////////////////////////////////////////////////////////
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_locEq;
		integrateExactDiffKernel.params =  &currentParams;
		currentParams.yObserv = yObs;
		currentParams.xObserv = xObs;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
		gsl_integration_qng(&integrateExactDiffKernel,0.00/*xGridG1[addGridPoints+1]*/, tailFactor*xGridG1[addGridPoints+numberGridPoints-2],0,1e-7,&Integral,&error,&n_evaluations);
		rightIntegral_locEq = Integral;

		locEqIntegral =  rightIntegral_locEq;
////////////////////////////////////////////////////////////////
		integrateExactDiffKernel.function = &exactDiffusionIntegrand_l1l2;
		integrateExactDiffKernel.params =  &current_l1l2Params;
		current_l1l2Params.yObserv = yObs;
		current_l1l2Params.xObserv = xObs;
		
		gsl_integration_qng(&integrateExactDiffKernel, 0., yGridG12[addGridPoints+numberGridPoints_y-2],0,1e-7,&l1l2Integral,&error,&n_evaluations);

		
//////////////////////////////////////////////////////////////////////
		Integral = 1.*conservationLawIntegral + 1.*locEqIntegral + l1l2Integral;

		equation[i+numberGridPoints] = -0.5*DcLDelta_DcLGamma*kappa*sigma-(pG1/(2.*PI))*Integral; // 1.0 here comes from assumption that d_12 = d_1S !
		
		
		integrand.at(i+numberGridPoints) = pG1*(1./(2.*PI))*Integral; // integrand is checked for -Delta/2 - behaviour.
		curvature.at(i+numberGridPoints) = kappa;
 
	
	//cout << " in L Delta Block: i, eq[i] " << i << " " << equation[i+numberGridPoints] << endl; 
	
}	
	
	//exit(1);
	
	//for (int i = 0; i < numberGridPoints+numberGridPoints_y-1; i++)
//	{
//		cout << " i, eq[i] " << i << " " << equation[i]  << endl;
//	}
	
	
	for(int i = 0; i< numberGridPoints+numberGridPoints_y-1; i++)
	{
		gsl_vector_set(equationVector,i,equation[i]);
	}
	gsl_integration_workspace_free(w);

	gsl_interp_accel_free(current_accelG1);
	gsl_interp_accel_free(current_accelG12);
	
	gsl_spline_free(current_splineG1);
	gsl_spline_free(current_splineG12);

	gsl_interp_accel_free(kappaRHP_accelG1);
	gsl_interp_accel_free(kappa_accelG12);
	gsl_spline_free(kappaRHP_splineG1);
	gsl_spline_free(kappa_splineG12);


	return GSL_SUCCESS;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
int set_jacobian_equationSystem_Meiron(const gsl_vector* x, void* params, gsl_matrix* jacobian)
{
	double d_x(1e-8);
	gsl_vector* equationVectorAt_x = gsl_vector_alloc(numberGridPoints+numberGridPoints_y-1);
	gsl_vector* equationVectorAt_xdx = gsl_vector_alloc(numberGridPoints+numberGridPoints_y-1);
	gsl_vector* xVector = gsl_vector_alloc(numberGridPoints+numberGridPoints_y-1);
	gsl_vector* xdxVector = gsl_vector_alloc(numberGridPoints+numberGridPoints_y-1);
 	gsl_vector_memcpy(xVector,x);
	gsl_vector_memcpy(xdxVector,x);
	double jacobian_ij(1e5);
	
	set_equationSystem_Meiron(xVector, params, equationVectorAt_x);

	for (int j=0;j<numberGridPoints+numberGridPoints_y-1;j++)
	{
		gsl_vector_memcpy(xdxVector,xVector);
		gsl_vector_set(xdxVector,j,gsl_vector_get(xVector,j)+d_x);
		set_equationSystem_Meiron(xdxVector, params, equationVectorAt_xdx);
		for (int i=0;i<numberGridPoints+numberGridPoints_y-1;i++)
		{
			
			jacobian_ij = ( gsl_vector_get(equationVectorAt_xdx,i) - gsl_vector_get(equationVectorAt_x,i) )/d_x;
			gsl_matrix_set(jacobian, i, j, jacobian_ij);
		}
	}

	gsl_vector_free(equationVectorAt_x);
	gsl_vector_free(equationVectorAt_xdx);
	gsl_vector_free(xVector);
	gsl_vector_free(xdxVector);


	return GSL_SUCCESS;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
int set_equationSystem_Meiron_fdf(const gsl_vector* x, void* params, gsl_vector* equationVector, gsl_matrix* jacobian)
{
	set_equationSystem_Meiron(x, params, equationVector);
////////////////////////////////////////////////////////////////////////////////////////////////
	set_jacobian_equationSystem_Meiron(x, params, jacobian);
	
	return GSL_SUCCESS;
}
//////////////////////////////////////////////////////////////////////////////
void printState (gsl_multiroot_fdfsolver* solver)
{
	double sum_error(0.);
	for (int i=0; i<numberGridPoints+numberGridPoints_y-1;i++)	
	{
		if (i>=0 && (i<numberGridPoints)) 
		{
			cout << " i, x, y(x) = " << i << " " << xGridG1[addGridPoints+i] << " " << gsl_vector_get(solver->x, i) << endl; 
		}
		else if (i>=(numberGridPoints))
		{
			cout << " i, y, x(y) = " << i << " " << yGridG12[addGridPoints-(numberGridPoints)+i] << " " << gsl_vector_get(solver->x, i) << endl;
		}
		sum_error += fabs(gsl_vector_get(solver->f, i));
	}
	cout << "iterations= " << iterations << " x_tip(sigma)= " << gsl_vector_get(solver->x, 0) << " y_tip(lambda)= " << gsl_vector_get(solver->x, numberGridPoints) << " sum_i |f_i(x)|= " << sum_error << endl;    


}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void setFileNames(string& pathName, string& shapeFileName)
{

	string lamb = DoubleToString(lambda);
	string d1SOR1 = DoubleToString(gammaD);
	//string delt = DoubleToString(angle_delta);
	string twoPhi = DoubleToString(angle_two_phi); 
	string PG1 = DoubleToString(pG1);
	//string rGap1S12 = DoubleToString(Dc1S_Dc12);
	//string rGap2S12 = DoubleToString(Dc2S_Dc12);
	string rGapDeltaGamma = DoubleToString(DcLDelta_DcLGamma);
	string n = IntToString(numberGridPoints+numberGridPoints_y-1);
	string tipS1 = DoubleToString(slopeG1); // Gamma 1 is right interface: x >0 , y<0
	string tipS2 = DoubleToString(slopeG2); // Gamma 2 is left interface: x<0,y<0
	string Ver = "peritectic";
	string BC = IntToString(BC_switch);
	string Disc = DoubleToString(initialDisc);
	string shape = "shape_";
	string valExpG1 = DoubleToString(expG1);
	
	//string goodConverge = IntToString(goodConvergeSwitch); // goodConvergeSwitch = 1 means good Convergence
//////////////////////////////////////////////////////////
/////////// change path  /////////////////////
	pathName = "peritecticRunsPrepareForPaper/";

	string configuration = Ver+"_n_"+n+"_disc_"+Disc+"_expG1_"+valExpG1+"_pG1_"+PG1+"_rG1S12_"+rGapDeltaGamma+"_2Phi_"+twoPhi+"_slopeG1_"+tipS1+"_slopeG2_"+tipS2+"_d1SOR1_"+d1SOR1+"_lambda_"+lamb+"_BC_"+BC+".dat";

/////////////////////////////////////////////////////////////////////////////////////////////////
	shapeFileName = shape+configuration;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
string IntToString ( int number )
{
  std::ostringstream oss;
  
  oss<< number;
 
  return oss.str();
}
//////////////////////////////////////////////////////////
string DoubleToString ( double doubleNumber )
{
  std::ostringstream oss;
 
  oss<< doubleNumber;
  
  return oss.str();
}
////////////////////////////////////////////////////////////
void displayInitialGuess(gsl_vector* x)
{
	cout << " display initial grids and guess " << endl;

	for (int i=0; i< (numberGridPoints+2*addGridPoints);i++)
	{
		cout << "xGridG1: " << xGridG1[i] << endl;
	}
	for (int i=0; i< numberGridPoints;i++)
	{
		cout << "(x,root_ini): " << gsl_vector_get(x,i) << endl;
	}
	for (int i=0; i< numberGridPoints_y-1;i++)
	{
		cout << " (y,root_ini): " << (yGridG12[i]) << "," << gsl_vector_get(x,i+(numberGridPoints));
	}
}
//////////////////////////////////////////////////////////
void updateSplines(const gsl_vector* x, gsl_spline* localSplineG1,  gsl_interp_accel* localAccelG1,  gsl_spline* kappaRHP_splineG1,  gsl_interp_accel* kappaRHP_accelG1,  gsl_spline* localSplineG12, gsl_interp_accel* localAccelG12, gsl_spline* kappa_splineG12, gsl_interp_accel* kappa_accelG12)
{
static int nrAccess(0);
++nrAccess;

	double interpolatedDataG1[numberGridPoints+ 2*addGridPoints],interpolatedDataG12[numberGridPoints_y+2*addGridPoints],interpolatedKappaRHP[numberGridPoints], interpolatedKappa_G12[numberGridPoints_y];
	
	//for (int i = 0; i<numberGridPoints+2*addGridPoints; i++)
//	{
//		cout << " i, prolonguedYGridG1[i] " << i << " " << prolonguedYGridG1[i] << endl;
//	}
	
	defineFittingParabola(x);

	//for (int i = 0; i<numberGridPoints+2*addGridPoints; i++)
//	{
//		cout << " i, prolonguedYGridG1[i] " << i << " " << prolonguedYGridG1[i] << endl;
//	}
	//exit(1);
	for(int i =0; i< numberGridPoints+2*addGridPoints;i++)
	{	
		interpolatedDataG1[i] = prolonguedYGridG1[i];
	}

	for(int i =0; i< numberGridPoints_y+2*addGridPoints;i++)
	{
		interpolatedDataG12[i] = prolonguedXGridG12[i];
	}
	
	gsl_spline_init(localSplineG1,xGridG1,interpolatedDataG1,numberGridPoints+2*addGridPoints);
	gsl_spline_init(localSplineG12,yGridG12,interpolatedDataG12,numberGridPoints_y+2*addGridPoints);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vector <double> s_vecRHP(numberGridPoints,0.);
	vector <double> x_vecRHP(numberGridPoints,0.);
	vector <double> y_vecRHP(numberGridPoints,0.);
	vector <double> nx_vecRHP(numberGridPoints,0.);
	vector <double> ny_vecRHP(numberGridPoints,0.);
	vector <double> curv_vecRHP(numberGridPoints,0.);
	vector <double> s_vecG12(numberGridPoints_y,0.);
	vector <double> x_vecG12(numberGridPoints_y,0.);
	vector <double> y_vecG12(numberGridPoints_y,0.);
	vector <double> nx_vecG12(numberGridPoints_y,0.);
	vector <double> ny_vecG12(numberGridPoints_y,0.);
	vector <double> curv_vecG12(numberGridPoints_y,0.); // before +1)/2

	for(int i =0; i< numberGridPoints;i++)
	{	
		x_vecRHP.at(i) = xGridG1[addGridPoints+i];
		
		y_vecRHP.at(i) = prolonguedYGridG1[addGridPoints+i];
		
	}
	
	for (int i = 0; i < numberGridPoints_y;i++)
	{
		x_vecG12.at(i) = prolonguedXGridG12[addGridPoints+i];
		
		y_vecG12.at(i) = yGridG12[addGridPoints+i];
	}

	calc_curv_HMK_G12(s_vecG12, y_vecG12, x_vecG12, ny_vecG12, nx_vecG12, curv_vecG12);
	curv_vecG12.at(0) = curv_vecG12.at(1);
	curv_vecG12.at(curv_vecG12.size()-1) = curv_vecG12.at(curv_vecG12.size()-2);
	
	//calc_curv_HMK_G12(s_vecG12, x_vecG12, y_vecG12, nx_vecG12, ny_vecG12, curv_vecG12);
 	calc_curv_HMK_RHP(s_vecRHP, x_vecRHP, y_vecRHP, nx_vecRHP, ny_vecRHP, curv_vecRHP);
	
	size_t xGridIndex_1(0), xGridIndex_2(0);
	for (vector<double>::iterator iterRHP = curv_vecRHP.begin(); iterRHP != curv_vecRHP.end(); ++iterRHP) 
	{
		interpolatedKappaRHP[xGridIndex_1] = *iterRHP;
		++xGridIndex_1;
	}
	for (vector<double>::iterator iterG12 = curv_vecG12.begin(); iterG12 != curv_vecG12.end(); ++iterG12) 
	{
 		interpolatedKappa_G12[xGridIndex_2] = (*iterG12);
		//cout << "i, interpolatedKappa_G12[i] , *iterG12  " << xGridIndex_2 << " " << interpolatedKappa_G12[xGridIndex_2] << " " << *iterG12 << endl; // here correct

		++xGridIndex_2;
	}
	
	
	interpolatedKappa_G12[numberGridPoints_y-1] = interpolatedKappa_G12[numberGridPoints_y-2];
	interpolatedKappa_G12[0] = interpolatedKappa_G12[1];	

	double shortXGridG1[numberGridPoints],shortYGridG12[numberGridPoints_y];
	for (int i=0; i< numberGridPoints;i++)
	{
		shortXGridG1[i] = xGridG1[addGridPoints+i];	
		
	}
	for (int i=0; i< numberGridPoints_y;i++)
	{
		shortYGridG12[i] = yGridG12[addGridPoints+i];	
	} 
	gsl_spline_init(kappa_splineG12,shortYGridG12,interpolatedKappa_G12,numberGridPoints_y);
	gsl_spline_init(kappaRHP_splineG1,shortXGridG1,interpolatedKappaRHP,numberGridPoints);
	
}
///////////////////////////////////////////////////////////////////////////
double exactDiffusionIntegrand_conservationLaw(double x, void* params)
{
	
	double xObs =((struct integralParameters*) params)-> xObserv;
	
	double KernelSource(10.),Kernel(10.),eta(10.),kappaLoc(10.)/*,Kernel_nGradGreensFunc(10.)*/;
	
	double yObs(1e5);
	double y(1e5),locDisc(10e-7);
	double dydx(0.),d2ydx2(0.);
	double y_x_m_locDisc(15.),y_x_p_locDisc(12.);
	double ratioDc(100.);

	
		yObs = ((struct integralParameters*) params)-> yObserv;

	
		if ( (x >= 0.) && (x<xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc));

			dydx = ( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x-locDisc,( ((struct integralParameters*) params)-> accelG1loc)) )/locDisc;
			
			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaRHP_splineG1),x,( ((struct integralParameters*) params)-> kappaRHP_accelG1));
			ratioDc = 1.;
		}
		else if ( (x >= xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = A_N_R_end+B_N_R_end*x+C_N_R_end*x*x;
			y_x_m_locDisc = A_N_R_end+B_N_R_end*(x-locDisc)+C_N_R_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_R_end+B_N_R_end*(x+locDisc)+C_N_R_end*(x+locDisc)*(x+locDisc);

			dydx = (y-y_x_m_locDisc)/locDisc;

			d2ydx2 = (y_x_p_locDisc-2.*y+y_x_m_locDisc)/(locDisc*locDisc);

			kappaLoc = - ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
			ratioDc = 1.;
		}
		else 
		{
			cout << " xObs= " << xObs << endl; 
			cout << " error in exactDiffusionIntegrand_conservationLaw c " << endl;
			exit(1);
		}

	eta = sqrt( (x-xObs)*(x-xObs) + (y-yObs)*(y-yObs) );
 
	if ((pG1*eta) >= 0.00000001)
	{
		if ((pG1*eta) >= 2.5)
		{
			KernelSource = sqrt(PI/(2.*pG1*eta))*2.*exp( -pG1*( (yObs-y)+ eta ))*ratioDc;
		}
		else 
		{
			KernelSource = 2.*exp(-(pG1*(yObs-y)))*gsl_sf_bessel_K0(pG1*eta)*ratioDc;
		}
		Kernel = KernelSource;
	}

	else 
	{
		Kernel = 0.;
	}
	return Kernel;
}
////////////////////////////////////////////////////////////////
double exactDiffusionIntegrand_l1l2(double y, void* params)
{
	
	double xObs =((struct integralParameters_l1l2Interface*) params)-> xObserv;
	double yObs = ((struct integralParameters_l1l2Interface*) params)-> yObserv;
	double sigmaLoc = ((struct integralParameters*) params)-> sigma;

	double KernelSource(10.),Kernel(10.),eta(10.), Kernel1sidedPart(1e5);
	double x(1e5),locDisc(10e-5),dxdy(1e5);
	double kappaLoc(1e5),ratioDc(1e5),locEq(1e5);
	
	yObs = ((struct integralParameters_l1l2Interface*) params)-> yObserv;	

	x = gsl_spline_eval(( ((struct integralParameters_l1l2Interface*) params)-> splineG12loc),y,( ((struct integralParameters_l1l2Interface*) params)-> accelG12loc)); // x(y) = x(-y)

	dxdy = -( gsl_spline_eval(( ((struct integralParameters_l1l2Interface*) params)-> splineG12loc),y-locDisc,( ((struct integralParameters_l1l2Interface*) params)-> accelG12loc)) - gsl_spline_eval(( ((struct integralParameters_l1l2Interface*) params)-> splineG12loc),y,( ((struct integralParameters_l1l2Interface*) params)-> accelG12loc)) )/locDisc; 

	eta = sqrt( (x-xObs)*(x-xObs) + (y-yObs)*(y-yObs) ); 
 
	kappaLoc = gsl_spline_eval(( ((struct integralParameters_l1l2Interface*) params)-> kappa_splineG12),y,( ((struct integralParameters_l1l2Interface*) params)-> kappa_accelG12));
	
	ratioDc = DcLDelta_DcLGamma;
	// for tests: just integrate kappa along branch for y = 0.5 x^2 initialization to check if curvature spline & integral do good
	locEq = -kappaLoc*sigmaLoc;// multiplied by ratioDc when returned value defined
	

	//cout << "x, y, dxdy, kappaLoc " << x << " " << y << " " << dxdy << " " << kappaLoc  << endl;
	//exit(1);
	
	
	if ((pG1*eta) >= 0.00000001)
	{
		if ((pG1*eta) >= 4.)
		{
			
			//Kernel1sidedPart = locEq*sqrt(PI/(2.*pG1*eta))*exp( -pG1*( (yObs-y)+eta ) )*( ((-dydx*(xObs-x)+(yObs-y))/eta)*(1.+(3./8.)*(1./(pG1*eta))) - 1. ); // would be copy of expression for y(x) interface solution
			Kernel1sidedPart = locEq*sqrt(PI/(2.*pG1*eta))*exp( -pG1*( (yObs-y)+eta ) )*( ((-1.*(xObs-x)+dxdy*(yObs-y))/eta)*(1.+(3./8.)*(1./(pG1*eta))) - dxdy*1. );
			
		}
		else
		{
			//Kernel1sidedPart = locEq*exp(-(pG1*(yObs-y)))*( ((-dydx*(xObs-x)+(yObs-y))/eta)*gsl_sf_bessel_K1(pG1*eta) - gsl_sf_bessel_K0(pG1*eta) ); // would be copy of expression for y(x) interface solution
			Kernel1sidedPart = locEq*exp(-(pG1*(yObs-y)))*( ( (-1.*(xObs-x)+dxdy*(yObs-y) )/eta)*gsl_sf_bessel_K1(pG1*eta) - dxdy*gsl_sf_bessel_K0(pG1*eta) );
			
		}
		
		
		
		if ((pG1*eta) >= 2.5) //gave nice Ivantsov relation matching, error < 1%
		{
			KernelSource = dxdy*sqrt(PI/(2.*pG1*eta))*2.*exp( -pG1*( (yObs-y)+ eta )); 		}
		else 
		{
			KernelSource = dxdy*2.*exp(-(pG1*(yObs-y)))*gsl_sf_bessel_K0(pG1*eta);
		}
		//Kernel = KernelSource;
		//Kernel = locEq*exp(-(pG1*(yObs-y)))*( ((-1.*(xObs-x)+dxdy*(yObs-y))/eta)*gsl_sf_bessel_K1(pG1*eta) - dxdy*gsl_sf_bessel_K0(pG1*eta) )  + dxdy*2.*exp(-(pG1*(yObs-y)))*gsl_sf_bessel_K0(pG1*eta);
		
		
		Kernel = (Kernel1sidedPart + KernelSource)*ratioDc;
	}

	else 
	{
		Kernel = 0.;
	}

	return Kernel;
	//return kappaLoc;
	
}
////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
double exactDiffusionIntegrand_locEq(double x, void* params)
{
	
	double sigmaLoc = ((struct integralParameters*) params)-> sigma;
	double xObs =((struct integralParameters*) params)-> xObserv;
	
	double Kernel1sidedPart(10.),Kernel(10.),eta(10.),locEq(10.),kappaLoc(10.);
	
	double yObs(1e5);
	double y(1e5),locDisc(10e-4);
	double dydx(0.),d2ydx2(0.);
	double y_x_m_locDisc(15.),y_x_p_locDisc(12.);
	double ratioDc;

	
	
		yObs = ((struct integralParameters*) params)-> yObserv;

		if ( (x >= 0.) && (x<xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc));

			dydx = ( gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x,( ((struct integralParameters*) params)-> accelG1loc)) - gsl_spline_eval(( ((struct integralParameters*) params)-> splineG1loc),x-locDisc,( ((struct integralParameters*) params)-> accelG1loc)) )/locDisc;
			
			kappaLoc = gsl_spline_eval(( ((struct integralParameters*) params)-> kappaRHP_splineG1),x,( ((struct integralParameters*) params)-> kappaRHP_accelG1));
			ratioDc = 1.;
		}
		else if ( (x >= xGrid[addGridPoints+numberGridPoints-1]) )
		{
			y = A_N_R_end+B_N_R_end*x+C_N_R_end*x*x;
			y_x_m_locDisc = A_N_R_end+B_N_R_end*(x-locDisc)+C_N_R_end*(x-locDisc)*(x-locDisc);
			y_x_p_locDisc = A_N_R_end+B_N_R_end*(x+locDisc)+C_N_R_end*(x+locDisc)*(x+locDisc);

			dydx = (y-y_x_m_locDisc)/locDisc;

			d2ydx2 = (y_x_p_locDisc-2.*y+y_x_m_locDisc)/(locDisc*locDisc);

			kappaLoc =  ( d2ydx2/pow((1.+dydx*dydx),1.5) );  // straightforward kappa, bad close to zero
			ratioDc = 1.;
		}
		else 
		{
			cout << " error in integration in testIvantsov relation c " << endl;
			exit(1);
		}
	

	eta = sqrt( (x-xObs)*(x-xObs) + (y-yObs)*(y-yObs) );
 	
 	
	locEq = 1.*(deltaG1-kappaLoc*sigmaLoc);
 	

	if (((pG1*eta) >= 0.00000001) && (fabs(x)>0.000001))
	{
		if ((pG1*eta) >= 4.)
		{
			
			Kernel1sidedPart = locEq*sqrt(PI/(2.*pG1*eta))*exp( -pG1*( (yObs-y)+eta ) )*( ((-dydx*(xObs-x)+(yObs-y))/eta)*(1.+(3./8.)*(1./(pG1*eta))) - 1. );


		}
		else
		{
			Kernel1sidedPart = locEq*exp(-(pG1*(yObs-y)))*( ((-dydx*(xObs-x)+(yObs-y))/eta)*gsl_sf_bessel_K1(pG1*eta) - gsl_sf_bessel_K0(pG1*eta) );

		}
		Kernel = Kernel1sidedPart;
	}

	else 
	{
		Kernel = 0.;
	}

	return Kernel;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void defineFittingParabola(const gsl_vector* x)
{

double discl1_R,discl2_R,discm1_R,discm2_R,discu1_R,discu2_R,y_k_R,y_k_plus_R,y_k_minus_R,x_k_minus_R,x_k_R,x_k_plus_R;
////////////////////////////// for right branch x>0 //////////////////////
	discl1_R = -Discretization.at(numberGridPoints+addGridPoints-3);
	discl2_R = -Discretization.at(numberGridPoints+addGridPoints-3)-Discretization.at(numberGridPoints+addGridPoints-2);
	discm1_R = Discretization.at(numberGridPoints+addGridPoints-3);
	discm2_R = -Discretization.at(numberGridPoints+addGridPoints-2);
	discu1_R = Discretization.at(numberGridPoints+addGridPoints-2)+Discretization.at(numberGridPoints+addGridPoints-3);
	discu2_R = Discretization.at(numberGridPoints+addGridPoints-2);
	y_k_minus_R = gsl_vector_get(x,numberGridPoints+1-4);
	y_k_R =gsl_vector_get(x,numberGridPoints+1-3);
	y_k_plus_R = gsl_vector_get(x,numberGridPoints+1-2);
	x_k_minus_R = xGrid[numberGridPoints+addGridPoints-4];
	x_k_R = xGrid[numberGridPoints+addGridPoints-3];
	x_k_plus_R = xGrid[numberGridPoints+addGridPoints-2];

	A_N_R.at((numberGridPoints-1)) = y_k_minus_R*x_k_R*x_k_plus_R/(discl1_R*discl2_R) + y_k_R*x_k_minus_R*x_k_plus_R/(discm1_R*discm2_R) + y_k_plus_R*x_k_minus_R*x_k_R/(discu1_R*discu2_R);
	B_N_R.at((numberGridPoints-1)) = - y_k_minus_R*(x_k_R+x_k_plus_R) /(discl1_R*discl2_R) -  y_k_R*(x_k_minus_R+x_k_plus_R)/(discm1_R*discm2_R) -  y_k_plus_R*(x_k_minus_R+x_k_R)/(discu1_R*discu2_R);
	C_N_R.at((numberGridPoints-1)) =  y_k_minus_R/(discl1_R*discl2_R) + y_k_R/(discm1_R*discm2_R) + y_k_plus_R/(discu1_R*discu2_R);

	A_N_R_end = A_N_R.at((numberGridPoints-1));
	B_N_R_end = B_N_R.at((numberGridPoints-1));
	C_N_R_end = C_N_R.at((numberGridPoints-1));

	///////////////////////////  for left branch x<0  /////////////////////////
	

	for (int sndGridIndex = 0; sndGridIndex < addGridPoints; ++sndGridIndex)
	{
		prolonguedYGridG1[sndGridIndex] = slopeG1*xGrid[sndGridIndex];
		prolonguedYGridG1[sndGridIndex+numberGridPoints+addGridPoints] = A_N_R_end+B_N_R_end*xGrid[sndGridIndex+numberGridPoints+addGridPoints]+C_N_R_end*xGrid[sndGridIndex+numberGridPoints+addGridPoints]*xGrid[sndGridIndex+numberGridPoints+addGridPoints];
		
		// rechter branch
	}

	for (int sndGridIndex = 0; sndGridIndex < addGridPoints; ++sndGridIndex)
	{
		prolonguedXGridG12[sndGridIndex] = slopeG12*yGridG12[sndGridIndex];
		prolonguedXGridG12[sndGridIndex+numberGridPoints_y+addGridPoints] = gsl_vector_get(x,numberGridPoints);	
		//cout << " extended parts of Grid XG12: " << sndGridIndex << " " << prolonguedXGridG12[sndGridIndex] << " " << prolonguedXGridG12[sndGridIndex+numberGridPoints_y+addGridPoints] << endl; 
		//cout << " extended parts of Grid YG1: " << sndGridIndex << " " << prolonguedYGridG1[sndGridIndex] << " " << prolonguedYGridG1[sndGridIndex+numberGridPoints+addGridPoints] << endl; 
		
		// l1l2-interface 
	}
	for (int sndGridIndex = 0; sndGridIndex < numberGridPoints; ++sndGridIndex)
	{
		prolonguedYGridG1[sndGridIndex+addGridPoints] = gsl_vector_get(x,sndGridIndex);
		//cout << " sndGridIndex+addGridPoints, prolonguedYGridG1[sndGridIndex+addGridPoints] " << sndGridIndex << " " << prolonguedYGridG1[sndGridIndex+addGridPoints] << endl;

	}
	
	for (int sndGridIndex = 1; sndGridIndex < numberGridPoints_y-1; ++sndGridIndex)
	{
		prolonguedXGridG12[sndGridIndex+addGridPoints] = gsl_vector_get(x,numberGridPoints+sndGridIndex);
		//cout << " sndGridIndex+addGridPoints, prolonguedXGridG12[sndGridIndex+addGridPoints] " << sndGridIndex << " " << prolonguedXGridG12[sndGridIndex+addGridPoints] << endl;

	}
	
	
	prolonguedXGridG12[addGridPoints+numberGridPoints_y-1] = gsl_vector_get(x,numberGridPoints); 
	prolonguedYGridG1[addGridPoints] = 0.;
	
	prolonguedXGridG12[addGridPoints] = 0.;
	
	
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void displayFittingParabolaAndDisc()
{
	cout << " show fitting parabolas " << endl;
	cout << " parameters: " << endl;
	cout << " A_N_R_end,B_N_R_end,C_N_R_end: " << A_N_R_end << " " << B_N_R_end << " " << C_N_R_end << endl;
	for (int i = 0; i< numberGridPoints+2*addGridPoints;i++)
	{
		cout << "(i,x,y_prol,y_Iv)G1 " << i << " " << xGridG1[i]  << "," << prolonguedYGridG1[i] << " " << -xGridG1[i]*xGridG1[i]*0.5 << endl;
	}
	for (int i = 0; i< (numberGridPoints_y)+2*addGridPoints;i++)
	{
		cout << "(i,y,x_prol)G12 " << i << " "<< yGridG12[i] << " , " << prolonguedXGridG12[i]  <<  endl;
	}
	cout << " ------------------------------------------------------- " << endl;

}
////////////////////////////////////////////////////////////////
void printChangingParams(gsl_vector* x_solution)
{
	cout << " coefficients of parabolic prolongation " << endl;
	cout << " right branch: ANRend,BNRend,CNRend: " << A_N_R_end << " " << B_N_R_end << " " << C_N_R_end << endl;
	cout << " coefficients of calculated shape last point of each branch/(-xObs*xObs/2) approx 1? " << endl;
	
	cout << " right side: y(xMax)/y_Iv(xMax)" <<  gsl_vector_get(x_solution,numberGridPoints)/(-0.5*(xGridG1[numberGridPoints+addGridPoints-1]*xGridG1[numberGridPoints+addGridPoints-1])) << endl; 
	cout << " xSolution(xMax)= " << gsl_vector_get(x_solution,numberGridPoints) << " xMax= " << xGridG1[numberGridPoints+addGridPoints-1] << endl;
}
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void calc_sep_RHP(vector<double>& s,vector<double>& x,vector<double>&
y,vector<double>& xm,vector<double>& ym,vector<double>& sp)
{
	double dx, dy, ds,x1,y1;
	int n=x.size();
	x1=-x[1];
	y1=y[1];
	ds=s[0];
	dx=x[0]-x1;
	dy=y[0]-y1;
	sp[0]=sqrt(dx*dx+dy*dy);
	xm[0]=dx/sp[0];
	ym[0]=dy/sp[0];


  for(int m=1;m<n;++m)
     {
      ds=s[m]-s[m-1];
      dx=x[m]-x[m-1];
      dy=y[m]-y[m-1];
      sp[m]=sqrt(dx*dx+dy*dy);
      xm[m]=dx/sp[m];
      ym[m]=dy/sp[m];
  
     }
}
////////////////////////////////////////////////////////////////////////////
void calc_sep_LHP(vector<double>& s,vector<double>& x,vector<double>&
y,vector<double>& xm,vector<double>& ym,vector<double>& sp)
{
	double dx, dy, ds,x1,y1;
	int n=x.size();
	x1=-x[1];
	y1=y[1];
	ds=s[0];
	dx=x[0]-x1;
	dy=y[0]-y1;
	sp[0]=sqrt(dx*dx+dy*dy);
	xm[0]=dx/sp[0];
	ym[0]=dy/sp[0];


  for(int m=1;m<n;++m)
     {
      ds=s[m]-s[m-1];
      dx=x[m]-x[m-1];
      dy=y[m]-y[m-1];
      sp[m]=sqrt(dx*dx+dy*dy);
      xm[m]=dx/sp[m];
      ym[m]=dy/sp[m];
  
     }

}

////////////////////////////////////////////////////////////////////////////////////
void calc_curv_HMK_RHP(vector<double>& s,vector<double>& x,vector<double>&
y, vector<double>& nx, vector<double>& ny, vector<double>& curv)
{
	
  int n=x.size();
  vector<double> xm(n),ym(n),sp(n);
  
  double bunsen;  
  double bunshi;
  double bbx, bby;
  double bunbo;
	
  calc_sep_RHP(s,x,y,xm,ym,sp);
	
  for(int m=0;m<n-1;m++)
     {
      bunsen=xm[m]*xm[m+1]+ym[m]*ym[m+1];
      bunshi=xm[m]*ym[m+1]-ym[m]*xm[m+1];
      bbx=-sp[m]*ym[m+1]-sp[m+1]*ym[m];
      bby=sp[m]*xm[m+1]+ sp[m+1]*xm[m];

      bunbo = sqrt(bbx*bbx + bby*bby);
 
      curv[m]=-2.*bunshi/bunbo;
 
     if ((bunbo==0))
     {
	cout << " RHP/m=" << m << " bbx,bby= " << bbx << "," << bby << endl; 
	cout << " bunsen, bunshi= " << bunsen << " " << bunshi << endl;
	cout << " m=0: sp,ym,xm " << sp[m] << " " << ym[m] << " " << xm[m] << endl;
	cout << " m+1: sp,ym,xm " << sp[m+1] << " " << ym[m+1] << " " << xm[m+1] << endl; 
     	cout << "nan curv. Warning!\n";
 
        exit(1);
     }
     }
 
 
}
////////////////////////////////////////////////////////////////////////////////////
void calc_curv_HMK_G12(vector<double>& s,vector<double>& x,vector<double>&
y, vector<double>& nx, vector<double>& ny, vector<double>& curv)
{
	
  int n=x.size();
  vector<double> xm(n),ym(n),sp(n);
  
  double bunsen;  
  double bunshi;
  double bbx, bby;
  double bunbo;
	
  calc_sep_LHP(s,x,y,xm,ym,sp);
	
  for(int m=0;m<n-1;m++)
     {
      bunsen=xm[m]*xm[m+1]+ym[m]*ym[m+1];
      bunshi=xm[m]*ym[m+1]-ym[m]*xm[m+1];
      bbx=-sp[m]*ym[m+1]-sp[m+1]*ym[m];
      bby=sp[m]*xm[m+1]+ sp[m+1]*xm[m];

      bunbo = sqrt(bbx*bbx + bby*bby);
 
      curv[m]=2.*bunshi/bunbo;
 
       if (bunbo==0)
        {
	cout << " LHP/m=" << m << " bbx,bby= " << bbx << "," << bby << endl; 
	cout << " bunsen, bunshi= " << bunsen << " " << bunshi << endl;
	cout << " m=0: sp,ym,xm " << sp[m] << " " << ym[m] << " " << xm[m] << endl;
	cout << " m+1: sp,ym,xm " << sp[m+1] << " " << ym[m+1] << " " << xm[m+1] << endl; 
        cout << "nan curv. Warning!\n";
 
        exit(1);
        }
     }
 
 
}
/////////////////////////////////////////////////////////////////////////////
float sgn(float val)
{
double outValue(10.);
    if (val<0)
    {
        outValue = -1.;
    }
    else if (val>0)
    {
        outValue= 1.;
    }
    else if (val==0)
    {
        outValue= 1.;
    }

	return outValue;
} 
////////////////////////////////////////////////////
// void saveResultsInFile(string& filename, gsl_vector* x)
// {
// 
// 
// }
/////////////////////////////////////////////////////////////////////////
void createDirectory()
{
	char commandString[255] = "mkdir ";
	char argumentString[255] = "peritecticRunsPrepareForPaper";
	
	strcat( commandString, argumentString);
	cout << commandString << endl;
	system( commandString );	
}
///////////////////////////////////////////////////////////////////////////

