/*
Util.h

TODO: All common headers, google styles says this is not good

Common data type definition and functions for vector manipulation.

*/

#ifndef ALL_H
///< .
#define ALL_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <iomanip>
#include <sstream>
#include <string>

//Eigen lib
#include <Eigen/Dense>
//Boost lib
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/odeint.hpp>

//#include <boost/math/interpolators/barycentric_rational.hpp> for not aproximate jacobian

#define DEBUG0 true
#define DEBUG1 false

/**********************************************************************************************//**
 * @typedef	std::vector<Doub> VecDoub_IO
 *
 * @brief	Defines an alias representing the vector Doub i/o
 **************************************************************************************************/
typedef double Doub; // default floating type useful if we need to change to long double
typedef std::vector< Doub > StateType;
typedef boost::numeric::ublas::vector< Doub > VectorBoost;
typedef boost::numeric::ublas::matrix< Doub > MatrixBoost;
typedef Eigen::VectorXd VectorEigen;
typedef Eigen::VectorXcd VectorcEigen; //Complex version
typedef Eigen::MatrixXd MatrixEigen;
typedef const std::vector<Doub> VecDoub_I; // Utility vector for input
typedef std::vector<Doub>  VecDoub_IO;// Utility vector for input & ouput
/**********************************************************************************************//**
 * @class	fixPoint, The class representing a fixpoint [x,y,z] and if info about
//convergence and stability
 *
 * @brief	A fix point.
 *
 * @author	Iron
 * @date	7/25/2018
 **************************************************************************************************/

class FixPoint
{
public:
	bool convergent;
	bool stability;
	VectorEigen solution;

	FixPoint(bool convergent_, bool stability_, VectorEigen solution_) :
		convergent(convergent_),
		stability(stability_),
		solution(solution_) {}
};

struct point {
	Doub x=0;
	Doub y=0;

	point(Doub x, Doub y) : x(x), y(y) {};

};
/*
*Basically a square 
*/
class cell2D {
public:
	point p1;
	point p2;
	point p3;
	point p4;

	cell2D(point p1, point p2, point p3, point p4) :
		p1(p1), p2(p2), p3(p3), p4(p4) {}

	friend std::ostream& operator << (std::ostream& out, const cell2D &cell) {
		out << cell.p1.x << "," << cell.p1.y << "\n";
		out << cell.p2.x << "," << cell.p2.y << "\n";
		out << cell.p3.x << "," << cell.p3.y << "\n";
		out << cell.p4.x << "," << cell.p4.y << "\n";
		return out;
	}

};
/* The container for the critical cells*/
typedef std::vector<cell2D> setCriticalCells;

/*
Converts to doubles to string and gives the concatenation. To use for
the index in a "multiindex" map:
std::map<std::string, int> g;
double x = 3.1416123254;
double y = x/10;
g.insert(std::pair<std::string,int> (strtokey(x,y),1 ) );
*/

std::string strtokey(double a, double b, int precision = 5) {
	std::stringstream stream;
	stream << std::fixed << std::setprecision(precision) << a;
	std::string as = stream.str();
	//https://stackoverflow.com/questions/20731/how-do-you-clear-a-stringstream-variable
	stream.str(std::string());
	stream << std::fixed << std::setprecision(precision) << b;
	std::string bs = stream.str();
	return as + bs;
}

/**********************************************************************************************//**
 * @fn	template <class T> void transpose(const T u, std::vector < StateType > & state )
 *
 * @brief	Transpose a matrix of made with std::vector's. Not being used in the current version.
 *
 * @author	Iron
 * @date	7/31/2018
 *
 * @tparam	T	Generic type parameter.
 * @param 		  	u	 	A T to process.
 * @param [in,out]	state	The state.
 **************************************************************************************************/

template <class T>
void transpose(const T u, std::vector < StateType > & state)
{
	//Copy the data as transpose
	//for some reason in windows we need i < u.size in ubuntu ins i <= u.size
	for (int i = 0; i < u.size(); i++)
		// STATE_SIZE before, posibly a source of bugs
		for (int j = 0; j < state.size(); j++)
			state[j][i] = u[i][j];
}

/**********************************************************************************************//**
 * @fn	MatrixEigen reshapeVectorToMatrix(const T x)
 *
 * @brief	Reshape a vector Size N^2 into a Matrix NxN,vector is assumed to come in row order
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	x	A T to process.
 *
 * @return	A MatrixEigen.
 **************************************************************************************************/
template <class T>
MatrixEigen reshapeVectorToMatrix(const T x)
{
	int N = sqrt(x.size());
	MatrixEigen A(N, N);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A(i, j) = x[i*N + j];
		}
	}

	return A;
}

//
// the matrix is stable()true.

/**********************************************************************************************//**
 * @fn	bool isStable(MatrixEigen A)
 *
 * @brief	Query if Matrix 'A' is stable. Calculate all eigenvalues of matrix A, if all of them
 * 			are less than 1, is stable.
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	A	A MatrixEigen to process.
 *
 * @return	True if stable, false if not.
 **************************************************************************************************/

bool isStable(MatrixEigen A)
{
	VectorcEigen eivals = A.eigenvalues();
	for (int i = 0; i < eivals.size(); ++i)
		if (std::abs(eivals[i]) >= 1)
			return false;

	return true;
};

/**********************************************************************************************//**
 * @fn	bool equal(Doub A, Doub B, Doub epsilon = 0.000005f)
 *
 * @brief	Cheks for aproximate  equality.
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	A	   	A Doub to process.
 * @param	B	   	A Doub to process.
 * @param	epsilon	(Optional) The epsilon.
 *
 * @return	True if it succeeds, false if it fails.
 **************************************************************************************************/

bool equal(Doub A, Doub B, Doub epsilon = 0.000005f)
{
	return (fabs(A - B) < epsilon);
}

/**********************************************************************************************//**
 * @fn	std::vector<Doub> toStdVectorD(const Vector3d v)
 *
 * @brief	Get a std vector from eigen Vector3d
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	v	A Vector3d to process.
 *
 * @return	V as a std::vector&lt;Doub&gt;
 **************************************************************************************************/

StateType toStateType(const VectorEigen v)
{
	std::vector<Doub> v2;
	v2.resize(v.size());
	VectorEigen::Map(&v2[0], v.size()) = v;
	return v2;
}

/**********************************************************************************************//**
 * @fn	void toStdVectorD(const Vector3d v, StateType &w)
 *
 * @brief	Copy a Vector3d into a std vector
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param 		  	v	A Vector3d to process.
 * @param [in,out]	w	A StateType to process.
 **************************************************************************************************/

void toStateType(const VectorEigen v, StateType &w)
{
	w.resize(v.size());
	VectorEigen::Map(&w[0], v.size()) = v;
}

/**********************************************************************************************//**
 * @fn	std::vector<Doub> toStdVectorD(VectorBoost v)
 *
 * @brief	Converts a vector drom Boost libs to a standard vector.
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	v	A VectorBoost to process.
 *
 * @return	V as a std::vector&lt;Doub&gt;
 **************************************************************************************************/

StateType toStateType(VectorBoost v)
{
	std::vector<Doub> w(v.size());
	std::copy(v.begin(), v.end(), w.begin());
	return w;
}

/**********************************************************************************************//**
 * @fn	VectorBoost toBoostVectorD(const StateType v)
 *
 * @brief	Converts a standard vector to a boost vector.
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	v	A StateType to process.
 *
 * @return	V as a VectorBoost.
 **************************************************************************************************/

VectorBoost toBoostVectorD(const StateType v)
{
	VectorBoost x(v.size());
	std::copy(v.begin(), v.end(), x.begin());
	return x;
}

/**********************************************************************************************//**
 * @fn	Vector3d toEigenVector(const StateType v)
 *
 * @brief	Converts a v to an eigen vector
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	v	A StateType to process.
 *
 * @return	V as a Vector3d.
 **************************************************************************************************/

VectorEigen toEigenVector(StateType v) {
	//Doub* ptr = &v[0];
	//Eigen::Map< Eigen::VectorEigen> v2(ptr, v.size());
	////Vector3d v2(v.data());
	//return v2;
	return VectorEigen::Map(v.data(), v.size());
}

/**********************************************************************************************//**
 * @fn	bool pointIsInSet(fixPoint p, std::vector<fixPoint> S)
 *
 * @brief	Check if the fixpoint is in the Set S .
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	p	A fixPoint to process.
 * @param	S	A std::vector&lt;fixPoint&gt; to process.
 *
 * @return	True if it succeeds, false if it fails.
 **************************************************************************************************/

bool pointIsInSet(FixPoint p, std::vector<FixPoint> S)
{
	for (FixPoint i : S) {
		//Check
		int j;
		//STATE_SIZE before
		for (j = 0; j < p.solution.size(); j++) {
			if (!equal(i.solution[j], p.solution[j]))
				break;	//Stop comparing, vectors ain't equal
		}
		//STATE_SIZE before
		if (j == p.solution.size()) // never break, then all were equal
			return true;
	}
	return false;
}

/**********************************************************************************************//**
 * @fn	bool pointHaveNegatives(fixPoint p)
 *
 * @brief	Check if the fix point have negatives values
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	p	A fixPoint to process.
 *
 * @return	True if it succeeds, false if it fails.
 **************************************************************************************************/

bool pointHaveNegatives(FixPoint p) {
	for (int i = 0; i < p.solution.size(); i++)
		if (p.solution[i] < 0)
			return true;

	return false;
}

/**********************************************************************************************//**
 * @class	LogAndStdout
 *
 * @brief	A log and stdout. Warper for writing both to the console and to a file.
 *
 *
 * @author	Iron
 * @date	7/25/2018
 **************************************************************************************************/

class LogAndStdout {
private:
	std::ofstream file;
public:

	/**********************************************************************************************//**
	 * @fn	TODO::LogAndStdout(const std::string fileName)
	 *
	 * @brief	Opens the file for writing and appending.
	 *
	 * @author	Iron
	 * @date	7/25/2018
	 *
	 * @param	fileName	Filename of the file.
	 **************************************************************************************************/

	LogAndStdout(const std::string fileName) {
		file.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);
	}

	template<class T>
	LogAndStdout& operator<<(const T data) {//TODO: this dont work with std::endl
		std::cout << data;
		file << data;

		/*TODO : This clone the object? if this is true, then this is bad */
		return *this;
	}
};

/**********************************************************************************************//**
 * @fn	template <class T> VectorEigenevalFunInLast(T &functionName, StateType initialCondition, Doub tau, Doub d)
 *
 * @brief	Evaluate F(x,y,z)  at the final time.
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @return	A VectorEigen.
 **************************************************************************************************/

template <class T>
VectorEigen evalFunInLast(T &functionName, StateType initialCondition, Doub tau, Doub d)
{
	int controlIndex = functionName.getControlIndex();
	initialCondition[controlIndex] = initialCondition[controlIndex] + d;
	int N = functionName.getSystemSize();
	std::vector<Doub> res(N);
	integrateSystem(functionName, initialCondition, res, 0, tau);

	//To create the VectorEigenfrom the std::vector
	//Doub* ptr = &res[0];
	//Eigen::Map<VectorEigen> v(ptr, N);

	return toEigenVector(res);
	//initialCondition has the last
}

/**
 *Has a jump of size d at the start 
 */
template <class T>
StateType evalFunInLastJump(T &functionName, StateType initialCondition, Doub t0, Doub tf,Doub d)
{
	int controlIndex = functionName.getControlIndex();
	initialCondition[controlIndex] = initialCondition[controlIndex] + d;
	int N = functionName.getSystemSize();
	std::vector<Doub> res(N);
	integrateSystem(functionName, initialCondition, res, t0, tf);

	//To create the VectorEigenfrom the std::vector
	//Doub* ptr = &res[0];
	//Eigen::Map<VectorEigen> v(ptr, N);

	return res;
	//initialCondition has the last
}



/**********************************************************************************************//**
 * @fn	std::vector< std::vector<T> > cartesianProduct(const std::vector<T> A, const std::vector<T> B)
 *
 * @brief	Create the cartesian set A x B
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @param	A	Set A.
 * @param	B	Set B.
 *
 * @return	A std::vector&lt;std::vector&lt;T&gt; &gt;
 **************************************************************************************************/
template<class T>
std::vector< std::vector<T> > cartesianProduct(const std::vector<T> A, const std::vector<T> B)
{
	std::vector<std::vector<T> > prod;
	for (int i = 0; i < A.size(); i++)
		for (int j = 0; j < B.size(); j++)
			prod.push_back(std::vector<T> {A[i], B[j]});

	return prod;
}

/**********************************************************************************************//**
 * @fn	template<class T> std::vector< std::vector<T> > cartesianProduct(const std::vector<std::vector<T>> A, const std::vector<T> B)
 *
 * @brief	Create the cartesian set A x B = (a x b) x B, where A is already a cartesian set.
 * 			Then to create the cartesian set of {x,y,z} x {x,y,z} x {x,y,z}. We have to : s = {x,
 * 			y,z};
 * 			std::vector&lt;std::vector&lt;Doub&gt; &gt; prod = cartesianProduct(s,s);
 * 			prod = cartesianProduct(prod,s);
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @tparam	T	Generic type parameter.
 * @param	A	A std::vector&lt;std::vector&lt;T&gt;&gt; to process.
 * @param	B	A std::vector&lt;T&gt; to process.
 *
 * @return	A std::vector&lt;std::vector&lt;T&gt; &gt;
 **************************************************************************************************/

template<class T>
std::vector< std::vector<T> > cartesianProduct(const std::vector<std::vector<T>> A, const std::vector<T> B)
{
	std::vector<std::vector<T> > prod;
	for (int i = 0; i < A.size(); i++)
		for (int j = 0; j < B.size(); j++) {
			std::vector<T> AB;
			AB.reserve(A.size() + B.size());
			AB.insert(AB.end(), A[i].begin(), A[i].end());
			AB.insert(AB.end(), B[j]);

			prod.push_back(AB);
		}
	return prod;
}

template<class T>
void printVector(std::vector<T> M)
{
	for (T j : M)
		std::cout << j << " ";
	std::cout << "\n";
}

/**********************************************************************************************//**
 * @fn	template<class T> void printVectorVector(std::vector< std::vector<T> > M)
 *
 * @brief	Print vector of vectors
 *
 * @author	Iron
 * @date	7/25/2018
 *
 * @tparam	T	Generic type parameter.
 * @param	M	A std::vector&lt;std::vector&lt;T&gt;&gt; to process.
 **************************************************************************************************/

template<class T>
void printVectorVector(std::vector< std::vector<T> > M)
{
	for (std::vector<T> i : M) {
		for (T j : i)
			std::cout << j << " ";
		std::cout << "\n";
	}
}

/**********************************************************************************************/ /**
 * @fn	template<class T> void nextPointSubdomain(std::vector<Doub> &b, Doub min, Doub max, Doub step)
 *
 * @brief	This function ge the next subdomain point saved in b. Simulates to have an N-Dimensional
 * 			for, i.e N-for's nested iterating i=min to max with increments of step
 **************************************************************************************************/

void nextPointSubdomain(std::vector<Doub> &b, Doub min, Doub max, Doub step)
{
	int i = 0;

	while (equal(b[i], max)) // do it with equals
	{
		b[i] = min;
		if (++i >= b.size())
			return; // prevents to go further of array size
	}
	b[i] += step;
}
/**********************************************************************************************/ /**
 *
 * @brief	This function ge the next subdomain point saved in b. Simulates to have an N-Dimensional
 * 			for, i.e N-for's nested iterating i=min to max with increments of step. This work for
 * 			diferent values for each index given by the vectors xmin and xmax. See notes p.18
 **************************************************************************************************/

template <class T>
void nextPointSubdomain(T &b, T xmin, T xmax)
{
	int i = 0;
	double step;
	while (equal(b[i], xmax[i])) // do it with equals
	{
		b[i] = xmin[i];
		if (++i >= b.size())
			return; // prevents to go further of array size
	}
	step = (xmax[i] - xmin[i]) / 2;
	b[i] += step;
}

#endif
