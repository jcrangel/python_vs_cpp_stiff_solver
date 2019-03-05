
#include "murinemodelv2.h"
#include "util.h"
#include <chrono>
#include <ctime>

using namespace boost::numeric::odeint;

struct push_back_state_and_time
{
	std::vector< StateType>& m_states;
	std::vector<Doub> & m_times;

	push_back_state_and_time(std::vector< StateType > &states, std::vector<Doub> &times)
		: m_states(states), m_times(times) { }

	void operator()(const StateType &x, Doub t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}

	void operator()(const VectorBoost &x, Doub t)
	{
		m_states.push_back(toStateType(x));
		m_times.push_back(t);
	}
};


template <class T>
void integrateSystem(T &system, StateType &initialCondition, Doub t0, Doub tf,
                    std::vector<StateType> &u, StateType &t) {
	typedef runge_kutta_dopri5<StateType> error_stepper_type;

	size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
		system, initialCondition, t0, tf, 0.001,push_back_state_and_time(u, t));
}

int main( int argc , char **argv )
{
	Doub aH, MuH, rH, KH, aC, MuC, rC, KC, KT, MuD, rI, MuIC, MuI;
	Doub rT, aT, eT, hT, aTBeta, eTBeta, rTBeta, MuBeta, aGammaC;
	Doub MuGamma, gMl, aMlGamma, eMlGamma, MuMl;
    aH = 1e-4;
	MuH = 0.005;
	rH = 10e-2;
	KH = 1;
	aC = 1e-4;
	MuC = 0.01925;
	rC = 0.00004e-2;
	KC = 1;
	KT = 1e12;
	MuD = 0.009625;
	rI = 1e-2;
	MuIC = 1e-7;
	MuI = 1e-2;
	rT = 0.002;
	aT = 0.1136;
	eT = 50;
	hT = 5.2e5;
	aTBeta = 0.69;
	eTBeta = 1e4;
	rTBeta = 5.57e-6;
	MuBeta = 6.93;
	aGammaC = 1.02e-4;
	MuGamma = 0.102;
	gMl = 1.44;
	aMlGamma = 2.89;
	eMlGamma = 3.38e5;
	MuMl = 0.0144;
	std::vector<double> PARAMETERS = { aH, MuH, rH, KH, aC, MuC, rC, KC, KT, MuD, rI, MuIC,
									MuI, rT, aT, eT, hT, aTBeta, eTBeta, rTBeta, MuBeta, aGammaC, MuGamma, gMl, aMlGamma, eMlGamma, MuMl };

    MurineModelv2 fun(PARAMETERS);
    StateType initialCondition = {6e4, 0, 0, 0, 0, 0, 0, 0};
    Doub t0 = 0.0 ;
    Doub delay = 232;
    Doub tf = 168 + delay;
    Doub ef=0.05;
    std::vector<StateType> u1,u2,u3,u4;
    StateType tt1,tt2,tt3,tt4;
	// Integration starts
	auto start = std::chrono::system_clock::now();
    integrateSystem(fun, initialCondition, t0,  tf, u1, tt1);
	//First inyection 
	initialCondition[3]+=1e6 * ef; // initialCondition holds the last value i.e. x(tf)
	t0 = tf ;
	tf = tf + 168;  
	integrateSystem(fun, initialCondition, t0,  tf, u2, tt2);
	//Second Inyection
	initialCondition[3]+=1e6 * ef;
	t0 = tf ;
	tf = tf + 168;  
	integrateSystem(fun, initialCondition, t0,  tf, u3, tt3);
	//Third inyection
	initialCondition[3]+=1e6 * ef;
	t0 = tf ;
	tf = 1400;  
	integrateSystem(fun, initialCondition, t0,  tf, u4, tt4);
	

	auto end = std::chrono::system_clock::now();



	std::ofstream out("points.txt",std::ios::trunc);

	for (int i = 0; i <= tt1.size() - 1; i++)
	{
		out<< tt1[i] <<"," << u1[i][0] << '\n';
	}
	for (int i = 0; i <= tt2.size() - 1; i++)
	{
		out<< tt2[i] <<"," << u2[i][0] << '\n';
	}
	for (int i = 0; i <= tt3.size() - 1; i++)
	{
		out<< tt3[i] <<"," << u3[i][0] << '\n';
	}
	for (int i = 0; i <= tt4.size() - 1; i++)
	{
		out<< tt4[i] <<"," << u4[i][0] << '\n';
	}

	std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout <<"Integration took :" << elapsed_seconds.count() <<"sec \n";
	// std::cin.get();
}
