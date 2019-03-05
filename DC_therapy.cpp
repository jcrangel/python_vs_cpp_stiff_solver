
#include "murinemodelv2.h"
#include "util.h"


//template <class T>
struct push_back_state_and_time2
{
    std::vector< std::vector<double> >& m_states;
    std::vector<double> & m_times;

    push_back_state_and_time2( std::vector< std::vector<double> > &states, std::vector<double> &times )
        : m_states( states ), m_times( times ) { }

    void operator()( const std::vector<double> &x, double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
	void operator()(const vector_type &x, double t)
	{
		//std::vector<double> v(x.size());
		//std::copy(x.begin(), x.end(), v.begin());
		m_states.push_back(toStdVectorD(x));
		m_times.push_back(t);
	}
};
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
void integrateSystem(T &system, StateType initialCondition, Doub t0, Doub tf,
                    std::vector<StateType> &u, StateType &t) {
	typedef runge_kutta_dopri5<StateType> error_stepper_type;

	size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
		system, initialCondition, t0, tf, 0.001,push_back_state_and_time(u, t));
}

int main( int argc , char **argv )
{
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
    t0 = 0.0 ;
    delay = 232;
    tf = 168 + delay;
    ef=0.05;
    std::vector<StateType> u;
    StateType tt;
    integrateSystem(fun, initialCondition, t0,  tf, u, tt)

    //]
    //clog << num_of_steps << end;
    for( size_t i=0; i<= tt.size() - 1 ; i++ )
    {
        cout<< tt[i] << '\t' << u[i][0] << '\t' << u[i][1] << '\n';
    }
    cout<< x(0) << " " << x(1)<< endl;


	cin.get();
}
