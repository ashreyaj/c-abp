#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <chrono>
#include <ctime>
#include <cstdlib>

using namespace std;

// ---- Helper functions ----
template<typename T>
vector<double> linspace(T start_in, T end_in, int num_in)
{

  vector<double> ls;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in+1);

  if (num == 1) {return ls;}
  if (num == 2) 
    {
      ls.push_back(start);
      return ls;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      ls.push_back(start + delta * i);
    }

  return ls;
}

template <typename T>
string to_string_with_precision(const T a_value, const int n = 6)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return std::move(out).str();
}

// ---- Parameters ----
const double PI = 3.14159265359;
const double kT = 1.0;

string outdir;
double Phi;
int N;
double sig;
double rcut;
double eps;
double frcut;
double al;
double vexp;
double v0;
double D0;
double taur;
double rcon;
double width;
double f;
double run;
double cycle;
double dt;
double t0;
double l0;

double rcut2;
double rin2;
double rout2;
double frcut2;
double circumference;
double tdiff;
double rdiff;
int cycle_steps;
int run_steps;

vector<double> x;
vector<double> y;
vector<double> phi;
vector<vector<double> > fwca; 
vector<int> neigh;

void readParameters(const char* ParameterFile)
{
	string line;
	string param;
	string value;
	ifstream IN(ParameterFile);
	while(getline(IN,line))
	{
		stringstream DATA(line);
		DATA >> param >> value;
		if (param=="outdir") {outdir=value;}
		else if (param=="phi") {Phi=stod(value);}
		else if (param=="sig") {sig=stod(value);}
		else if (param=="rcut") {rcut=stod(value);}
		else if (param=="eps") {eps=stod(value);}
		else if (param=="frcut") {frcut=stod(value);}
		else if (param=="al") {al=stod(value);}
		else if (param=="v0") {vexp=stod(value);}
		else if (param=="D0") {D0=stod(value);}
		else if (param=="taur") {taur=stod(value);}
		else if (param=="rcon") {rcon=stod(value);}
		else if (param=="width") {width=stod(value);}
		else if (param=="f") {f=stod(value);}
		else if (param=="run") {run=stod(value);}
		else if (param=="cycle") {cycle=stod(value);}
		else if (param=="dt") {dt=stod(value);}
		else if (param=="t0") {t0=stod(value);}
		else if (param=="l0") {l0=stod(value);}
	}
	IN.close();
	
	// derived parameters
	v0 = vexp*1e-6*t0/l0;
	rcut2 = rcut*rcut;
	frcut2 = frcut*frcut;
	rin2 = rcon*rcon;
	rout2 = (rcon+width)*(rcon+width);
	circumference = 2*PI*rcon;
	tdiff = sqrt(3*2*D0*dt);
	rdiff = sqrt(3*2*dt/taur);
	run_steps = static_cast<unsigned long int>(run/dt);
	cycle_steps = static_cast<unsigned long int>(cycle/dt);
}

// --- Initialize particles ----
void initialize(int N, double rcon)
{
	default_random_engine generator;
	uniform_real_distribution<double> theta_distribution(0.0,2.0*PI);

	double N_full = static_cast<int>(1.0*circumference/sig);
	if (N>N_full) // overpacking
	{
		std::vector<double> ang = linspace(0.0, 2*PI, N_full);
		for (int i=0;i<N_full;i++)
		{
			x.push_back(rcon*cos(ang[i]));
			y.push_back(rcon*sin(ang[i]));
			phi.push_back(theta_distribution(generator));
		}
		int N_rem = N - N_full;
		std::vector<double> ang2 = linspace(0.0, 2*PI, N_rem);
		for (int i=0;i<N_rem;i++)
		{
			x.push_back((rcon+sig)*cos(ang2[i]));
			y.push_back((rcon+sig)*sin(ang2[i]));
			phi.push_back(theta_distribution(generator));
		}
	}

	else
	{
		std::vector<double> ang = linspace(0.0, 2*PI, N);
		for (int i=0;i<N;i++)
		{
			x.push_back(rcon*cos(ang[i]));
			y.push_back(rcon*sin(ang[i]));
			phi.push_back(theta_distribution(generator));
		}
	}
}

// ---- Confinement ----
pair<double, double> confine(double xp, double yp)
{
	double fx = 0;
	double fy = 0;
	double theta;
	if((xp*xp + yp*yp) > rout2)
	{
		theta = atan2(yp,xp);
		fx -= f*cos(theta);
		fy -= f*sin(theta);
	}
	else if ((xp*xp + yp*yp) < rin2)
	{
		theta = atan2(yp,xp);
		fx += f*cos(theta);
		fy += f*sin(theta);
	}
	return make_pair(fx,fy);
}

// ---- WCA potential and neighbours ----
void wca_neigh()
{
	for(int i=0;i<N;i++)
	{
		fwca[i][0]=0.0; // reset forces and neighbours to 0
		fwca[i][1]=0.0;
		neigh[i] = 0;
	} 
	vector<double> dr(2, 0.0);
	for(int i=0;i<N;i++)
	{
		for (int j=0;j<i;j++)
		{
			// WCA
			dr[0] = x[i]-x[j];
			dr[1] = y[i]-y[j];

			double r2 = (dr[0]*dr[0]) + (dr[1]*dr[1]);
			if (r2<rcut2)
			{
				double r6i = 1.0 /(r2 * r2 * r2);
				double fij = 48 * eps * (r6i * r6i - 0.5 * r6i);
				for (int k=0;k<2;k++)
				{
					fwca[i][k] += dr[k] * (fij / r2);
					fwca[j][k] -= dr[k] * (fij / r2);
				}
				
			}

			// neighbours
			if (r2<frcut2)
			{
				neigh[i] += 1;
				neigh[j] += 1;
			}
		}
	}
}

// ---- Integrate ----
void step()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	uniform_real_distribution<double> uniform(-1.0,1.0);

	wca_neigh();
	for (int i=0;i<N;i++)
	{	
		pair<double,double> fcon = confine(x[i],y[i]);
		vector<double> fc = fwca[i];
		double Dfact = 1.0 - al * neigh[i] / (4*(frcut2-0.25));
		
		if (Dfact<0)
		{
			cout << "Negative diffusion coefficient. alpha is too large." << endl;
			exit(0);
		}

		double dx = v0*cos(phi[i])*dt + D0*fcon.first*dt/kT + D0*Dfact*fc[0]*dt/kT + tdiff*sqrt(Dfact)*uniform(gen);
		double dy = v0*sin(phi[i])*dt + D0*fcon.second*dt/kT + D0*Dfact*fc[1]*dt/kT + tdiff*sqrt(Dfact)*uniform(gen);

		if (dx>10 || dy>10)
		{
			cout << "Particle moved a lot in one time step." << endl;
			exit(0);
		}

		x[i] = x[i] + dx;
		y[i] = y[i] + dy;
		phi[i] = phi[i] + rdiff*uniform(gen);
	}
}

// ---- Write data ----
void write(ofstream& OUT)
{
	for(int j=0;j<N;j++)
	{
		OUT << x[j] << "\t" << y[j] << "\t" << phi[j] << endl;
	}
}

// ---- Simulate ----
int main(int argc, char *argv[])
{
	readParameters(argv[1]);
	N = static_cast<int>(Phi*circumference/sig);
	fwca.resize(N, std::vector<double>(2, 0.0));
	neigh.resize(N);
	
	srand(unsigned(time(0)));
	initialize(N,rcon);
	cout << outdir << endl;
	ofstream OUT(outdir+"/"+"v"+to_string_with_precision(vexp,1)+"phi"+to_string_with_precision(Phi,1)+"f"+to_string_with_precision(f,1)+"al"+to_string_with_precision(al,1)+"_0.dat");
	write(OUT);
	OUT.close();

	cout << cycle_steps << "\t" << run_steps << endl;	

	// reach steady state
	int steady_steps = static_cast<unsigned long int>(2/dt);
	for(int i=0;i<steady_steps;i++)
	{
		step();
	}
	
	// record data
	int j=0;
	for(int i=0;i<run_steps;i++)
	{
		if((i%cycle_steps)==0)
		{
			ofstream OUT(outdir+"/"+"v"+to_string_with_precision(vexp,1)+"phi"+to_string_with_precision(Phi,1)+"f"+to_string_with_precision(f,1)+"al"+to_string_with_precision(al,1)+"_"+to_string(j)+".dat");
			write(OUT);
			OUT.close();
			j+=1;
		}
		step();
	}
}
