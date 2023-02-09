

#include "droplet_diffusion.h"
#include "config_file.h"
#include "radius_frequency.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>


using namespace std;


struct FunctorC0Out {
  FunctorC0Out(double a, double b, double c, double d, double x0)
    : a(a), b(b), c(c), d(d), x0(x0) {}

  double operator() (double x, double y, double z, double t) {
    return a + b * x + c * tanh( (x - x0) / d ) ; 
  }

  double a, b, c, d, x0;
};


class RadiusDistribution {
 public:
  RadiusDistribution(double Rmin, double Rmax)
    : Rmin(Rmin), Rmax(Rmax) {}

  double operator() (double random_uniform_number) {
    return Rmin + (Rmax - Rmin) * random_uniform_number;
  }

  double Rmin, Rmax;
};

double get_mean(const vector<double>& radii)
{
  double mean = 0.0;
  for (unsigned int i = 0; i < radii.size(); ++i) {
    mean += radii[i];
  }
  return mean / radii.size(); 
}


int main()
{

  // read configuration file
  Config params("input.txt");
 

  // set parameter values 
  double D = params.get_parameter<double>("D");
  double Dd1 = params.get_parameter<double>("Dd1");
  double l_gamma = params.get_parameter<double>("l_gamma");
  double c0_in = params.get_parameter<double>("c0_in");
  double c_out = params.get_parameter<double>("c_out");
  double Rco = params.get_parameter<double>("Rco");
  unsigned int number_of_droplets =
      params.get_parameter<unsigned int>("number_of_droplets");
  double Lx = params.get_parameter<double>("Lx");
  double Ly = params.get_parameter<double>("Ly");
  double Lz = params.get_parameter<double>("Lz");
  unsigned int Nx =params.get_parameter<unsigned int>("Nx");
  unsigned int Ny =params.get_parameter<unsigned int>("Ny");
  unsigned int Nz =params.get_parameter<unsigned int>("Nz");
  double dt = params.get_parameter<double>("dt");
  unsigned int seed =params.get_parameter<unsigned int>("seed");
  double integration_time =
      params.get_parameter<double>("integration_time");
  double save_every =
      params.get_parameter<double>("save_every");


  double c0_out_a = params.get_parameter<double>("a");
  double c0_out_b = params.get_parameter<double>("b");
  double c0_out_c = params.get_parameter<double>("c");
  double c0_out_d = params.get_parameter<double>("d");
  double c0_out_x0 = params.get_parameter<double>("x0");

  // functor for the equilibrium concentration outide of the drop
  FunctorC0Out f_c0_out(c0_out_a, c0_out_b, c0_out_c,
                        c0_out_d, c0_out_x0);


  // simulation object
  DropletDiffusion<FunctorC0Out>
        droplet_diffusion(D, Dd1, l_gamma, c0_in, c_out,
            f_c0_out, Rco, Lx, Ly, Lz, dt, Nx, Ny, Nz, seed);


  double epsilon = (c_out / c0_out_a) - 1;
  double Rcrit = l_gamma / epsilon;
  Rco = Rcrit / 10.0;
  RadiusDistribution radius_distribution(0.9*Rcrit, 1.5 * Rcrit);
  cout << radius_distribution.Rmin << endl;
  cout << radius_distribution.Rmax << endl;
  cout << Rco << endl;
  cout << Rcrit << endl;

  droplet_diffusion.InitializeDroplets(number_of_droplets, radius_distribution);


  cout << "Concentration integration time step: "
       << droplet_diffusion.GetMaxTimeStepConcentration() << endl;
  
  cout << "Number of droplets start: " 
       << droplet_diffusion.GetNumberOfDroplets()  << endl;

  cout << "Total concentration start: " 
       << droplet_diffusion.GetTotalConcentration() << endl;

  unsigned int ti = 0;
  droplet_diffusion.SaveDroplets(
      "data/drops_" + to_string(ti) + ".dat");
  droplet_diffusion.SaveConcentration(
      "data/c_" + to_string(ti) + ".dat");

  vector<double> mean_radii;
  vector<unsigned int> current_number_of_droplets;
  while (droplet_diffusion.GetTime() < integration_time) {
    // integrate time
    droplet_diffusion.TimeEvolve(save_every);
    ti++;
    mean_radii.push_back(get_mean( droplet_diffusion.GetRadii()));
    current_number_of_droplets.push_back(droplet_diffusion.GetRadii().size());
    cout << droplet_diffusion.GetTime() << "\t"<< integration_time << endl;
    droplet_diffusion.SaveDroplets("data/drops_"+to_string(ti)+".dat");
  }

  cout << "Total concentration end: "
       << droplet_diffusion.GetTotalConcentration() << endl;
  cout << "Number of droplets end: " 
       << droplet_diffusion.GetNumberOfDroplets()  << endl;

  // save final droplets and concentration
  droplet_diffusion.SaveDroplets("data/drops.dat");
  droplet_diffusion.SaveConcentration("data/c.dat");

  ofstream mean_out("mean_radius.dat");
  for (unsigned int i = 0; i < mean_radii.size(); ++i) {
 
      mean_out << (i + 1) * save_every << "\t"
               << mean_radii[i] << "\t"
               << current_number_of_droplets[i] << "\n";
  }
  mean_out.close();

  return 0;
}
