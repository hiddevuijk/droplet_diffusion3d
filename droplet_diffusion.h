#ifndef GUARD_DROPLET_DIFFUSION_H
#define GUARD_DROPLET_DIFFUSION_H

#include "concentration.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include <boost/random.hpp>


/*
  The structure "Droplet" contains all info
  (x, y position and radius) about a single droplet.
  The radius can only be changed by SetR and ChangeR,
  be cause the object keeps track of the last change in
  the radius and area.
  A radius of R = -1, means the droplet has dissolved.
*/
struct Droplet {
 public:
  Droplet() {};
  Droplet(double x, double y, double z, double R, double Rco)
    : x(x), y(y), z(z), R(R), Rco(Rco) {} 

  // position
  double x, y, z;

  double GetR() const { return R; }
  double GetRco() const { return Rco;}

  // Set droplet radius, after this the change in the radius
  // (dR) is set to zero.
  void SetR(double newR) {
    R = newR;
    dR = 0;
  }

  void SetRco(double newRco) { Rco = newRco; }

  // change the radius by delta_R
  void ChangeRadius(double delta_R) {
      dR = delta_R;
      R += delta_R;
  }

  // current area of the droplet
  double Volume() const { return 4 * M_PI * R * R * R / 3;}

  // last change in the area of the droplet
  // if R < Rco the droplet dissolves,
  // and the previous Area is returned
  double VolumeChange() const
  {
    if (R > Rco) {
      return 4 * M_PI * ( R * R * R - (R-dR)*(R-dR)*(R-dR) ) / 3;
    } else {
      // droplet has disolved
      return -4 * M_PI * (R - dR) * (R - dR) * (R - dR) / 3;
    }
  }

 private:

  // radius
  double R;

  // cut-off radius
  // droplet dissolves if R < Rco
  double Rco;

  // last change in the radius
  // (needed to calculate the last change in the area)
  double dR;
};



/*
 The DropletDiffusion class is the main simulation class.

 Funct is a class that gives the c0_out concentration
 it is specified by the user, and should have a function call
 operator overloaded that takes three double arguments,
 the x,y,z coordinates and the time,
 and returns the c0_out concentration at that place and time.
 For example,
  Funct F {

    double operator() (double x, double y, double z, double t) {
      return {double};
    }

  };

*/
template<class Funct>
class DropletDiffusion {
 public:
  DropletDiffusion(
    double D,     // diffusion constant of the concentration
    double Dd1,   // diffusion const. of a drop with unit radius
    double gamma, // surface tension
    double c0_in, // equilibiurm concentration inside the droplet
    double c_out, // initial uniform concentration outside
    Funct c0_out, // Functor for the equilibrium outside concentration
    double Rco,   // cut of radius, drop dissolves if radius < Rco
    double Lx,    // x dimension of the box
    double Ly,    // same for y
    double Lz,    // same for z
    double dt,    // integration time step, 
                  // the concentration field is integrated with 
                  // time step min(dt, dt_von_neumann / 10),
                  // where dt_von_neumann is the time step
                  // corresponding to marginal stability
    unsigned int Nx, // number of grid points in the x dimension
                     // for the concentration
    unsigned int Ny,    // same for y
    unsigned int Nz,    // same for z
    unsigned int seed); // seed for the random number generator

  // integrate t -> t + delta_t in steps of dt_
  void TimeEvolve(double delta_t);

  double GetTime() const { return t_; }
  unsigned int GetNumberOfDroplets() const;

  // return the total number of concentration particles
  // This should always return the same amount.
  double GetTotalConcentration() const;

  std::vector<double> GetRadii() const;

  // save current droplet positions and radii to file {name}
  void SaveDroplets(std::string name) const;
  // save current concentration field to file {name}
  void SaveConcentration(std::string name) const;

  // returns max time step used to integrate the concentration.
  double GetMaxTimeStepConcentration() const
    { return concentration_.GetMaxTimeStep(); }

  // Initialize N droplets randomly and uniformly in space
  // RadiusDistribution is a functor
  // such that radius_distribution(double random_uniform_number)
  // samples the distribution of the initial radii of the droplets
  template<class RadiusDistribution> 
  void InitializeDroplets(unsigned int N, RadiusDistribution radius_distribution);

 private:

  ////////////////////////////////////////
  /// Imutable members variables
  ////////////////////////////////////////

  // diffusion constant of the of the concentration field
  double D_;
  // Diffusion constant of droplet with unit radius
  double Dd1_;
  // equilibrium concentration inside the droplet
  double c0_in_;

  // capillary length
  double l_gamma_;
  // function pointer to c0, the equilibrium concentration outside
  Funct c0_out_;
  // cutt-off radius. Droplet dissolves if R < Rco
  double Rco_;
  // box dimensions
  double Lx_, Ly_, Lz_;
  // integration time step
  double dt_;

  ////////////////////////////////////////
  /// Mutable members variables
  ////////////////////////////////////////
 
  // current time 
  double t_;

  // concentration field (see concentration.h)
  Concentration concentration_;
   
  // random number generator for diffusion of droplets
  const boost::normal_distribution<double> normal_distribution_;
  boost::mt19937 random_number_generator_;
  boost::variate_generator<boost::mt19937&,
         boost::normal_distribution<double> > random_normal_distribution_;

  // list with all the droplets in the system
  // length can change as droplets dissolve
  std::vector<Droplet> droplets_;   

  // number of droplets
  unsigned int number_of_droplets_;

  ////////////////////////////////////////
  /// private member functions
  ////////////////////////////////////////

   // make a time step of dt
   void MakeTimeStep( double dt);
};

template<class Funct>
DropletDiffusion<Funct>::DropletDiffusion(
    double D,
    double Dd1,
    double l_gamma,
    double c0_in,
    double c_out,
    Funct c0_out,
    double Rco,
    double Lx,
    double Ly,
    double Lz,
    double dt,
    unsigned int Nx,
    unsigned int Ny,
    unsigned int Nz,
    unsigned int seed)
  : D_(D),
    Dd1_(Dd1),
    c0_in_(c0_in),
    l_gamma_(l_gamma),
    c0_out_(c0_out),
    Rco_(Rco),
    Lx_(Lx),
    Ly_(Ly), 
    Lz_(Lz), 
    dt_(dt),
    t_(0),
    concentration_(c_out, D, Lx, Ly, Lz, Nx, Ny, Nz),
    normal_distribution_(0.0,1.0),
    random_number_generator_(seed),
    random_normal_distribution_(random_number_generator_, normal_distribution_)
{}


template<class Funct>
void DropletDiffusion<Funct>::MakeTimeStep(double dt)
{
  // evolve droplet positions and radii
  double alpha, epsilon;
  for (unsigned int i = 0; i < number_of_droplets_; ++i) {

    double R = droplets_[i].GetR();
    if (R < 0) continue; // droplet has dissolved
    double &x = droplets_[i].x;
    double &y = droplets_[i].y;
    double &z = droplets_[i].z;

    // move droplet
    x += sqrt(2 * dt * Dd1_ / R) * random_normal_distribution_();
    y += sqrt(2 * dt * Dd1_ / R) * random_normal_distribution_();
    z += sqrt(2 * dt * Dd1_ / R) * random_normal_distribution_();

    // reflecting b.c. for the droplet center.
    // reflect droplet back in the box
    if (x > Lx_) {
      x = 2 * Lx_ - x;
    } else if (x < 0) {
      x = -x;
    }

    if (y > Ly_) {
      y = 2 * Ly_ - y;
    } else  if (y < 0) {
      y = -y;
    }

    if (z > Lz_) {
      z = 2 * Lz_ - z;
    } else  if (z < 0) {
      z = -z;
    }


    // evolve droptlet radius
    alpha = D_ * c0_out_(x,y,z,t_) / c0_in_;
    epsilon = -1 + concentration_.GetConcentration(x,y,z) /
              c0_out_(x,y,z,t_);
    double dR = epsilon - l_gamma_ / R;
    dR *= alpha / R;
    droplets_[i].ChangeRadius(dt * dR);
  }

  // change concentration due to change in droplet radii
  // if radius < Rco, remove droplet (set radius = -1)
  for (unsigned int i = 0; i < number_of_droplets_; ++i) {
    if (droplets_[i].GetR() < 0) continue; //droplet has dissolved

    //Transfer the change in the concentration particles inside
    //the droplet due to the area change in the bin corresponding
    //to the center of the droplet
    concentration_.Source(
                 droplets_[i].x, droplets_[i].y, droplets_[i].z,
                 -1 * droplets_[i].VolumeChange() * c0_in_);

    // if the droplet has dissolved, set radius to -1
    if (droplets_[i].GetR() < Rco_) droplets_[i].SetR(-1); 
  }

  // Evolve the concentratoin field in time
  concentration_.TimeEvolve(dt);

  t_ += dt;
}

// Time evolve the system t->t+delta_t, 
// in steps dt_ or smaller
template<class Funct>
void DropletDiffusion<Funct>::TimeEvolve(double delta_t) 
{
  // integrate in steps of dt_
  while (delta_t > dt_) {
    MakeTimeStep(dt_);
    delta_t -= dt_;
  }

  // integrate remaining time
  MakeTimeStep(delta_t);
}


template<class Funct>
void DropletDiffusion<Funct>::SaveDroplets(std::string name) const
{
  std::ofstream out;
  out.open(name);
  out << std::setprecision(16);
  for (unsigned int i = 0; i < number_of_droplets_; ++i) {
    out << droplets_[i].x << "\t";
    out << droplets_[i].y << "\t";
    out << droplets_[i].z << "\t";
    out << droplets_[i].GetR() << "\n";
  }

  out.close();
}

template<class Funct>
void DropletDiffusion<Funct>::SaveConcentration(std::string name) const 
{
  concentration_.Save(name);
}

template<class Funct>
double DropletDiffusion<Funct>::GetTotalConcentration() const
{
  double c_tot = 0.0;

  // Add concentration from droplets 
  for (unsigned int i = 0; i < number_of_droplets_; ++i) {
    if (droplets_[i].GetR() > 0 )
      c_tot += droplets_[i].Volume() * c0_in_;
  }

  // add concentration from outside the droplets
  c_tot += concentration_.GetAverageConcentration() *Lx_*Ly_*Lz_;

  return c_tot;
}

template<class Funct> 
unsigned int DropletDiffusion<Funct>::GetNumberOfDroplets() const
{
  unsigned int Ndroplets = 0;
  for (unsigned int i = 0; i < number_of_droplets_; ++i) {
    if (droplets_[i].GetR() > Rco_ ) Ndroplets++;
  }

  return Ndroplets;
}


template<class Funct> template<class RadiusDistribution> 
void DropletDiffusion<Funct>::InitializeDroplets(
              unsigned int N, RadiusDistribution radius_distribution )
{
  number_of_droplets_ = N;
  droplets_ = std::vector<Droplet>(number_of_droplets_);
  
  boost::uniform_real<double> udist(0,1);
  boost::variate_generator<boost::mt19937&,
      boost::uniform_real<double> >
        rudist(random_number_generator_, udist);

  // initialize droplets with random positions in the box
  for (unsigned int i = 0; i < N; ++i) {
    droplets_[i].x = rudist() * Lx_;
    droplets_[i].y = rudist() * Ly_;
    droplets_[i].z = rudist() * Lz_;

    // initialize droplet radius
    droplets_[i].SetR( radius_distribution( rudist() ) );
    droplets_[i].SetRco(Rco_);
  }

}

template<class Funct>
std::vector<double> DropletDiffusion<Funct>::GetRadii() const
{
  std::vector<double> radii;
  for (unsigned int i = 0; i < droplets_.size(); i++) {
    if (droplets_[i].GetR() > Rco_) {
      radii.push_back( droplets_[i].GetR());
    }

  }
  
 return radii; 
    
}
#endif
