#ifndef GUARD_CONCENTRATION_H
#define GUARD_CONCENTRATION_H

#include <vector>
#include <string>
#include <fstream>


typedef std::vector<std::vector<double> > VVec;
typedef std::vector<std::vector<std::vector<double> > > VVVec;

class Concentration
{
 public:
  Concentration(double c_init, double D, double Lx, double Ly, double Lz, unsigned int Nx, unsigned int Ny, unsigned int Nz);

  // Source of dn particles in bin centered ad x,y,z
  void Source(double x, double y, double z, double dn);
 
  // evolve time by delta_t
  void TimeEvolve(double delta_t); 

  // save concentration field to file named {name}
  void Save(std::string name) const;

  // return the concentration in the bin centered at x,y,z
  double GetConcentration(double x, double y, double z) const;

  double GetAverageConcentration() const;

  double GetMaxTimeStep() const { return dt_max_;}
 private:

  // diffusion constant
  double D_;
  // box size
  double Lx_, Ly_, Lz_;
  // number of bins in the x and y directions
  double Nx_, Ny_, Nz_;

  // bin dimensions
  double dx_, dy_, dz_;

  
  VVVec concentration1_;
  VVVec concentration2_;

  // pointer to the current concentration  (either concentration1_ or 
  // concentratoin2_)
  VVVec *current_concentration_ptr_;

  // pointer to concentration array that current_concentration_ptr
  // does not point to
  // used as a temporary container
  VVVec *temp_concentration_ptr_;

  // set to 0.2 time the Von Nuemann stability
  double dt_max_;

  // make a single time step of length dt
  void MakeTimeStep(double dt);
};

Concentration::Concentration(double c_init, double D, double Lx, double Ly,double Lz, unsigned int Nx, unsigned int Ny, unsigned int Nz)
  : D_(D), Lx_(Lx), Ly_(Ly), Lz_(Lz), Nx_(Nx), Ny_(Ny), Nz_(Nz), dx_(Lx/Nx), dy_(Ly/Ny), dz_(Lz/Nz),
    concentration1_(Nx, VVec(Ny, std::vector<double>(Nz,c_init))),
    concentration2_(Nx, VVec(Ny, std::vector<double>(Nz,c_init))),
    current_concentration_ptr_(&concentration1_),
    temp_concentration_ptr_(&concentration2_)
{

  // set time step for the integration to
  // 2/10 times the Von Neumann stability limit
  double dmin;
  if (dx_ < dy_ and dx_ < dz_) {
    dmin = dx_;
  } else if (dy_ < dz_) {
    dmin = dy_; 
  } else {
    dmin = dz_;
  }
  dt_max_ = 0.1 * dmin * dmin / D_;

}


void Concentration::TimeEvolve(double delta_t)
{
  while (delta_t > dt_max_) {
    MakeTimeStep(dt_max_);
    delta_t -= dt_max_;
  }
  MakeTimeStep(delta_t);
}

void Concentration::MakeTimeStep(double dt)
{

  double Dx = dt * D_ / (dx_ * dx_);
  double Dy = dt * D_ / (dy_ * dy_);
  double Dz = dt * D_ / (dz_ * dz_);
 
  for (unsigned int ix = 0; ix < Nx_; ++ix) {
  for (unsigned int iy = 0; iy < Ny_; ++iy) {
  for (unsigned int iz = 0; iz < Nz_; ++iz) {
    (*temp_concentration_ptr_)[ix][iy][iz] = (*current_concentration_ptr_)[ix][iy][iz];

    // x derivative
    if (ix == 0) {
      (*temp_concentration_ptr_)[ix][iy][iz] += Dx * (
          (*current_concentration_ptr_)[ix+1][iy][iz]
          - (*current_concentration_ptr_)[ix][iy][iz]);
    }else if (ix == Nx_-1) {
      (*temp_concentration_ptr_)[ix][iy][iz] += Dx * (
          (*current_concentration_ptr_)[ix-1][iy][iz]
          - (*current_concentration_ptr_)[ix][iy][iz]);
    } else {
      (*temp_concentration_ptr_)[ix][iy][iz] += Dx * (
          (*current_concentration_ptr_)[ix+1][iy][iz]
          + (*current_concentration_ptr_)[ix-1][iy][iz]
          - 2 * (*current_concentration_ptr_)[ix][iy][iz]);
    }

    // y derivative
    if (iy == 0) {
      (*temp_concentration_ptr_)[ix][iy][iz] += Dy * (
          (*current_concentration_ptr_)[ix][iy+1][iz]
          - (*current_concentration_ptr_)[ix][iy][iz]);
    } else if (iy == Ny_ - 1) {
      (*temp_concentration_ptr_)[ix][iy][iz] += Dy * (
          + (*current_concentration_ptr_)[ix][iy-1][iz]
          - (*current_concentration_ptr_)[ix][iy][iz]);
    } else {
      (*temp_concentration_ptr_)[ix][iy][iz] += Dy * (
          (*current_concentration_ptr_)[ix][iy+1][iz]
          + (*current_concentration_ptr_)[ix][iy-1][iz]
          - 2 * (*current_concentration_ptr_)[ix][iy][iz]);
    }
    
    // z derivative
    if (iz == 0) {
      (*temp_concentration_ptr_)[ix][iy][iz] += Dz * (
          (*current_concentration_ptr_)[ix][iy][iz+1]
          - (*current_concentration_ptr_)[ix][iy][iz]);
    } else if (iz == Nz_ - 1) {
      (*temp_concentration_ptr_)[ix][iy][iz] += Dz * (
          + (*current_concentration_ptr_)[ix][iy][iz-1]
          - (*current_concentration_ptr_)[ix][iy][iz]);
    } else {
      (*temp_concentration_ptr_)[ix][iy][iz] += Dz * (
          (*current_concentration_ptr_)[ix][iy][iz+1]
          + (*current_concentration_ptr_)[ix][iy][iz-1]
          - 2 * (*current_concentration_ptr_)[ix][iy][iz]);
    }


  }}} // end loop over xi, yi, zi 

  // swap pointers
  VVVec *temp_ptr;
  temp_ptr = temp_concentration_ptr_;
  temp_concentration_ptr_ = current_concentration_ptr_;
  current_concentration_ptr_ = temp_ptr;
}

void Concentration::Source(double x, double y, double z, double dn)
{
  // CHECK FOR UNPHYSICAL VALUES
  unsigned int xi = x / dx_; 
  unsigned int yi = y / dy_; 
  unsigned int zi = z / dz_; 
  (*current_concentration_ptr_)[xi][yi][zi] += dn/(dx_*dy_*dz_);  
}

double Concentration::GetConcentration(double x, double y, double z) const
{
  // CHECK FOR UNPHYSICAL VALUES
  unsigned int xi = x / dx_; 
  unsigned int yi = y / dy_; 
  unsigned int zi = z / dz_; 
  return (*current_concentration_ptr_)[xi][yi][zi];
}

void Concentration::Save(std::string name) const 
{
  // TO DO
}

double Concentration::GetAverageConcentration() const
{

  double c_average = 0.0;
  for (unsigned int ix = 0; ix < Nx_; ++ix) {
  for (unsigned int iy = 0; iy < Ny_; ++iy) {
  for (unsigned int iz = 0; iz < Nz_; ++iz) {
    c_average += (*current_concentration_ptr_)[ix][iy][iz] * dx_ * dy_ * dz_;
  }}}

  return c_average / (Lx_ * Ly_ * Lz_);
}
#endif
