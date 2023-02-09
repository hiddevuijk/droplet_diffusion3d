#ifndef GUARD_RADIUS_FREQUENCY_H
#define GUARD_RADIUS_FREQUENCY_H

#include <vector>

class RadiusFrequency {
 public:
  RadiusFrequency(double dr, unsigned int N)
    : dr_(dr), N_(N),
      NumberOfSamples_(0), NumberOfSamplesOutOfRange_(0),
      radius_frequency_(N,0.0)
    {}

  void Sample(const std::vector<double>& radii);
 private:
  double dr_;
  unsigned int N_;

  long unsigned int NumberOfSamples_;
  long unsigned int NumberOfSamplesOutOfRange_; 

  std::vector<double> radius_frequency_;
};


void RadiusFrequency::Sample(const std::vector<double>& radii)
{
  unsigned int index; 
  for (unsigned int i = 0; i < radii.size(); ++i) {
    index = std::floor(radii[i] / dr_);
    if (index < N_) {
      radius_frequency_[index] += 1;
      NumberOfSamples_ += 1;
    } else {
      NumberOfSamplesOutOfRange_ += 1;
    }
  }

}

#endif
