#include <iostream>
#include <map>
#include <vector>
#include <array>
#include <any>
#include <functional>
#include <memory>


double lj_cutoff = 2.5;
double lj_cutoff2 = lj_cutoff * lj_cutoff;
double lj_potential_at_cutoff;

/*! \brief Evaluate the Lennard-Jones potential associated with a specific particle separation.
 *
 * \param [in]  r2
 *                   Square of the distance between two particles.
 */
double lj_potential(double r2) {
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    return 4.0 * (inv_r6 * inv_r6 - inv_r6);
}

/*! \brief Evaluate the Lennard-Jones potential associated with a specific particle separation, with a cutoff.
 *
 * \param [in]  r2
 *                   Square of the distance between two particles.
 */
double lj_potential_with_cutoff(double r2) {
    if (r2 < lj_cutoff2) {
        return lj_potential(r2) - lj_potential_at_cutoff;
    }
    else {
        return 0.0;
    }
}

/*! \brief Evaluate the Lennard-Jones force for a specific particle separation.
 *
 * \param [in]  r2
 *                   Square of the distance between two particles.
 */
double lj_force(double r2) {
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    return 24.0 * inv_r2 * (2.0 * inv_r6 * inv_r6 - inv_r6);
}

/*! \brief Evaluate the Lennard-Jones force for a specific particle separation, with a cutoff.
 *
 * \param [in]  r2
 *                   Square of the distance between two particles.
 */
double lj_force_with_cutoff(double r2) {
    if (r2 < lj_cutoff2) {
        return lj_force(r2);
    }
    else {
        return 0.0;
    }
}

/*! \brief Function to evaluate all the forces using a Lennard-Jones potential
 *
 * \param [in]  nparticles
 *                   Number of particles in the system
 * \param [out] potential_energy
 *                   Total potential_energy of the system
 * \param [in]  box_size
 *                   Length of the simulation cell
 * \param [in]  positions
 *                   Position of the nuclei
 * \param [out] forces
 *                   Forces on the nuclei
 */
void evaluate_lj_forces(
        const int &nparticles,
        double &potential_energy,
        const double &box_size,
        const std::vector<std::array<double, 3>> &positions,
        std::vector<std::array<double, 3>> &forces) {

  for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
    for (int jparticle = 0; jparticle < nparticles; ++jparticle) {
      if ( iparticle != jparticle ) { // Only compute the interactions between different particles
        double dx = positions[iparticle][0] - positions[jparticle][0];
        if (dx > 0.5 * box_size) dx -= box_size;
        if (dx < -0.5 * box_size) dx += box_size;

        double dy = positions[iparticle][1] - positions[jparticle][1];
        if (dy > 0.5 * box_size) dy -= box_size;
        if (dy < -0.5 * box_size) dy += box_size;

        double dz = positions[iparticle][2] - positions[jparticle][2];
        if (dz > 0.5 * box_size) dz -= box_size;
        if (dz < -0.5 * box_size) dz += box_size;

        double r2 = (dx * dx) + (dy * dy) + (dz * dz);

        double f = lj_force_with_cutoff(r2);

        forces[iparticle][0] += f * dx;
        forces[iparticle][1] += f * dy;
        forces[iparticle][2] += f * dz;

        potential_energy += 0.5 * lj_potential_with_cutoff(r2);
      }
    }
  }

}

/*! \brief Initialization function for the plugin.
 *
 * \param [in]  state
 *                   Map with all state accessible to the plugin
 */
extern "C"
void initialize(
       std::map<std::string, std::shared_ptr<std::any>> &state) {

  // Determine the Lennard-Jones potential at the cutoff
  lj_potential_at_cutoff = lj_potential(lj_cutoff2);

}

/*! \brief Function to execute the plugin.
 *
 * \param [in]  state
 *                   Map with all state accessible to the plugin
 */
extern "C"
void evaluate_forces(
       std::map<std::string, std::shared_ptr<std::any>> &state) {

  // You need to call evaluate_lj_forces, but first you need to extract all its required arguments from "state"
  /*
  evaluate_lj_forces(nparticles, 
                     potential_energy,
                     box_size,
                     positions,
                     forces);
  */
}
