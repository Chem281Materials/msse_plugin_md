#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <dlfcn.h>

class MDSimulation {
  public:
    MDSimulation(double box_size_in, int nparticles_in);
    void run(int nsteps, double dt);
  private:
    double box_size;               // Length of each side of the periodic simulation cell, which is cubic.
    double potential_energy;
    double kinetic_energy;
    int nparticles;                // Number of particles in the simulation
    std::vector<std::array<double, 3>> positions;  // Position of the particles
    std::vector<std::array<double, 3>> velocities; // Velocities of the particles
    std::vector<std::array<double, 3>> forces;     // Forces on the particles
};

/*! \brief Initialize a molecular dynamics simulation
 *
 * \param [in]  box_size_in
 *                   Length of each side of the periodic simulation cell, which is cubic.
 * \param [in]  nparticles_in
 *                   Number of particles in the simulation.
 */
MDSimulation::MDSimulation(double box_size_in, int nparticles_in) {
    box_size = box_size_in;
    nparticles = nparticles_in;

    // Initialize the particles on a rough grid
    int particles_per_side = std::ceil( std::pow(nparticles, 1.0/3.0) );
    double particle_spacing = box_size / (particles_per_side + 1);
    for (int iparticle = 0; iparticle < nparticles; iparticle++) {
        int ix = iparticle % particles_per_side;
        int iy = (iparticle / particles_per_side) % particles_per_side;
        int iz = iparticle / (particles_per_side * particles_per_side);
        positions.push_back({
            particle_spacing * ix + ( 0.5 * particle_spacing ),
            particle_spacing * iy + ( 0.5 * particle_spacing ),
            particle_spacing * iz + ( 0.5 * particle_spacing )
            });
    }

    // Initialize the velocities randomly
    for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
        /* The random number generators here use the particle index as the seed.
           This isn't something you would normally do, but it is quite helpful
           in this case for the purpose of ensuring that the velocities are 
           reproducible with respect to parallelization. */
        std::mt19937 gen(iparticle);
        std::uniform_real_distribution<double> random_vel(-0.5, 0.5);
        velocities.push_back({random_vel(gen), random_vel(gen), random_vel(gen)});
    }

    // Initialize the forces
    for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
        forces.push_back({0.0, 0.0, 0.0});
    }

    // Call the plugin initialize function here
}

/*! \brief Run a molecular dynamics simulation.
 *
 * \param [in]  nsteps
 *                   Number of time integration steps to perform.
 * \param [in]  dt
 *                   Size of the timestep (reduced Lennard-Jones units).
 */
void MDSimulation::run(int nsteps, double dt) {

    // Main simulation loop
    for (int istep = 0; istep < nsteps; ++istep) {

        // Update the particle velocities and positions
        for (int iparticle = 0; iparticle < nparticles; ++iparticle) {

            // Update the positions
            positions[iparticle][0] += velocities[iparticle][0] * dt;
            positions[iparticle][1] += velocities[iparticle][1] * dt;
            positions[iparticle][2] += velocities[iparticle][2] * dt;

            // Apply periodic boundary conditions; ensure that particles outside the box wrap to the other side
            for (int idimension = 0; idimension < 3; ++idimension) {
                if (positions[iparticle][idimension] < 0.0) positions[iparticle][idimension] += box_size;
                if (positions[iparticle][idimension] >= box_size) positions[iparticle][idimension] -= box_size;
            }

        }

        // Zero the energy and forces
        potential_energy = 0.0;
        kinetic_energy = 0.0;

        for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
            forces[iparticle] = {0.0, 0.0, 0.0};
        }
        
        // Call the plugin evaluate_forces function here

        // Compute the kinetic energy
        for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
            kinetic_energy += 0.5 * velocities[iparticle][0] * velocities[iparticle][0];
            kinetic_energy += 0.5 * velocities[iparticle][1] * velocities[iparticle][1];
            kinetic_energy += 0.5 * velocities[iparticle][2] * velocities[iparticle][2];
        }

        // Update the particle velocities
        for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
            velocities[iparticle][0] += forces[iparticle][0] * dt;
            velocities[iparticle][1] += forces[iparticle][1] * dt;
            velocities[iparticle][2] += forces[iparticle][2] * dt;
        }

        // Print output
        std::cout << "Iteration " << istep << std::endl;
        std::cout << "    Potential Energy: " << potential_energy << std::endl;
        std::cout << "    Kinetic Energy:   " << kinetic_energy << std::endl;
        std::cout << "    Total Energy:     " << potential_energy + kinetic_energy << std::endl << std::endl;
    }

    std::cout << "Simulation completed." << std::endl;
}

int main(int argc, char** argv) {
    MDSimulation mysimulation(20.0, 1000);
    mysimulation.run(100, 0.005);
    return 0;
}