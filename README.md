I have used the code to obtain all the results for my thesis and a future paper. The code creates the Hamiltonian matrix for either a 1D or quasi-1D tight binding model. The user defines the model, and the Hamiltonian is created according to the definition. The matrix is diagonalized to find the energys(eigenvalues) and states(eiegenvectors).

In interacting systems, the code performs iterations, so that in each iteration the Hamiltonian is updated according to the Hartree-Fock equations, using the eigenstates.

The code is written in C and uses the Matlab engine API for all calculations related to Eigenvalue decompositions.

optional functions:

- interactions

- magnetic field perpendicular to the ring's/ladder's plain - possible in various configurations

- next-nearest neighbor hopping

- impurities

- Cosine modulated on site potential - Harper model (and therefore the user has to define amplitude and frequency)

- Calculations of persistent currents and current operators

- Entanglement entropy

Depending on the user's definiiton, the code saves .mat files to a directory. The user may save the following observables:

- Eigenvalues (energy)

- Probability amplitude (wave function amplitude)

- Persistent currents

- Current operator

- Entanglement Entropy

The code can run at one instance through a range of:

- Fillings

- Magnetic field amplitudes

- Impurities
