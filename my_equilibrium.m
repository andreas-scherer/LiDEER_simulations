function rho=my_equilibrium(spin_system)


% Set the assumptions
spin_system=assume(spin_system,'labframe');

% Get the Hamiltonian
[I,Q]=hamiltonian(spin_system,'left');


rho=equilibrium(spin_system,I,Q,[0 0 0]);

end