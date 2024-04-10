function E = general_ising_energy(sigmas,J)
    E = 0.5 * sigmas * J * sigmas';
end
