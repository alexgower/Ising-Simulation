function E = cube_ising_energy(sigma_0, W)
    d = length(size(sigma_0));
    Sigma_0 = get_V(sigma_0);
    E = double(sum(sigma_0.*sum(Sigma_0.*W,d+1),'all')/2);
end