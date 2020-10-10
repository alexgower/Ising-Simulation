function Ediff = multiple(vars,falgo,sz,flist,fRBM,runs,T)

algo = get_suffix(fRBM,falgo);
[Wlist,Esol] = ensemble(sz,flist,fRBM,runs);
Elist = zeros(1,runs);

parfor run = 1:runs
W = Wlist{run}; 
if strcmp(algo,'SA')
Elist(run) = SA(vars,Esol(run),W,fRBM,T,Inf,falgo,[1 0 0]);
elseif strcmp(algo,'PT')
Elist(run) = PT(vars,Esol(run),W,fRBM,T,Inf,[],[1 0 0]);
elseif strcmp(algo,'ICM')
Elist(run) = PTI(vars,Esol(run),W,fRBM,T,Inf,[],[1 0 0]);
elseif strcmp(algo,'mem')
Elist(run) = mem(vars,Esol(run),W,fRBM,T,Inf,[],[1 0 0]);
end
end
Ediff = Esol-Elist;

end