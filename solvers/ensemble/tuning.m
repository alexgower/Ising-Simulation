function var = tuning(vars0,varlb,varub,indices,falgo,sz,flist,fRBM,steps,runs,T,conf,monitor,fSA,fAPX,gap,oim)

ll = length(vars0); % = number of variables in vars0
varc = vars0(setdiff(1:ll,indices));  % values of variables to keep constant
varv = vars0(indices); % initial values of variables to tune

[falgo,~,flist,fname] = get_suffix(falgo,sz,flist,fRBM);
fname = strcat('tune',fname,'.mat');
quiet = monitor(1); fsave = monitor(3);

if fAPX
    fun = @(var) approx(myinsert(varc,var,indices),falgo,sz,flist,fRBM,runs,T,conf,gap);
elseif gap >= 0
    fun = @(var) mean(multiple(myinsert(varc,var,indices),falgo,sz,flist,fRBM,runs,T,conf,oim));
else
    fun = @(var) log(mean(multiple(myinsert(varc,var,indices),falgo,sz,flist,fRBM,runs,T,conf,oim))/prod(sz) - gap);
end

if fSA
    options = optimoptions(@simulannealbnd,'display','iter','maxiter',steps);
    var = simulannealbnd(fun, varv, varlb, varub, options);
else
    options = optimset('display','iter','maxiter',steps);
    var = fminsearch(fun, varv, options);
end

if ~quiet
    fprintf(num2str(var)); fprintf('\n');
end
if fsave
    save(fname,'var'); 
end

end