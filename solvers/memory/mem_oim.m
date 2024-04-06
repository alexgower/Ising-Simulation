function [Ebest,tt,step,conf,state] = mem_oim(vars,falgo,Esol,W,fRBM,T,conf,monitor)

% Monitor = [quiet, record, save] boolean (0/1 values)
% Extract quiet and record boolean variables from monitor
quiet = monitor(1); record = monitor(2);

% fp = if true then plaquette based memory, if false then bond based memory
fp = falgo(1)==2.5; % DOUBLE-CHECK why is this 2.5?
% fd = boolean value = 0/1 = single/double precision for W later
fd = falgo(2);


%% Non-RBM case
if ~fRBM

    % sz = (n,m,k,12) etc for 3d, this sets sz = (n,m,k) for 3d etc
    sz = size(W); sz = sz(1:end-1); 

    % d = dimension of lattice, N = number of spins in lattice
    d = length(sz); N = prod(sz);

    % Extract memory dynamics variables from vars
    alpha = vars(1); beta = vars(2); gamma = vars(3); delta = vars(4); zeta = vars(5);
    xini = vars(6); 
    t0 = floor(2^vars(7)); % t0 = 2^vars(7) rounded down to nearest integer
    dtlist = 2.^[-5 vars(8)]; % dtlist = [2^-5, 2^vars(8)]


%% RBM case - MAYBE later
else

    sz = size(W); n = sz(1); m = sz(2); W = W/max(abs(W(:)));
    d = Inf; N = sum(sz);
    alpha = vars(1); beta = vars(2); gamma = vars(3); delta = vars(4); zeta = vars(5);
    xini = vars(6); t0 = floor(2^vars(7)); dtlist = 2.^[-5 vars(8)];

end


% Convert W to single precision if fd = 0 (i.e. not double precision)
if ~fd
    W = single(W);
end

% dd is an effective dimension parameter used for plaquette based memory - just equals 12 for 3D and 4 for 2D (i.e. nearest neighbours nubmer)
dd = 8*d-12;

% Probably just validation that memory initial value xini and restart time t0 are positive - CHECK
if xini < 0 || t0 <= 0
    Ebest = -Esol; tt = Inf; state = 0;
    return;
end


dt0 = 2^-5; % DELETE?

% Number of restarts = ceil(total simulation time/time per restart)
nr = ceil(T/t0); 
% Total simulation time = number of restarts * time per restart (which has been made an integer by rounding down to nearest integer earlier)
T = t0*nr;



%% INITIALISE SPIN CONFIGURATION

% Create checkboard instance in case we use plaquette based memory
% Shade 0 will create a checkerboard that aligns with the U lattice in which the subproblems are defined
% Shade 1 will create a checkerboard that aligns offset to the U lattice in which the subproblems are defined
% The point of this is rather than having symmetric X,Y (double counting) bond memory matrices where X_ij = X_ji etc
% we instead only keep one of X_ij or X_ji etc (this wont affect soft spin constraint dynamics as are summed over anwyway)
% and then either align it with the U lattice or offset it to the U lattice (i.e. checkerboard) from which the cubic subproblems are defined
% to show that 'knowing' this bond structure in advance makes no difference
check = checkerboard(sz,0);

% If no configuration is provided, initialize v, X, Y with random values
if isempty(conf)
    disp("Initialising random configuration")
    % Non-RBM case
    if ~fRBM

    % Initialise spins at all vectors v in the lattice with random values -1 or 1
    % NOW CHANGE THIS TO BETWEEN -PI AND PI via phi = arccos(sigma) FOR OIM CHECK
    % v = acos(double(-1+2*round(rand(sz,'single')))); % SAME AS PEI CASE BUT WITH PHASES
    % v = acos(double(-1+2*round(rand(sz,'single')))) +0.1*(-1+2*rand(sz,'single')); % + SOME NOISE
    v = acos(-1+2*rand(sz,'single')); % CHECK MAYBE A GRADIENT AT 0/PI ERROR THING


        % Bond based memory case
        if ~fp

            % Set all x values to xini, and all y values to 1
            % Note X and Y are (n,m,k,6) etc for 3d i.e. 6 bond values for each spin in a 3d cubic lattice so DOUBLE COUNTING
            X = xini*ones([sz 2*d]);
            Y = ones([sz 2*d]);


        % Plaquette based memory
        else
            
            % i.e. set all the 6 nearest neighbour x values to 0 if they are in the 'white' part of the checkerboard and xini if not
            X = xini*ones([sz dd]).*check;
            % i.e. set all the 6 nearest neighbour y values to 1 if they are in the 'white' part of the checkerboard and 1 if not
            Y = check;

            % Make sure zeta parameter only subtracts from the black part of the checkerboard
            zeta = zeta*check;
        end

    % RBM case - MAYBE - understand what this does later
    else

        % v = double(-1 + 2*round(rand([1 n+m],'single')));
        % X = xini*ones([n m],'single');
        % X = double(X.*logical(W)); 
        % Y = double(logical(W));

    end


% If a configuration IS provided, assign v, X, Y to the values in the provided configuration
else
    v = conf{1}; X = conf{2}; Y = conf{3};
end

% Display initial conf
if ~quiet
    fprintf('Initial conf: \n');
    disp(v);
end

% If not double precision, convert v, X, Y to single precision
if ~fd
    v = single(v); X = single(X); Y = single(Y);
end








%% MONITORING SECTION 
Ebest = 0; flag = 0;

%% INITIAL RECORDING SECTION
if record
    % Create a state struct to store record
    state = struct;

    %% Energy recording

    ttlist = unique(round(geoseries(1,(T),10*round(log2(T))))); % Create list of total times at which to record
    rrecs = length(ttlist); % Get total number of recordings to be made (rr as over all restarts)

    state.Et = zeros([1 rrecs]); % Create empty array to store energy (density) difference values at each recording time
    state.Eb = zeros([1 rrecs]); % Create empty array to store best energy (density) values at each recording time


    %% Bond cluster recording

    rlist = unique(round(geoseries(1,nr,10*round(log2(nr))))); % Create list of restarts at which to record bond cluster values
    nrl = length(rlist); % Get total number of restarts in which to record
    state.tlist = rlist*t0; % Get list of restart start times for restarts in which recordings will take place

    % Create empty array to store bond cluster values at each recording time
    % N generates a dimension from 1:N of all possible bond cluster sizes bg
    % nrl generates a dimension from 1:nrl for each restart
    state.bclus = zeros([1 N nrl]); 

    %% Alex Extra Recording
    state.best_conf = zeros(sz);
    state.confs = zeros([sz rrecs]);
    state.confs_tlist = zeros([1 rrecs]);
else
    state = 0; % If no recording, set state to 0
end


%%% MAIN SIMULATION SECTION
% tt = total time over all restarts, step = total number of steps over all restarts
tt = 0; step = 0; 

% Loop over restarts
for r = 1:nr
    
    %% MEMORY FLIP (SW)
    % On a restart, do a SW cluster update based on clusters formed from rounded X (bond memory) values
    % Only for non-RBM case
    if ~fRBM
        % Only if at least one restart has already been done
        if r > 1

            % Bond-based memory case
            if ~fp
                % Get list of spins to flip from SW algorithm (hence 1 in get_bclus call)
                % i.e. finds clusters based on active bonds (via rounding of X to 0 or 1) and flips them with probability 0.5
                [list,~,~] = get_bclus(get_ww(X),1); 

            % Plaquette-based memory case 
            else
                [list,~,~] = get_bclus(wp_to_w(X,check,0),1); 
            end

            % Flip the spins in the list
            % NOW WE WANT TO FLIP V NOT PHI SO WE WANT COS(PHI) --> -COS(PHI) FOR OIM
            if list
                v(list) = acos(-cos(v(list)));
                % v(list) = -v(list);
            end

        end
    end

    % t = time over current restart
    t = 0; 

    % Only allow time within restart time t0
    while t < t0

        %%% GRADIENT
        % Get rounded Ising energy E,
        % C matrix value (affects memory), 
        % G matrix value (affects continuised spins) 
        % from local gradient descent on Hamiltonian
        [E,C,G] = get_L_oim(v,X,W,fp,check,alpha,beta,fRBM);

        % Update Ebest 
        % and flag (i.e. decide whether to break out of the loop for this restart after this time step or not)
        E_best_previous = Ebest;
        [Ebest,flag] = breakout(E,Ebest,Esol,record);

        % If the energy has changed, update the best configuration
        if E_best_previous ~= Ebest
            state.best_conf = v;
        end

        % Adjust stepsize dt inversely proportional to the maximum absolute value of G (i.e. how fast continuised spin values v are changing)
        % With cutoffs at 2^-5 and 2^vars(8)
        % TODO SHOULD MODIFY THIS BIT LATER
        dt = bound( 1/max(abs(G(:))) , dtlist);

        %%% UPDATE

        % Update continuised spins v but keep them within the range -1 to 1
        % BUT MAYBE THIS DOESN'T MATTER FOR PHASES AS WILL USE COS LATER ANYWAY FOR OIM
        v = v+dt*G;
        % v = min(max( v+dt*G ,-1),1);

        % Non-RBM case
        if ~fRBM
            % Bond-based memory case
            if ~fp 
                % Update memory values X but keep them within the range 0.01 to 1
                X = min(max( X+dt*(gamma*C-Y) ,0.01),1);
                % Update long term memory values Y but keep them within the range 1 to 10
                Y = min(max( Y+dt*(delta*X-zeta), 1),10);

            % Plaquette-based memory case
            else
                % The 'checks' just force that the X and Y values are only non-zero in the black part of the checkerboard
                X = min(max( X+dt*(gamma*C-Y) ,0.01*check),check);
                Y = min(max( Y+dt*(delta*sum(X,d+1)-zeta), check),dd*check);
            end

        % RBM case - has same memory and long term memory updates as non-RBM case but no plaquette based memory option
        else
            X = min(max( X+dt*(gamma*C-Y) ,0.01),1);
            Y = min(max( Y+dt*(delta*X-zeta), 1),10);
        end

        %%% RECORD
        if record
            % If time step is changing into a new integer value and this value is in the totaltime record list ttlist then record
            if (floor(t+dt)>floor(t)) && ismember(floor(tt+dt),ttlist) 
                rrec = find(ttlist == floor(tt+dt),1,'first'); % Record index value (first index in ttlist equal to floor(tt+dt))
                state.Et(rrec) = (Esol-E)/N; % Record energy (density) difference at this time
                state.Eb(rrec) = -Ebest/N; % Record best energy (density) at this time

                state.confs_tlist(rrec) = tt; % Record time at this time
                state.confs(:,:,:,rrec) = v; % Record configuration at this time 
            end

            % If the restart number is in the record list rlist of restart values in which to record cluster values
            if ismember(r,rlist)
                rec = find(rlist == r,1,'first'); % 'Other' record index value (first index in rlist equal to r)

                % tlist = [tlist t]; vt = [vt v(1)]; xt = [xt X(1)]; - DELETE?
                
                if floor(t+dt) > floor(t) % If time step is changing into a new integer value

                    if ~fp % Bond-based memory case
                    
                        % Note that get_ww(X) converts a (n,m,k,6) array of bond memory values to a (n,m,k,3) array of bond memory values 
                        % (for 3D etc) (so get_bclus() doesn't double count bonds)
                        % get_bclus() returns bond cluster information (for active bonds formed from +random then rounded bond memory X values)
                        % bg = list of bond cluster sizes, b = corresponding list of number of bond clusters of this size
                        [~,b,bg] = get_bclus(get_ww(X),0); 
                    
                    else % Plaquette-based memory case
                    
                        % Again check just forces that the re-organised X matrix is only non-zero in the black part of the checkerboard
                        [~,b,bg] = get_bclus(wp_to_w(X,check,0),0);   
                    
                    end

                    % So long as there are non-trivial bond clusters
                    if sum(b)
                        % Add b (number of bond clusters of each size) to the bond cluster size values for this restart (rec)
                        % Note that we add another sampled b value whenever we enter a new integer timestep - CHECK why?
                        state.bclus(1,bg,rec) = state.bclus(1,bg,rec) + b;
                    end

                end

                % if t >= t0
                % state.vspec(:,:,rec) = spec(vt,tlist,nt,dt0);
                % state.xspec(:,:,rec) = spec(xt,tlist,nt,dt0);
                % end

            end
        end

        %% TIME AND STEP UPDATE
        % Increase time within restart by dt, set total time over all restarts to tt, increase step count by one
        t = t+dt; tt = (r-1)*t0+t; step = step+1;

        %% PRINT
        if ~quiet mydot(t,t0,dt,1); end
        if flag break; end

    end

    %% More printing
    if ~quiet
        fprintf(strcat( '\n','t = ',num2str(floor(tt)),',  dE = ',num2str(Esol-Ebest),'\n' ));
    end

    if flag break; end

end

% Set current conf to v, X, Y
conf = {v,X,Y};

end