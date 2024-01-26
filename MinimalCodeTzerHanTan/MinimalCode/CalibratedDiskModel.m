function [] = CalibratedDiskModel
% Solve ODE system for disk model with effective hydrodynamic interactions

close all
rng('shuffle');

% Video of swimmer dynamics?
mk_vid = 1;

% Plot solution;
plot_sol = 1;

% Path to save all data
sv_file = 'Example/';

% Simulation ID
simID = 'sym';

% Simulation time for ODE model
simtimes = linspace(0,1000,100);

N = 5;

% Size of periodic box (in units of embryo size)
L = 100;

% Periodic domain?
per_dom = 0;

% Compressing surface potential?
Surf_Pot = 1;
WellLength = 30; % Length scale of the well potential

% Stokeslet strength (must be >0, sets the attraction strength between disks)
% Fg = 219 + 20*randn([N,1]);  % Units: [Embryo-Radius^2]/second
Fg = 219 + 70*randn([N,1]);  % Units: [Embryo-Radius^2]/second
% Fg = [150; 250];

% Maximum interaction distance for attractive Stokeslet interaction
RFg_int = 3.8; % 2*sqrt(2) is the second next nearest neighbour in hexagonal grid

% Strength of rotational near-field interactions of neighbouring particles
% Free spinning calibration
f0 = -0.06;
tau0 = 0.12;

% Minimal distance of disk boundaries from which near-field interactions start
Rnf_int = 0.5;

% Single disk angular frequency (= time-scale)
omega0 = 0.05*2 * pi * ( 0.72 + 0.17 * randn([N,1]) );

% Flow interactions between cells 
% = 0: each cell will only interact with its image
% = 1: each cell interacts with all other cells and with its image
Flowinteract = 1;

% Lateral steric repulsion
Sterinteract = 1;

% Spin-Frequency near-field interactions to slow down nearby spinners?
NFinteract = 1;

% Far-field attraction from embryos with up to two neighbours 
% (otherwirse only near-field interactions)
SelectFarField = 1;

% Modulation of intrinsic torques through presence of nearby embryos
ModOmega0 = 1;
N0_damping = 80;

%%%%%% GRAVITY (=Stokeslet direction) (DON'T CHANGE IN DISK MODEL) %%%%%%%
grav_vec = [0,0,-1];
grav_vec = grav_vec/vecnorm(grav_vec);

% Distance of flow singularity below the surface
h = 1; 

% Strength of steric repulsion
Frep_str = 3200 * ( 1 + 0.4 * (rand([N,1]) - 0.5 ) );  % For 1/r^12 repulsive potential (CORRECT ONE)

%%%%%%%%%%%%%%%%%%%%%%% SET INITIAL POSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% (1) RANDOM INITIAL CONDITIONS (SIMPLE) %%%%%%%
% % % First guess
% Pos_init = 2 * L * ( rand([N,2]) - 0.5 );
% 
% % % Eliminate steric overlaps (Brute force, takes longer at high density)
% test_cl = 0; 
% while test_cl == 0
%     dx = Pos_init(:,1) - Pos_init(:,1)';
%     dy = Pos_init(:,2) - Pos_init(:,2)';    
%     distall = sqrt(dx.^2+dy.^2);
%     distall = distall + 100*eye(N);
%     [row, ~] = ind2sub(size(dx),find(distall<2.1));
%     row = unique(row);
%     if isempty(row)
%         % All particles sufficiently far away
%         test_cl = 1; 
%     else
%         % Remove all too close particles and add new ones
%         Pos_init(row,:) = [];
%         Pos_init = [Pos_init; 2 * L * ( rand([length(row),2]) - 0.5 )];
%     end       
% end
% 
% Pos_init = [Pos_init(:,1)'; Pos_init(:,2)'];
% Pos_init = Pos_init(:);

%%%%%%%% (2) RANDOM INITIAL CONDITIONS (PAPER) %%%%%%%
%%%%%% (Starts with particles far outside) %%%%%%%
% phi_part = 2*pi*rand([1,N]); % Random angles
% R_part = 0.8 * ( 350 + 200*sqrt(rand([1,N])).^(1/6) ); % Random radii
% Pos_init = [R_part.*cos(phi_part); R_part.*sin(phi_part)];
% 
% % Apply random stretch factors
% for i = 1:1
%     Pos_init = (1 + 0.4*randn([1,N])).*Pos_init;
% end
% 
% % Move to center of domain
% Pos_init = (Pos_init - mean(Pos_init,2));
% 
% % Arrange for ODE solver
% Pos_init = Pos_init(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% (3) INITIALIZE WITH DATA FROM A PREVIOUS RUN %%%%%
% % Load previous simulation data
% PrevParams = load('Example/Parameters.mat');
% PrevData = load('Example/Simdata.mat');
% 
% % Position of previous data set at last time point
% Pos_prev = PrevData.Cent_coords(end-PrevParams.N+1:end,:)';
% 
% % Position of previous data set at first time point
% % Pos_prev = PrevData.Cent_coords(1:PrevParams.N,:)'; 
% 
% % Arrange for ODE solver
% Pos_init = Pos_prev(:);
% 
% % Parameters in case of microscopic variability (pick whichever from PrevParams)
% Fg = PrevParams.Fg;
% omega0 = PrevParams.omega0;
% Frep_str = PrevParams.Frep_str;
% WellLength = PrevParams.WellLength;
% 
% % Clear memory
% clear('PrevParams');
% clear('PrevData');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% (4) A FEW SMALL GROUPS INITIAL POSITIONS FOR TESTING %%%%%%
% Pos_init = [-1.2, 0; 1.2, 0]'; % N = 2
% Pos_init = [-1.4, 0; 1.4, 0; 0, 2]'; % N = 3
% Pos_init = [-1.4, 0.0; 1.4, -0; 0, 2; 0, -2]'; % N = 4
Pos_init = [-1.4, 0.0; 1.4, -0; 0, 2; 0, -2; 2.3, -2]'; % N = 5

% Arrange for ODE solver
Pos_init = Pos_init(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% END SET INITIAL POSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Singularity implementation  %%%%%%%%%%%%%%%%%%%%%%%
% Functions require input of position vectors [rx1, ry1, rz1; rx2, ry2, rz2; ...]
% and similarly for the parameterizing Stokeslet unit vector f

%%%%%%%%%%%%%%%%%%% STOKESLET flow %%%%%%%%%%%%%%%%%%%%
u_st = @(f,r) (1/(8*pi)) * ( f./repmat(vecnorm(r,2,2),[1,3]) ...
            + r.*repmat(sum(f.*r,2),[1,3])./repmat(vecnorm(r,2,2).^3,[1,3]) );                 

%%%%%%%%%%%%%%%%% Steric repulsion of nearby embryos %%%%%%%%%%%%%%%%%%%%%
% 1/r^12 potential repulsion between midpoints -> 1/r^13 force 
% (written \vec{r}/r^14 because of vec(r) in nominator!!)
Frep = @(r) 12*r./repmat(vecnorm(r,2,2).^14,[1,3]);

%%%%%%%%%%%%%%% Radial dependence of near-field forces %%%%%%%%%%%%%%%%%%%
feta = @(r) log(Rnf_int./abs(r - 2));

%%%%%%%%%%%%%% Radial dependence of near-field torques %%%%%%%%%%%%%%%%%%%
taueta = @(r) log(Rnf_int./abs(r - 2));

% Transformations of image vector and pseudo-vector orientations
vec_img = @(e) [e(:,1), e(:,2), -e(:,3)]; % Stokeslet, force- and source-dipole

function dydt = EmbryoDynamics(t,y)
    % ODE function of embryo dynamics for axisymmetric embryos: 
    % Input vector y contains for each particle 2D position and 3D orientation
    % [x1; y1; x2; y2; ... xN; yN]
    
    % Extract and format positions as needed for flow functions
    Pos3D = [reshape(y(1:2*N),[2,N])', -h*ones([N,1])];
    
    if Surf_Pot == 1 
        % Cylindircal coordinates of the particle positions
        R_pos = vecnorm(Pos3D(:,1:2),2,2);
        Phi_pos = pi + atan2(-Pos3D(:,2),-Pos3D(:,1));
    end
       
    %%% Signed distance matrices r_i - r^0_i where flows from 
    %%% singularities placed at r^0_i are evaluated at r_i    
    dist_x = Pos3D(:,1)' - Pos3D(:,1);
    dist_y = Pos3D(:,2)' - Pos3D(:,2);
    
    % Stokeslet force orientation (DON'T CHANGE FOR DISK MODEL)
    grav = repmat(grav_vec,[N,1]);
    grav = grav./vecnorm(grav,2,2);
    
    % Fixed global orientation of gravity
    fst = grav;   
    fst_img = vec_img(fst);            
    
    %%%%%%%%%%%%% Determine angular frequency of each particle %%%%%%%%%%%%%
    % Compute neighbourhood matrix for near-field interactions
    rij = sqrt(dist_x.^2 + dist_y.^2); % Distance matrix
    NH_matrix = sparse( rij < (2 + Rnf_int) ); % Adjacency matrix of particles within near-field interaction distance
    
    % Determine connected near-field neighbourhood components
    % p: Linear permuted list of all particles
    % r: Array of indices in p from which on a new connected component starts
    [p,~,r,~] = dmperm(NH_matrix); % Has to include diagonal zeros to work
    
    % Remove diagonals for further steps
    NH_matrix = round(NH_matrix - speye(N));
    
    % Build and solve for each connected component the linear system 
    % that determines angular frequencies        
    omega_all = omega0';
    
    if SelectFarField == 1
        % Assume all particles participate in far-field interactions
        idx_FF = true([N,1]);        
    end
    
    %%%%%%%%%%%%%%%%% ANGULAR FREQUENCY CALCULATION %%%%%%%%%%%%%%%%%%%%%
    if NFinteract == 1
        for j = 1:length(r) % Loop through connected components        
            idx_num = [];
            if j < length(r)
                if r(j+1) - r(j) > 1 % If at least two elements in connected component
                    % Extract numeric indices of all disks in this component
                    idx_num = sort( p(r(j):(r(j+1)-1)) );
                end
            else
                % Seperate criterion for at least two elements in cc for last entry
                if N - r(end) > 0
                    % Extract numeric indices of all disks in this component
                    idx_num = sort( p(r(end):end) );
                end
            end

            if ~isempty(idx_num)
                % Number of disks in current connected component
                nrd_cc = size(idx_num,2);

                % Sorted index vector needed to fill linear system matrix
                idx_lin = 1:nrd_cc;

                % Linear matrix of the torque balance for given connected component
                M = zeros(nrd_cc);                                
                
                % Loop through those disks to build linear system
                for l = 1:nrd_cc
                    % Current disk
                    curr_disk = idx_num(l);

                    % Near-field interaction neighbours for current disk
                    nh_vec = ( NH_matrix(curr_disk,:) > 0 );

                    % Pair-wise distances of disks within interaction distance
                    rij_curr = rij(curr_disk,nh_vec);                    
                    
                    % Torque interactions strengths for those distances
                    tau = tau0*taueta(rij_curr);

                    % Fill linear system matrix row for given particle
                    M(l,l) = 1 + sum(tau);                        
                    M(l,idx_lin(ismember(idx_num,find(nh_vec)))) = tau;
                end
                
                if ModOmega0 == 1
                    % Renormalize intrinsic rotation frequencies                    
                    omega0_M = omega0(idx_num)/(1 + (nrd_cc/N0_damping)^2);                    
                end

                %%%%%%%% IS THE OMEGA FUNCTION FOR SUBSEQUENT PLOTTING CORRECT??? %%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Solve for the angular frequencies in this connected component                
                omega = M\omega0_M;

                % Add into array of all angular frequencies according to disk IDs
                omega_all(idx_num) = omega;
                
                if SelectFarField == 1
                    % If these particle belong to group with more than 3 particles
                    % remove those indices from the far-field interaction list
                    if size(idx_num,2) > 3
                        idx_FF(idx_num) = 0;
                    end
                end
            end
        end       
    end
    %%%%%%%%%%%%%%%%% END ANGULAR FREQUENCY CALCULATION %%%%%%%%%%%%%%%%%%%%%
    
    if Surf_Pot == 1 
        % Emulate centering effect of well curvature
        upot_x = - WellLength.^(-2) * (R_pos) .* cos(Phi_pos);
        upot_y = - WellLength.^(-2) * (R_pos) .* sin(Phi_pos);        
    end    
    
    %%%% Sum up all flow contributions that affect a given disk %%%%
    % Cross-product vectors for near-field force interactions 
    RotForce_dir_x = dist_y./rij; % Near-field rotation force direction x
    RotForce_dir_y = -dist_x./rij; % Near-field rotation Force direction y
    
    u_e = zeros([N,2]); % Initiate translational velocities array              

    for j = 1:N % Compute velocity for each particle
        if Flowinteract == 1 % Full flow interactions
            % Logical index for particles whose lateral interaction is taken into account
            idx_lat = true(N,1);        
            idx_lat(j) = 0;
                            
            if (SelectFarField == 1) && (idx_FF(j) == 1)
                % If particle is in far-field group (at most one neighbour)
                % it will interact with all particles
                % idx_far = idx_lat;
                % idx_lat = idx_far;
            else
                % Find all particles further than a distance away 
                % and exclude them from Stokeslet interactions
                idx_far = (rij(:,j)>RFg_int);
                idx_lat(idx_far) = 0;
            end
            
            % Logical index for particles whose image interaction is taken into account            
            idx_img = idx_lat; % Only images that participate in lateral interactions
        else % Each cell only interacts with its image
            idx_lat = false(N,1); % No lateral interactions
            idx_img = false(N,1);
            idx_img(j) = 1; % Cell interacts only with its image
        end
        
        if Sterinteract == 1
            idx_ster = idx_lat; % Same neighbourhood as Stokeslet interaction
        end               
        
        %%%% Because all particles are at the same distance below %%%%%%%
        %%%% the surface only image flow interaction play a role %%%%%%%%
        r_curr = [dist_x(idx_img,j), dist_y(idx_img,j), -2*h*ones(sum(idx_img),1)];
        
        %%%% Collect distances for steric interactions %%%%%
        if Sterinteract == 1 
            % Same neighbourhood as Stokeslet interaction
            r_curr_ster = r_curr;
            r_curr_ster(:,3) = 0;
        end        
        
        % Prepare according arrays of vectors parameterizing flow singularities        
        fst_curr = fst_img(idx_img,:);
        
        % Collect all attractive Stokeslet flow interactions (no additional weighting)
        u_star =  0.5 * sum((Fg(j) + Fg(idx_img)) .* u_st(fst_curr,r_curr),1);  % SYMMETRIZED                    
        % u_star =  sum(Fg(idx_img) .* u_st(fst_curr,r_curr),1);                 % UNSYMMETRIZED        

        u_e(j,:) = u_star(1:2); % Only vx and vy are relevant

        % Steric repulsion only laterally between embryos
        if Sterinteract == 1
            u_rep = 0.5 * sum((Frep_str(j) + Frep_str(idx_ster)) .* Frep(r_curr_ster),1);
            u_e(j,:) = u_e(j,:) + u_rep(1:2);            
        end  
        
        %%%%%%% Contributions from transverse force interactions %%%%%%%
        idx_neighb = (NH_matrix(j,:) > 0); % Neighbour-indices (no diagonals!)
        
        if sum(idx_neighb,2) ~= 0
            nr_neighb = full(sum(idx_neighb)); % Number of neighbours
            feta_curr = feta(rij(j,idx_neighb));
            
            % Build omega array that can be used for rotation force summation
            omega_full = repmat(omega_all(j),[1,nr_neighb]) + omega_all(idx_neighb);
            
            % Sum up all transverse force contributions for given particles
            u_e(j,1) = u_e(j,1) + f0*sum(feta_curr.*omega_full.*RotForce_dir_x(j,idx_neighb),2);
            u_e(j,2) = u_e(j,2) + f0*sum(feta_curr.*omega_full.*RotForce_dir_y(j,idx_neighb),2);                       
        end
    end

    % Fill RHS output vector of ODE system
    ue_transp = u_e';
    dydt = ue_transp(:);

    % If well curvature effect is included
    if Surf_Pot == 1
        u_pot = [upot_x'; upot_y'];
        dydt = dydt + u_pot(:);
    end    
end %dydt

%%%%%%%%%%%%%%%% Solve ODEs %%%%%%%%%%%%%%%%%%%%
y_init = Pos_init(:);
opts = odeset('RelTol',1e-5,'AbsTol',1e-5,'Vectorized','off','OutputFcn',@odetpbar);
[t,y] = ode113(@EmbryoDynamics,simtimes,y_init,opts);
if per_dom == 1
    y = mod(y,L);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save parameters to workspace to folder
save([sv_file,'Parameters_',simID,'.mat']);

% Save position and radius data
Cent_coords = y(:,1:(2*N))';
Cent_coords = Cent_coords(:);
Cent_coords = reshape(Cent_coords',[2,N*length(simtimes)])';

if plot_sol == 1
    % Extract and plot all tracks
    figure(1);
    scatter(y(:,1),y(:,2),'filled');
    hold on
    for i = 2:N
        scatter(y(:,2*i-1),y(:,2*i),'filled');    
    end
    % Red squares to label starting points of trajectories
    scatter(y(1,1:2:(2*N-1)),y(1,2:2:(2*N)),40,'r','square','filled');
    axis equal
    if per_dom == 1
        axis([0,L,0,L]);
    else
        axis([-L,L,-L,L]);
    end

    % Plot the time courses of all orientation z-components
    f2 = figure(2);
    f2.InvertHardcopy = 'off'; % Ensures figure is saved with given background color
    set(gcf,'color','w'); % WHITE BACKGROUND FOR FIGURES

    drawnow;

    if mk_vid == 1
        v = VideoWriter(strcat(sv_file,simID),'MPEG-4');
        v.Quality = 50;
        v.FrameRate = 15;
        open(v)        
    end            
end

Omega_all = [];

for k = 1:length(t)   
    Pos_x = y(k,1:2:(2*N-1));
    Pos_y = y(k,2:2:(2*N));

    % Calculate angular spinning frequency from current position
    omega_curr = ComputeOmega(Pos_x,Pos_y);
    Omega_all = [Omega_all, omega_curr];
    
    if plot_sol == 1
        if k == 1 % First loop assign axis                
            set(0,'CurrentFigure',f2);
            sc = scatter(y(k,1:2:(2*N-1)),y(k,2:2:(2*N)),3e3/L,omega_curr'/(2*pi),'filled',...
                    'MarkerEdgeColor','k');
            hold on
            set(gca,'FontSize',15);

            axis equal
            if per_dom == 1
                axis([0,L,0,L]);
            else
                axis([-L,L,-L,L]);                
            end
            box on
            
            caxis([0,0.7]);            
            colorbar
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.55, 0.5, 0.7]);
        else % Following loops just update axis content
            sc.XData = y(k,1:2:(2*N-1));
            sc.YData = y(k,2:2:(2*N));
            sc.CData = omega_curr'/(2*pi);            
            caxis([0,0.7]);                    
        end

        drawnow;

        if (mk_vid == 1)&&(mod(k-1,1)==0)
            set(0,'CurrentFigure',f2);
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end
end

save([sv_file,'Simdata_',simID,'.mat'],'Cent_coords','Omega_all');
    
if mk_vid == 1
    close(v);
end
    
    function outp = ComputeOmega(x,y)
        %%% Signed distance matrices r_i - r^0_i where flows from 
        %%% singularities placed at r^0_i are evaluated at r_i
        dist_x = x' - x;
        dist_y = y' - y;

        if per_dom == 1
            % Correct distances that where computed across the boundary
            dist_x = dist_x + ( (dist_x<-L/2) - (dist_x>L/2) ) * L;
            dist_y = dist_y + ( (dist_y<-L/2) - (dist_y>L/2) ) * L;                   
        end

        %%%%%%%%%%%%% Determine angular frequency of each particle %%%%%%%%%%%%%
        % Compute neighbourhood matrix for near-field interactions
        rij = sqrt(dist_x.^2 + dist_y.^2); % Distance matrix
        NH_matrix = sparse( rij < 2 + Rnf_int ); % Adjacency matrix of particles within near-field interaction distance

        % Determine connected components
        % p: Linear permuted list of all particles
        % r: Array of indices in p from which on a new connected component starts
        [p,~,r,~] = dmperm(NH_matrix); % Has to include diagonal zeros to work

        % Remove diagonals for further steps
        NH_matrix = round(NH_matrix - speye(N));

        % Build and solve for each connected component the linear system that
        % determines angular frequencies            
        outp = omega0';

        if NFinteract == 1
            for j = 1:length(r) % Loop through connected components        
                idx_num = [];
                if j < length(r)
                    if r(j+1) - r(j) > 1 % If at least to elements in connected component
                        % Extract numeric indices of all disks in this component
                        idx_num = sort( p(r(j):(r(j+1)-1)) );
                    end
                else
                    % Seperate criterion for at least two elements in cc for last entry
                    if N - r(end) > 0
                        % Extract numeric indices of all disks in this component
                        idx_num = sort( p(r(end):end) );
                    end
                end

                if ~isempty(idx_num)
                    % Number of disks in current connected component
                    nrd_cc = size(idx_num,2);

                    % Sorted index vector needed to fill linear system matrix
                    idx_lin = 1:nrd_cc;

                    % Linear matrix of the torque balance for given connected component
                    M = zeros(nrd_cc);                    

                    % Loop through those disks to build linear system
                    for l = 1:nrd_cc
                        % Current disk
                        curr_disk = idx_num(l);

                        % Near-field interaction neighbours for current disk
                        nh_vec = ( NH_matrix(curr_disk,:) > 0 );

                        % Pair-wise distances of disks within interaction distance
                        rij_curr = rij(curr_disk,nh_vec);

                        % Torque interactions strengths for those distances
                        tau = tau0*taueta(rij_curr);

                        % Fill linear system matrix row for given particle   
                        M(l,l) = 1 + sum(tau);
                        M(l,idx_lin(ismember(idx_num,find(nh_vec)))) = tau;
                    end

                    if ModOmega0 == 1
                        % Renormalize intrinsic rotation frequencies                        
                        omega0_M = omega0(idx_num)/(1+(nrd_cc/N0_damping)^2);
                    end

                    % Solve for the angular frequencies in this connected component                
                    omega = M\omega0_M;

                    % Add into array of all angular frequencies according to disk IDs
                    outp(idx_num) = omega;
                end
            end       
        end
    end % END function to compute OMEGA

end %END main

