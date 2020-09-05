%% constants
kB = 1.3806504e-23; %Boltzmann constant
R = 8.314472;
mHe = 6.6465e-27; %Helium mass
dref = 2.33e-10; %Helium reference diameter
Tref = 273; %reference temperature for above diameter
omega = 0.66; %viscosity-temperature power law exponent
VSS_param = 1; %1 for VHS model, VSS model not implemented
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


%% inputs to be written by the user
T0 = ; %stagnation temperature
v_input = sqrt(5*kB*T0/mHe); %input stream velocity in x-direction, default set to terminal velocity
Tskim = ; %skimmer surface temperature
n_input = ; %input number density
max_iter = ;
geometryFile = ; %stl file for skimmer geometry mesh, save in units of meters

%simulation frame size, x is flow direction
xmin = ;
xmax = ;
ymin = ; 
ymax = ;
zmin = ;
zmax = ;

%cell setup
%number of cells in each direction
numCellX = ;
numCellY = ;
numCellZ = ;
N_est = ; %estimate of how many simulated molecules in one cell at steady state




%% initialize

n_iter = 0;
%calculate input stream properties
stream_mfp = (T0/Tref)^(omega-0.5)/sqrt(2)/n_input/pi/dref^2; %stream mean free path with VHS model
if (xmax-xmin)/numCellX > stream_mfp 
    error('Too few cells in x direction')
elseif (ymax-ymin)/numCellY > stream_mfp
    error('Too few cells in y direction')
elseif (ymax-ymin)/numCellY > stream_mfp
    error('Too few cells in z direction')
elseif numCellX*numCellY*numCellZ*N_est > 1e8 %an array of size above roughly 1e8*3 will be too large for matlab's preference, program designed to house ~27 molecules in a collision cell at steady state
        error('Too many molecules for matlab to process')
end

mex -setup C++;
mex generateVx.cpp;
mex moveMolecules.cpp;
mex collision.cpp;
mex sampleCells.cpp;

rng('shuffle');
simulation_time = 0; %simulation time
X = []; %stores positions of particles
V = []; %stores velocities of particles
index_sample = []; %stores indices of particles to each sampling cell
mct = NaN(numCellX*numCellY*numCellZ,1); %stores local mean collision time in each collision cell
mcs = zeros(numCellX*numCellY*numCellZ,1); %stores mean collision separation in each collision cell
mfp = NaN(numCellX*numCellY*numCellZ,1); %stores local mean free path in collision cell
N_sample = zeros(numCellX*numCellY*numCellZ,1); %stores number of particles in each collision cell
Vcell = zeros(numCellX*numCellY*numCellZ,3); %stores average velocity in each sampling cell
Tcell = zeros(numCellX*numCellY*numCellZ,3); %stores translational temperature in each sampling cell along three directions
sample_cell_volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/numCellX/numCellY/numCellZ;

%mesh for sampling cells
Xsample = linspace(xmin,xmax,numCellX+1);
Ysample = linspace(ymin,ymax,numCellY+1);
Zsample = linspace(zmin,zmax,numCellZ+1);
%indexing functions, note that this marks points at the max boundary as
%numCell+1, this needs to be corrected in the indexing part of the routine
index_sample_X = @(a)find(Xsample <= a, 1, 'last');
index_sample_Y = @(b)find(Ysample <= b, 1, 'last');
index_sample_Z = @(c)find(Zsample <= c, 1, 'last');

%these are for plotting purposes of the cell properties
centerX = (Xsample(1:end-1)+ Xsample(2:end))/2;
centerX = reshape(repmat(centerX,numCellY*numCellZ,1),numCellX*numCellY*numCellZ,1);
centerY = (Ysample(1:end-1)+ Ysample(2:end))/2;
centerY = repmat(reshape(repmat(centerY,numCellZ,1),numCellY*numCellZ,1),numCellX,1);
centerZ = (Zsample(1:end-1)+ Zsample(2:end))/2;
centerZ = repmat(centerZ',numCellX*numCellY,1);
center = [centerX centerY centerZ];


%read in skimmer geometry, note that Nskim has to point into the flow
%region
[Fskim, Vskim, Nskim] = stlread(geometryFile);
wall_thickness = 1e-4;
wall_position = -4.88e-3;
% Fskim = [1 2 3; 2 4 3; 1 5 2; 5 6 2; 6 5 8; 5 7 8];
% Vskim = [wall_position 0 zmin; wall_position 0 zmax; wall_position ymax zmin; wall_position ymax zmax; wall_position+wall_thickness 0 zmin; wall_position+wall_thickness 0 zmax; wall_position+wall_thickness ymax zmin; wall_position+wall_thickness ymax zmax];
% Nskim = [-1 0 0; -1 0 0; 0 -1 0; 0 -1 0; 1 0 0; 1 0 0];
% 


Cm = sqrt(2*kB*T0/mHe);
s = v_input/Cm; %speed ratio at input
stream_mct = stream_mfp/sqrt(8*kB*T0/pi/mHe); %stream mean collision time
stream_traversal = (xmax-xmin)/numCellX/v_input; %time to traverse a cell in x-direction
sigma_ref = pi*dref^2*(2*kB*Tref*2/mHe)^(omega-0.5)/gamma(2.5-omega); %a constant to be passed in the VHS collision model, not really the collisional cross section
sgmax = ones(numCellX*numCellY*numCellZ,1)*sigma_ref*(mHe/4/kB/T0)^(omega-1)*gamma(3.5-2*omega)/gamma(2.5-omega); 
%stores estimated maximum expected value of sigma*g in collision cell,
%initialized to average value at stagnation pressure

%% routine

while n_iter<max_iter
    time_step = min([stream_mct/5; stream_traversal/2; mct/5]); %default value set as 1/5 of mean collision time or 1/2 of cell traversal time
    Fn = n_input*Cm/2/sqrt(pi)*(exp(-s^2)+sqrt(pi)*s*(1+erf(s)))*(ymax-ymin)*(zmax-zmin)*time_step; %total influx of particles
    Wp = n_input*(xmax-xmin)*(ymax-ymin)*(zmax-zmin)/numCellX/numCellY/numCellZ/N_est; %such that number of element is the same as number of simulated molecules in steady state 
    N_entry = floor(Fn/Wp+rand); %number of entry molecules
    
    %initialize entry molecules
    entryX = NaN(N_entry,3);
    entryV = NaN(N_entry,3);
    entryX(:,1) = xmin;
    entryX(:,2) = ymin + (ymax-ymin)*rand(N_entry,1);
    entryX(:,3) = zmin + (zmax-zmin)*rand(N_entry,1);
    entryV(:,2:3) = Cm*sqrt(-log(rand(N_entry,2))).*sin(2*pi*rand(N_entry,2));
    entryV(:,1) = generateVx(T0,v_input,entryV(:,1));
    X = [entryX;X];
    V = [entryV;V];
    
    %move all molecules while applying boundary conditions
    [X, V] = moveMolecules(Tskim, time_step, X, V, Fskim, Vskim, Nskim);
    ind_exclude = X(:,1) > xmax | X(:,1) < xmin | X(:,2) > ymax | X(:,2) < ymin | X(:,3) > zmax | X(:,3) < zmin;
    X(ind_exclude,:) = [];
    V(ind_exclude,:) = [];
    
    %index all molecules to the corresponding cell
    indx_sample = arrayfun(index_sample_X,X(:,1));
    indx_sample(indx_sample == numCellX+1) = numCellX;
    indy_sample = arrayfun(index_sample_Y,X(:,2));
    indy_sample(indy_sample == numCellY+1) = numCellY;
    indz_sample = arrayfun(index_sample_Z,X(:,3));
    indz_sample(indz_sample == numCellZ+1) = numCellZ;
    index_sample = (indx_sample-1)*numCellY*numCellZ + (indy_sample-1)*numCellZ + indz_sample;
    
    
    %update cell number counts and sort the data by collision cell
    [ind_occupied, ia, ic] = unique(index_sample);
    occupation_counts = accumarray(ic, 1);
    N_sample(ind_occupied) = occupation_counts;
    partitions = cumsum(occupation_counts);
    partitions = [0;partitions]; %change for interfacing with C++
    center_occupied = center(ind_occupied);
    [index_sample, order] = sort(index_sample);
    X = X(order,:);
    V = V(order,:);
    
    
    %collision routine with VHS model
    [V, mct, mcs, sgmax] = collision(sigma_ref,Wp,time_step,sample_cell_volume,X,V,partitions,sgmax,mct,mcs);
    
    %sample the cell properties
    mfp = sqrt(sum(Tcell,2)*8*kB/3/pi/mHe).*mct; %use mean thermal speed and mct to estimate mfp
    [Vcell(ind_occupied,:), Tcell(ind_occupied,:)] = sampleCells(V,partitions,Vcell(ind_occupied,:),Tcell(ind_occupied,:));
    
    
    %if want to plot cell properties, each of them can be anchored to the
    %cell center position variable center or one could scatter plot X
    %directly
    end
    n_iter = n_iter + 1;
end
