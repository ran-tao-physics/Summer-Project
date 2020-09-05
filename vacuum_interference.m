%% constants and program initialization
id = fopen('DS2VD.DAT','w');
kB = 1.3806504e-23;
mHe = 6.6465e-27;

%% inputs to be written by the user
T0 = 37; %stagnation temperature
vref = sqrt(5*kB*T0/mHe); %uniform velocity for input
Tsurr = 293; %temperature of skimmer surface
dskimmer = 0.0005;
n_input = 2e21; %input number density from the left

% simulation boundaries
xmin = -0.015;
xmax = 0.035;
ymin = 0;
ymax = 0.025;

%assigned memory in MB
memory = 500;


%the vertices, note that for each segment the flow is on the right side
%when you look from starting point to ending point if the surface
V = [0.024948 0.01145+dskimmer/2-0.0001; 0 0.00011+dskimmer/2-0.0001; 0 0.0001+dskimmer/2-0.0001; 0.024948+8e-5 0.0114+dskimmer/2-0.0001+5.4106e-5];

%Number of sampling intervals along the skimmer and mount, default arc
%geometry needs 5 inputs
sample = [100 300 30 300 100];

%whether to use default arc+mount geometry or the polygons given by user-defined vertices and sampling
%intervals
geom_default = true;

%whether surface defined is closed or open if user chooses non-default
%geometry
is_closed = false;




%% skimmer geometry generation, could be replaced by reading in data from outside file
%radius and center of first arc
R1 = 83.507697/1000;
cx1 = -21.613382/1000;
cy1 = 80.772241/1000+dskimmer/2-0.0001;
%radius and center of second arc
R2 = 81.5278855/1000;
cx2 = -20.685868/1000;
cy2 = 78.959945/1000+dskimmer/2-0.0001;
%generate surface for default arc geometry
thetai1 = atan((cy1-V(1,2))/(V(1,1)-cx1));
thetaf1 = atan((cy1-V(2,2))/(V(2,1)-cx1));
theta1 = linspace(thetai1, thetaf1, sample(2)+1);
arcx1 = cx1 + R1*cos(theta1);
arcy1 = cy1 - R1*sin(theta1);
arc1 = [arcx1; arcy1];
arc1(:,1) = V(1,:)';
arc1(:,end) = V(2,:)';
thetai2 = atan((cy2-V(3,2))/(V(3,1)-cx2));
thetaf2 = atan((cy2-V(4,2))/(V(4,1)-cx2));
theta2 = linspace(thetai2, thetaf2, sample(4)+1);
arcx2 = cx2 + R2*cos(theta2);
arcy2 = cy2 - R2*sin(theta2);
arc2 = [arcx2; arcy2];
arc2(:,1) = V(3,:)';
arc2(:,end) = V(4,:)';


%% write in simulation parameters
sim_cond1 = [4 5 memory 1 1]; 
frame_cond = [xmin xmax ymin ymax]; 
sim_cond2 = [1 1 1 0 0 0];
fprintf(id, '%d\n', sim_cond1);
fprintf(id, '%f\n', frame_cond); 
fprintf(id, '%d\n', sim_cond2);
fprintf(id, '%e\n', 2.33e-10); 
gas_cond = [273 0.66 1]; %Helium VHS parameters
fprintf(id, '%d\n%.2f\n%d\n', gas_cond);
fprintf(id, '%e\n', 6.6465e-27);
fprintf(id, 'Helium\n');
fprintf(id, '%d\n', [0 0]);
background_cond = [n_input T0 vref 1];
fprintf(id, '%d\n', background_cond);




%% write in surface data
surf_cond = [1 sum(sample)+1 sum(sample)+1 1 sum(sample) 0];
fprintf(id,'%d\n', surf_cond);
if geom_default
    fprintf(id,'%d\n', sample(2)+sample(4)+3);
    fprintf(id,'%d\n',1);
    fprintf(id,'%f\n',V(1,1));
    fprintf(id,'%f\n',ymax);
    fprintf(id,'%f\n',V(1,1));
    fprintf(id,'%f\n',V(1,2));
    fprintf(id,'%d\n',sample(1));
    for i=1:sample(2)
        fprintf(id,'%d\n',1);
        fprintf(id,'%.8f\n',arc1(1,i+1));
        fprintf(id,'%.8f\n',arc1(2,i+1));
        fprintf(id,'%d\n',1);
    end
    fprintf(id,'%d\n',1);
    fprintf(id,'%f\n',V(3,1));
    fprintf(id,'%f\n',V(3,2));
    fprintf(id,'%d\n',sample(3));
    for j=1:sample(4)
        fprintf(id,'%d\n',1);
        fprintf(id,'%.8f\n',arc2(1,j+1));
        fprintf(id,'%.8f\n',arc2(2,j+1));
        fprintf(id,'%d\n',1);
    end
    fprintf(id,'%d\n',1);
    fprintf(id,'%f\n',V(4,1));
    fprintf(id,'%f\n',ymax);
    fprintf(id,'%d\n',sample(5));
else
    fprintf(id, '%d\n',length(sample));
    for k = 1:length(V)
        fprintf(id,'%d\n',1);
        if k == 1
            fprintf(id,'%f\n',V(k,1));
            fprintf(id,'%f\n',V(k,2));
        end
        if k == length(V) && is_closed
            fprintf(id,'%f\n',V(1,1));
            fprintf(id,'%f\n',V(1,2));
        else
            fprintf(id,'%f\n',V(k+1,1));
            fprintf(id,'%f\n',V(k+1,2));
        end
        fprintf(id,'%d\n',sample(k));
    end
end

group_cond = [1 sum(sample) 1 Tsurr 0 0 0 0 0 0];
%[0 10 1 0.1 0 0 0 0 0 1];
fprintf(id, '%d\n', group_cond);


%% write in customary options for input
sim_cond = [3 4 1 4 2 0 0 1 1 10000 1 10 1 5 2 1 1 1 1 1]; 
fprintf(id, '%d\n', sim_cond);
fclose('all');