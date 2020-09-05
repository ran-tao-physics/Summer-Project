id = fopen('DS2VD.DAT','w');
kB = 1.3806504e-23;
mHe = 6.6465e-27;


%% inputs to be written by the user
T0 = 37;
vref = sqrt(5*kB*T0/mHe);
Tsurr = 293;
dskimmer = 5e-4;
n_input = 1e20;
memory = 300; %assigned memory

xmin = -dskimmer*20/2;
xmax = dskimmer*20/2;
ymin = 0;
ymax = dskimmer*20/2;

%Vertices of the skimmer tip, the geometry is assumed to be open since we
%are only concerned with the tip here I defined the two functions to
%calculate the coords of an arc in the default geometry
V = [];
%Number of sampling intervals along the skimmer surface 
sample = [300 30 300];

%whether to use default arc+mount geometry or the polygons given by user-defined vertices and sampling
%intervals
geom_default = true;



%% default geometry
if geom_default
    %radius and center of first arc
    R1 = 83.507697/1000;
    cx1 = -21.613382/1000;
    cy1 = 80.772241/1000+dskimmer/2-0.0001;
    %radius and center of second arc
    R2 = 81.5278855/1000;
    cx2 = -20.685868/1000;
    cy2 = 78.959945/1000+dskimmer/2-0.0001;
    %the vertices
    y1calc = @(x) cy1-sqrt(R1^2-(x-cx1)^2);
    y2calc = @(x) cy2-sqrt(R2^2-(x-cx2)^2);
    V = [xmax y1calc(xmax); 0 0.00001+dskimmer/2; 0 dskimmer/2; xmax y2calc(xmax)];


    %generate surface
    thetai1 = atan((cy1-V(1,2))/(V(1,1)-cx1));
    thetaf1 = atan((cy1-V(2,2))/(V(2,1)-cx1));
    theta1 = linspace(thetai1, thetaf1, sample(1)+1);
    arcx1 = cx1 + R1*cos(theta1);
    arcy1 = cy1 - R1*sin(theta1);
    arc1 = [arcx1; arcy1];
    arc1(:,1) = V(1,:)';
    arc1(:,sample(1)+1) = V(2,:)';
    thetai2 = atan((cy2-V(3,2))/(V(3,1)-cx2));
    thetaf2 = atan((cy2-V(4,2))/(V(4,1)-cx2));
    theta2 = linspace(thetai2, thetaf2, sample(3)+1);
    arcx2 = cx2 + R2*cos(theta2);
    arcy2 = cy2 - R2*sin(theta2);
    arc2 = [arcx2; arcy2];
    arc2(:,1) = V(3,:)';
    arc2(:,end) = V(4,:)';
end

%% write in simulation parameters
sim_cond1 = [4 5 memory 1 1]; 
frame_cond = [xmin xmax ymin ymax]; 
sim_cond2 = [1 1 1 0 0 0];%n_input*kB*Tsurr/Pvacuum
fprintf(id, '%d\n', sim_cond1);
fprintf(id, '%f\n', frame_cond); %depends on input size
fprintf(id, '%d\n', sim_cond2);
fprintf(id, '%e\n', 2.33e-10); 
gas_cond = [273 0.66 1]; %Helium
fprintf(id, '%d\n%.2f\n%d\n', gas_cond);
fprintf(id, '%e\n', 6.6465e-27);
fprintf(id, 'Helium\n');
fprintf(id, '%d\n', [0 0]);
background_cond = [n_input T0 vref 1];
fprintf(id, '%d\n', background_cond);



%%skimmer data
surf_cond = [1 sum(sample)+1 sum(sample)+1 1 sum(sample) 0];
fprintf(id,'%d\n', surf_cond);
if geom_default
    fprintf(id,'%d\n', sample(1)+sample(3)+1);
    for i=1:sample(1)
        fprintf(id,'%d\n',1);
        if i == 1
            fprintf(id,'%.12f\n',arc1(1,i));
            fprintf(id,'%.12f\n',arc1(2,i));
        end
        fprintf(id,'%.12f\n',arc1(1,i+1));
        fprintf(id,'%.12f\n',arc1(2,i+1));
        fprintf(id,'%d\n',1);
    end
    fprintf(id,'%d\n',1);
    fprintf(id,'%f\n',V(3,1));
    fprintf(id,'%f\n',V(3,2));
    fprintf(id,'%d\n',sample(2));
    for j=1:sample(3)
        fprintf(id,'%d\n',1);
        fprintf(id,'%.12f\n',arc2(1,j+1));
        fprintf(id,'%.12f\n',arc2(2,j+1));
        fprintf(id,'%d\n',1);
    end
else
    fprintf(id,'%d\n', length(sample));
    for k=1:length(sample)
        fprintf(id,'%d\n',1);
        if k == 1
            fprintf(id,'%f\n',V(1,1));
            fprintf(id,'%f\n',V(1,2));
        end
        fprintf(id,'%f\n',V(k+1,1));
        fprintf(id,'%f\n',V(k+1,2));
        fprintf(id,'%d\n',sample(k));
    end
end

group_cond = [1 sum(sample) 1 Tsurr 0 0 0 0 0 0];
fprintf(id, '%d\n', group_cond);


sim_cond = [3 4 1 4 2 0 0 1 1 10000 1 10 1 5 2 1 1 1 1 1]; %final boolean is nearest neighbor or not
fprintf(id, '%d\n', sim_cond);
fclose('all');