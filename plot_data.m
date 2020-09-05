kB = 1.3806504e-23;
mHe = 6.6465e-27;
T0 = 37; %stagnation temperature
vref = sqrt(5*kB*T0/mHe);
dskimmer = 0.0005;
n_input = 5e20;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');



file = ;%.dat file produced by DS2V
data = prepare_data(file);
nonzero = data(data(:,1)~=0,1);
epsilon = min(abs(nonzero)); 
epsilon = 2*epsilon; % tolerance around point zero
    
sample_range = data(data(:,2) < dskimmer/2 + epsilon,:);
sample_range = sample_range(sample_range(:,1) > -epsilon,:); %construct roughly where we are sampling and getting the unique z values to plot over
Z = unique(sample_range(:,1));
Z = sort(Z);
R = linspace(0, dskimmer/2, 200);
[zinterp, rinterp] = meshgrid(Z,R);
n = scatteredInterpolant(data(:,1:2),data(:,3),'nearest');
v = scatteredInterpolant(data(:,1:2),data(:,5),'nearest');
I = n(zinterp, rinterp).*v(zinterp, rinterp);
N = pi*R(1)*R(1)*I(1,:)+trapz(R(2:end), 2*pi*R(2:end)'.*I(2:end,:), 1);
%can plot axial data if needed
axial = sample_range(sample_range(:,2) < epsilon,:);
axial = sortrows(axial,1);


plot(Z,N/(n_input*pi*dskimmer^2/4*vref))
xlabel('Distance along axis / m')
ylabel('Normalized Entry Flux')