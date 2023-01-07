clc;
clear;
close all;

%% geometry information
ftable=readtable('ModelInfo.xlsx','sheet','GeometryInfo');
Geometry_parameter = table2struct(ftable,'ToScalar',true);
x_start_stage = Geometry_parameter.Coord_Start_m;
x_end_stage = Geometry_parameter.Coord_End_m;
thickness_stage = 10^-3*Geometry_parameter.Thickness_mm;
diameter_stage = Geometry_parameter.Diameter_m;

%% lumped mass information
ftable=readtable('ModelInfo.xlsx','sheet','MassInfo');
Mass_parameter = table2struct(ftable,'ToScalar',true);
coordinate_node = Mass_parameter.Coord_X_m;
num_node = length(coordinate_node); % number of elements
load_mass = Mass_parameter.LoadMass_kg;
fuel_mass = Mass_parameter.FuelMass_kg;

%% fuel consumption rate
ftable=readtable('ModelInfo.xlsx','sheet','FuelConsumptionRate');
FuelConsumptionRate = table2struct(ftable,'ToScalar',true);
fuel_rate = FuelConsumptionRate.FuelConsumptionRate_kg_s_1;

%% Material information
ftable=readtable('ModelInfo.xlsx','sheet','MaterialInfo');
Material_parameter = table2struct(ftable,'ToScalar',true);
Young_modulus = Material_parameter.Modulus_Mpa*10^5;
density = Material_parameter.Density_kg_m_3;
Possion_ratio = Material_parameter.PossionRatio;

%%
stiff_matrix_element = zeros(4,4,num_node-1); % 
mass_matrix_element = zeros(4,4,num_node-1);
for ii = 1:num_node-1
    L_element = coordinate_node(ii+1)-coordinate_node(ii);
    x_stage = [x_start_stage,x_end_stage];
    temp = bsxfun(@minus,x_stage,coordinate_node(ii:ii+1,1)');
    ind = find((temp(:,1).*temp(:,2))<=0);
    ii_in_stage = ind(1);
    outer_radius = diameter_stage(ii_in_stage)*0.5;
    inner_radius = outer_radius-thickness_stage(ii_in_stage);
    [Me,Ke] = DynamicMatrix(L_element,density,Young_modulus,inner_radius,outer_radius,'pipe');
    stiff_matrix_element(:,:,ii) = Ke;
    mass_matrix_element(:,:,ii) = Me;
end

%% assemble the element matrix
M_str = zeros(2*num_node,2*num_node); % mass matrix only considering the structural mass
K = zeros(2*num_node,2*num_node); % stiffness matrix
for ii = 1:num_node-1
    ind = (2*ii-1):(2*ii+2);
    M_str(ind,ind) = M_str(ind,ind)+mass_matrix_element(:,:,ii);
    K(ind,ind) = K(ind,ind)+stiff_matrix_element(:,:,ii);
end

%% propotional/Rayleigh damping matrix
C = 7*M_str+7e-5*K;

%% baseline values of modal frequencies and mode shape vectors
ind_lumped_mass = 1:2:2*num_node; % DOFs with lumped mass

ts = 72; % maximal time to calculate the baseline values
fs_FRF_max = 60; % maximal frequency to calculate the frequency response function
fs_FRF_min = 2; % minimal frequency to calculate the frequency response function

t = 0:2:ts; % time instants to calculate the baseline values and FRFs
N_mode = 4; % number of modes
TrueFreq = zeros(N_mode,length(t)); % baseline frequencies
TrueDamp = zeros(N_mode,length(t)); % baseline damping ratios
TrueMode = cell(N_mode,1); % baseline mode shape vectors
No = size(M_str,1); 
faxis = linspace(fs_FRF_min,fs_FRF_max,1000); % frequency lines to calculate FRFs
FRF_matrix = zeros(No,No,length(faxis),length(t));
for ii = 1:length(t)
    tii = t(ii);
    fuel_mass_remain = FuelMassCalculation(fuel_mass,tii,fuel_rate); % the remained fuel mass of each node at time tii
    lumped_mass = fuel_mass_remain+load_mass; % total lumped mass at time tii
    temp = zeros(2*num_node,1);
    temp(ind_lumped_mass) = lumped_mass;
    M = M_str+diag(temp);
    [modes,freq]=struct_eigMK(M,K);
    TrueFreq(:,ii) = freq(3:3+N_mode-1,1); % the first and second frequencies are due to rigid body motion
    for jj = 1:N_mode
        TrueMode{jj,1}(:,ii) = modes(1:2:end,3+jj-1); % mode shape vectors of translational DOFs
    end
    [modes,freq,damp]=struct_eigMCK(M,C,K);
    TrueDamp(:,ii) = damp(41:41+N_mode-1,1)*100; % damping ratios in percentage
    
    for jj = 1:length(faxis)
        ff = faxis(jj)*2*pi*1i;
        dynamic_stiff = M*ff^2+C*ff+K;
        FRF_matrix(:,:,jj,ii) = dynamic_stiff\eye(No); % FRF matrix
    end
end
% save('Baseline.mat','t','TrueFreq','TrueDamp','TrueMode');

%% plot baseline values of modal frequencies and damping ratios
color = [255 0 0;0 0 255;255 30 255;85 151 85;0 0 0;0 0 128]/255;%
% set(groot,'DefaultAxesColorOrder',color,'DefaultAxesLineStyleOrder','o-|s-|d-|>')
figure
plot(t,TrueFreq(1,:),'r','Marker','^','MarkerFaceColor','r','Markersize',4);
hold on
plot(t,TrueFreq(2,:),'b','Marker','d','MarkerFaceColor','b','Markersize',4);
hold on
plot(t,TrueFreq(3,:),'color',color(3,:),'Marker','s','MarkerFaceColor',color(3,:),'Markersize',5);
hold on
plot(t,TrueFreq(4,:),'color',color(4,:),'Marker','o','MarkerFaceColor',color(4,:),'Markersize',4);
ylabel('Frequency/Hz')
xlabel('Time/s');
legend('1st mode','2nd mode','3rd mode','4th mode');
set(gca,'Fontsize',14,'Linewidth',1);
xlim([t(1),t(end)])

figure
plot(t,TrueDamp(1,:),'r','Marker','^','MarkerFaceColor','r','Markersize',4);
hold on
plot(t,TrueDamp(2,:),'b','Marker','d','MarkerFaceColor','b','Markersize',4);
hold on
plot(t,TrueDamp(3,:),'color',color(3,:),'Marker','s','MarkerFaceColor',color(3,:),'Markersize',5);
hold on
plot(t,TrueDamp(4,:),'color',color(4,:),'Marker','o','MarkerFaceColor',color(4,:),'Markersize',4);
ylabel('Damping ratio/%')
xlabel('Time/s');
legend('1st mode','2nd mode','3rd mode','4th mode');
set(gca,'Fontsize',14,'Linewidth',1);
xlim([t(1),t(end)])

% plot frozen-time FRF function
H11 = squeeze(FRF_matrix(1,1,:,:));
figure
for tii = 1:length(t)
    scatter3(ones(length(faxis),1)*t(tii),faxis,abs(H11(:,tii)),2,log10(abs(H11(:,tii))),'v');
    hold on
end
set(gca,'ZScale','log')
% surf(t,faxis,10*log10(abs(H11)));
zlabel('Amp./m·N^-^1');
ylabel('Frequency/Hz','rotation',-1)
xlabel('Time/s','rotation',85);
set(gca,'Fontsize',14,'Linewidth',1);
colormap('jet')
% shading interp
xlim([0,t(end)])
ylim([0,faxis(end)])
view([93,30])

%% Newmark-beta method to compute dynamic responses
beta = 0.5;  % parameters of Newmark-Beta method
gama = 0.25;

dt0 = 1/4000; % time integration step
ts = 75; % time duration of responses
t0 = 0:dt0:ts;
Nt0 = length(t0);
No = size(M,1);
fs = 128; % final sampling frequency of the response
t = linspace(0,ts,size(dis,1)); % time axis of the resampled response

% parameters of SDOF system to generate colored noise excitation
m_sdof = 1;  % mass of SDOF system
excitation_freq_start = 35;
excitation_freq_end = 45;
k_sdof = 2*pi*linspace(excitation_freq_start,excitation_freq_end,Nt0); % predominant frequency of non-white excitation changes linearly
extra_length = 50;  % extra length to eliminate the initial effect of SDOF system
k_sdof = [2*pi*excitation_freq_start*ones(1,extra_length),k_sdof];
k_sdof = k_sdof.^2;  % stiffness of SDOF system
c_sdof = 0.005*2*(k_sdof(1)*m_sdof)^0.5;  % damping of SDOF system

% add the SDOF parameters to the baseline values of modal parameters
N_baseline = size(TrueFreq,2);
temp_k = 2*pi*linspace(excitation_freq_start,excitation_freq_end,N_baseline);
temp_k = temp_k.^2;  % stiffness of SDOF system
temp_c = 0.005*2*(temp_k(1)*m_sdof)^0.5;  % damping of SDOF system
for tii = 1:N_baseline
    freq_sdof_tii = (temp_k(tii)/m_sdof)^0.5/2/pi;
    damp_sdof_tii = temp_c/(2*sqrt(temp_k(tii)*m_sdof));
    TrueFreq(N_mode+1,tii) = freq_sdof_tii;
    TrueDamp(N_mode+1,tii) = 100*damp_sdof_tii;
    TrueMode{N_mode+1,1} = NaN*zeros(size(TrueMode{N_mode,1}));
end
save('Baseline.mat','t','TrueFreq','TrueDamp','TrueMode');

% constants in the Newmark-beta method
alpha_0  = 1/gama/dt0^2; 
alpha_1  = beta/gama/dt0;
alpha_2  = 1/gama/dt0;
alpha_3  = 1/2/gama-1;
alpha_4  = beta/gama-1;
alpha_5  = dt0/2*(alpha_4-1);
alpha_6  = dt0*(1-beta);
alpha_7  = beta*dt0;

MCn = 200; % number of Monte Carlo simulations
for MCii = 1:MCn
    u = zeros(No,Nt0); % input force
    u(1:2:end,:) = randn(No/2,Nt0)*1000; % only translational DOFs are excited
    
    % non-white excitation
    u_sdof = randn(1,Nt0+extra_length)*2e8;  % white noise excitation of SDOF system
    dis_sdof = zeros(1,Nt0+extra_length);
    vel_sdof = zeros(1,Nt0+extra_length);
    acc_sdof = zeros(1,Nt0+extra_length);
    acc_sdof(:,1) = m_sdof\(u_sdof(:,1)-k_sdof(1)*dis_sdof(:,1)-c_sdof*vel_sdof(:,1));
    for ii = 1:Nt0+extra_length-1
        kt = k_sdof(ii);
        K_e = kt+alpha_0*m_sdof+alpha_1*c_sdof;
        
        F_e = u_sdof(:,ii+1)+m_sdof*(alpha_0*dis_sdof(:,ii)+alpha_2*vel_sdof(:,ii)+alpha_3*acc_sdof(:,ii))+...
            c_sdof*(alpha_1*dis_sdof(:,ii)+alpha_4*vel_sdof(:,ii)+alpha_5*acc_sdof(:,ii));
        dis_sdof(:,ii+1) = K_e\F_e;
        acc_sdof(:,ii+1) = alpha_0*(dis_sdof(:,ii+1)-dis_sdof(:,ii))-alpha_2*vel_sdof(:,ii)-alpha_3*acc_sdof(:,ii);
        vel_sdof(:,ii+1) = vel_sdof(:,ii)+alpha_6*acc_sdof(:,ii)+alpha_7*acc_sdof(:,ii+1);
    end
    u(1,:) = dis_sdof(1,extra_length+1:end);  % only the first translational DOF is excited with non-white excitation
    
    % responses computation by Newmark-beta method
    dis = zeros(No,Nt0);
    vel = zeros(No,Nt0);
    acc = zeros(No,Nt0);
    dis(:,1) = 0;  % initial values of displacement
    vel(:,1) = 0;  % initial values of velocity
    acc(:,1) = M\(u(:,1)-K*dis(:,1)-C*vel(:,1));  % initial values of acceleration
    
    for k = 1:Nt0-1
        tk = t0(k);
        fuel_mass_remain = FuelMassCalculation(fuel_mass,tk,fuel_rate);
        lumped_mass = fuel_mass_remain+load_mass;
        temp = zeros(2*num_node,1);
        temp(ind_lumped_mass) = lumped_mass;
        M = M_str+diag(temp);
        K_e = K+alpha_0*M+alpha_1*C;
        
        F_e = u(:,k+1)+M*(alpha_0*dis(:,k)+alpha_2*vel(:,k)+alpha_3*acc(:,k))+...
            C*(alpha_1*dis(:,k)+alpha_4*vel(:,k)+alpha_5*acc(:,k));
        dis(:,k+1) = K_e\F_e;
        acc(:,k+1) = alpha_0*(dis(:,k+1)-dis(:,k))-alpha_2*vel(:,k)-alpha_3*acc(:,k);
        vel(:,k+1) = vel(:,k)+alpha_6*acc(:,k)+alpha_7*acc(:,k+1);
    end
    
    % resample the responses
    dis = resample(dis',fs,1/dt0);
    dis = dis(:,1:2:end); % only consider the translational displacements
    
    savedir = strcat('.\Response');
    if ~exist(savedir,'dir')
        mkdir(savedir);
    end
    datafile = strcat(savedir,'\Dis_MC',num2str(MCii),'.mat');
    save(datafile,'t','dis');
end
