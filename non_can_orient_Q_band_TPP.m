clear all;
clc;
constants;


spin_system.disable={'trajlevel', 'krylov'};
spin_system.output='hush';

spin_system.enable={'caching'};


bas.formalism='sphten-liouv';
bas.approximation='none';

g_1=2.00687*eye(3); g_2=2*eye(3);
B=1.2097;     %T
D=1.159e9;  %Hz
E=-0.238e9; %Hz
parameters.Px=0.33;                           
parameters.Py=0.41;                            
parameters.Pz=0.26;
r=2.2e-9;   %m
T=30;       %K

[spin_system_t2, spin_system_d, spin_system_dt, spin_system_t_liouv, spin_system_t_hilb, spin_system_d1]=...
    create_spin_systems_lideer(g_1, g_2, B, r, D, E, T, spin_system, bas);

R=sphten2zeeman(spin_system_t_liouv);
R=inv(R);

Ep=operator(spin_system_t2,'L+','E3')+operator(spin_system_t2,'L+','E');
parameters.Ex_obs=(Ep+Ep')/2; parameters.Ey_obs=(Ep-Ep')/2i;

Ep_pump=operator(spin_system_t2,'L+','E');
parameters.Ex_pump=(Ep_pump+Ep_pump')/2; parameters.Ey_pump=(Ep_pump-Ep_pump')/2i;

parameters.coil_prob=state(spin_system_t2,{'L-'},{2})+state(spin_system_t2,{'L-'},{1});


parameters.rho_eq=my_equilibrium(spin_system_d1);
parameters.grid_t='rep_3ang_12800pts';  
parameters.grid_d='Grid_one_orient';  
parameters.grid_dt='rep_1ang_200pts_hem'; 
parameters.counts=0;

parameters.method='expm';
parameters.pump_frq=-33.979e9;
parameters.carrier_frq=-33.25e9; 
parameters.obs_frq=parameters.carrier_frq;
parameters.bandwidth=100e6;
parameters.zero_fill=100000;

% DEER pulse timings
parameters.p1_p2_gap=0.1e-6;  
parameters.p2_p4_gap=2.1e-6;  
parameters.nsteps=250;
parameters.stepsize=8e-9;
parameters.int_stepsize=1e-9;
parameters.int_nsteps=50;


% Pulse parameter first pulse, observer pulse, pi/2
parameters.pulse_frq(1)=parameters.obs_frq;
parameters.pulse_dur(1)=10e-9;
parameters.pulse_pwr(1)=pi/2/parameters.pulse_dur(1)/sqrt(2);
parameters.pulse_rnk(1)=2;

% Pulse parameter second pulse, observer pulse
parameters.pulse_frq(2)=parameters.obs_frq;
parameters.pulse_dur(2)=20e-9;
parameters.pulse_pwr(2)=pi/parameters.pulse_dur(2)/sqrt(2);
parameters.pulse_rnk(2)=2;

% Pulse parameter third pulse, pump pulse
parameters.pulse_frq(3)=parameters.pump_frq;
parameters.pulse_dur(3)=10e-9;
parameters.pulse_pwr(3)=pi/parameters.pulse_dur(3);
parameters.pulse_rnk(3)=2;


% Pulse parameter fourth pulse, observer pulse
parameters.pulse_frq(4)=parameters.obs_frq;
parameters.pulse_dur(4)=20e-9;
parameters.pulse_pwr(4)=pi/parameters.pulse_dur(4)/sqrt(2);
parameters.pulse_rnk(4)=2;



parameters.tau_start=0.15e-6+8e-9;  
parameters.tau=linspace(parameters.tau_start,parameters.tau_start+parameters.stepsize*parameters.nsteps,parameters.nsteps+1);
parameters.echo_max=4.244e-6;

parameters.t0=0;
parameters.t1=parameters.pulse_dur(1);
parameters.t2=parameters.t1+parameters.p1_p2_gap;
parameters.t3=parameters.t2+parameters.pulse_dur(2);
parameters.t4=parameters.tau(1);
parameters.t5=parameters.tau(end);
parameters.t6=parameters.t3+parameters.p2_p4_gap;
parameters.t7=parameters.t6+parameters.pulse_dur(4);
parameters.t8=parameters.echo_max-parameters.int_nsteps*parameters.int_stepsize/2;
parameters.t9=parameters.t8+parameters.int_nsteps*parameters.int_stepsize;
parameters.int_time_axis=linspace(parameters.t8,parameters.t9,parameters.int_nsteps+1);

parameters.pulse_phi(1)=0;
parameters.pulse_phi(2)=2*pi*parameters.obs_frq*parameters.t2;
parameters.pulse_phi(3)=2*pi*parameters.obs_frq*parameters.t4;
parameters.pulse_phi(4)=2*pi*parameters.obs_frq*parameters.t6;


parameters.transient=false;


%%
tic
[deer,parameters]=powder_lideer_soft_pulse(spin_system_t_hilb,spin_system_t2,spin_system_d,...
    spin_system_dt,@DEER_obs_soft_pump_soft,parameters,'labframe',...
    @(Px,Py,Pz,beta,gamma,H)spin_polarized_one_spin(Px,Py,Pz,beta,gamma,R,H)); toc
simulation_time=toc;
%%

[deer_pc,phi]=phase_corr(deer);

figure(1); hold all;
plot(parameters.tau,real(deer_pc));
plot(parameters.tau,imag(deer_pc));
%%
M=[parameters.tau(:)*1e9,imag(deer(:)),real(deer(:))];         
writematrix(M, 'non_can_orient_Q_band_TPP.txt');

save('non_can_orient_Q_band_TPP.mat','spin_system','parameters','D','E','simulation_time','phi','r') 

