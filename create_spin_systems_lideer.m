function [spin_system_t2, spin_system_d, spin_system_dt, spin_system_t_liouv, spin_system_t_hilb, spin_system_d1]=...
    create_spin_systems_lideer(g_1, g_2, B, r, D, E, T, spin_system, bas)

mu_0 = 1.25663706e-6; %vacuumm permeability
hplanck = 6.626070040e-34; % Js
mu_e=927.4009994e-26; %J/T
g_e=-2.00231930436182;

%%
spin_system_t2.isotopes={'E','E3'};
spin_system_t2.magnet=B;
inter_t2.zeeman.matrix={zeros(3,3),g_2};
inter_t2.temperature=T; %K
inter_t2.coupling.matrix=cell(2,2);
inter_t2.coupling.matrix{2,2}=[-1/3*D+E,0,0;0,-1/3*D-E,0;0,0,2/3*D];

spin_system_t2.disable= spin_system.disable;
spin_system_t2.disable= spin_system.disable;
spin_system_t2.output= spin_system.output;
spin_system_t2.enable= spin_system.enable;
if isfield(spin_system, 'scratch')
    spin_system_t2.scratch= spin_system.scratch;
end

spin_system_t2=create(spin_system_t2,inter_t2);
spin_system_t2=basis(spin_system_t2,bas);

%%

spin_system_d.isotopes={'E','E3'};
spin_system_d.magnet=B;
inter_d.zeeman.matrix={g_1,zeros(3,3)};
inter_d.temperature=T; %K

spin_system_d.disable= spin_system.disable;
spin_system_d.disable= spin_system.disable;
spin_system_d.output= spin_system.output;
spin_system_d.enable= spin_system.enable;
if isfield(spin_system, 'scratch')
    spin_system_d.scratch= spin_system.scratch;
end

spin_system_d=create(spin_system_d,inter_d);
spin_system_d=basis(spin_system_d,bas);



%%

D_dipolar=mu_0/4/pi*mu_e^2*g_e^2/r^3/hplanck;

spin_system_dt.isotopes={'E','E3'};
spin_system_dt.magnet=B;
inter_dt.zeeman.matrix={zeros(3,3), zeros(3,3)};
inter_dt.coupling.matrix=cell(2,2);
dipolar=D_dipolar*[-1,0,0;0,-1,0;0,0,2];
inter_dt.coupling.matrix{1,2}=dipolar;
%inter_dt.coordinates={[0 0 0] ...
%                       [r*1e10 0 0]};
inter_dt.temperature=T; %K

spin_system_dt.disable= spin_system.disable;
spin_system_dt.disable= spin_system.disable;
spin_system_dt.output= spin_system.output;
spin_system_dt.enable= spin_system.enable;
if isfield(spin_system, 'scratch')
    spin_system_dt.scratch= spin_system.scratch;
end

spin_system_dt=create(spin_system_dt,inter_dt);
spin_system_dt=basis(spin_system_dt,bas);



%%
spin_system_t_liouv.isotopes={'E3'};
spin_system_t_liouv.magnet=B;
inter_t_liouv.zeeman.matrix={g_2};
inter_t_liouv.temperature=T; %K
inter_t_liouv.coupling.matrix=cell(1,1);
inter_t_liouv.coupling.matrix{1,1}=[-1/3*D+E,0,0;0,-1/3*D-E,0;0,0,2/3*D];

spin_system_t_liouv.disable= spin_system.disable;
spin_system_t_liouv.disable= spin_system.disable;
spin_system_t_liouv.output= spin_system.output;
spin_system_t_liouv.enable= spin_system.enable;
if isfield(spin_system, 'scratch')
    spin_system_t_liouv.scratch= spin_system.scratch;
end

spin_system_t_liouv=create(spin_system_t_liouv,inter_t_liouv);
spin_system_t_liouv=basis(spin_system_t_liouv,bas);


%%

spin_system_t_hilb.isotopes={'E3'};
spin_system_t_hilb.magnet=B;
inter_t_hilb.zeeman.matrix={g_2};
inter_t_hilb.temperature=T; %K
inter_t_hilb.coupling.matrix=cell(1,1);
inter_t_hilb.coupling.matrix{1,1}=[-1/3*D+E,0,0;0,-1/3*D-E,0;0,0,2/3*D];

spin_system_t_hilb.disable= spin_system.disable;
spin_system_t_hilb.disable= spin_system.disable;
spin_system_t_hilb.output= spin_system.output;
spin_system_t_hilb.enable= spin_system.enable;
if isfield(spin_system, 'scratch')
    spin_system_t_hilb.scratch= spin_system.scratch;
end
bas2=bas;
bas2.formalism='zeeman-hilb';
spin_system_t_hilb=create(spin_system_t_hilb,inter_t_hilb);
spin_system_t_hilb=basis(spin_system_t_hilb,bas2);


%%

spin_system_d1.isotopes={'E'};
spin_system_d1.magnet=B;
inter_d1.zeeman.matrix={g_1};
inter_d1.temperature=T; %K

spin_system_d1.disable= spin_system.disable;
spin_system_d1.disable= spin_system.disable;
spin_system_d1.output= spin_system.output;
spin_system_d1.enable= spin_system.enable;
if isfield(spin_system, 'scratch')
    spin_system_d1.scratch= spin_system.scratch;
end

spin_system_d1=create(spin_system_d1,inter_d1);
spin_system_d1=basis(spin_system_d1,bas);




end