% This is a modification of the powder function of spinach 2.6.5625
% to handle the powder average over the three grids for the doublet,
% triplet and the dipolar grid

function [answer,parameters]=powder_lideer_soft_pulse(spin_system_t_hilb, ...
    spin_system_t2,spin_system_d,spin_system_dt, pulse_sequence,parameters, ...
    assumptions,rho_polarized)

if ~isfield(parameters,'sum_up')
    parameters.sum_up=1;
end


% Set the assumptions
spin_system_t_hilb=assume(spin_system_t_hilb,assumptions);
spin_system_t2=assume(spin_system_t2,assumptions);
spin_system_d=assume(spin_system_d,assumptions);
spin_system_dt=assume(spin_system_dt,assumptions);

% Get the Hamiltonians
[I_t_hilb,Q_t_hilb]=hamiltonian(spin_system_t_hilb);
[I_t2,Q_t2]=hamiltonian(spin_system_t2);
[I_d,Q_d]=hamiltonian(spin_system_d);
[I_dt,Q_dt]=hamiltonian(spin_system_dt);




% Get the averaging grid as a structure
folder_grids='grids\';

sph_grid_d=load(fullfile(folder_grids,parameters.grid_d),'alphas','betas','gammas','weights');
sph_grid_t=load(fullfile(folder_grids,parameters.grid_t),'alphas','betas','gammas','weights');
sph_grid_dt=load(fullfile(folder_grids,parameters.grid_dt),'alphas','betas','gammas','weights');


% Assign local variables


alphas_t=sph_grid_t.alphas; betas_t=sph_grid_t.betas;
gammas_t=sph_grid_t.gammas; weights_t=sph_grid_t.weights;




% Preallocate answer array
ans_array=cell(numel(weights_t),1);

% Run serially if needed
if isfield(parameters,'serial')&&...
        parameters.serial

    % Serial execution
    nworkers=0;

    % Inform the user
    report(spin_system_t2,'WARNING: parallelisation turned off by the user.');

else

    % Parallel execution
    nworkers=min([poolsize numel(weights_t)]);

end




% MDCS diagnostics
spin_system=spin_system_t2;
parallel_profiler_start;


count_array=zeros(1,length(weights_t));

% Powder averaging loop
%parfor (n=1:numel(weights_t),nworkers)
for n=1:numel(weights_t)

    % Localise the parameter array

    sph_grid_d_loc=sph_grid_d;
    alphas_d=sph_grid_d_loc.alphas; betas_d=sph_grid_d_loc.betas;
    gammas_d=sph_grid_d_loc.gammas; weights_d=sph_grid_d_loc.weights;

    sph_grid_dt_loc=sph_grid_dt;
    alphas_dt=sph_grid_dt_loc.alphas; betas_dt=sph_grid_dt_loc.betas;
    gammas_dt=sph_grid_dt_loc.gammas; weights_dt=sph_grid_dt_loc.weights;

    len_d_grid=length(weights_d);
    len_dt_grid=length(weights_dt);
    localpar=parameters;
    pulse_seq_loc=pulse_sequence;
    rho_polarized_loc=rho_polarized;

    % Get the full Hamiltonian of the triplet at the current triplet orientation
    H_t_hilb=I_t_hilb+orientation(Q_t_hilb,[alphas_t(n) betas_t(n) gammas_t(n)]); H_t_hilb=(H_t_hilb+H_t_hilb')/2;

    % calculate the eigenvalues to see whether the transition in the triplet
    % are close to the excitation bandwidth of the observer pulses
    levels=eig(H_t_hilb)/2/pi;

    transitions=abs([levels(1)-levels(2), levels(2)-levels(3), levels(3)-levels(1)]);

    % The simulation for the current orientation is only carried out if at
    % least one of the three transitions lies reasonalby within the 
    % excitation bandwidth of the observer pulses, otherwise this
    % orientation is skipped

    if((abs(transitions(1))>abs(parameters.carrier_frq)-parameters.bandwidth && abs(transitions(1))<abs(parameters.carrier_frq)+parameters.bandwidth) || ...
            (abs(transitions(2))>abs(parameters.carrier_frq)-parameters.bandwidth && abs(transitions(2))<abs(parameters.carrier_frq)+parameters.bandwidth) || ...
            (abs(transitions(3))>abs(parameters.carrier_frq)-parameters.bandwidth && abs(transitions(3))<abs(parameters.carrier_frq)+parameters.bandwidth) )



        count_array(n)=1;
        H_t2=I_t2+orientation(Q_t2,[alphas_t(n) betas_t(n) gammas_t(n)]); H_t2=(H_t2+H_t2')/2;



        % Create the triplet density operator
        rho_triplet=rho_polarized_loc(parameters.Px,parameters.Py,parameters.Pz,...
            betas_t(n), gammas_t(n), H_t_hilb);

        localpar.rho0=kron(parameters.rho_eq,rho_triplet);



        for k=1:len_d_grid



            % Get the full Hamiltonian at the current orientation
            H_d=I_d+orientation(Q_d,[alphas_d(k) betas_d(k) ...
                gammas_d(k)]); H_d=(H_d+H_d')/2;

            for j=1:len_dt_grid

                % Get the full Hamiltonian at the current orientation
                H_dt=I_dt+orientation(Q_dt,[alphas_dt(j) betas_dt(j) ...
                    gammas_dt(j)]); H_dt=(H_dt+H_dt')/2;

                H=H_dt+H_d+H_t2;

                % Run the simulation )
                ans_array{n,k,j}=pulse_seq_loc(spin_system,localpar,H,zeros(size(H)),zeros(size(H)));

            end
        end


    else
        for k=1:len_d_grid
            for j=1:len_dt_grid
                ans_array{n,k,j}=0;
            end
        end
    end


end

parameters.counts=sum(count_array);

% Decide the return array
if parameters.sum_up

    % Return weighted sum
    answer=0*ans_array{1,1,1};
    for n=1:length(weights_t)
        for k=1:length(sph_grid_d.weights)
            for j=1:length(sph_grid_dt.weights)
                answer=answer+weights_t(n)*sph_grid_d.weights(k)*sph_grid_dt.weights(j)*ans_array{n,k,j};
            end
        end
    end

    % Inform the user
    report(spin_system,'returning powder averaged pulse sequence output...');

else

    % Return components and weights
    answer.components=ans_array;
    answer.weights=sph_grid_t.weights;
    % Inform the user
    report(spin_system,'returning pulse sequence outputs at each orientation...');

end



end