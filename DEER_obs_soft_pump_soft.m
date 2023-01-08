% <http://spindynamics.org/wiki/index.php?title=deer_4p_soft_deer.m>

function result=DEER_obs_soft_pump_soft(spin_system,parameters,H,R,K)

% Compose Liouvillian
L=H+1i*R+1i*K;

% First pulse
rho=shaped_pulse_af(spin_system,L,parameters.Ex_obs,parameters.Ey_obs,parameters.rho0,parameters.pulse_frq(1),parameters.pulse_pwr(1),...
    parameters.pulse_dur(1),parameters.pulse_phi(1),...
    parameters.pulse_rnk(1),parameters.method);

% Evolution
rho=evolution(spin_system,L,[],rho,parameters.t2-parameters.t1,1,'final');

% Second pulse
rho=shaped_pulse_af(spin_system,L,parameters.Ex_obs,parameters.Ey_obs,rho,parameters.pulse_frq(2),parameters.pulse_pwr(2),...
    parameters.pulse_dur(2),parameters.pulse_phi(2),...
    parameters.pulse_rnk(2),parameters.method);


% Evolution
rho=evolution(spin_system,L,[],rho,parameters.t4-parameters.t3,1,'final');
rho_stack=evolution(spin_system,L,[],rho,parameters.stepsize,parameters.nsteps,'trajectory');

% Third pulse
rho_stack=shaped_pulse_af(spin_system,L,parameters.Ex_pump,parameters.Ey_pump,rho_stack,parameters.pulse_frq(3),parameters.pulse_pwr(3),...
    parameters.pulse_dur(3),parameters.pulse_phi(3),...
    parameters.pulse_rnk(3),parameters.method);

% Evolution
rho_stack(:,end:-1:1)=evolution(spin_system,L,[],rho_stack(:,end:-1:1),parameters.stepsize,parameters.nsteps,'refocus');
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.t6-parameters.t5-parameters.pulse_dur(3),1,'final');

% Fourth pulse
rho_stack=shaped_pulse_af(spin_system,L,parameters.Ex_obs,parameters.Ey_obs,rho_stack,parameters.pulse_frq(4),parameters.pulse_pwr(4),...
    parameters.pulse_dur(4),parameters.pulse_phi(4),...
    parameters.pulse_rnk(4),parameters.method);

rho_stack=evolution(spin_system,L,[],rho_stack,parameters.t8-parameters.t7,1,'final');

echo=zeros(parameters.int_nsteps+1,size(rho_stack,2));


for k=1:size(rho_stack,2)
    
    rho_stack_echo=evolution(spin_system,L,[],rho_stack(:,k),parameters.int_stepsize,parameters.int_nsteps,'trajectory');
    
    echo(:,k)=parameters.coil_prob'*rho_stack_echo/norm(parameters.coil_prob,2);
    
    echo(:,k)=echo(:,k).*transpose(exp(1i*parameters.int_time_axis*parameters.carrier_frq*2*pi));
    
    if isfield(parameters,'bandwidth') && isfield(parameters,'zero_fill') && parameters.bandwidth ~=Inf
        [~,tmp] = bandpass_filter(parameters.int_time_axis,echo(:,k),parameters.bandwidth,parameters.zero_fill);
        echo(:,k)=tmp;
    end
end

% Echo integration
if parameters.transient==false     % integrate over echo
    result=sum(echo,1);
elseif parameters.transient==true  % return the echo
    result=echo;
end


end

