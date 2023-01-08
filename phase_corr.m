function [y,phi]=phase_corr(y)

%this function does a phase correction of the signal y by minimizing the
%the imaginary part


[phi,~]=fminsearch(@rms_phi,0,[],y);


if sum(real(y*exp(1i*phi)))<0
    phi=phi+pi; 
end

y=y*exp(1i*phi);

end

function rms=rms_phi(phi,tr)
% r.m.s. of imaginary part after phase correction
itr=imag(tr*exp(1i*phi));
rms=sqrt(sum(itr.*itr));
end



