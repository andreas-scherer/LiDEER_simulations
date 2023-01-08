function rho=spin_polarized_one_spin(Px,Py,Pz,beta,gamma,R, varargin)


% R=transformiert zeeman in spherical tensor operator
if length(varargin)==2
    H=varargin{1}+orientation(varargin{2}, [0, beta, gamma]); 
elseif length(varargin)==1
    H=varargin{1};
end

H=(H+H')/2;
   

rho=[Px,0,0;...
    0,Py,0;...
    0,0,Pz];


V1=wigner(1,0,beta,gamma)*[-1/sqrt(2), 1i/sqrt(2), 0; 0,0,1; 1/sqrt(2), 1i/sqrt(2),0];

[V2,~]=eigs(H);

%transform to high field basis
rho=V1*rho*V1';
%transform to H basis and remove coherences
rho=V2'*rho*V2;
rho=diag(diag(rho));

%transform to high field basis
rho=V2*rho*V2';

rho=rho(:);
rho=R*rho;

end