% EPG Matrix implementation of the Extended Phase Graph (EPG) algortihm.
%   H = EPG(N,TAU,R1,R2,ALPHA) calculate the echo amplitudes for a CPMG
%   sequence with N echoes and echo spacing of (2*TAU). The parameters are
%   R1=1/T1 (scalar) and R2=1./T2 (scalar or vector) and ALPHA the flip 
%   angle. The matrix H contains one column for each element in R2
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
%
function [ H ] = epg( n, tau, R1, R2vec, alpha )


nRates = length(R2vec);
tau=tau/2;

H = zeros(n,nRates);

% RF mixing matrix
%
T0 = [cos(alpha/2)^2, sin(alpha/2)^2, sin(alpha); ...
    sin(alpha/2)^2, cos(alpha/2)^2, -sin(alpha); ...
    -0.5*sin(alpha), 0.5*sin(alpha), cos(alpha)];

TArray = cell(1,n);
TArray(:) = {sparse(T0)};
T = blkdiag(1,TArray{:});

% Selection matrix to move all traverse states up one coherence level
%
S = sparse(3*n+1,3*n+1);
S(1,3)=1;
S(2,1)=1;
S(3,6)=1;
S(4,4)=1;
for o=2:n
    offset1=( (o-1) - 1)*3 + 2;
    offset2=( (o+1) - 1)*3 + 3;
    if offset1<=(3*n+1)
    S(3*o-1,offset1)=1;  % F_k <- F_{k-1}
    end
    if offset2<=(3*n+1)
    S(3*o,offset2)=1;  % F_-k <- F_{-k-1}
    end
    S(3*o+1,3*o+1)=1;  % Z_order
end
    
for iRate=1:nRates

    % Relaxation matrix
    R2=R2vec(iRate);
    R0 = diag([exp(-tau*R2),exp(-tau*R2),exp(-tau*R1)]);

    RArray = cell(1,n);
    RArray(:) = {sparse(R0)};
    R = blkdiag(exp(-tau*R2),RArray{:});

    % Precession and relaxation matrix
    P = (R*S);

    % Matrix representing the inter-echo duration
    E = (P*T*P);
    
    % Recursively apply matrix to get echo amplitudes
    %
    x = zeros(size(R,1),1);
    x(1)=1;
    for iEcho=1:n
        x=E*x;
        H(iEcho,iRate) = x(1);
    end

end

end

