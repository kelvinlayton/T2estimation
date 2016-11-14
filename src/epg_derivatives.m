% EPG_DERIVATIVES Derivatives of the Extended Phase Graph (EPG) signal.
%   [Jr Ja] = EPG_DERIVATIVES(N,TAU,R1,R2,ALPHA) calculate the derivatives
%   of the EPG signal with respect to rate (Jr) and flip angle (Ja). The
%   signal is from a CPMG sequence with N echoes and echo spacing of 
%   (2*TAU). The parameters are R1=1/T1 (scalar) and R2=1./T2 (scalar or 
%   vector) and ALPHA the flip angle. The matrices Jr and Ja contains one 
%   column for each element in R2.
%
% The following implementation uses sparse matrices for efficient
% compuation.
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
%
function [ Jr, Ja ] = epg_derivatives( n, tau, R1, R2vec, alpha )

nRates = length(R2vec);
tau=tau/2;

Jr = zeros(n,nRates);
Ja = zeros(n,nRates);

% RF mixing matrix, T
T0 = [cos(alpha/2)^2, sin(alpha/2)^2, sin(alpha); ...
    sin(alpha/2)^2, cos(alpha/2)^2, -sin(alpha); ...
    -0.5*sin(alpha), 0.5*sin(alpha), cos(alpha)];

% Derivative of T with respect to alpha
T0d = [-0.5*sin(alpha), 0.5*sin(alpha), cos(alpha); ...
    0.5*sin(alpha), -0.5*sin(alpha), -cos(alpha); ...
    -0.5*cos(alpha), 0.5*cos(alpha), -sin(alpha)];

% Build mixing matrix for all states
TArray = cell(1,n);
TArray(:) = {sparse(T0)};
T = blkdiag(1,TArray{:});
TArray(:) = {sparse(T0d)};
Td = blkdiag(0,TArray{:});
    
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
    S(3*o+1,3*o+1)=1;  %Z_order
end

% Loop over different relaxation rates
%
for iRate=1:nRates

    % Relaxation matrix
    R2 = R2vec(iRate);
    R0 = diag([exp(-tau*R2),exp(-tau*R2),exp(-tau*R1)]);
    R0d = diag(-tau*exp(-tau*R2)*[1 1 0]);  % derivative w.r.t R2 

    % Build relaxation matrix (and derivates) for all states
    RArray = cell(1,n);
    RArray(:) = {sparse(R0)};
    R = blkdiag(exp(-tau*R2),RArray{:});
    RArray(:) = {sparse(R0d)};
    Rd = blkdiag(-tau*exp(-tau*R2),RArray{:});

    % Precession and relaxation matrix (and derivative)
    P = (R*S);
    Pd = (Rd*S);

    % Matrix representing the inter-echo duration
    E = P*T*P;
    EdR0 = Pd*T*P + P*T*Pd;
    EdR=EdR0;
    EdA0 = P*Td*P;
    EdA=EdA0;

    % Vector of all states
    x = zeros(size(R,1),1);
    x(1)=1;
    
    % Deriviatives of state vector
    xdashR = EdR*x;
    xdashA = EdA*x;
    
    % First echo
    x = E*x;
    Jr(1,iRate) = xdashR(1);
    Ja(1,iRate) = xdashA(1);

    % Subsequent echoes
    for i=2:n
        % Calculate derivatives using the product rule
        xdashR = EdR0*x + E*xdashR;
        xdashA = EdA0*x + E*xdashA;
        
        Jr(i,iRate) = xdashR(1);
        Ja(i,iRate) = xdashA(1);
        
        % Calculate new state for next interation
        x=E*x;

    end

end

