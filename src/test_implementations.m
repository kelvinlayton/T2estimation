% TEST_IMPLEMENTATIONS Test different implementations of EPG algorithms
%    TEST_IMPLEMENTATIONS  This script tests two different implementations
%    of the EPG algorithm and the EPG derivative algorithm. The first
%    implementation of both algorithms is implemented in MATLAB and uses 
%    sparse matrices. The second implementation is a MEX C++ file and uses
%    the inherent structure to efficiently execute operations.

%% EPG Algorithm

% Simulation setup
%
tau=10e-3;
T1=1;
T2=[100e-3 20e-3]';

alpha=160*pi/180;
n=32;
R2=1./T2;
R1=1/T1;
nRuns = 10;

% Run algorithms and save execution time
%
tic
for i=1:nRuns
    aMat=epg(n,tau,R1,R2,alpha);
end
timeMat=toc;
tic
for i=1:nRuns
    aMex=epg_mex(n,tau,T1,R2,alpha);
end
timeMex=toc;

% Calculate difference
%
diff=norm(aMat-aMex);
strResults = {'FAILED','PASSED'};

fprintf('EPG MATLAB, %d runs: %f seconds\n',nRuns,timeMat);
fprintf('EPG MEX,    %d runs: %f seconds\n',nRuns,timeMex);
fprintf('Numerical difference: %g (%s)\n\n',diff, strResults{(diff<1e-10)+1});

%% EPG derivative algorithm

% Simulation setup
%
tau=10e-3;
T1=1000e-3;
T2=[100e-3 20e-3]';

alpha=160*pi/180;
n=32;
R2=1./T2;
R1=1/T1;
nRuns = 10;

% Run algorithms and save execution time
%
tic
for i=1:nRuns
    [dMat,dMat2]=epg_derivatives(n,tau,R1,R2,alpha);
end
timeMat=toc;
tic
for i=1:nRuns
    [dMex,dMex2]=epg_derivatives_mex(n,tau,R1,R2,alpha);
end
timeMex=toc;

% Calculate difference
%
diff=norm([dMat,dMat2]-[dMex dMex2]);
strResults = {'FAILED','PASSED'};

fprintf('EPG derivative MATLAB, %d runs: %f seconds\n',nRuns,timeMat);
fprintf('EPG derivative MEX,    %d runs: %f seconds\n',nRuns,timeMex);
fprintf('Numerical difference: %g (%s)\n\n',diff, strResults{(diff<1e-10)+1});
