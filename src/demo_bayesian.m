% Demonstration of the Bayesian estimation algorithm for multi-component 
% T2 estimation.
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
%

% Generate synthetic data
%
w = [0.3 0.7]';
T2 = [20e-3 100e-3]';
T1 = inf;
R1 = 1./T1; R2 = 1./T2;
nComponents = length(w);
distWidth = 7e-3;
phi = pi/8;
alpha = 160*pi/180;

% Sequence parmaeters
%
nEchoes = 32;
tau = 10e-3;
SNR = 100; 
sigma2 = (1./SNR)^2;
opts.tau = tau;
opts.sigma2 = sigma2;

% Calculate signal from a Gaussian mixture
%
T2grid = linspace(0e-3,200e-3,1000)';
dist = w(1)*normpdf(T2grid,T2(1),distWidth) + ...
      w(2)*normpdf(T2grid,T2(2),distWidth);
dist=dist./sum(dist);

s = epg_mex(nEchoes,tau,R1,1./T2grid,alpha)*dist*exp(1j*phi);

% Define prior (in terms of R1 and R2)
%
opts.priorMean = [0.5; 0.5; 60; 7; pi/2; 0];
opts.priorCov = diag([10^2*ones(nComponents,1); 1000^2*ones(nComponents,1); 100^2*ones(2,1)]);

%nComponents=3;
%opts.priorMean = [0.3; 0.3; 0.3; 60; 30; 5; pi/2; 0];
%opts.priorCov = diag([100*ones(nComponents,1); 1000^2*ones(nComponents,1); 50^2*ones(2,1)]);


% Define correction schedule
%
gammas = logspace(-6,1,30);
gammas = gammas ./ sum(gammas);
opts.correctionSchedule = gammas;

% Allocate output arrays
%
nTrials=50;
T2Hat = zeros(nComponents,nTrials);
wHat = zeros(nComponents,nTrials);

% Repeat for different noise realisations
%
for iTrial=1:nTrials
    
    % Add noise to signal
    %
    y = s + sqrt(sigma2)*randn(size(s));

    % Estimate distribution components
    %
    [T2Hat(:,iTrial), wHat(:,iTrial)] = bayesian_estimate(y,opts);
    
end

% Plot estimates
%
plotScale = 0.7./max(dist);
hold on;
plot(1e3*T2grid,plotScale*dist,'k','LineWidth',2);
plot(1e3*T2Hat(1,:),wHat(1,:),'x');
plot(1e3*T2Hat(2,:),wHat(2,:),'rx');
%plot(1e3*T2Hat(3,:),wHat(3,:),'gx');
stem(1e3*mean(T2Hat'),mean(wHat'),'Color',0.5*[1 1 1],'LineWidth',2);

xlabel('T2 (ms)');
ylabel('weight, w');
legend('True Distribution','Estimates Fast Mode',...
    'Estimates Slow Mode','Mean');

%print -dpng -r300 'example_bayesian.png'