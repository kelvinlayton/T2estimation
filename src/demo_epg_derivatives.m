% Demonstration of the recursive algorithm the computation of EPG signal 
% derivatives
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
%

T2 = 80e-3;
T1 = 1;
R1 = 1./T1; R2 = 1./T2;

% Sequence parmaeters
%
nEchoes = 32;
tau = 10e-3;

% Calculate the echo amplitudes as a function of flip angle
%
nFlips = 100;
flipAngles = linspace(60*pi/180,pi,nFlips)';

s = zeros(nEchoes,nFlips);
for iFlip = 1:nFlips   
    s(:,iFlip) = epg(nEchoes,tau,R1,R2,flipAngles(iFlip));
end

% Calculate derivate of amplitude with respect to flip angle and calculate
% tangent curve at a point using the epg derivative function
%
[ds_dT2 ds_dFlip] = epg_derivatives(nEchoes,tau,R1,R2,flipAngles(30));
tangent = ds_dFlip(2)*flipAngles + s(2,30) - ds_dFlip(2)*flipAngles(30) ;

% Plot the first two echo amplitudes and the tanget curve as functions of 
% flip angle
%
l(1:2)=plot(180/pi*flipAngles,s(1:2,:)); hold on;
l(3)=plot(180/pi*flipAngles,tangent,'r--');
l(4)=plot(180/pi*flipAngles(30),s(2,30),'ro');
set(l,'LineWidth',2);

legend('First echo','Second echo','Linear approx at \alpha=95 degrees')
ylabel('Amplitude');
xlabel('Flip angle, \alpha, (degrees)');
xlim([60,180]);

%print -dpng -r300 'example_epg_derivatives.png'