% Radioactive Decay of Radon by using Monte-Carlo-Simulation
% Version 1.0
% Date: 14.12.2017


clc;
clear;

N01       = input('Enter total number of nuclides: ');
N02       = 0; %starting number of nuclides of radon
dt        = input('Enter time steps in days: '); %time step [days] 

lambda1   = 1.187*10^(-6); %decay constant radium over days
lambda2   = 0.182; %decay constant radon over days

p1        = lambda1*dt; %probability per time step of radium decay
p2        = lambda2*dt; %probability per time step of radon decay

n_monteCarloRadium  = zeros(10000,1); 
n_monteCarloRadon   = zeros(10000,1); 


tmax=10000*dt;
t=(0:dt:tmax)';

n_monteCarloRadium(1) = N01;
n_monteCarloRadon(1)  = N02;
N1 = n_monteCarloRadium(1);
N2 = n_monteCarloRadon(1);

t_start = cputime;
% Monte-Carlo-Simulation Radium
for i=2:10001
 
  n1 = N1;
  n2 = N2;
%  printf(num2str(n2))
  for j=0:n1-1 
    r1 = rand();
    if(r1<p1)  % compare r with probability density function, here p = const.: is there a decay?
      N1 = N1-1;
      N2 = N2+1;
%        printf(num2str(N2))
    end
  end
  
  n_monteCarloRadium(i) = N1;
  
  for j=0:n2-1;
    r2 = rand();
    if(r2<p2) % compare r with probability density function, here p = const.: is there a decay?
      N2 = N2-1;
 %     printf('%i',N2)
    end
  end
  
  n_monteCarloRadon(i) = N2;
  
end
t_end = cputime - t_start;
% printf('\nSimulation time (min):\t %4d\n',t_end/60);
fprintf('\nSimulation time (min):\t %4d\n',t_end/60);

% numerical solution by solving ODEs
[t,n_numeric]   = ode45(@RadioactiveDecayMCRadonmodel_MATLAB,t,[N01,N02]); %%starting conditions for ra216,rn222
% analytical solution (s. also: http://www.physik.uni-bielefeld.de/~borghini/Teaching/Kernphysik/11_24.pdf, p. 62)
n_analytic  = N02*exp(-lambda2*t) + lambda1/(lambda2-lambda1)*N01*(exp(-lambda1*t) - exp(-lambda2*t));

% plotting
figure(1)
title(strcat ('Parameters used:  <N01>    <dt> = <',num2str(N01),'>   <',num2str(dt),'>'));
plot(t,n_analytic,'r',t,n_numeric(:,2),'b', t,n_monteCarloRadon,'k');
xlabel('time');
ylabel('number of nuclides');
legend('analytical solution','numerical solution','Monte Carlo Simulation');

figure(2)
title(strcat ('Parameters used:  <N01>    <dt> = <',num2str(N01),'>   <',num2str(dt),'>'));
subplot(3,1,1); plot(t,n_analytic,'r');
xlabel('time');
ylabel('number of nuclides');
legend('analytical solution');

subplot(3,1,2); plot(t,n_numeric(:,2),'b');
xlabel('time');
ylabel('number of nuclides');
legend('numerical solution');

subplot(3,1,3); plot(t,n_monteCarloRadon,'k');
xlabel('time');
ylabel('number of nuclides');
legend('Monte Carlo Simulation');
