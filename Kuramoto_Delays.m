function Phases = Kuramoto_Delays(C,D,f,K,dt,tmax,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to run simulations of coupled systems using the
% KURAMOTO MODEL OF COUPLED OSCILLATORS WITH TIME DELAYS
%
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% ============= MODIFIED VERSION!!! ================
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%
% - Each node is represented by a Phase Oscillator 
% - with Natural Frequency f
% - with couplig factor K
% - coupled according to a  connectivity matrix C
% - with a delay matrix D
% - with time resolution dt, stopping at tmax
%
% Optional argument: random generator seed (for reproducibility)
%  
% All units in Seconds, Meters, Hz
%
% Simulation code written by Joana Cabral, 
% January 2019, joana.cabral@psych.ox.ac.uk
% 
% Modified by Matteo Mancini,
% April 2021, ingmatteomancini@gmail.com
%
%%%%%%%%%%%%%%%%%%%%

if ~isempty(varargin)
    rng(varargin{1});
end

% Normalize parameters
N=size(C,1);   % Number of units
Omegas=2*pi*f*ones(N,1)*dt; % Frequency of all units in radians/second
kC=K*C*dt; % Scale matrix C with 'K' and 'dt' to avoid doing it at each step

Delays=floor(D/dt);
Delays(C==0)=0;

Max_History=max(Delays(:))+1;
% Revert the Delays matrix such that it contains the index of the History 
% that we need to retrieve at each dt
Delays=Max_History-Delays;


%% Initialization

% History of Phases is needed at dt resolution for as long as the longest 
% delay. The system is initialized all desinchronized and uncoupled
% oscillating at their natural frequencies

Phases_History = 2*pi*rand(N,1)+Omegas.*ones(N,Max_History).*(1:Max_History);
Phases_History = mod(Phases_History,2*pi); % all values between 0 and 2pi
    
% This History matrix will be continuously updated (sliding history).
% figure; plot(sin(Phases_History'))   % To see the History

% Simulated activity will be saved only at dt_save resolution
Phases=zeros(N,tmax/dt); 
sumz=zeros(N,1);

%% Run simulations
disp(['Now running for K=' num2str(K)])

tic % to count the time of each run
nt=0;

for t=dt:dt:tmax  % We only start saving after t_prev
    
    Phase_Now=Phases_History(:,end); % The last collumn of Z is 'now'
    
    % Input from coupled units
    for n=1:N
        sumzn=0; % Intitalize total coupling received into node n
        for p=1:N
            if kC(n,p) % Calculate only input from connected nodes (faster)             
                sumzn = sumzn + kC(n,p) * sin(Phases_History(p,Delays(n,p))-Phase_Now(n));
            end  
        end
        sumz(n)=sumzn;
    end
    
    Phases_History(:,1:end-1)=Phases_History(:,2:end);
    
    % Update the end of History with evolution in the last time step           
    Phases_History(:,end)= Phase_Now + Omegas + sumz;
    
    nt=nt+1;
    Phases(:,nt)=Phases_History(:,end);
end

disp(['Finished, lasted ' num2str(round(toc)) ' secs for real ' num2str(tmax) ' seconds']);
disp(['Simu speed ' num2str(round(toc/tmax)) ' Realseconds/SimuSecond']);

