%% This script computes the potential function for the only two games in which is
% is possible to visualize (2player coordination and and anti-coordination). It is
% possible to compute this analytically in these cases so really this piece
% of code is worthless but it confirms to us that the potential function
% works the way we expect it to. 

%Adjacency matrix (This is the only one that works)
W=[0,1;1,0];

% Payoff matrix 
A=[1,0;0,1]; %...................Coordination
%A=[0,1;1,0]; %...................Anti-Coordination

%Spatial mesh for generating surface plot
N = 80;
h = 1/(N-1);
Omegap1 = h.*(0:N-1);
Omegap2 = h.*(0:N-1);
[OmegaP1,OmegaP2]=meshgrid(Omegap1,Omegap2);


%function to compute surface
potential = @(x1,x2)TW(x1,x2,A);

%compute surface
Wsurface = potential(OmegaP1,OmegaP2);

%Plot 3d figure
figure
mesh(OmegaP1,OmegaP2,Wsurface)
title('Two player coordination potential')
xlabel('u_1^1');
ylabel('u_2^1');
zlabel('Total Fitness');
 


function totalfitness = TW(x1,x2,A)
%% TW: computing the potential at each point, for a coordination or 
% anticoordination game is equivalent to computing the total payoff
%
% x1: Strategy profile of player 1
% x2: Strategy profile of player 2
% A: Payoff matrix
%
% Each players pairwise payoff is determined by the bilienar form <x1,Ax2>
    totalfitness = 2.*(x1.*(A(1,1).*x2+A(1,2).*(1-x2))+(1-x1).*(A(2,1).*x2+A(2,2).*(1-x2)));
end
