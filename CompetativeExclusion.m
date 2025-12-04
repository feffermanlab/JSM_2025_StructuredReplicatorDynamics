%% This script shows the evolution of strategy profiles
% given an initial strategy profile, adjacency matrix, and payoff matrix 
% in the manner described the structured replicator dynamics from the 
% manuscript associated to this repository.
% It has been made obsolete because CompetativeExclusionTer.m does it much
% better

%Initial condition 
u0=[0.75,0.75,0.75,0.75;
    0.24,0.01,0.01,0.24;
    0.01,0.24,0.24,0.01]

%number of players
n=4
%number of strategies
m=3
%adjacency matrix (nXn matrix)
W=[0,1,1,1;
   1,0,1,1;
   1,1,0,1;
   1,1,1,0]

%payoff matrix (mXm matrix)
%A=eye(m) %.................................Coordination
A=[0,1,-1;-1,0,1;1,-1,0]; %.................RPS
%A=[-1,-3,-3;-3,-1,-3;-3,-3,-1]; %..........Shifted Coordination

%Start and end times of the simulation
tspan = [0,20]

%Solves the ODE system described in the manuscript
[t,u] = ode45(@(t,u)f(t,u,n,m,W,A),tspan,u0);

%Creates a figure showing each player's strategy profile as it evolves in
%time. It is written to show four players currently but can be easily
%modified to show any number of players. 
tiledlayout(2,2)
nexttile
h=plot(t,u(:,1), t,u(:,2), t,u(:,3));
set(h(1),'Color','#0000a4')
set(h(2),'Color','#bc272d')
set(h(3),'Color','#e9c716')
title('Vertex 1')
ylim([0,1])
xlabel('time')
ylabel('Strategic Proportion')

nexttile
h=plot(t,u(:,4), t,u(:,5), t,u(:,6));
set(h(1),'Color','#0000a4')
set(h(2),'Color','#bc272d')
set(h(3),'Color','#e9c716')
title('Vertex 2')
ylim([0,1])
xlabel('time')
ylabel('Strategic Proportion')

nexttile
h=plot(t,u(:,7), t,u(:,8), t,u(:,9));
set(h(1),'Color','#0000a4')
set(h(2),'Color','#bc272d')
set(h(3),'Color','#e9c716')
title('Vertex 3')
ylim([0,1])
xlabel('time')
ylabel('Strategic Proportion')

nexttile
h=plot(t,u(:,10), t,u(:,11), t,u(:,12));
set(h(1),'Color','#0000a4')
set(h(2),'Color','#bc272d')
set(h(3),'Color','#e9c716')
title('Vertex 4')
ylim([0,1])
xlabel('time')
ylabel('Strategic Proportion')




%ODE fucntion
function dudt =f(t,u,n,m,W,A)
%% f this function computes the time derivative of u 
% t:        time (used only by ODE solver)
% u:        present value of u (nm dimensional vector)
% n:        number of players
% m:        number of strategies
% W:        Adjacency matrix (nXn matrix)
% A:        Payoff matrix (mXm matrix)

%%This function computes the time derivative of u as a function of u
    if(n*m~=length(u))
        error('State is not the right shape')
    end
    %by default assume the payoff matrix is the coordination game
    if(nargin<6)
        A=eye(m);
    end
    
    %reshape vector into an mXn matrix 
    du=zeros(m,n);
    umat =reshape(u,m,[]);
    w=0;
    for v = 1:n
        %uv is player v's strategies
        uv = umat(:,v);
        %gv is the cumulative strategy in player v's neighborhood
        temp = W(v,:).*umat;
        ugv = sum(temp,2);
        
        %ODE definied in manuscript
        du(:,v)=uv.*(A*ugv-(uv'*A*ugv)*ones(m,1));
    end
    dudt=du(:);
end


    
