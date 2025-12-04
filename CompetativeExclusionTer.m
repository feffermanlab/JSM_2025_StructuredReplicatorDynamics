%% This script shows, in a tenary plot, the evolution of strategy profiles
% given an initial strategy profile, adjacency matrix, and payoff matrix 
% in the manner described the structured replicator dynamics from the 
% manuscript associated to this repository.
% 
% This script requires the ternary plots repository: https://github.com/lynch730/ternary_plots
% Written by lynch730 based off of Ulrich Theune's original package:
% https://www.mathworks.com/matlabcentral/fileexchange/7210-ternary-plots
% 
% Fork from that repository, read the documentation, 
% and add ternary_plots to the current path using
add_ternary_paths


%Initial strategy profile
u0=[0.2,0.98,0.01,0.2;
    0.6,0.01,0.98,0.3;
    0.2,0.01,0.01,0.5]
%Number of players given by the number of columns in u0
n=length(u0(1,:));
%Number of strategies given by the number of rows in u0
m=length(u0(:,1));

%Adjacency matrix must be an nXn sysmmetric matrix
W=[0,1,0,1;
   1,0,1,0;
   0,1,0,1;
   1,0,1,0]

%Payoff matrix must be an mxm matrix
A=-eye(m) %................................anti-coordination
%A=eye(m) %................................coordination
%A=[0,1,-2;-2,0,1;1,-2,0]; %...............Skew RPS
%A=[0,1,-1,0;0,0,1,-1;-1,0,0,1;1,-1,0,0] %.4x4 RPS
%A=[-1,-3,-3;-3,-1,-3;-3,-3,-1]; %.........Shifted coordination
%A=[-1,-2,0; 0,-1,-2; -4,0,-1]; %..........Other

%Start and End times for solution
tspan = [0,20]

%Solve the IVP with u0 initial data on tpsan
%At a time slice u(t,:) is an nxm dimensional vector 
%l is the total number of time slices in the solution
[time,u] = ode45(@(t,u)f(t,u,n,m,W,A),tspan,u0);
l = length(u(:,1));


%For potential games, the total payoff is an inverse energy (nondecreasing)
%for the ODE system. This routine computes the energy over the solution
energy = zeros(length(time),1);
for t = 1:length(time)
    %get the solution at a time slice and reshape it into a mXn matrix 
    ut = u(t,:);
    umat =reshape(ut,m,[]);
    w=0;
    for v = 1:n
        %uv is the strategy profile of player v
        uv = umat(:,v);
        %ugv is the cumulative strategies in the neighborhood of v
        temp = W(v,:).*umat;
        ugv = sum(temp,2);
        
        %compute player v's payoff with the bilinear form
        wv=uv'*A*ugv;
        %add to total payoff
        w= w+wv;
    end
    energy(t)=w;
end

%disp(energy)


%The following code generates the ternery plots. It is currently written
%for 4 players and an energy plot at the bottom. This can be changed to
%accomidate as many players as is necessary and to include or remove the
%energy plot easily

%labels for the strategies
vgen = {'titlelabels',{'Rock','Paper','Scissors'}}; %........RPS
%vgen  = {'titlelabels', {'A','B','C'}}; %...................Arbitrary


t=tiledlayout(3,2);
t.TileSpacing = "compact";
t.Padding = "compact";
nexttile

handle_base  = ternary_axes(vgen); 
% Create line plot Plot
var = {'linewidth',2,'color','b'};
dataplots(1).obj = ternary_plot3(ternary_axes_limits(), 'l', u(:,1), 'b', u(:,2),[],var{:});
var = {'MarkerFaceColor','g','MarkerEdgeColor','g'};
dataplots(1).obj = ternary_scatter3(ternary_axes_limits(), 'l', u(1,1), 'b', u(1,2),[],'none',var{:});
var = {'MarkerFaceColor','r','MarkerEdgeColor','r'};
dataplots(1).obj = ternary_scatter3(ternary_axes_limits(), 'l', u(l,1), 'b', u(l,2),[],'none',var{:});    
% Add Title
title('Vertex 1          ','FontSize',14,HorizontalAlignment='right')

nexttile

handle_base  = ternary_axes(vgen); 
% Create line plot Plot
var = {'linewidth',2,'color','b'};
dataplots(1).obj = ternary_plot3(ternary_axes_limits(), 'l', u(:,4), 'b', u(:,5),[],var{:});
var = {'MarkerFaceColor','g','MarkerEdgeColor','g'};
dataplots(1).obj = ternary_scatter3(ternary_axes_limits(), 'l', u(1,4), 'b', u(1,5),[],'none',var{:});
var = {'MarkerFaceColor','r','MarkerEdgeColor','r'};
dataplots(1).obj = ternary_scatter3(ternary_axes_limits(), 'l', u(l,4), 'b', u(l,5),[],'none',var{:});    
% Add Title
title('Vertex 2          ','FontSize',14,HorizontalAlignment='right')

nexttile

handle_base  = ternary_axes(vgen); 
% Create line plot Plot
var = {'linewidth',2,'color','b'};
dataplots(1).obj = ternary_plot3(ternary_axes_limits(), 'l', u(:,7), 'b', u(:,8),[],var{:});
var = {'MarkerFaceColor','g','MarkerEdgeColor','g'};
dataplots(1).obj = ternary_scatter3(ternary_axes_limits(), 'l', u(1,7), 'b', u(1,8),[],'none',var{:});
var = {'MarkerFaceColor','r','MarkerEdgeColor','r'};
dataplots(1).obj = ternary_scatter3(ternary_axes_limits(), 'l', u(l,7), 'b', u(l,8),[],'none',var{:});    
% Add Title
title('Vertex 3          ','FontSize',14,HorizontalAlignment='right')

nexttile

handle_base  = ternary_axes(vgen); 
% Create line plot Plot
var = {'linewidth',2,'color','b'};
dataplots(1).obj = ternary_plot3(ternary_axes_limits(), 'l', u(:,10), 'b', u(:,11),[],var{:});
var = {'MarkerFaceColor','g','MarkerEdgeColor','g'};
dataplots(1).obj = ternary_scatter3(ternary_axes_limits(), 'l', u(1,10), 'b', u(1,11),[],'none',var{:});
var = {'MarkerFaceColor','r','MarkerEdgeColor','r'};
dataplots(1).obj = ternary_scatter3(ternary_axes_limits(), 'l', u(l,10), 'b', u(l,11),[],'none',var{:});    
% Add Title
title('Vertex 4          ','FontSize',14,HorizontalAlignment='right')

nexttile([1,2])
plot(time,energy)


%ODE fucntion
function dudt =f(t,u,n,m,W,A)
%%f this function computes the time derivative of u 
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
