# JSM_2025_StructuredReplicatorDynamics
Code for the analysis and visualization of replicator dynamics with explicit relational structure in continuous time for both continuous and discrete player spaces. As it exists now, this code can be used to investigate the dynamics of replicator dynamics where each individual uses a mixed strategy and local information is pereseved throughout.  

The ODE and Nonlocal systems are described in a manuscript that will be linked here once the preprint is up: 

In <code>StructuredReplicator.m</code> and <code>StructuredReplicatorTer.m</code> solutions to the problem in a discrete player space are solved and visualized. <code>StructuredReplicatorTer.m</code> does the visualizaion through ternery plots while <code>StructuredReplicator.m</code> visualizes the time axis explicitly.
<code>continuousspacereplicator1D.m</code> Solves the nonlocal IVP in one spatial dimension and visualized it as a .fig animation. 

<code>connectionstrength.m</code> and <code>edgesensitivity.m</code> are used to investigate how edge removal, addition, or weight perturbation can change the stability of different strategy profiles.  The former can be used to understand the connection between two modules but it is made obsolete by the latter which can be used for general adjacency matricies

<code>potentialfunction.m</code> is a short script meant to visualize the potential function of the two player coordination or anticoordination game from game theoretic principles although it can be easily computed analytically and plotted far easier. 
