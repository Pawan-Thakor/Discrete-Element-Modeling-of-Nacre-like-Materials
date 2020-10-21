% MATLAB codes for Finite Element Analysis 
clear all 
clear
format long g

%Structure p
p=struct();

%% Defining Problem

% Element_Nodes: connections at elements 
p.Element_Nodes=[1 2;2 3];
% Num_Nodes: number of nodes 
p.Num_Nodes=3;
% Element_Stiffness: spring stiffness values
p.Element_Stiffness=[2e6 2e6 2e6];
% Num_Elements: number of Elements 
p.Num_Elements=size(p.Element_Nodes,1); 

%% Initializing Matrix 

% Displacement: Displacement vector
p.Displacement=zeros(p.Num_Nodes,1);
% Force: Force vector 
p.Force=zeros(p.Num_Nodes,1);
% Stiffness: Stiffness Matrix 
p.Stiffness=zeros(p.Num_Nodes);

%% Given Boundary Conditions

% %  Boundary condition for Homogeneous Solution
% % Prescribed_Dof: Value of displacement is given 
 p.Prescribed_Dof=[1;3];
 p.Force(2)=100;


%  Boundary condition for Non-Homogeneous Solution
% Prescribed_Dof: Value of displacement is given 
% p.Prescribed_Dof=[1;3];
% Displacement: Applied at node 3 
% p.Displacement(3)=0.06;

%% Solving Spring Problem 

% Calculating Global Stiffness Matrix
p.Stiffness=Global_Stiff_Matrix(p);
% Calculating Displacement
p.Displacement=Displacement(p,p.Stiffness,p.Displacement,p.Force);
% Call to Function to Display output
p.Force = p.Stiffness*p.Displacement;
Output(p);


