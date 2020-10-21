% MATLAB codes for Finite Element Analysis of Truss Member + Spring Member 
clear;
clear all;

%Structure p
p=struct();

%% Defining Problem

Nodes_1= xlsread('Input_file', 'B:C');
Nodes_angle_1= xlsread('Input_file', 'Y:Z');
Element_Nodes_1= xlsread('Input_file', 'E:F');
Spring_1 = xlsread('Input_file', 'L:O');
Element_E_1= xlsread('Input_file', 'H:H');
Element_Area_1= xlsread('Input_file', 'J:J');
Prescribed_Dof_1= xlsread('Input_file', 'Q:S');
Force_1= xlsread('Input_file', 'U:W');

%% Substituting Values in Variables

if Spring_1(1,1)==0
    % % Values for Spring Element
    % Spring_Element: Connection of Spring Element
    p.Spring_Element = [];
    % Spring_Stiffness: Value of Stiffness of Spring
    p.Spring_Stiffness = [];
    % Spring_angle: Value of angle of Spring
    p.Spring_angle = [];
else
    % % Values for Spring Element
    % Spring_Element: Connection of Spring Element
    p.Spring_Element = Spring_1(:,1:2);
    % Spring_Stiffness: Value of Stiffness of Spring
    p.Spring_Stiffness = Spring_1(:,3);
    % Spring_angle: Value of angle of Spring
    p.Spring_angle = Spring_1(:,4);
end
    
% % Values for each Elements
% Element_Nodes: connections at elements 
p.Element_Nodes = cat(1,Element_Nodes_1,p.Spring_Element);
% Num_Elements: number of Elements 
p.Num_Elements = size(Element_Nodes_1,1) + size(p.Spring_Element,1);
% Element_Stiffness: Stiffness values of Elements
p.Element_E = Element_E_1;
% Element_Area: Area values of Elements
p.Element_Area = Element_Area_1;

% % Values for each Nodes 
% Nodes_coordinate: position of nodes 
p.Node_coordinate = Nodes_1;
% Num_Nodes: number of nodes 
p.Num_Nodes = size(Nodes_1,1) + size(p.Spring_Element,1);
% Node Angle
p.Node_angle = Nodes_angle_1(:,2);

% BCs: Boundary Conditions
p.BCs=Prescribed_Dof_1;

% Load: Load Conditions
p.Load=Force_1;

% NDof: Node degree of freedom
p.NDof = 2;

% GDof: total number of degrees of freedom
p.GDof=p.NDof*p.Num_Nodes;

%% Initializing Matrix 

% Displacement: Displacement vector
p.Displacement=zeros(p.GDof,1);
% Force: Force vector 
p.Force=zeros(p.GDof,1);
% Stiffness: Stiffness Matrix 
p.Stiffness=zeros(p.GDof);
% Stress: Stress Matrix 
p.Stress=zeros((p.Num_Elements-size(p.Spring_Element,1)),1);

% Given Boundary Conditions

% % Boundary condition for Truss Problem
% % Prescribed_Dof: Value of displacement is given 
 p.Prescribed_Dof=[];
 for ii = 1:size(Prescribed_Dof_1,1)
     thisdof = p.NDof*(Prescribed_Dof_1(ii,1)-1) + Prescribed_Dof_1(ii,2);
     p.Prescribed_Dof = [p.Prescribed_Dof thisdof]; 
     p.Displacement(thisdof) = Prescribed_Dof_1(ii,3); 
 end
 
% % Force: Application of Force
 for ii = 1:size(Force_1,1)
     p.Force(p.NDof*(Force_1(ii,1)-1) + Force_1(ii,2)) = Force_1(ii,3); 
 end
 
% This Transfromation Matrix is used for Inclined Boundary Condition (w.r.t Global x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining Transformation Matrix at Nodes 
T = zeros(p.GDof);       % Transformation Matrix for Nodes
for n=1:p.Num_Nodes
    a=p.Node_angle(n);
    temp=[cosd(a) sind(a); -sind(a) cosd(a)];
    T(p.NDof*(n-1)+1:p.NDof*(n-1)+p.NDof, p.NDof*(n-1)+1:p.NDof*(n-1)+p.NDof) = T(p.NDof*(n-1)+1:p.NDof*(n-1)+p.NDof, p.NDof*(n-1)+1:p.NDof*(n-1)+p.NDof) + temp;
end

p.Displacement = T*p.Displacement;
p.Force = T*p.Force;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Solving Truss Problem 

% Calculating Global Stiffness Matrix
p.Stiffness = Truss_Global_Stiff_Matrix(p);
p.Stiffness = T*p.Stiffness*T';
% Calculating Displacement
p.Displacement = Displacement(p,p.Stiffness,p.Displacement,p.Force);
% Calculating Reaction Forces
p.Force = p.Stiffness*p.Displacement;
% Calculating Stress Matrix
p.Stress = Stress_2D_Truss(p);
% Call to Function to Display output
Output(p);

K = p.Stiffness;
% Plotting the Diplacement plot
if size(p.Spring_Element,1)==0
Plot_Displacement(p);
end
