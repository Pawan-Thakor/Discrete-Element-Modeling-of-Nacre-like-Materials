% MATLAB codes for Finite Element Analysis
% It consist of an Algorithm which consist Displacement Increment
clear all 
clear

format long g

%Structure p
p=struct();

%% Defining Problem

% Properties for Shear Stress-Strain Curve
% Properties for Shear Stress-Strain Curve
G=100*10^6;                    %Shear Modulus
H_1=6*10^6;                      %Hardening Value
tau_y=3.5*10^6;                %Yeild Shear Stress value
gamma_u=2;                     %Max Allowed Plastic Strain 
gamma_p_s=1.5;                 %Plastic Strain at start of Softening region

tau_s = tau_y + H_1*gamma_p_s;               %Max Shear Stress Point
G_1 = H_1*G/(G+H_1);                         %Slope of Hardening in Shear stress-strain curve 
G_2 = -tau_s/(gamma_u-(gamma_p_s+tau_s/G));  %Slope of Softening in Shear stress-strain curve
H_2 = G*G_2/(G-G_2);                         %Softening Value

% Parameters of the Tablets
rho=4;                         %Tablet Ratio
k=0.3;                         %Ovterlap Ratio
ti=1/9;                        %Interface Thickness
t=1;                           %Tablet Thickness
phi=t/(t+ti);                  %Tablet Volume concentration

%% Defining Coordinates 

% Node coordinates
p.Node_coordinate=[0 0; rho*k*t 0; rho*t 0 ];
% Num_Nodes: number of nodes 
p.Num_Nodes=size(p.Node_coordinate,1);
% Element_Nodes: connections at elements 
p.Element_Nodes=[1 2; 2 3];
% Num_Elements: number of Elements 
p.Num_Elements=size(p.Element_Nodes,1); 


Element_Length=zeros(p.Num_Elements,1);
for e=1:p.Num_Elements
    %Nodes of each Element 
    el_node = p.Element_Nodes(e, 1:2);
    %x-coordinate of Node
    node_xx = p.Node_coordinate(el_node);
    % Evaluation of Element Length
    Element_Length(e) = abs((node_xx(1,1)-node_xx(1,2)));          
end


% Element_Length: Length of each spring Element
p.Element_Length = Element_Length;
% Element_Area: Area of Element
p.Element_Area = ti*Element_Length;
% Element_Stiffness: truss stiffness values
p.Element_G=G*ones([1,p.Num_Elements]);
% NDof: Degree of Freedom at each node
p.NDof = 1;
% GDof: total number of degrees of freedom
p.GDof=p.Num_Nodes;


%% Given Boundary Conditions

% Boundary condition for Homogeneous Solution
% Prescribed_Dof: Value of displacement is given 
 p.Prescribed_Dof=[1; 3];
 
%% Initializing Matrix 

% Displacement: Displacement vector
p.Displacement=zeros(p.GDof,1);
% Displacment Increment:
Displacement_increment=zeros(p.GDof,1);
% Force: Force vector 
Force_increment=zeros(p.GDof,1);
% Stiffness: Stiffness Matrix 
p.Stiffness=zeros(p.GDof,p.GDof);

% Number of Displacements Steps
n=3000;

% u: Displacement for all Steps
u = zeros(p.Num_Nodes,n+1);
% Force: Force for all Steps
Force = zeros(p.Num_Nodes,n+1);
% Strain: for each Element for all Steps
Shear_Strain = zeros(p.Num_Elements,n+1);
% Stress: for each Element for all Steps
Shear_Stress = zeros(p.Num_Elements,n+1);

% Yeild Stress: at each Step
Shear_sigma_yeild = zeros(p.Num_Elements,1);
% Plastic Strain: for all Steps
Shear_Plastic_Strain = zeros(p.Num_Elements,n+1);
% Change in Plastic Strain at each Step
Shear_Delta_plastic_strain = zeros(p.Num_Elements,1);

% Active Dof: Nodes at which Force is supposed to be zero
% As Displacement Increment is considered 
p.active_Dof = setdiff([1:p.GDof]', [p.Prescribed_Dof]);

%% Solving Spring Problem 

% Displacement to Shear Strain Transformation Matrix
% (Displacement is converted to Shear Strain)
B = zeros(p.Num_Elements,p.Num_Nodes);       
for e=1:p.Num_Elements
    % Nodes of each Element 
    el_node = p.Element_Nodes(e, 1:2);
    %x-coordinate of Node
    node_xx = p.Node_coordinate(el_node);
    
    % Sign for Finding position of Node
    % (If Node 2 is ahead of Node 1 then Le is positive)
    % (Otherwise it is negative and accordingly B matrix would be made)
    Le =(node_xx(1,2)-node_xx(1,1));
    
    % Transformation term
    temp=[-sign(Le)/ti sign(Le)/ti];   
    
    B(e, el_node(1)) = B(e, el_node(1)) + temp(1);
    B(e, el_node(2)) = B(e, el_node(2)) + temp(2);
end


%Shear Stress to Force Transformation at each Node
% (Stress is converted to Force)
B_T = zeros(p.Num_Elements,p.Num_Nodes);       
for e=1:p.Num_Elements
    %Nodes of each Element
    el_node = p.Element_Nodes(e, 1:2);
    %x-coordinate of Node
    node_xx = p.Node_coordinate(el_node);
    % (If Node 2 is ahead of Node 1 then Le is positive)
    % (Otherwise it is negative and accordingly B_T matrix would be made)
    Le = (node_xx(1,2)-node_xx(1,1));

    % Transformation term
    temp=[-sign(Le)*p.Element_Area(e)  sign(Le)*p.Element_Area(e)]; 
    
    B_T(e, el_node(1)) = B_T(e, el_node(1)) + temp(1);
    B_T(e, el_node(2)) = B_T(e, el_node(2)) + temp(2);
end

% Loop Solving for Non-Linear Elastic Analysis

% Value for Accuracy
tol = 1.0e-8;
i=2;
while i<n+1
    
    % Displacement increment at each Step
    Displacement_increment(3) = 0.001;
    
    % Force Evaluation due to Displacement Increment
    Stiffness = Global_Stiff_Matrix(p,p.Element_G);
    u_increment = Displacement(p,Stiffness,Displacement_increment,Force_increment);
    u(:,i) = u(:,i-1) + u_increment(:,1);
    
    %Displacement change after each step 
    Delta_u = u(:,i)-u(:,i-1);
    %Displacement change after each iteration
    delta_u = zeros(p.GDof,1);
    
    % Convergence loop
    conv = 10; 
    iter = 0;
    
    while conv > tol 

        % Evaluating new Global Stiffness Matrix
        K_t = Global_Stiff_Matrix(p,p.Element_G);
        
        % Displacement increment
        Delta_u(:,1) = Delta_u(:,1) + delta_u(:,1);
        
        % Strain Increment
        % B=Disp_Strain(p,u(:,i));
        Shear_Delta_Strain(:,1) = B*Delta_u(:,1);
        
        % I3: Elements in Softening Region
        I3 = find(Shear_Plastic_Strain(:,i-1)>gamma_p_s);
        % I8: Elements in Hardening or Elasstic Region
        I8 = setdiff(1:p.Num_Elements, I3);
        
        % Trial Stress Evaluation
        tr_sigma(:,1) = Shear_Stress(:,i-1) + G*Shear_Delta_Strain(:,1);
        
        % Yeild Stress Evaluation
        % For Hardening Region
        Shear_sigma_yeild(I8,i) = tau_y + H_1*Shear_Plastic_Strain(I8,i-1); 
       
        % For Softening Region
        Shear_sigma_yeild(I3,i) = tau_y + H_1*gamma_p_s + H_2*(Shear_Plastic_Strain(I3,i-1)-gamma_p_s);
  
        % Evaluation of Change in Plastic Strain 
        % (Evaluating Delta Plastic Strain)
        fr = abs(tr_sigma(:,1))-Shear_sigma_yeild(:,i);
        I0 = find(fr<=0);   
        I1 = setdiff(1:p.Num_Elements, I0);
        
        Common_Elements = intersect(I0,I3);
        % I0: Elements in Elsstic/Unloading region
        I0 = setxor(I0,Common_Elements);
        Common_Elements = intersect(I1,I3);
        % I1: Elements in Hardening region
        I1 = setxor(I1,Common_Elements);
        
        Shear_Delta_plastic_strain(I8,1) = fr(I8)/(G + H_1);
        Shear_Delta_plastic_strain(I3,1) = fr(I3)/(G + H_2);
        Shear_Delta_plastic_strain(Shear_Delta_plastic_strain<0)=0;
        
        % Evaluation of Plastic Strain
        Shear_Plastic_Strain(:,i) = Shear_Plastic_Strain(:,i-1) + Shear_Delta_plastic_strain(:,1);
        
        % Evaluation of Stress
        Shear_Stress(I0,i) = tr_sigma(I0,1);
        Shear_Stress(I1,i) = tr_sigma(I1,1) - G*sign(tr_sigma(I1,1)).*Shear_Delta_plastic_strain(I1,1);
        Shear_Stress(I3,i) = tr_sigma(I3,1) - G*sign(tr_sigma(I3,1)).*Shear_Delta_plastic_strain(I3,1);
        % Below commented line could be used in Explicit method
        % Shear_Stress(I3,i) = Shear_Stress(I3,i-1) + G_2*Shear_Delta_Strain(I3,1);
        
        % Updating Stiffness values According to different regions
        p.Element_G(I0) = G;
        p.Element_G(I1) = G_1;
        p.Element_G(I3) = G_2;
        
        % Ressidual Evaluation
        R = -B_T'*Shear_Stress(:,i);

        % convergence is checked at active Dof 
        conv = (norm(R(p.active_Dof)))^2;       
        
        % Displacement increment in each Iteration at Active Dof
        delta_u(p.active_Dof,1) = K_t(p.active_Dof,p.active_Dof)\R(p.active_Dof);
        
        iter = iter + 1; 
    end  
    
    % Displacement Evaluation
    u(:,i) = u(:,i-1) + Delta_u;
    % Strain Evaluation
    Shear_Strain(:,i) = B*u(:,i);
    
    % The Crack has been taken place  if the below conditions are satisfied 
    if any(Shear_Stress(:,i) < 0) %&& any(Strain(:,i) >= epsilon_p_max)
        disp('Crack is happenend in the Structure')
        % Delete Elements Presents after Crack Formation
        Shear_Stress(:,i+1:n+1)=[];
        Shear_Strain(:,i+1:n+1)=[];
        u(:,i+1:n+1)=[];        %Displacement at which crack happens
        Force(:,i+1:n+1)=[];    %Force at the time of crack
        break;
    end   
    i=i+1;
end


%%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Formation of Stress-Strain Curve for whole Structure
%%%

% Evaluating Strain values at each step for whole Tablet
% (Thus divided by Tablet Length = rho*t)
Strain = (u(3,:)-u(1,:))./(rho*t);

% Tensile Stress Evaluation
%(Both the Elements will have same value of Stress )
%(ThusEvaluating Stress using one of the Element)
Stress = phi*k*rho*Shear_Stress(1,:); 

% Plotting Stress-Strain Curve
figure(1)
hold on
plot(Strain,Stress/10^6,'y','LineWidth',1.5)
axis([0 inf 0 inf])
xlabel('Tensile Strain');
ylabel('Tensile Stress (MPa)');
title('Tensile Stress-Strain Curve');
grid on;
grid minor;
