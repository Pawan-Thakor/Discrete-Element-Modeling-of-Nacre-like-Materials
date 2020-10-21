% MATLAB codes for Finite Element Analysis
% It consist of an Algorithm which consist Displacement Increment
clear all 
clear

format long g

%Structure p
p=struct();

%% Defining Problem

% Properties for Shear Stress-Strain Curve
G=100*10^6;                    %Shear Modulus
H=6*10^6;                      %Hardening Value
tau_y=3.5*10^6;                %Yeild Shear Stress value
gamma_u=2;                     %Max Allowed Plastic Strain 
gamma_p_s=1.5;                 %Plastic Strain at start of Softening region
tau_s = tau_y + H*gamma_p_s;   %Max Shear Stress Point

% Parameters of the Tablets
rho=4;                         %Tablet Ratio
k=0.3;                         %Ovterlap Ratio
ti=1/9;                        %Interface Thickness
t=1;                           %Tablet Thickness
phi=t/(t+ti);                  %Tablet Volume concentration

% Properties of Stress-Strain curve
E = (phi^2/(1-phi))*(k*(1-k)*rho^2*G);
E_T_1_1 = (phi^2/(1-phi))*(k*rho^2*G)/(G/H + 1/(1-k));
E_T_1_2 = (phi^2/(1-phi))*(k*(1-k)*rho^2*G)/(G/H+ 1);
E_T_2 = (phi^2/(1-phi))*(k*(1-k)*rho^2*G)/(((G/tau_s)*(1-k))*(gamma_p_s-gamma_u) + 1);
H_1_1 = E*E_T_1_1/(E-E_T_1_1);
H_1_2 = E*E_T_1_2/(E-E_T_1_2); 
H_2 = E*E_T_2/(E-E_T_2);
sigma_y_1 = tau_y*phi*k*rho;
sigma_y_2 = tau_y*phi*(1-k)*rho; 

if (H*gamma_p_s/tau_y) > (1-2*k)/k
    gamma_2 = (k/(1-k))*(gamma_p_s - ((1-2*k)/k)*tau_y/H);
else
    gamma_2 = 0;
end 
epsilon_p_u = (1/rho)*((1-phi)/phi)*(gamma_p_s + gamma_2); 

%% Defining Coordinates 

% Node coordinates
p.Node_coordinate=[0 0; rho*k*t 0; rho*t 0 ];
% Num_Nodes: number of nodes 
p.Num_Nodes=size(p.Node_coordinate,1);

% Element_Nodes: connections at elements 
p.Element_Nodes=[1 2; 2 3];
% Num_Elements: number of Elements 
p.Num_Elements=size(p.Element_Nodes,1); 
% Element_Area: Area of Element
p.Element_Area=t*1*ones([1,p.Num_Elements]);

% Element_Stiffness: truss stiffness values
p.Element_E=E*ones([1,p.Num_Elements]);
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
% Force: Force vector 
p.Force=zeros(p.GDof,1);
% Stiffness: Stiffness Matrix 
p.Stiffness=zeros(p.GDof,p.GDof);

% Number of Displacements Steps
n=3000;

% u: Displacement for all Steps
u = zeros(p.Num_Nodes,n+1);
% Force: Force for all Steps
Force = zeros(p.Num_Nodes,n+1);
% Strain: for each Element for all Steps
Strain = zeros(p.Num_Elements,n+1);
% Stress: for each Element for all Steps
Stress = zeros(p.Num_Elements,n+1);

% Yeild Stress: at each Step
sigma_yeild = zeros(p.Num_Elements,1);
% Plastic Strain: for all Steps
Plastic_Strain = zeros(p.Num_Elements,n+1);
% Change in Plastic Strain at each Step
Delta_plastic_strain = zeros(p.Num_Elements,1);

% Active Dof: Nodes at which Force is supposed to be zero
% As Displacement Increment is considered 
active_Dof = setdiff([1:p.GDof]', [p.Prescribed_Dof]);

%% Solving Spring Problem 

% Displacement to Strain Transformation Matrix
B = zeros(p.Num_Elements,p.Num_Nodes);       
for e=1:p.Num_Elements
    %Nodes of each Element 
    el_node = p.Element_Nodes(e, 1:2);
    %x-coordinate of Node
    node_xx = p.Node_coordinate(el_node);
    % Evaluation of Transformation matrix
    Le = abs((node_xx(1,1)-node_xx(1,2)));          
    temp=[-1/Le 1/Le];         
    B(e, e:e+1) = B(e, e:e+1) + temp;
end

%Stress to Force Transformation at each Node
B_T = zeros(p.Num_Elements,p.Num_Nodes);       
for e=1:p.Num_Elements
    %Nodes of each Element 
    el_node = p.Element_Nodes(e, 1:2);
    % Evaluation of Transformation matrix         
    temp=[-p.Element_Area(e) p.Element_Area(e)];         
    B_T(e, e:e+1) = B_T(e, e:e+1) + temp;
end

% Loop Solving for Non-Linear Elastic Analysis

% Value for Accuracy
tol = 1.0e-8;
for i=2:n+1
    
    % Displacement increment at each Step 
    p.Displacement(3) = 0.0005*(i-1);
    
    % Force Evaluation due to Displacement Increment
    p.Stiffness = Global_Stiff_Matrix(p);
    u(:,i) = Displacement(p,p.Stiffness,p.Displacement,p.Force);
    Force(:,i) = p.Stiffness*u(:,i);
    
    %Displacement change after each step 
    Delta_u = u(:,i)-u(:,i-1);
    %Displacement change after each iteration
    delta_u = zeros(p.GDof,1);
    
    % Convergence loop
    conv = 10; 
    iter = 0;
    
    while conv > tol  
        
        % Displacement increment
        Delta_u = Delta_u + delta_u;
        
        % Strain Increment
        Delta_Strain = B*Delta_u;
        
        I3 = find(Plastic_Strain(:,i-1)>=epsilon_p_u);
        I8 = setdiff([1:p.Num_Elements]', [I3]);
        
            if (H*gamma_p_s/tau_y) > (1-2*k)/k
                I4 = find(Plastic_Strain(:,i-1)>=(sigma_y_2-sigma_y_1)/H_1_1);
            else
                I4 = [];
            end 
   
        % Trial Stress Evaluation
        tr_sigma = Stress(:,i-1) + E*Delta_Strain(:,1);
        
        % Yeild Stress Evaluation
        % For First Hardening Region
        sigma_yeild(I8) = sigma_y_1 + H_1_1*Plastic_Strain(I8,i-1);
        % For Second Hardening Region
        sigma_yeild(I4) = sigma_y_2 + H_1_2*(Plastic_Strain(I4,i-1)-(sigma_y_2-sigma_y_1)/H_1_1);       
        % For Softening Region
        if isempty(I4)
            sigma_yeild(I3) = sigma_y_1 + H_1_1*epsilon_p_u + H_2*(Plastic_Strain(I3,i-1)-epsilon_p_u);
        else
            sigma_yeild(I3) = sigma_y_2 + H_1_2*(epsilon_p_u-(sigma_y_2-sigma_y_1)/H_1_1) + H_2*(Plastic_Strain(I3,i-1)-epsilon_p_u);
        end
        
        % Evaluation of Change in Plastic Strain
        fr = abs(tr_sigma)-sigma_yeild;
        I = find(fr<=0);   
        I1 = find(fr>0);

        Delta_plastic_strain(I8,1) = fr(I8)/(E+H_1_1);
        Delta_plastic_strain(I4,1) = fr(I4)/(E+H_1_2);
        Delta_plastic_strain(I3,1) = fr(I3)/(E+H_2);
        Delta_plastic_strain(Delta_plastic_strain<0)=0;
        
        % Evaluation of Stress and Plastic Strain
        Plastic_Strain(:,i) = Plastic_Strain(:,i-1) + Delta_plastic_strain(:,1);
        Stress(I,i) = tr_sigma(I);
        Stress(I1,i) = tr_sigma(I1)-sign(tr_sigma(I1))*E.*Delta_plastic_strain(I1,1);
        
        % Updating Stiffness values According to different regions
        p.Element_E(I) = E;
        p.Element_E(I1) = E_T_1_1;
        p.Element_E(I4) = E_T_1_2;
        p.Element_E(I3) = E_T_2;
        
        % Evaluating new Global Stiffness Matrix
        K_t = Global_Stiff_Matrix(p);
        % Ressidual Evaluation
        R = Force(:,i) - B_T'*Stress(:,i);
        
        % convergence is checked at active Dof as they hav zero
        % displacement
        conv = norm(R(active_Dof))^2;
        
        % Displacement increment in each Iteration at Active Dof
        delta_u(active_Dof) = K_t(active_Dof,active_Dof)\R(active_Dof);

        iter = iter + 1;
    end  
    
    % Displacement Evaluation
    u(:,i) = u(:,i-1) + Delta_u;
    % Strain Evaluation
    Strain(:,i) = B*u(:,i);
    
    % The Crack has been taken place  if the below conditions are satisfied 
    if any(Stress(:,i) <= 0) %&& any(Strain(:,i) >= epsilon_p_max)
        disp('Crack is happenend in the Structure')
        % Delete Elements Presents after Crack Formation
        Stress(:,i+1:n+1)=[];
        Strain(:,i+1:n+1)=[];
        u(:,i+1:n+1)=[];        %Displacement at which crack happens
        Force(:,i+1:n+1)=[];    %Force at the time of crack
        break;
    end    
end

figure(1)
plot(Strain(1,:),Stress(1,:)/10^6,'y','LineWidth',1.5)
hold on
axis([0 inf 0 inf])
xlabel('Tensile Strain');
ylabel('Tensile Stress (MPa)');
title('Tensile Stress-Strain Curve');
grid on;
grid minor;
