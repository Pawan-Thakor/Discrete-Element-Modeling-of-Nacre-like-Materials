 MATLAB codes for Finite Element Analysis
% It consist of an Algorithm which consist Displacement Increment
clear
clear all
close all
format long g

%Structure p
p=struct();

%%% Probem Defination

% Number of Tablets in x-direction
Nx = 3;
% Number of Tablets in y-direction
Ny = 4;
% Finding Even Number in Ny array 
j = 1:1:Ny;
iseven_2=rem(j,2);             %Odd number=1, Even Number=0

% Properties for Shear Stress-Strain Curve of each spring
G=100*10^6;                    %Shear Modulus
H_1=6*10^6;                    %Hardening Value
tau_y=3.5*10^6;                %Yeild Shear Stress value
gamma_u=2;                     %Max Allowed Plastic Strain 
gamma_p_s=1.5;                 %Plastic Strain at start of Softening region

tau_s = tau_y + H_1*gamma_p_s;               %Max Shear Stress Point
G_1 = H_1*G/(G+H_1);                         %Slope of Hardening in Shear stress-strain curve 
G_2 = -tau_s/(gamma_u-(gamma_p_s+tau_s/G));  %Slope of Softening in Shear stress-strain curve
H_2 = G*G_2/(G-G_2);                         %Softening Value

% Parameters of the Tablets
rho_mean=4;                    %Mean Tablet Ratio
k_mean=0.5;                    %Mean Ovterlap Ratio
ti=1/9;                        %Interface Thickness
t=1;                           %Tablet Thickness
phi=t/(t+ti);                  %Tablet Volume concentration

ratio = 0.0;                   %Delta_rho/rho_mean (Used for Statistical Variation)
delta_rho = ratio*rho_mean;    

% d: offset to the tablet given to each row  
d = iseven_2*(1-k_mean)*rho_mean;

% Used for Columnar Layer
% d: For columnar Tablets (put inside the loop)
% d = rho_mean*iseven_2.*rand(Ny,1)';

% Total Number of Tablets
Number_of_Tablets = Nx*Ny + sum(iseven_2(:) == 0);
% Total Number of Non-Linear springsn
Number_of_Springs = Nx*(Ny-1)*2;
Number_of_Springs_in_row = Nx*2;

% Total Number of Tablets in which Spring is Present
% (Last row of Tablet does not contain spring)
% (Note: Even row has one more tablet as compared to odd row)
iseven_1 = iseven_2;
iseven_1(end) = [];
Number_of_Spring_Tablet = Nx*(Ny-1) + sum(iseven_1(:) == 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Properties of Tablets when no Statistical Variation is Applied

% Positions of corner of Tablets when no statistical  variation is applied
x_temp_corner = zeros(Nx+3,Ny);
q=1;
for i=1:1:Nx+3
    for j=1:1:Ny 
        x_temp_corner(i,j)= d(j) + (i-1)*rho_mean*t;
        q=q+1;
    end
end

% rho: Aspect Ratio of each Tablet when no statistical varition is applied
% (Even row has one more Talet as compared to odd row)
rho_temp = zeros(Nx+2,Ny);
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            rho_temp(i,j) = (x_temp_corner(i+1,j)-x_temp_corner(i,j))/(t);
        end 
    else
        for i=2:1:Nx+1
            rho_temp(i,j) = (x_temp_corner(i+1,j)-x_temp_corner(i,j))/(t); 
        end   
    end
end

% x_spring_temp: Spring coordinates when no statistical variation is
% applied
x_spring_temp = zeros(Number_of_Springs_in_row+1,Ny-1);
e=1;
for j=1:1:Ny-1 
    e=1;
    if iseven_2(j)==0
        for i=2:1:Nx+2
            if i==Nx+2
                x_spring_temp(e,j) = x_temp_corner(i,j);
            else
                x_spring_temp(e,j) = x_temp_corner(i,j);
                e=e+1;
                x_spring_temp(e,j) = x_temp_corner(i,j+1);
            end
            e=e+1;
        end 
    else
        for i=2:1:Nx+1
            if i==Nx+1
                x_spring_temp(e,j) = x_temp_corner(i,j+1);
                e=e+1;
                x_spring_temp(e,j) = x_temp_corner(i,j);
                e=e+1;
                x_spring_temp(e,j) = x_temp_corner(i+1,j+1);
            else
                x_spring_temp(e,j) = x_temp_corner(i,j+1);
                e=e+1;
                x_spring_temp(e,j) = x_temp_corner(i,j);       
            end
            e=e+1;               
        end   
    end
end

% Minimum Spring Length: x_spring_minimum when no statistical variaiton is
% applied
min_length = zeros(Number_of_Springs,1);
n=1;
for j=1:1:Ny-1 
    if iseven_2(j)==0
        for i=1:1:(Number_of_Springs_in_row)
            min_length(n) = x_spring_temp(i+1,j) - x_spring_temp(i,j);
            n=n+1;
        end 
    else
        for i=1:1:(Number_of_Springs_in_row)
            min_length(n) = x_spring_temp(i+1,j) - x_spring_temp(i,j);
            n=n+1;
        end   
    end
end
Length =  min(min_length)/2; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Properties of Tablets when Statistical Variation is Applied

% Statistical variation Appliedat the corner positions of the Tablets
if ratio==0
    z=zeros(Nx+3,Ny);
else
    stats2=zeros(3,Ny);
    pd = makedist('Normal','mu',0,'sigma',2*delta_rho);
    pd = truncate(pd,-Length/t,Length/t);
    z = random(pd,(Nx+3)*Ny,1);
end


% Positions of Tablet Nodes with Statistical variation
x_corner = zeros(Nx+3,Ny);
y_corner = zeros(Nx+3,Ny);
q=1;
for i=1:1:Nx+3
    for j=1:1:Ny 
        x_corner(i,j)= d(j) + (i-1)*rho_mean*t + z(q)*t;
        y_corner(i,j)= -t*(j-1)+t/2;
        q=q+1;
    end
end

% rho: Aspect Ratio of each Tablet
% (Even row has one more Talet as compared to odd row)
rho = zeros(Nx+2,Ny);
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            rho(i,j) = (x_corner(i+1,j)-x_corner(i,j))/(t);
        end 
    else
        for i=2:1:Nx+1
            rho(i,j) = (x_corner(i+1,j)-x_corner(i,j))/(t); 
        end   
    end
end

% k: overlap Ratio of Each Interface
% (This overlap ration is not defined as given in the paper)
% (It is defined on the basis of the x_corner position)
k = zeros(Number_of_Spring_Tablet,2);
n=1;
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2       
            if i==2
                k(n,2) = (x_corner(i+1,j)-x_corner(i,j+1))/(rho(i,j)*t); 
                n=n+1;               
            elseif i==Nx+2
                k(n,1) = (x_corner(i,j+1)-x_corner(i,j))/(rho(i,j)*t); 
                n=n+1;                   
            else
                k(n,1) = (x_corner(i,j+1)-x_corner(i,j))/(rho(i,j)*t); 
                k(n,2) = (x_corner(i+1,j)-x_corner(i,j+1))/(rho(i,j)*t);
                n=n+1; 
            end
        end
    else
        for i=2:1:Nx+1
            k(n,1) = (x_corner(i+1,j+1)-x_corner(i,j))/(rho(i,j)*t); 
            k(n,2) = (x_corner(i+1,j)-x_corner(i+1,j+1))/(rho(i,j)*t);
            n=n+1;
        end        
    end
end


% Node_coordinate: Coordinates of each Tablet
Node_coordinate = zeros(Number_of_Tablets,2);
n=1;
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            if j==Ny
                Node_coordinate(n,:) = [x_corner(i,j-1) (y_corner(i,j)-t/2)];
            else
                Node_coordinate(n,:) = [x_corner(i,j+1) (y_corner(i,j)-t/2)];
            end
            n=n+1;
        end 
    else
        for i=2:1:Nx+1
            if j==Ny
                Node_coordinate(n,:) = [x_corner(i+1,j-1) (y_corner(i,j)-t/2)];
            else
                Node_coordinate(n,:) = [x_corner(i+1,j+1) (y_corner(i,j)-t/2)];
            end            
            n=n+1;
        end   
    end
end


% Node_Number: number corresponding to each Tablet Node for Defining
%              Element Connections
% (Each Node is given a number row-wise till all the rows are complete)
% (It is provided as Spring connections could be described)
Node_Number = zeros(Nx+2,Ny);
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            if j==1
                Node_Number(i,j) = (i-1);
            else      
                Node_Number(i,j) = (i-1) + (Node_Number(Nx+1,j-1));
            end
        end 
    else
        for i=2:1:Nx+1
            if j==1
                Node_Number(i,j) = (i-1);
            else      
                Node_Number(i,j) = (i-1) + (Node_Number(Nx+2,j-1));
            end
        end   
    end
end

% Element_Nodes: Nodes cooresponding to each spring Element
Element_Nodes = zeros(Number_of_Springs,2);
% Element_Length: Interface Length for each Spring
% (Length of Each spring is determined using Node number, k, rho, t)
Element_Length = zeros(Number_of_Springs,1);
n=1;
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2            
            if i==2
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i,j+1)];  
                Element_Length(n,1) = k(Node_Number(i,j),2)*rho(i,j)*t;
                n=n+1;
            elseif i==Nx+2
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i-1,j+1)];
                Element_Length(n,1) = k(Node_Number(i,j),1)*rho(i,j)*t;
                n=n+1;
            else
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i-1,j+1)];
                Element_Length(n,1) = k(Node_Number(i,j),1)*rho(i,j)*t;
                n=n+1;
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i,j+1)];
                Element_Length(n,1) = k(Node_Number(i,j),2)*rho(i,j)*t;
                n=n+1;
            end        
        end
    else
        for i=2:1:Nx+1
            Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i,j+1)];
            Element_Length(n,1) = k(Node_Number(i,j),1)*rho(i,j)*t;
            n=n+1;
            Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i+1,j+1)];
            Element_Length(n,1) = k(Node_Number(i,j),2)*rho(i,j)*t;
            n=n+1;
        end   
    end
end

% Prescribed_Dof: Nodes at which Displacement Boundary Condition is defined
% (Nodes which are either fixed or increment is provided)
% (Disp_Dof: Nodes at which Displacment is applied)
n=1;
m=1;
for j=1:1:Ny
    if iseven_2(j)==0
        Prescribed_Dof_disp(n,1) = Node_Number(2,j);
        n=n+1;
        Prescribed_Dof_disp(n,1) = Node_Number(Nx+2,j);
        n=n+1;
        Disp_Dof(m,1) = Node_Number(Nx+2,j);
        m=m+1;
    end
end

% Prescribed_Dof: Nodes at which Displacement Boundary Condition is defined
% (Nodes which are either fixed or increment is provided)
% (Force_Dof: Nodes at which Force is applied)
n=1;
m=1;
for j=1:1:Ny
    if iseven_2(j)==0
        Prescribed_Dof_force(n,1) = Node_Number(2,j);
        n=n+1;
        Force_Dof(m,1) = Node_Number(Nx+2,j);
        m=m+1;
    end
end

%% Parameters for Stress-Strain Curve

% % Properties of Stress-Strain curve
% % (These properties are obtained from Table 1 in abid17 paper)
% % (It is defined for two spring which is converted to single) 
% % (spring using relations in the equations given in the paper above)
% p.E = (phi^2/(1-phi))*(Element_Length.^2*G/t^2);
% p.E_T_1 = (phi^2/(1-phi))*(Element_Length.^2*G/t^2)/(G/H + 1);
% p.E_T_2 = (phi^2/(1-phi))*(Element_Length.^2*G/t^2)/((G/tau_s)*(gamma_p_s-gamma_u) + 1);
% p.H_1 = (p.E).*(p.E_T_1)./(p.E-p.E_T_1);
% p.H_2 = (p.E).*(p.E_T_2)./(p.E-p.E_T_2);
% p.sigma_y_1 = tau_y*phi*Element_Length/t;
% p.epsilon_p_u = (1./(Element_Length/t))*((1-phi)/phi)*(gamma_p_s); 
% 

% Properties for Shear Stress-Strain Curve of each spring
p.G = G;              
p.G_T_1 = G_1;                        
p.G_T_2 = G_2; 
p.H_1 = H_1; 
p.H_2 = H_2; 
p.tau_y = tau_y;
p.gamma_p_s=gamma_p_s;
p.gamma_u=gamma_u;   
p.t=t;
p.phi=phi;


%% Initiaize variables for the code

% Node coordinates
p.Node_coordinate=Node_coordinate;
% Num_Nodes: number of nodes 
p.Num_Nodes=size(p.Node_coordinate,1);

% Element_Nodes: connections at elements 
p.Element_Nodes=Element_Nodes;
% Num_Elements: number of Elements 
p.Num_Elements=size(p.Element_Nodes,1); 
% Element_Length: Length of each spring Element
p.Element_Length = Element_Length;
% Element_Area: Area of each spring Element
p.Element_Area = Element_Length*ti;


% NDof: Degree of Freedom at each node
p.NDof = 1;
% GDof: total number of degrees of freedom
p.GDof=p.Num_Nodes;

%% Given Boundary Conditions
% Boundary condition for Homogeneous Solution

%%%Boundary condition for Displacment Increment
% Prescribed_Dof: Value of displacement is given 
p.Prescribed_Dof_disp = Prescribed_Dof_disp;
% Displacement Applied at Dof
p.Displacement_Dof = Disp_Dof;
% Active Dof: Nodes at which Force is supposed to be zero
p.active_Dof_disp = setdiff([1:p.GDof]', [p.Prescribed_Dof_disp]);

%%%Boundary Condition for Force Increment
% Force_Dof: Value of Force is given
p.Force_Dof = Force_Dof;
% p.Prescribed_Dof_force
p.Prescribed_Dof_force = Prescribed_Dof_force;
% Active Dof force: At Force Boundary Condtion
p.active_Dof_force = setdiff([1:p.GDof]', [p.Prescribed_Dof_force]);

%%%Right side nodes are always fixed
% p.Fixed_Dof: Nodes which are fixed
p.Fixed_Dof = setdiff(p.Prescribed_Dof_disp, p.Displacement_Dof);


%% Initializing Matrix 

% Number of Arc Length Increment Steps
n=10000;
% Displacment Increment:
Displacement_increment=zeros(p.GDof,1);

% u: Displacement for all Steps
u = zeros(p.Num_Nodes,n+1);
% Force: Force for all Steps
Force = zeros(p.Num_Nodes,n+1);
% Strain: for each Element for all Steps
Shear_Strain = zeros(p.Num_Elements,n+1);
% Stress: for each Element for all Steps
Shear_Stress = zeros(p.Num_Elements,n+1);
% Plastic Strain: for all Steps
Shear_Plastic_Strain = zeros(p.Num_Elements,n+1);
% Lambda: Force Addition
lambda = zeros(1,n+1);

%Stiffness value at each iteration
Element_G_array = zeros(p.Num_Elements,n+1);
Element_G_array(:,1) = p.G;


%% Solving Spring Problem 

% Displacement to Shear Strain Transformation Matrix
% (Displacement is converted to Strain)
p.B = zeros(p.Num_Elements,p.Num_Nodes);       
for e=1:p.Num_Elements
    % Nodes of each Element 
    el_node = p.Element_Nodes(e, 1:2);
    %x-coordinate of Node
    node_xx = p.Node_coordinate(el_node);
    
    % Sign for Finding position of Node
    % (If Node 2 is ahead of Node 1 then Le is positive)
    % (Otherwise it is negative and accordingly B matrix would be made)
    Le =(node_xx(1,2)-node_xx(1,1));
    
    % Length of each Element
    temp=[-sign(Le)/ti sign(Le)/ti];   
    
    p.B(e, el_node(1)) = p.B(e, el_node(1)) + temp(1);
    p.B(e, el_node(2)) = p.B(e, el_node(2)) + temp(2);
end

%Shear Stress to Force Transformation at each Node
% (Stress is converted to Force)
p.B_T = zeros(p.Num_Elements,p.Num_Nodes);       
for e=1:p.Num_Elements
    %Nodes of each Element
    el_node = p.Element_Nodes(e, 1:2);
    %x-coordinate of Node
    node_xx = p.Node_coordinate(el_node);
    % (If Node 2 is ahead of Node 1 then Le is positive)
    % (Otherwise it is negative and accordingly B_T matrix would be made)
    Le = (node_xx(1,2)-node_xx(1,1));

    % Area for each Tablet
    temp=[-sign(Le)*p.Element_Area(e)  sign(Le)*p.Element_Area(e)]; 
    
    p.B_T(e, el_node(1)) = p.B_T(e, el_node(1)) + temp(1);
    p.B_T(e, el_node(2)) = p.B_T(e, el_node(2)) + temp(2);
end


%Tensile Stress to Force Transformation at each Node
% (Stress is converted to Force)
p.B_T_1 = zeros(p.Num_Elements,p.Num_Nodes);       
for e=1:p.Num_Elements
    %Nodes of each Element
    el_node = p.Element_Nodes(e, 1:2);
    %x-coordinate of Node
    node_xx = p.Node_coordinate(el_node);
    % (If Node 2 is ahead of Node 1 then Le is positive)
    % (Otherwise it is negative and accordingly B_T matrix would be made)
    Le = (node_xx(1,2)-node_xx(1,1));

    % Area for each Tablet
    temp=[-sign(Le)*t  sign(Le)*t]; 
    
    p.B_T_1(e, el_node(1)) = p.B_T_1(e, el_node(1)) + temp(1);
    p.B_T_1(e, el_node(2)) = p.B_T_1(e, el_node(2)) + temp(2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop Solving for Non-Linear Elasto-Plastic Analysis

% Cell array of colours 
p.C = {[0 0 0],[0.8 0.8 0.8],[1 1 0],[1 0 0]}; 
% Black:    'c'         : Elastic Region
% Grey:     'g'         : Unloading Region
% Yellow:   'y'         : Strain Hardening Region
% Red:      'r'         : Strain Softening Region
Colour=cell(p.Num_Elements,n+1);
Colour(:,1) = p.C(1);
z=1;  

%%% Method==1 Newtons Method || Method==2 Arc Length Method
Method = 1;
i=2;
tol = 1.0e-6;

%%% Newtons Method
Method_control=2;
%Method_control==1 Dispalcement Increment
d_increment = 0.001;
%Method_control==2 Force Increment
f_increment = 10000;
max_iter = 1000;
iter=0;

%%% Arc Length Method 
arcL = 1000;        % specified arc-length to follow solution curve
% Evaluating Force unit vector for Arc length Method
p.F = zeros(p.Num_Nodes,1);
p.F(p.Force_Dof)=1;

if Method_control==1
    increment=d_increment;
elseif Method_control==2
    increment=f_increment;
end

while i<n+1
    
    if Method==1
        u_0 = u(:,i-1);
        [u(:,i),lambda(i),Force(:,i),Shear_Stress(:,i),Shear_Plastic_Strain(:,i),Element_G_array(:,i),Colour(:,i)] = Newton(increment,u_0,Force(:,i-1),p,Shear_Stress(:,i-1),Shear_Plastic_Strain(:,i-1),Element_G_array(:,i-1),Colour(:,i-1),tol,max_iter,Method_control);
        Shear_Strain(:,i) = p.B*u(:,i); 
    elseif Method==2
        u_0 = u(:,i-1);          
        lambda_0 = lambda(i-1); 
        [u(:,i),lambda(i),Force(:,i),Shear_Stress(:,i),Shear_Plastic_Strain(:,i),Element_G_array(:,i),Colour(:,i)] = arclength(arcL,u_0,lambda_0,p,Shear_Stress(:,i-1),Shear_Plastic_Strain(:,i-1),Element_G_array(:,i-1),Colour(:,i-1),tol,max_iter);  
        Shear_Strain(:,i) = p.B*u(:,i);
        iter=iter+1;
    end
%    Stress(Stress<0)=0;
%     I3 = find(Plastic_Strain(:,i)>=p.epsilon_p_u);  
%     if ((isempty(I3)==0))&&(iter==0)
%         Method=2;
%         i=i-1;
%     end 
    
    % The Crack has been taken place  if the below conditions are satisfied 
    if any(Shear_Stress(:,i) < 0) %&& any(Strain(:,i) >= p.epsilon_p_max)
        disp('Crack is happenend in the Structure')
        % Delete Elements Presents after Crack Formation
        Shear_Stress(:,i+1:n+1)=[];
        Shear_Strain(:,i+1:n+1)=[];
        u(:,i+1:n+1)=[];        %Displacement at which crack happens
        Force(:,i+1:n+1)=[];    %Force at the time of crack
        Colour(:,i+1:n+1)=[];
        Shear_Plastic_Strain(:,i+1:n+1)=[];
        lambda(i+1:n+1)=[];
        Element_G_array(:,i+1:n+1)=[];
        break;
    end 
    i=i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Formation of Stress-Strain Curve for whole Structure
%%%

% Evaluating Total length for even rows to evaluate strain
Combined_L = zeros(size(p.Displacement_Dof,1),1);
q=1;
for j=1:1:Ny
    if iseven_2(j)==0
            Combined_L(q,1) = x_corner(Nx+3,j)-x_corner(2,j);
            q=q+1;
    end
end

% Evaluating Strain values at each step for whole structure
Strain_combined = (u(p.Displacement_Dof,:)-u(p.Fixed_Dof,:))./(Combined_L);
% (As multiple rows are present thus average Strain is taken which would be obtained in multiple rows)
if size(Strain_combined,1)==1
    Strain_combined_1 = Strain_combined;
else
    Strain_combined_1 = sum(Strain_combined)/size(Strain_combined,1);
end


% Tensile Stress Evaluation
Stress = zeros(size(Shear_Stress,1),size(Shear_Stress,2));
for i=1:1:size(Shear_Stress,2)
    Stress(:,i) = phi*p.Element_Length.*Shear_Stress(:,i)/t; 
end

% Force Evaluation
Force = zeros(p.Num_Nodes,size(Shear_Stress,2));
for i=1:1:size(Stress,2) 
    Force(:,i) = p.B_T_1'*Stress(:,i);
end

% Evaluating total thickness of Tablets to evaluate Stress
total_t=0;
for j=1:1:Ny
    if j==1
        total_t = total_t + t/2;
    elseif j==Ny
        total_t = total_t + t/2;
    else
        total_t = total_t + t;
    end
end

% Evaluating Stress values at each step for whole structure
Stress_combined = Force(p.Displacement_Dof,:)/(total_t*1);
% (As multiple rows are present thus all Stress values are added)
if size(Stress_combined,1)==1
    Stress_combined_1 = Stress_combined;
else
    Stress_combined_1 = sum(Stress_combined);
end


% Plotting combined Stress-Strain for whole structure
h = figure;
plot(1)

figure(1)
plot(Strain_combined_1,Stress_combined_1/10^6,'k','LineWidth',1.5)
hold on
axis([0 inf 0 inf])
xlabel('Tensile Strain');
ylabel('Tensile Stress (MPa)');
title('Tensile Stress-Strain Curve');
grid on;
grid minor;

% % Figure for differen springs in the Structure
% figure(2)
% plot(Strain(1,:),Stress(1,:),'r','LineWidth',1.5)
% hold on
% axis([0 inf 0 inf])
% xlabel('Strain');
% ylabel('Stress');
% title('Stress-Strain Curve for both the Elements');
% grid on;
% grid minor;
% 
% figure(3)
% plot(Strain(2,:),Stress(2,:),'b','LineWidth',1.5)
% hold on
% axis([0 inf 0 inf])
% xlabel('Strain');
% ylabel('Stress');
% title('Stress-Strain Curve for both the Elements');
% grid on;
% grid minor;
% 
% figure(4)
% plot(Strain(3,:),Stress(3,:),'b','LineWidth',1.5)
% hold on
% axis([0 inf 0 inf])
% xlabel('Strain');
% ylabel('Stress');
% title('Stress-Strain Curve for both the Elements');
% grid on;
% grid minor;
% 
% figure(5)
% plot(Strain(8,:),Stress(8,:),'b','LineWidth',1.5)
% hold on
% axis([0 inf 0 inf])
% xlabel('Strain');
% ylabel('Stress');
% title('Stress-Strain Curve for both the Elements');
% grid on;
% grid minor;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tablet Formation of the Structure

% Rectangle formation with the help of corner position of Tablets
figure(10)
for j=1:1:Ny 
    if iseven_2(j)==0
        for i=2:1:Nx+2
            hold on;
            if i==2
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor','b','EdgeColor','k','LineWidth',2);
            else
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor','b','EdgeColor','k','LineWidth',2);
            end
        end 
    else
        for i=2:1:Nx+1
            hold on;
            if i==2
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor',[0 0 1],'EdgeColor','k','LineWidth',2);
            else             
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor',[0 0 1],'EdgeColor','k','LineWidth',2);
            end
        end   
    end
end

% Colour lines at the last increment stp are formed
% (Colour Lines can be formed with the help of Colour information stored at
% each step)
n=1;
Col = size(Colour,2);
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2 
            if i==2
                x_1=[x_corner(i,j)+rho(i,j)*t-Element_Length(n) , x_corner(i,j)+rho(i,j)*t];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
            elseif i==Nx+2
                x_1=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
            else
                x_1=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
                x_1=[x_corner(i,j)+Element_Length(n-1) , x_corner(i,j)+Element_Length(n-1)+Element_Length(n)];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
            end        
        end
    else
        for i=2:1:Nx+1
            x_1=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
            y_1=[y_corner(i,j) , y_corner(i,j)];
            line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
            n=n+1;
            x_1=[x_corner(i,j)+Element_Length(n-1) , x_corner(i,j)+Element_Length(n-1)+Element_Length(n)];
            y_1=[y_corner(i,j) , y_corner(i,j)];
            line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
            n=n+1;
        end   
    end
end

set(gca, 'Visible', 'off');
