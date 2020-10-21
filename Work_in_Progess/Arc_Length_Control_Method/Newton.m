% Arc-length method
function[u,lambda,F,S_S,S_PS,G_array,C] = Newton(increment,u_0,Force,p,Shear_Stress,Shear_Plastic_Strain,Element_G,Colour,tol,max_iter, Method_control)
    
if Method_control==1
    
    % Displacment Increment:
    Displacement_increment=zeros(p.GDof,1); 
    % Force Increment:
    Force_increment=zeros(p.GDof,1);   
    
    % Displacement increment at each Step
    Displacement_increment(p.Displacement_Dof) = increment;
    
    %As it is not a force incrment thus lambda=0
    lambda=0;
    
    % Stiffness Evaluation
    Stiffness = Global_Stiff_Matrix(p,Element_G);
    
    % u_increment = p.Displacement;
    u_increment = Displacement(p,Stiffness,Displacement_increment,Force_increment,p.Prescribed_Dof_disp);
    u = u_0 + u_increment(:,1); 
    
    % Force Evaluation due to Displacement Increment
    F = Force + Stiffness*u_increment(:,1);
    
    % Displacement change after each step 
    Delta_u(:,1) = u-u_0;
    %Displacement change after each iteration
    delta_u(:,1) = zeros(p.GDof,1);  

    % Convergence loop
    conv = 0; 
    iter = 1;
    
    while conv == 0 && iter < max_iter 
        
        % Evaluating new Global Stiffness Matrix
        K_t = Global_Stiff_Matrix(p,Element_G);
        
        % Displacement increment
        Delta_u(:,1) = Delta_u(:,1) + delta_u(:,1);
        
        % Strain Increment
        Delta_Strain(:,1) = p.B*Delta_u(:,1);
           
        % Determination of Stress and Plastic Strain
        [S_S, S_PS, Element_G, C] = State_Determination(Shear_Stress, Shear_Plastic_Strain, Delta_Strain, Element_G, Colour, p); 
        
        % Ressidual Evaluation
        R = 0-p.B_T_1'*S_S;

        % convergence is checked at active Dof 
        convergence = (norm(R(p.active_Dof_disp)))^2; 
        
        if convergence < tol
            conv = 1;           % converged
        end
        
        % Displacement increment in each Iteration at Active Dof
        delta_u(p.active_Dof_disp,1) = K_t(p.active_Dof_disp,p.active_Dof_disp)\R(p.active_Dof_disp);
        
        iter = iter + 1; 
    end 
    G_array=Element_G;
    
elseif Method_control==2
    
    % Displacment Increment:
    Displacement_increment=zeros(p.GDof,1);
    % Force Increment:
    Force_increment=zeros(p.GDof,1);  
    
    % Displacement increment at each Step
    Force_increment(p.Force_Dof) = increment;
    
    % Force Evaluation due to Displacement Increment
    F_ext = Force + Force_increment;
    lambda = max(F_ext);
    
    Stiffness = Global_Stiff_Matrix(p,Element_G);
    
    %Displacement change after each step      
    Delta_u(:,1) = Displacement(p,Stiffness,Displacement_increment,Force_increment,p.Prescribed_Dof_force);
    %Displacement change after each iteration
    delta_u(:,1) = zeros(p.GDof,1);  

    % Convergence loop
    conv = 0; 
    iter = 1;
    
    while conv == 0 && iter < max_iter 
        
        % Evaluating new Global Stiffness Matrix
        K_t = Global_Stiff_Matrix(p,Element_G);
        
        % Displacement increment
        Delta_u(:,1) = Delta_u(:,1) + delta_u(:,1);
        
        % Strain Increment
        Delta_Strain(:,1) = p.B*Delta_u(:,1);
           
        % Determination of Stress and Plastic Strain
        [S_S, S_PS, Element_G, C] = State_Determination(Shear_Stress, Shear_Plastic_Strain, Delta_Strain, Element_G, Colour, p); 
        
        Stress = p.phi*p.Element_Length.*S_S/p.t; 
        Internal_Force = p.B_T_1'*Stress;
        % Ressidual Evaluation
        R = F_ext-p.B_T'*S_S;
        
        % convergence is checked at active Dof 
        convergence = (norm(R(p.active_Dof_force)))^2; 
        
        if convergence < tol
            conv = 1;           % converged
        end
        
        % Displacement increment in each Iteration at Active Dof
        delta_u(p.active_Dof_force,1) = K_t(p.active_Dof_force,p.active_Dof_force)\R(p.active_Dof_force);
        
        iter = iter + 1; 
    end  
    u = u_0 + Delta_u;
    F = p.B_T'*S_S;
    G_array = Element_G;
else
    Displ("Incorrect Method Control Number");
end
end
