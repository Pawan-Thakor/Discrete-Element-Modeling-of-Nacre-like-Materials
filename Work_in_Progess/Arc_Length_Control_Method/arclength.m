% Arc-length method
function[u,lambda,Force,S,PS,Element_E,C] = arclength(arcL,u_0,lambda_0,p,Stress,Plastic_Strain,E_array,Colour,tol,max_iter)

    [row,~] = size(u_0);
    du_g = zeros(p.GDof,1);

    % arcL     -user specified arc-length
    % u_0      -guessed or initial state
    % lambda_0 -guessed or initial parameter
    [K_0,I_F,F,~,S,PS,~,~] = sys_eq(u_0,u_0,lambda_0,p,Stress,Plastic_Strain,E_array,Colour);  %output system data for point u_0
    
    del_lambda_g = sign(det(K_0(p.active_Dof_force,p.active_Dof_force)));
%     del_lambda_g = 1;
%     I3 = find(PS(:,1)>p.epsilon_p_u);  
%     if ((isempty(I3)==0))
%         del_lambda_g = -1;
%     end 
     
    du_g(p.active_Dof_force) = K_0(p.active_Dof_force,p.active_Dof_force)\(del_lambda_g*F(p.active_Dof_force));    % du_g = inv(K_0)*del_lambda_g*F;
    arcL_g                   = sqrt(del_lambda_g^2 + du_g'*du_g);
    
    del_lambda_0 = (arcL/arcL_g)*del_lambda_g;
    du_0         = (arcL/arcL_g)*du_g;
    
    lambda_i     = lambda_0 + del_lambda_0; % lambda at end of arc-length  
    u_i          = u_0 + du_0;              % u at end of arc-length
    r_i          = [du_0; del_lambda_0];    % store vector defining arc-length
    L_squared    = r_i'*r_i;                % arc-length magnitude

    du_I = zeros(p.GDof,1);
    du_II = zeros(p.GDof,1);

    %%%% begin iterations to determine equilibrium point
    conv  = 0;                               % convergence criteria 
    iter  = 1;                               % starting iteration count
    count = 0;
    
    while conv == 0 && iter < max_iter
        
        r_p = r_i;                 % previous vector for arc-length
        
        [K_i,~,F,R_i,~,~,~,~] = sys_eq(u_i,u_0,lambda_i,p,S,PS,E_array,Colour);   % output system data for point u_i
        
        du_I(p.active_Dof_force)  = K_i(p.active_Dof_force,p.active_Dof_force)\F(p.active_Dof_force);               % du_I  = inv(K_i)*F;
        du_II(p.active_Dof_force) = K_i(p.active_Dof_force,p.active_Dof_force)\R_i(p.active_Dof_force);             % du_II = inv(K_i)*R_i;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            c(1)         = 1 + du_I'*du_I;
            c(2)         = 2*(r_p(row + 1) + r_p(1:row)'*du_I + du_I'*du_II);
            c(3)         = 2*r_p(1:row)'*du_II + du_II'*du_II;
            del_lambda_t = roots(c);
    
            if isreal(del_lambda_t) == 0
                del_lambda_t = real(del_lambda_t);% use with caution
                fprintf('Complex roots for iteration %g\n',iter)
            end
        
            du_t1        = du_II + del_lambda_t(1)*du_I;
            du_t2        = du_II + del_lambda_t(2)*du_I;
            r_t1         = r_p + [du_t1; del_lambda_t(1)];
            r_t2         = r_p + [du_t2; del_lambda_t(2)];
            cos_theta1   = (r_p'*r_t1)/L_squared;
            cos_theta2   = (r_p'*r_t2)/L_squared;
        
            if cos_theta1 > cos_theta2
                dui         = du_t1;
                del_lambdai = del_lambda_t(1);
            else
                dui         = du_t2;
                del_lambdai = del_lambda_t(2);
            end
        
            du_i         = dui;
            del_lambda_i = del_lambdai;
            
        u_i      = u_i + du_i;
        lambda_i = lambda_i + del_lambda_i;            
        r_i  = r_p + [du_i; del_lambda_i]; % updated vector

        [~,I_F,~,R_i,S_i,PS_i,C_i,Element_E_i] = sys_eq(u_i,u_0,lambda_i,p,S,PS,E_array,Colour);
        convergence = (norm(R_i(p.active_Dof_force)))^2; 
        
        if convergence < tol
            conv = 1;           % converged
        end
        iter = iter + 1;
            
        if(iter==max_iter)&&(count==0)
            max_iter=max_iter+1000;
            count=count+1;
        end
    end 
    
    if(iter>=max_iter)
        disp("Convegence did not happen");
    end                                    
    u      = u_i;               % new equilibrium displacement
    lambda = lambda_i;          % new equilibrium load factor
    Force  = I_F;
    S      = S_i;
    PS     = PS_i;
    C      = C_i;
    Element_E = Element_E_i;
end
