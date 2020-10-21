function[K,Internal_Force,F,Residual,Shear_Stress,Shear_Plastic_Strain,Colour,Element_G] = sys_eq(u,u_0,lambda,p,S_S,S_PS,G,C) % Define system of nonlinear equations
% System can be general or based on finite element model. K is tangent
% stiffness or Jacobian matrix. This is the slope of tangent hyperplane at
% displacement u. fu is system value at state u. F is a reference vector 
% for scalar lambda.

Delta_u = u-u_0;
Shear_Delta_Strain(:,1) = p.B*Delta_u(:,1);

[Shear_Stress,Shear_Plastic_Strain,Element_G,Colour] = State_Determination(S_S, S_PS, Shear_Delta_Strain(:,1),G,C, p);

Stress(:,1) = phi*p.Element_Length.*Shear_Stress(:,1)/p.t; 
Internal_Force = p.B_T_1'*Stress(:,1);

K = Global_Stiff_Matrix(p,Element_E);
F = p.F;
Residual = lambda*p.F-Internal_Force;

end