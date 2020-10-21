% Evalustion of Displacements at each node

function Displacements = Displacement(p,Stiffness,Displacement,Force)
% Active Dof = subtract Prescribed Dof from Set of all Dof 
active_Dof=setdiff([1:p.Num_Nodes]', [p.Prescribed_Dof]);

% Modify forces at active nodes to account for homogenous/non-homogenous BC
Force_Modified=Force(active_Dof)-Stiffness(active_Dof,p.Prescribed_Dof)*Displacement(p.Prescribed_Dof);

%Solution U = Inverse of modified[K] x F for Active Dof
U=Stiffness(active_Dof,active_Dof)\Force_Modified;  

%Substituting Displacement at Active_Dof
Displacement(active_Dof)=U;

%Returning Displacement obtained at all Dof
Displacements =Displacement;
end