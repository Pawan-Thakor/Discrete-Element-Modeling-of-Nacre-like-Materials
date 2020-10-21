% Evaluation of Global Stiffness matrix

function Global_Matrix = Global_Stiff_Matrix(p)
 Stiffness=zeros(p.Num_Nodes);
   for e=1:p.Num_Elements 
     %Element_Stiff_Matrix: Element stiffness Matrix  
     Element_Stiff_Matrix=[1 -1;-1 1]*p.Element_Stiffness(e);    
     %Element_Dof: Nodes for the element   
     Element_Dof=p.Element_Nodes(e,:) ; 
     %Stiffness: Assembly of Element Stiffness Matrix
     Stiffness(Element_Dof,Element_Dof) = Stiffness(Element_Dof,Element_Dof) + Element_Stiff_Matrix; 
   end
   
 %Returning Global Stiffness Matrix  
 Global_Matrix = Stiffness; 
end
