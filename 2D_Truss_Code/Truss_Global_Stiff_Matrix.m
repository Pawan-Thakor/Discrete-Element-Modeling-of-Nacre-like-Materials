% Evaluation of Global Stiffness matrix

function Global_Matrix = Truss_Global_Stiff_Matrix(p)
 Stiffness=zeros(p.GDof);
  
   for e=1:p.Num_Elements        
     el_node = p.Element_Nodes(e, 1:2); 
     
     if e<=(p.Num_Elements-size(p.Spring_Element,1))
         node_xy = p.Node_coordinate(el_node, :);
        
         %Element_Stiff_Matrix: Element stiffness Matrix  
         Element_Stiff_Matrix = Truss_Element_Stiffness_2D(node_xy, p.Element_E(e), p.Element_Area(e));
         
     else
         q=e-(p.Num_Elements-size(p.Spring_Element,1));
         Element_Stiff_Matrix = Spring_Element_Stiffness_2D(p.Spring_Stiffness(q), p.Spring_angle(q));
     end
     
     %Element_Dof: Nodes for the element   
     temp = p.NDof*(el_node(1)-1)+1:p.NDof*(el_node(1)-1)+p.NDof; 
     Element_Dof = [temp p.NDof*(el_node(2)-1)+1:p.NDof*(el_node(2)-1)+p.NDof];
     
     %Stiffness: Assembly of Element Stiffness Matrix
     Stiffness(Element_Dof,Element_Dof) = Stiffness(Element_Dof,Element_Dof) + Element_Stiff_Matrix; 
   end
   
 %Returning Global Stiffness Matrix  
 Global_Matrix = Stiffness;
end