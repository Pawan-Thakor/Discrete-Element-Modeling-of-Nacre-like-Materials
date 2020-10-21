% Evaluation of Stresses at each node

function Stress = Stress_2D_Truss(p)
   for e=1:p.Num_Elements-size(p.Spring_Element,1) 
       
     el_node = p.Element_Nodes(e, 1:2);
     node_xy = p.Node_coordinate(el_node, :);
 
     %Element_Dof: Nodes for the element   
     temp = p.NDof*(el_node(1)-1)+1:p.NDof*(el_node(1)-1)+p.NDof; 
      
     % Evaluation of Transformation matrix
     Le = norm([(node_xy(1,1)-node_xy(2,1)) (node_xy(1,2)-node_xy(2,2))]);
     l = (node_xy(2,1)-node_xy(1,1))/Le;
     m = (node_xy(2,2)-node_xy(1,2))/Le;
     T = [-l -m l m];         %Transformation Matrix 
     
     %Element_Dof: Nodes for the element   
     temp = p.NDof*(el_node(1)-1)+1:p.NDof*(el_node(1)-1)+p.NDof; 
     Element_Dof = [temp p.NDof*(el_node(2)-1)+1:p.NDof*(el_node(2)-1)+p.NDof]; 
     
     %Stress Matrix Evaluation
     Stress(e,1)=p.Element_E(e)/Le*T*p.Displacement(Element_Dof); 
   end
  
end

