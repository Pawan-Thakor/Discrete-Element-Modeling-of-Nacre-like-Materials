%Output of Displacements and Reactions in Tabular form

function Output(p)

nn = 1:p.GDof; 
%Displacements 
disp('Displacements') 
[nn' p.Displacement]           %Display displacements at each node

%Reactions forces at all nodes  
disp('Reactions') 
[nn' p.Force]                        %Display Forces at each node

%Stiffness Matrix  
disp('K') 
[p.Stiffness]

%Stress for each Element
mm = 1:(p.Num_Elements-size(p.Spring_Element,1));   
disp('Stress') 
[mm' p.Stress]                 %Display Stress for Each element

end