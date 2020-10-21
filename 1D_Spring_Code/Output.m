%Output of Displacements and Reactions in Tabular form

function Output(p)
%Displacements 
disp('Displacements') 
nn = 1:p.Num_Nodes; 
[nn' p.Displacement]            %Display displacements at each node

%Reactions 
F = p.Stiffness*p.Displacement;  %F = [K] x U   
p.Reaction = F; 
disp('Reactions') 
[nn' p.Reaction]     %Display Forces at prescribed nodes
end