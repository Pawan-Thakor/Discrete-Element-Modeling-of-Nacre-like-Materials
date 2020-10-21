% Evaluation of Truss Element Stifffness matrix

function Element_Stiff_Matrix = Element_Stiffness(node_xx, E, A)
% This function evaluates 2x2 element stiffness Matrix
% This Matrix must be in GLOBAL Coordinates

% Evaluation of Transformation matrix
Le = abs((node_xx(1,1)-node_xx(1,2)));

% Element Stiffness Matrix
Element_Stiff_Matrix = (E*A/Le)*[1 -1 ; -1 1];
end


