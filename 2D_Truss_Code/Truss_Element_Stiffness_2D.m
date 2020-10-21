% Evaluation of Truss Element Stifffness matrix

function Element_Stiff_Matrix = Truss_Element_Stiffness_2D(node_xy, E, A)
% This function evaluates 4x4 element stiffness Matrix
% This Matrix must be in GLOBAL Coordinates

% Evaluation of Transformation matrix
Le = norm([(node_xy(1,1)-node_xy(2,1)) (node_xy(1,2)-node_xy(2,2))]);
l = (node_xy(2,1)-node_xy(1,1))/Le;
m = (node_xy(2,2)-node_xy(1,2))/Le;
T = [l m 0 0 ; 0 0 l m];         %Transformation Matrix

% Element Stiffness Matrix
Ke = (E*A/Le)*[1 -1 ; -1 1];

%Returning Truss Element Stiffness Matrix in Global Coordinate
Element_Stiff_Matrix = T'*Ke*T;
end
