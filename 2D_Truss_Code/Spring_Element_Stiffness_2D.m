% Evaluation of Truss Element Stifffness matrix

function Element_Stiff_Matrix = Spring_Element_Stiffness_2D(k,a)
% This function evaluates 4x4 element stiffness Matrix
% This Matrix must be in GLOBAL Coordinates

% Evaluation of Transformation matrix
l = cosd(a);
m = sind(a);
T = [l m 0 0 ; 0 0 l m];         %Transformation Matrix

% Element Stiffness Matrix
Ke = (k)*[1 -1 ; -1 1];

%Returning Truss Element Stiffness Matrix in Global Coordinate
Element_Stiff_Matrix = T'*Ke*T;
end
