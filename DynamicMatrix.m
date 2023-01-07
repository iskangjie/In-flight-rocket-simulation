function [Me,Ke] = DynamicMatrix(Le,rho,Young,x,y,varargin)
% Return the mass and stiffness matrix of an element
% Le: length of the element
% rho: density of the material
% Young: Young's modulus
% x,y: size of the cross section of the element
% The type of the cross section could be rectanular or pipe, and
% rectangular as default
if isempty(varargin) || strcmp(varargin{1},'rectangular')
    width = x;
    height = y;
    area_cross = width*height;
    inertial_moment = 1/12*width*height^3;
elseif strcmp(varargin{1},'pipe')
    inner_radius = x;
    outer_radius = y;
    area_cross = pi*(outer_radius^2-inner_radius^2);
    inertial_moment = pi/4*(outer_radius^4-inner_radius^4);
end
Me = [156     22*Le   54      -13*Le
    22*Le  4*Le^2  13*Le  -3*Le^2
    54      13*Le   156     -22*Le
    -13*Le  -3*Le^2 -22*Le 4*Le^2]*rho*area_cross*Le/420;

Ke = [12     6*Le    -12     6*Le
    6*Le  4*Le^2  -6*Le  2*Le^2
    -12    -6*Le     12    -6*Le
    6*Le  2*Le^2  -6*Le  4*Le^2]*Young*inertial_moment/Le^3;
end