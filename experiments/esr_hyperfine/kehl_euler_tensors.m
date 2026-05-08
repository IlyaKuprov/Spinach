% This function calculates the euler angles (alpha,beta,gamma) to rotate
% from the v-principal axes frame into the g-principal axes frame
% Since from DFT calculations usually the g-values are ordered
% g1<g2<g3 and usually g1>g2>g3 is used in simulations the order of the
% values is switched and the calulated angles refer turn into the g1>g2>g3
% frame.
%
% input parameters:
% g: the g-tensor
% v: the tensor to be transformed into the g-frame
%
% December 2023 A. Kehl (akehl@gwdg.de)
%


function [alpha,beta,gamma]=kehl_euler_tensors(g,v)

    % Check consistency
    grumble(g,v);
v1=v(:,1);
v2=v(:,2);
v3=v(:,3);

d11=cos(atan2(norm(cross(-g(:,3),v1)),dot(-g(:,3),v1)));
d21=cos(atan2(norm(cross(g(:,2),v1)),dot(g(:,2),v1)));
d31=cos(atan2(norm(cross(g(:,1),v1)),dot(g(:,1),v1)));

d12=cos(atan2(norm(cross(-g(:,3),v2)),dot(-g(:,3),v2)));
d22=cos(atan2(norm(cross(g(:,2),v2)),dot(g(:,2),v2)));
d32=cos(atan2(norm(cross(g(:,1),v2)),dot(g(:,1),v2)));

d13=cos(atan2(norm(cross(-g(:,3),v3)),dot(-g(:,3),v3)));
d23=cos(atan2(norm(cross(g(:,2),v3)),dot(g(:,2),v3)));
d33=cos(atan2(norm(cross(g(:,1),v3)),dot(g(:,1),v3)));

DCM=[d11,d12,d13;d21,d22,d23;d31,d32,d33];

gamma=atand(-d23/d13);
beta=acosd(d33);
alpha=atand(d32/d31);





end

function grumble(g,v)
if ~isnumeric(g)
    error('g must be numeric.');
end
if ~isnumeric(v)
    error('v must be numeric.');
end
end

