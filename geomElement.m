function [diameter,area,centroid] = geomElement(verts)
% geomElement computes diameter, area and centroid of an element
%
% SYNOPSIS: [diameter,area,centroid] = geomElement(verts)
%
% INPUT: verts:	coordinates of the vertices of the element
%
% OUTPUT: diameter: diameter of the element
%         centroid: centroid of the element
%         area:     area of the element
%

% AUTHOR: Juan G. Calvo and collaborators, 2021

% compute diameter (maximum distance between two vertices)
allPairs  = combnk(1:size(verts,1),2);
diffVerts = verts(allPairs(:,1),:)-verts(allPairs(:,2),:);
diameter  = sqrt(max(sum(diffVerts.*diffVerts,2)));
% formulas for area and centroid
x  = verts(:,1); y  = verts(:,2);
x1 = x([2:end,1]);  y1 = y([2:end,1]);  % cyclic ordering
area_sum = x.*y1 - x1.*y;
area     = sum(area_sum)/2;
centroid = sum([(x+x1), (y+y1)].*area_sum)/(6*area);
end