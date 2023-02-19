function [x, y] = getCircPoints(xCenter, yCenter, numPoints,radius)
    theta = linspace(0,2*pi,numPoints);
    x = radius * cos(theta) + xCenter;
    y = radius * sin(theta) + yCenter;

