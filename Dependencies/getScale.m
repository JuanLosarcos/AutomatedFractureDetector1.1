%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: getScale
% by: Juan Manuel Losarcos (UIC)
%
% Description:
%   Allows the user to interactively draw a reference line over an object 
%   of known length in the image. Computes and returns the distance 
%   between the endpoints in pixels.
%
% Input:
%   - (none) — user interaction via GUI
%
% Output:
%   - pixelDist: distance in pixels between the two points of the line
%   - hLine: handle to the drawn line (useful for deletion or updates)
%
% Dependencies:
%   - Requires MATLAB Image Processing Toolbox for drawline()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pixelDist, hLine] = getScale()

    disp(' ');
    disp('> Click and drag to draw a line over the reference object');

    % Draw the line ROI in red
    hLine = drawline('Color','r');

    % Extract endpoints
    x1 = hLine.Position (1,1);
    x2 = hLine.Position (2,1);
    y1 = hLine.Position (1,2);
    y2 = hLine.Position (2,2);

    % Compute Euclidean distance in pixels
    pixelDist = ((x2-x1)^2+(y2-y1)^2)^.5;
end
