%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: getNorthCorrection
% By: Juan Manuel Losarcos (UIC)
%
% Description:
%   Lets the user draw a reference line on the image, prompts for that
%   fracture’s true azimuth, and returns the angular offset required to
%   rotate all MATLAB-measured bearings so that 0° aligns with geographic
%   north. Also returns the ROI handle for optional deletion/redraw.
%
% Inputs  : (interactive) – mouse line ROI and numeric azimuth
% Outputs : correctionAngle [deg], roi (handle to drawn line)
%
% Dependency:
%   – MATLAB Image Processing Toolbox  (drawline)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [correctionAngle, roi] = getNorthCorrection()

  %-- Step 1: Prompt and draw --
  disp(' ');
  disp('> Click and drag to draw a line along a fracture of known azimuth.');
  roi = drawline('Color','g');            % requires Image Processing Toolbox

  % Extract endpoints

    x1 = roi.Position (1,1);
    x2 = roi.Position (2,1);
    y1 = roi.Position (1,2);
    y2 = roi.Position (2,2);
    
   
  %-- Step 2: Compute the “raw” photo angle relative to MATLAB’s built-in north (0–360° CCW from +X) --
  photoAngle = atan2d( y2 - y1, x2 - x1 );    % returns –180..+180
  photoAngle = mod(photoAngle, 360);         % map to 0..360

  %-- Step 3: Ask for the real azimuth of that fracture (0–360°) --
  referenceAz = mod(input('Enter the true azimuth of the drawn fracture [0–360]: '),360);

  % Helper for shortest signed difference
  angDiff = @(a,b) mod(a - b + 180,360) - 180;

  %-- Step 4: Compute the correction to apply to all MATLAB angles
  correctionAngle = angDiff(referenceAz, photoAngle);

  fprintf('Apply a correction of %+0.2f° to MATLAB angles to align with true north.\n', ...
          correctionAngle);
end