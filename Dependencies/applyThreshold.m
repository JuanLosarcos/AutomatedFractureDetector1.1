%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: applyThreshold
% Modified by: Juan Manuel Losarcos (UIC), based on Prabhakaran et al. (2019)
%
% Description:
%   Applies a user-defined global ridge-intensity threshold to a sigmoid-
%   transformed ridge-strength image. Pixels below the threshold are set to
%   0 and pixels above the threshold are set to 1, producing a binarized
%   ridge map. The function also displays the result for visual QC.
%
% Input:
%   - C_Ridges_norm_sigmoid: ridge ensemble image after normalization and
%     sigmoid enhancement (expected range ~[0,1])
%
% Output:
%   - C_Ridges_norm_thresh: binarized ridge image (uint8, values 0 or 255)
%   - C_Ridges_norm_sigmoid: unchanged input image returned for convenience
%   - threshold: numeric threshold selected by the user
%
% Dependencies:
%   - MATLAB Image Processing Toolbox: im2uint8, imshow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C_Ridges_norm_thresh, C_Ridges_norm_sigmoid, threshold] = applyThreshold(C_Ridges_norm_sigmoid)
    disp(' ');
    disp('> Ridge-intensity Global Threshold:');
    disp('Enter a threshold value to segment the detected ridges');
    disp(['This step filters isolated pixels and disconnected regions,' newline ...
          'removing false positives and yielding a more accurate fracture map.']);
    disp([newline '             NOTE             ' newline '---------------------------------']);
    disp('Threshold: Higher values reduce noise by raising the cutoff for ridge detection');
    disp('-------------> Single numeric value required');
    disp('-------------> Examples: [0.4 0.5 0.6 0.7], suggested: 0.58');
    disp(' ');
    
    % Prompt for threshold
    C_Ridges_norm_thresh = C_Ridges_norm_sigmoid;
    threshold = input('Threshold = ');

    % Apply threshold to normalize ridge intensities

    C_Ridges_norm_thresh(C_Ridges_norm_thresh < threshold) = 0;
    C_Ridges_norm_thresh(C_Ridges_norm_thresh > threshold) = 1;

    % Convert to uint8 for saving/display
    C_Ridges_norm_thresh = im2uint8(C_Ridges_norm_thresh);

    % Display the binarized result
    figure(3);
    imshow(C_Ridges_norm_thresh);
    title(['Normalized Ridge Intensity Threshold: ' num2str(threshold)]);
end