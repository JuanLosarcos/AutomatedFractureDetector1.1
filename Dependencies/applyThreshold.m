function [C_Ridges_norm_thresh, C_Ridges_norm_sigmoid, threshold] = applyThreshold(C_Ridges_norm_sigmoid)
% applyThreshold  Prompt for and apply an intensity threshold to segment ridges
%   [binaryRidges, ridgesNormSigmoid, threshold] = applyThreshold(ridgesNormSigmoid)
%   1) Prompts the user to enter a threshold value.
%   2) Filters out isolated pixels and false positives.
%   3) Returns a binary image of detected ridges.

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
    C_Ridges_norm_thresh(C_Ridges_norm_thresh >= threshold) = 1;

    % Convert to uint8 for saving/display
    C_Ridges_norm_thresh = im2uint8(C_Ridges_norm_thresh);

    % Display the binarized result
    figure(3);
    imshow(C_Ridges_norm_thresh);
    title(['Normalized Ridge Intensity Threshold: ' num2str(threshold)]);
end