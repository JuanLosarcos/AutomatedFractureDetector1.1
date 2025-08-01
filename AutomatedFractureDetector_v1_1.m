% SPDX-License-Identifier: GPL-3.0-only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: AutomatedFractureDetector_v1_1
% Authors: Losarcos, J. M.; Bernardi, M. I.; and Hryb, D. E.
% Corresponding author: Juan Manuel Losarcos, University of Illinois Chicago (jlosa@uic.edu)
% Copyright (c) 2025 Losarcos, J. M.; Bernardi, M. I.; and Hryb, D. E.
% Portions Copyright (c) 2019 Rahul Prabhakaran, TU Delft
% Based on and modified from Prabhakaran et al. (2019) – see README for details.
%
% Description:
% This MATLAB pipeline automates the detection and extraction of rock
% fracture traces from outcrop/satellite images. It combines and modifies four
% original scripts from Prabhakaran et al. (2019) into a single, user-
% friendly workflow, allowing researchers to run each step independently
% or sequentially.
%
% Dependencies:
%   - CoSHREM Toolbox by Rafael Reisenhofer
%       -> CSHRMgetContRidgeSystem.m
%       -> CSHRMsheardec.m
%       -> CSHRMgetConeOris.m
%       -> CSHRMgetContShearlet.m
%       -> SLdshear.m
%       -> CSHRMgetRidges.m
%       -> CSHRMmapOrientationsToAngles.m
%       -> SLpadArray.m
%       -> yapuls.m
%   - MATLAB Image Processing Toolbox
%   - NEW* functions
%       -> parameters.m
%       -> applyThreshold.m
%       -> Ridge_Detection.m
%       -> getScale.m
%       -> getNorthCorrection.m
%
% License: GPL-3.0-only (GNU General Public License v3.0 only)
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU GPL v3.0 as published by the Free Software Foundation.
% See the LICENSE file for full text. Distributed WITHOUT ANY WARRANTY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFIED SCRIPT 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  This code extracts fracture traces from images of fractured rock. The
%  extraction is performed using the complex shearlet ridge measure. The
%  ridge measure is computed on an ensemble of shearlet parameters to
%  obtain a ridge ensemble. A highly probable ridge realization is obtained
%  from the ridge ensemble using a threshold ridge strength. The ridges
%  are converted to fractures using image processing steps such as Otsu
%  thresholding, skeletonization and polyline fitting (Prabhakaran et al., 
%  2019). This version replaces hardcoded paths with user input and adds 
% interactive image selection and scale definition for improved usability.
%--------------------------------------------------------------------------

clc
clear
close all
format compact

% Ask for the base folder where all outputs will be saved
outfolder1 = input(...
    ['> Enter the full path of the folder where results will be saved, in quotes\n' ...
     '  e.g. ''C:\\Users\\JuanLosarcos\\Desktop\\FractureDetector\\Run1''\n' ...
     'Folder = ']);
outfolder0 = char(outfolder1 + "\");

% Specific folder where the Ridge Esembles will be saved 
output_folder = char(outfolder0 + "P_ridges\");
mkdir(output_folder)

% Select a single rock-fracture image to analyze
disp(' ');
disp('> Select a rock-fracture image to analyze');
[InFileListShort, pathname] = uigetfile('*.jpg;*.tif;*.png','Image file (*.jpg,*.tif,*.png)', 'Select your image','MultiSelect','on');
    
if not(iscell(InFileListShort))
    InFileListShort = {InFileListShort};
end
    
replicatePath = repmat(cellstr(pathname),size(InFileListShort));
InFileList = strcat(replicatePath,InFileListShort);    
          
clear pathname replicatePath;
funPath = fileparts(which('W125_G63_SPO2_SL3_AL0.5_OCT_3.5_MC00.png'));
addpath(genpath(funPath));
    
disp(['You selected ', num2str(length(InFileList)) , ' imágenes']);
disp (InFileListShort);

% Draw a reference line to set the scale
for k=1:length(InFileList)   
    imageIN = imread(InFileList{k});
    imshow (imageIN);
end

% Allow up to 10 attempts to draw the reference line
for attempt = 1:10
    [pixelDist, hLine] = getScale();  
    redo = input('Redraw scale line? [Yes/No]: ','s');
    if any(strcmpi(redo,{'No','N'}))
        close;
        break;
    end
    delete(hLine);
end

% Ask user for the real-world length corresponding to that line
scaleLength = input(['\nEnter the real-world length of the drawn line\n' ...
                     'numeric only: ']);
scaleUnits = input('Enter the units for that length (e.g., cm, m, mm): ', 's');

%% Specify Complex Shearlet Input Parameters    
% Specify parameters for Complex Shearlet and ridge extraction (parameters function), then perform primary fracture detection (ridge_detection) 

for i=1:50
    close all
    % Prompt for and retrieve parameter values
    
    [no_Shearlets, RidgeExtract_Combs, onlyPositiveOrNegativeRidges, total_combs,...
    waveletEffSuppFraction_Combs, scalesPerOctave_Combs, shearLevel_Combs,...
    alpha_Combs, minContrast_Combs, offset_Combs] = parameters(outfolder0);
    
    chosenParms.WaveletEffSuppFraction  = waveletEffSuppFraction_Combs;
    chosenParms.ScalesPerOctave         = scalesPerOctave_Combs;
    chosenParms.ShearLevel              = shearLevel_Combs;
    chosenParms.Alpha                   = alpha_Combs;
    chosenParms.MinContrast             = minContrast_Combs;
    chosenParms.Offset                  = offset_Combs;

    % Run the edge‐detection stage with those settings
    [C_Edges_norm] = Ridge_Detection(InFileList,InFileListShort,no_Shearlets,RidgeExtract_Combs,onlyPositiveOrNegativeRidges,total_combs,outfolder0,output_folder);
    
    % Show the normalized edge image
    imshow(C_Edges_norm)
    % Ask whether to tweak parameters and rerun
   
    answer = input('Based on this image, would you like to adjust the parameters? [Yes/No]: ','s');
    
    fprintf('\n');
    if any(strcmpi(answer, {'No','N'}))
        break
    end
end
   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFIED SCRIPT 2 %%%%%%%%%%%%%%%%%%%%%%
%
%  This script converts ridge ensemble images into binary images based on
%  a ridge threshold. The ridge ensemble image is first normalized from
%  grayscale (0-255) to a range from (0-1). A non-linear sigmoid function
%  is then applied on the normalized image to improve contrast, so that
%  picking the ridge intensity threshold is easier. Once the threshold is 
%  set, all ridge pixels above the threshold are set to 1 and the remaining 
%  pixels to zero. The user compares the thresholded ridges with the 
%  original image to assess quality. If some fractures are missed, the
%  shearlet parameters may need adjustment (Prabhakaran et al., 2019).
%
%  Summary of Operations:  
%   ridge ensemble -> normalized ridge ensemble -> sigmoided image -> binarized ridge
%
%  New functionality: the script includes an interactive loop for threshold
%  fine-tuning (up to 50 iterations) and automatic saving of outputs in organized folders.
%
%  Uses MATLAB Image Processing Toolbox functions.
%--------------------------------------------------------------------------

% Create a folder to save the binarized image
outfolder = char(output_folder + "Binarized\");
mkdir(outfolder)

% Read the normalized edge strength image from the previous step
C_Ridges = C_Edges_norm;

% Normalize to [0,1] and display
C_Ridges_norm = mat2gray(C_Ridges);
% figure(1)
% f=imshow(C_Ridges_norm);
% title('Normalized Ridge Strength');

% Apply a sigmoid nonlinearity to enhance strong ridges
[m,n]=size(C_Ridges_norm);
C_Ridges_norm_sigmoid = C_Ridges_norm;
for i=1:m
  for j=1:n  
    if   C_Ridges_norm(i,j)~=0
      C_Ridges_norm_sigmoid (i,j) = 1 / (1 + exp((-1)*C_Ridges_norm(i,j)));
    end  
  end
  i;
end
% figure(2)
% imshow(C_Ridges_norm_sigmoid)
% title('Simoid Nonlinearity')


% Loop up to 50 times to allow the user to fine-tune the binarization threshold
for attempt = 1:50

    % Apply a threshold interactively and return the binarized image & chosen value
    
    [C_Ridges_norm_thresh, C_Ridges_norm_sigmoid, threshold] = applyThreshold(C_Ridges_norm_sigmoid);
    chosenParms.Threshold                   = threshold;
    
    % Ask if the user wants to adjust the threshold again
    answer = input('Based on this image, would you like to adjust the threshold? [Yes/No]: ', 's');
    fprintf('\n');
    if any(strcmpi(answer, {'No','N'}))
        close all
        break
    end
end

% Save the final binarized image
disp ('')
disp ('Saving image after applying threshold')
outfolder = char(output_folder + "Thresholding\");
mkdir(outfolder)
imwrite(C_Ridges_norm_thresh,strcat(outfolder,'threslholding',num2str(threshold),'.jpg'));

% After selecting the threshold, this segment creates and saves the binarized image
% Define and create the output folder for the binarized result
outfolder = char(output_folder + "Binarized\");

% Save the binarized image with the same base name as the input
imwrite(C_Ridges_norm_thresh,strcat(outfolder,'Bin_',InFileListShort{1}));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCRIPT 3 (NO MOD) %%%%%%%%%%%%%%%%%%%%%%
% This script performs the post-processing on the binarized ridges. The
% operations performed are segmentation using Otsu thresholding,
% skeletonization and polyline fitting
%
% binary ridge image -> segmentation -> skeletonization -> poyline fitting
%
% Calls functions from Geom2D Toolbox by David Legland, CoSHREM Toolbox
% by Rafael Reisenhofer and MATLAB Image Processing Toolbox
%       -> CSHRMgetOverlay.m (from CoSHREM Toolbox)
%       -> polynomialCurveSetFit.m (from Geom2D Toolbox)
%       -> drawPolynomialCurve.m (from Geom2D Toolbox)
%       -> minDistancePoints.m (from Geom2D Toolbox)
%       -> parametrize.m (from Geom2D Toolbox)
%       -> polynomialCurveFit.m (from Geom2D Toolbox)
%       -> polynomialCurvePoint.m (from Geom2D Toolbox)
%       -> polynomialCurveSetFit.m (from Geom2D Toolbox)
%
% (Prabhakaran et al., 2019).
%-------------------------------------------------------------------

% Create directories for post-processing outputs
output_folder= char(outfolder + "Segmented_Ridges\");
mkdir(output_folder)
output_folder1 = char(outfolder + "Segmented_Ridges_Overlay\");
mkdir(output_folder1);
output_folder2 = char(outfolder + "Skeletons\");
mkdir(output_folder2)
output_folder3 = char(outfolder + "Fitted_Curves\");
mkdir(output_folder3)

%%  Loop that performs the post-processing steps

for m = 1:size(InFileList,2)
    tic   
    % counter
    disp(' ');
    disp(['Performing Otsu Thresholding Segmentation (fixed parameters) for the image  ' num2str(m) ' of ' num2str(size(InFileList, 2))]);
   
    
    % read the binarized ridge image
    imageIN =C_Ridges_norm_thresh;
    %imageIN = C_Ridges_norm_sigmoid_thresh;
    % read the source image

    imageIN_1 = imread(InFileList{m});
    imageIN_1 = rgb2gray(imageIN_1); 
    imageIN_4 = imageIN;
  
    % estimating background using morphological opening
    background = imopen(imageIN,strel('disk',5));
    %imshow(background)
    
    % subtracting the background image from the original image
    imageIN_2 = imageIN - background;
    %imshow(imageIN_2)
    
    % increasing the image contrast
    imageIN_3 = imadjust(imageIN_2);
    %imshow(imageIN_3)

    bw = imbinarize(imageIN_3);
    bw = bwareaopen(bw, 10); 
    
    % identifying objects within the image
    cc = bwconncomp(bw, 4);
    
    % examining one object
    % grain = false(size(bw));
    % grain(cc.PixelIdxList{1049}) = true;
    
    % creating labels
    labeled = labelmatrix(cc);
    whos labeled;
    RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
    
    % removing isolated clusters
    Pixel_List = cc.PixelIdxList';
    for j=1:length(Pixel_List)
     Pixel_List_length(j,1)=length(Pixel_List{j,1});
    end
    imageIN_4 = labeled;
    %imshow(imageIN_4)

    Isolation_Threshold = 10;
    for j=1:length(Pixel_List)
     if Pixel_List_length(j,1)< Isolation_Threshold
      imageIN_4(imageIN_4==j)=0;
     end
    end
    %imshow(imageIN_4)
    
    j=1;
    tic
    for k=1:length(Pixel_List)
     [r,c]=find(imageIN_4==k);
     if isempty(r)~=1  
      Cluster_List{j,1}=[r c];
      j=j+1;
     end
    end
    toc
    
%     set(gca,'Ydir','normal')
%     for i=1:length(Cluster_List)
%        set(gca,'Ydir','reverse')
%        scatter(Cluster_List{i,1}(:,1),Cluster_List{i,1}(:,2),'Filled','b');
%        hold on
%     end    
    
   %  converting to image
   imageIN_4(imageIN_4>0)=255;
   % imageIN_4=mat2gray(imageIN_4);
   imageIN_1=double(255 * mat2gray(imageIN_1));
   imageIN_4=double(255 * mat2gray(imageIN_4));
   
   %  writing the segmented ridges and its overlay to the output folders
   disp ('')
   disp ('Saving image with segmented ridges')
   imwrite(imageIN_4, strcat(output_folder,'Segmented_Ridges',InFileListShort{1}));

   % creating overlay of segmented ridges on the source image and saving it
   disp ('Saving image with segmented ridges overlaid on the original image')
   overlay_plus_segmented_ridges=CSHRMgetOverlay(imageIN_1,imageIN_4);
   imwrite(overlay_plus_segmented_ridges, strcat(output_folder1,'Segmented_Ridges_Overlay',InFileListShort{1}));
  
   % Skeletonizing the segmented ridges  and
   % calculating branch points and end points for each cluster
   skelImg   = bwmorph(imageIN_4, 'thin', 'inf');
   branchImg = bwmorph(skelImg, 'branchpoints');
   endImg    = bwmorph(skelImg, 'endpoints');

   [row, column] = find(endImg);
   endPts        = [row column];

   [row, column] = find(branchImg);
   branchPts     = [row column];
   cNumBranchPoints = length(branchPts);
   
   %  writing the skeletonized ridges to the output folders
   disp(['Writing the Skeletons for Image ' num2str(m) ' out of ' num2str(size(InFileList, 2))]);  
   skelImg2 = mat2gray(skelImg);
   imwrite(~skelImg2, strcat(output_folder2,'Skeletons',InFileListShort{1}));
   
   % concatenating the clusters into one large matrix
    Clusters=Cluster_List{1,1};
    for i=2:length(Cluster_List)
     A= Cluster_List{i,1};
     Clusters=[Clusters;A];
    end

   % finding the branch points and end points associated with each cluster
   % and storing them in two additional columns
   % clusters which have no branches and end points are identified and
   % stored in no_branches and no_ends arrays
   
   y=1;
   z=1;
   for i=1:length(Cluster_List)
    Cluster_List{i,2}=intersect(branchPts,Cluster_List{i,1},'rows') ;
    Cluster_List{i,3}=intersect(endPts,Cluster_List{i,1},'rows') ;
     if isempty(Cluster_List{i,2})==1
      no_branches(y,1) = i;
      y=y+1;
     end
     if isempty(Cluster_List{i,3})==1
      no_ends(z,1) = i;
      z=z+1;
     end   
   end

   % removing clusters that have no end points. It seems these clusters
   % are very close to existing clusters and escape the bwmorph function
   if exist('Nodes','var')==1
    length_Cluster_List = length(Cluster_List);
    idx = find(no_ends==length_Cluster_List); 
    Cluster_List(no_ends(1:idx),:)=[];
   end 
   
   % fitting curves through the clusters using function from David Legland.
   % compute coeffs of each individual branch and returns a matrix of
   % labels for each fitted curve. I use this function for the time being.
   % It does not use the endPts and branchPts calculated using the
   % bwmorph call but calculates end points and branch points using bwlabel,
   % bwconncomp and regionprops functions.
   skelImg = imrotate(skelImg,180);
   [coeffs, curve_matrix] = polynomialCurveSetFit(skelImg, 5); 
     
   %writing the polynomial fit to a table and stored in the output folder
   if isempty(coeffs)~=1
     
    % Obtaining the polynomial points for each fitted curve using coeffs for
     %each curve
    for n = 1:length(coeffs)
     Poly_Points{n,1}= drawPolynomialCurve([0 1], coeffs{n});
    end
     
     
%     figure(9); imshow(~skelImg); hold on;
%         for i = 1:length(coeffs)
%            Poly_Points{i,1}= drawPolynomialCurve([0 1], coeffs{i});
%    
%             set(hc, 'linewidth', 2, 'color', 'r');
%         end

    % converting cell array to a table
    for i=1:length(Poly_Points)
      Poly_Points_Table_Header {i}= (['Polyline_' ,num2str(i)]);
    end
    Poly_Points_Table=cell2table(Poly_Points);
       
    OutFileName = InFileListShort{1};
    OutFileName = OutFileName(1:length(OutFileName)-4);
    OutFileName = strcat(output_folder3,OutFileName);
    disp(['Writing Polyline Points for Image ' num2str(m) ' out of ' num2str(size(InFileList, 2))]);
    save(OutFileName,'Poly_Points_Table');
   
   else
      disp(['Image ' num2str(i) ' is empty. No Polylines to write' ]); 
   writetable(Poly_Points_Table,OutFileName);
   end
   
   toc
end


% Load the original image
    image = imread(InFileList{1});
    [rows, columns, numberOfColorChannels] = size(image);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFIED SCRIPT 4 %%%%%%%%%%%%%%%%%%%%%%
%  
%   Final post-processing in image space (no georeferencing): rescale
%   lengths from pixel→units (getScale), flip/rotate polylines, correct
%   azimuths to geographic north (getNorthCorrection), plot/save QC figures,
%   and compute/save summary stats.
%
% Diff vs. original:
%   – Georeferencing/shapefile export removed (commented).
%   – Added interactive north correction + real-world length scaling.
%   – Produces rose plots (raw/corrected/length-weighted), length histograms,
%     fracture overlays, and an Excel summary.
%
% Dependencies:
%   MATLAB Image Processing Toolbox; Geom2D; getScale.m; getNorthCorrection.m
%   (optional: DouglasPeucker).
%
% Outputs:
%   LengthComparison.jpg, RoseComparison.jpg, Fractures.jpg,
%   SummaryStatistics.xlsx
%-------------------------------------------------------------------

% creates a list of files which have polylines...some tiles of the original
% orthomosaic (on the boundary) do not have no detected features
% creates a cell array of the fractures detected in each tile

% an aspect ratio can be applied here to adjust the polylines
aspect_ratio = 1/1;
for i=1:length(InFileList)
   imread(InFileList{i});   
   [table_size,~] = size(Poly_Points_Table);
    for j=1:table_size
     Poly_Points_Table.Poly_Points{j,1}(:,1)=Poly_Points_Table.Poly_Points{j,1}(:,1).*aspect_ratio;
     Poly_Points_Table.Poly_Points{j,1}(:,2)=Poly_Points_Table.Poly_Points{j,1}(:,2)./aspect_ratio;       
    end
    PolyLines{i,1} = table2cell(Poly_Points_Table);

end

% %% Creating an N1N2 matrix
% m=length(PolyLines{1,1});
% 
% % creating row indices matrix
% N1N2=zeros(length(PolyLines),2);
% N1N2(1,1)=1;
% N1N2(1,2)=length(PolyLines{1,1});
% for i=2:length(PolyLines)
%   N1N2(i,1)=  N1N2(i-1,2)+1;
%   N1N2(i,2)=  length(PolyLines{i,1})+N1N2(i-1,2);   
% end 

% %% Georeferencing the PolyLines
% PolyLines_Georeferenced = PolyLines;
% 
% for i=1:length(PolyLines)
% 
%     R = worldfileread(InFileList2{i}, 'planar', [938 1062]);
% 
%  % calculating the corner points of the orthotile (in lat-long )
%  Lr_x = R.XWorldLimits(1);
%  Lr_y = R.YWorldLimits(1);
% 
%  % converting the pixel values into lat-long and then setting the origin
%  % based on the corner values
%  for j=1:length(PolyLines{i,1})
%    PolyLines_Georeferenced{i,1}{j,1}(:,1)= PolyLines_Georeferenced{i,1}{j,1}(:,1).*...
%        R.CellExtentInWorldX + Lr_x;
% 
% 
%    PolyLines_Georeferenced{i,1}{j,1}(:,2)= PolyLines_Georeferenced{i,1}{j,1}(:,2).*...;   
%        R.CellExtentInWorldY + Lr_y;
% 
%  end  
%  clearvars  Lr_x Lr_y R
%  i
% end
% 
% 
% %% Rotating the Georeferenced PolyLines
% PolyLines_Georeferenced_Rotated = PolyLines_Georeferenced;
% j=1;
% for i=1:length(PolyLines_Georeferenced) 
%     %[~,R] = geotiffread(InFileList2{i});  
%     R = worldfileread(InFileList2{i}, 'planar', [938 1062]);
%     for k=1:length(PolyLines_Georeferenced{j,1})  
% %     C=PolyLines_Georeferenced{j,1};  
% 
%       % the following step flips the shapefile about the x-axis
%         %PolyLines_Georeferenced_Rotated{j,1}{k,1}(:,2) = PolyLines_Georeferenced{j,1}{k,1}(:,2)*-1 + R.YWorldLimits(1,1) + R.YWorldLimits(1,2) ; 
% 
%       % the following step flips the shapefile about the y-axis
%         PolyLines_Georeferenced_Rotated{j,1}{k,1}(:,1) = PolyLines_Georeferenced{j,1}{k,1}(:,1)*-1 + R.XWorldLimits(1,1) + R.XWorldLimits(1,2) ; 
%     end
%      j=j+1;          
%      disp(i)
% end

%% Rotating non-Georeferenced Polylines
PolyLines_Rotated = PolyLines;
j=1;
for i=1:length(PolyLines) 
    
    for k=1:length(PolyLines{j,1})  
%     C=PolyLines_Georeferenced{j,1};  

      % rotating about the x-axis
      PolyLines_Rotated{j,1}{k,1}(:,2) = PolyLines{j,1}{k,1}(:,2) *-1 + rows;
      
      % rotating about the y-axis
       PolyLines_Rotated{j,1}{k,1}(:,1) = PolyLines{j,1}{k,1}(:,1) *-1 + columns;
      
      % translating down the y-axis
      %PolyLines_Rotated{j,1}{k,1}(:,2) = PolyLines{j,1}{k,1}(:,2) - 1000;
    end  
     j=j+1;          
     %disp(i)
end

% %% Creating the ShapeFile Structure (not georeferenced)
% j=1;
% %P = PolyLines_Rotated;
% P = PolyLines;
% 
% for i=1:length(P) 
% 
%   C=P{i,1};
%   if i>1
%     m=length(C)+length(P{i-1,1});
%   end
%   z=1;
%   for k=N1N2(i,1):N1N2(i,2)  
%     [PolyLines_Shape(k).Tile_ID] = i;
% %     [PolyLines_Shape(k).Ortho_ID] = P{i,2};
% %      [PolyLines_Shape(k).BoundingBox] = [R{P{i,2},1}.LongitudeLimits(1,1) ...
% %                                         R{P{i,2},1}.LatitudeLimits(1,1); ...
% %                                         R{P{i,2},1}.LongitudeLimits(1,2) ...
% %                                         R{P{i,2},1}.LatitudeLimits(1,2)];
%     [PolyLines_Shape(k).Geometry] = 'PolyLine';
%     [PolyLines_Shape(k).Polyline_ID] = k;
%     [PolyLines_Shape(k).Polyline_in_Ortho_ID] = z;
%     % if C is a cell array
%     % [PolyLines_Shape(k).X] =  C{k-j+1,1}(:,1);
%     % if C is not a cell array
%     [PolyLines_Shape(k).X] =  C(k,1);
%     % if C is a cell array
%     % [PolyLines_Shape(k).Y] =  C{k-j+1,1}(:,2);   
%     % if C is not a cell array
%     [PolyLines_Shape(k).Y] =  C(k,2); 
%     z=z+1;
%   end 
%   j=length(PolyLines_Shape)+1;
% 
%   clearvars C; clearvars z;
% end 
% 
% %% Creating the Shapefile Structure (not georeferenced, but rotated)
% j=1;
% %P = PolyLines;
% %P = PolyLines_Rotated;
% P = P_Simplified;
% 
% for i=1:length(P) 
% 
%   C=P{i,1};
%   if i>1
%     m=length(C)+length(P{i-1,1});
%   end
%   z=1;
%   for k=N1N2(i,1):N1N2(i,2)  
%     [PolyLines_Shape(k).Tile_ID] = i;
% %     [PolyLines_Shape(k).Ortho_ID] = PolyLines{i,2};
% %      [PolyLines_Shape(k).BoundingBox] = [R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,1) ...
% %                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,1); ...
% %                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,2) ...
% %                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,2)];
%     [PolyLines_Shape(k).Geometry] = 'PolyLine';
%     [PolyLines_Shape(k).Polyline_ID] = k;
%     [PolyLines_Shape(k).Polyline_in_Ortho_ID] = z;
%     % if C is a cell array
%      [PolyLines_Shape(k).X] =  C{k-j+1,1}(:,1);
%     % if C is not a cell array
%     %[PolyLines_Shape(k).X] =  C(k,1);
%     % if C is a cell array
%      [PolyLines_Shape(k).Y] =  C{k-j+1,1}(:,2);   
%     % if C is not a cell array
%     %[PolyLines_Shape(k).Y] =  C(k,2); 
%     z=z+1;
%   end 
%   j=length(PolyLines_Shape)+1;
% 
%   clearvars C; clearvars z;
% end 
% 
% 
% 
% %% Creating the ShapeFile Structure (georeferenced and rotated)
% j=1;
% %P = PolyLines_Georeferenced_Rotated;
% %P = PolyLines_Georeferenced;
% P = P_Simplified;
% for i=1:length(P) 
% 
%   C=P{i,1};
%   if i>1
%     m=length(C)+length(P{i-1,1});
%   end
%   z=1;
%   for k=N1N2(i,1):N1N2(i,2)  
%     [PolyLines_Shape(k).Tile_ID] = i;
% %     [PolyLines_Shape(k).Ortho_ID] = PolyLines{i,2};
% %      [PolyLines_Shape(k).BoundingBox] = [R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,1) ...
% %                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,1); ...
% %                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,2) ...
% %                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,2)];
%     [PolyLines_Shape(k).Geometry] = 'PolyLine';
%     [PolyLines_Shape(k).Polyline_ID] = k;
%     [PolyLines_Shape(k).Polyline_in_Ortho_ID] = z;
%     % if C is a cell array
%      [PolyLines_Shape(k).X] =  C{k-j+1,1}(:,1);
%     % if C is not a cell array
%     %[PolyLines_Shape(k).X] =  C(k,1);
%     % if C is a cell array
%     [PolyLines_Shape(k).Y] =  C{k-j+1,1}(:,2);   
%     % if C is not a cell array
% %     [PolyLines_Shape(k).Y] =  C(k,2); 
%     z=z+1;
%   end 
%   j=length(PolyLines_Shape)+1;
%   disp(i)
%   clearvars C; clearvars z;
% end 
% 
% %% Creating ShapeFile Structures and writing shapefile for each image separately
% 
% P = P_Simplified;
% outfolder='D:\Github_Test\';
% 
% % j=1;
% for i=1:length(P)
%   tic  
%   filename=InFileListShort{i};
%   filename=filename(1:length(filename)-4);
%   outfilename=strcat(filename,'.shp');
% 
%   C=P{i,1};
% %   if i>1
% %     m=length(C)+length(PolyLines_Georeferenced_Rotated{i-1,1});
% %   end
%   z=1;
%   for k=1:length(C)  
%     [PolyLines_Shape(k).Tile_ID] = i;
% %     [PolyLines_Shape(k).Ortho_ID] = PolyLines{i,2};
% %      [PolyLines_Shape(k).BoundingBox] = [R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,1) ...
% %                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,1); ...
% %                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,2) ...
% %                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,2)];
%     [PolyLines_Shape(k).Geometry] = 'PolyLine';
%     [PolyLines_Shape(k).Polyline_ID] = k;
%     [PolyLines_Shape(k).Polyline_in_Ortho_ID] = z;
%     % if C is a cell array
%      %[PolyLines_Shape(k).X] =  C{k-j+1,1}(:,1);
%     % if C is not a cell array
%     [PolyLines_Shape(k).X] =  C{k,1}(:,1);
%     % if C is a cell array
%     %[PolyLines_Shape(k).Y] =  C{k-j+1,1}(:,2);   
%      %if C is not a cell array
%      [PolyLines_Shape(k).Y] =  C{k,1}(:,2); 
%     z=z+1;
%   end 
%   j=length(PolyLines_Shape)+1;
%   shapewrite(PolyLines_Shape,strcat(outfolder,outfilename));
%   disp(i)
%   clearvars C  PolyLines_Shape;
%   toc
% end 


% %% performing line simplification on the fitted curves using Douglas-Pecker
% % algorithm
% %epsilon = 1E-10;
% epsilon = 1E-02;
% 
% %P_Simplified = PolyLines_Georeferenced_Rotated;
% %P_Simplified = PolyLines_Georeferenced;
% %P_Simplified = PolyLines_Rotated;
% P_Simplified = PolyLines_Rotated;
% for i=1:length(P_Simplified)
%     tic
%  for j=1:length(P_Simplified{i,1})
%     Points = P_Simplified{i,1}{j,1}';
%     Simplified_Points = DouglasPeucker(Points,epsilon);
%     P_Simplified{i,1}{j,1} = Simplified_Points';
%  end
%  %disp(i)
%  toc
% end
% 
% % %% Writing shapefile
% % tic
% % %shapewrite(PolyLines_Shape,'D:\PhD\Automatic_Detection\Core_Fractures\Shape_Files\xz769_Rotated_Simplified.shp');
% % shapewrite(PolyLines_Shape,'D:\PhD\Automatic_Detection\Probabilistic_Edges\Ortho_1_107\Ortho_107_New\Shapefiles_Georef\Tile_107_Georef.shp');
% toc

%% Classifying Polylines by length

Poly_Length = zeros(length(PolyLines{1,1}),1);  % Preallocate

for i = 1:length(PolyLines{1,1})
    coords = PolyLines{1,1}{i,1};  % Get the [x y] coordinates of the i-th polyline
    polyline_length = 0;
    
    for j = 1:(size(coords,1) - 1)
        segment = [coords(j,1), coords(j,2), coords(j+1,1), coords(j+1,2)];
        polyline_length = polyline_length + Lengths2D(segment);
    end

    Poly_Length(i,1) = polyline_length;
end

% Modification of fracture length values in pixels to scaled values 
% based on the value entered by the user in the getScale function
Length_real= (Poly_Length (:,1) * scaleLength)/ pixelDist;

%% Plots

fLen = figure('Name','Length Check','NumberTitle','off');
subplot(1,2,1)
histogram(Poly_Length)


title('Raw Lengths (px)')
ylabel('Count')
xlabel('Length (px)')

subplot(1,2,2)
histogram(Length_real);
title(['Scaled Lengths (' scaleUnits ')'])
ylabel('Count')
xlabel(['Length (' scaleUnits ')'])


LenGraph = char(outfolder + "LengthComparison.jpg");
saveas (fLen,LenGraph)


%% Polylines azimuth correction

% % 1) Load & display the first image
% for k=1:length(InFileList)   
%     imageIN = imread(InFileList{k});
%     imshow (imageIN);
% end

hFigNorth = figure('Name','Draw Scale Line','NumberTitle','off');
imshow(imread(InFileList{1}));

% 2) Let user redraw up to 10× and compute deltaNorth
for attempt = 1:10
    [correctionAngle, roi] = getNorthCorrection();  
    resp = input('Redraw reference line? [Yes/No]: ','s');
    if startsWith(lower(resp),'n')
        close(hFigNorth);
        break;
    end
    delete(roi);
end

% 3) Compute and correct every raw MATLAB azimuth
nLines = numel(PolyLines_Rotated{1,1});
rawAz = zeros(nLines,1);
correctedAz = zeros(nLines,1);

for k = 1:nLines
    pts   = PolyLines_Rotated{1,1}{k,1};
    dx  = pts(end,1) - pts(1,1);
    dy  = pts(end,2) - pts(1,2);
    % raw photo angle in MATLAB coords (0 = +X/right, CCW +)           
    rawAz(k) = mod(atan2d(dy, dx),360);

  % Apply correction:
    correctedAz(k) = mod(rawAz(k) + correctionAngle, 360);
   
end

fprintf('Correcting azimuths...');

%% Build weigthed polar diagrams

% fold 0–180° and mirror to 0–360 for raw
Orient_raw   = mod(rawAz,180);
az_full_raw  = [Orient_raw; Orient_raw + 180];
% azq_raw      = round(az_full_raw/5)*5;

% fold 0–180° and mirror to 0–360 for corrected
Orient_corr  = mod(correctedAz,180);
az_full_corr = [Orient_corr; Orient_corr + 180];

binEdges = deg2rad(0:5:360);
% % Weighted: 
lengths_full = [ Length_real; Length_real ];

% 1) figure out which bin each angle falls into
binIdx = discretize(mod(deg2rad(az_full_corr),2*pi), binEdges);

% 2) sum up the lengths_full for each bin, ignoring possible NaNs
counts_w = accumarray(binIdx(~isnan(binIdx)), lengths_full(~isnan(binIdx)), [numel(binEdges)-1, 1]);

%% Plot
% set up a 1×3 polar layout
fRose = figure('Name','Rose Check','NumberTitle','off');
tiledlayout(1,3);

% raw rose
nexttile
h1 = polarhistogram(deg2rad(az_full_raw), 'BinEdges', binEdges);
h1.DisplayStyle       = 'bar';
h1.EdgeColor          = 'k';
h1.LineWidth          = 0.8;
ax = gca;
ax.ThetaDir          = 'clockwise';
ax.ThetaZeroLocation = 'top';
title('Raw Azimuths');

% corrected rose
nexttile
h2 = polarhistogram(deg2rad(az_full_corr), 'BinEdges', binEdges);
h2.DisplayStyle       = 'bar';
h2.EdgeColor          = 'k';
h2.LineWidth          = 0.8;
ax = gca;
ax.ThetaDir          = 'clockwise';
ax.ThetaZeroLocation = 'top';
title('Corrected Azimuths');

% weighted rose
nexttile
h3 = polarhistogram( 'BinEdges', binEdges, ...
                    'BinCounts', counts_w );
h3.DisplayStyle       = 'bar';
h3.EdgeColor          = 'k';
h3.LineWidth          = 0.8;
ax = gca;
ax.ThetaDir           = 'clockwise';
ax.ThetaZeroLocation  = 'top';
title('Weighted by Length');

RoseGraph = char(outfolder + "RoseComparison.jpg");
saveas (fRose,RoseGraph)

%% Displaying all polylines 
figure(4)
imshow(image); 
hold on;
for i=1:length(PolyLines_Rotated{1,1})
    scatter(PolyLines_Rotated{1,1}{i,1}(:,1),PolyLines_Rotated{1,1}{i,1}(:,2),'.')
    grid on
    pbaspect([1 1 1])
    hold on   
end  
title('Lineas de Fracturas')

drawnow;
fFract = figure(4);
FractGraph = char(outfolder + "Fractures.jpg");
saveas (fFract,FractGraph)


%% ---------------------- SUMMARY STATISTICS ----------------------

% Number of fractures
numFractures = numel(Length_real);

% Length stats (in your chosen units)
meanLen = mean(Length_real);
medLen  = median(Length_real);
stdLen  = std(Length_real);
minLen  = min(Length_real);
maxLen  = max(Length_real);

% Build corrected‐azimuth counts in 5° bins

binEdges180 = deg2rad(0:5:180);
centers180 = rad2deg(binEdges180(1:end-1) + diff(binEdges180)/2);
counts180 = histcounts(deg2rad(Orient_corr), 'BinEdges', binEdges);

% 3) assign each fracture to one of those bins
binIdxW    = discretize(deg2rad(Orient_corr), binEdges180);

% 4) sum lengths in each bin
counts_wW   = accumarray( binIdxW(~isnan(binIdxW)), ...
                        Length_real(~isnan(binIdxW)), ...
                        [numel(binEdges180)-1,1] );



% Pick top 5 bins by count
K = 10;
topOrientsW = zeros(K,1);
[sortedW, idxW] = sort(counts_wW,'descend');
topCountsW  = sortedW(1:K);
topBinsW    = idxW(1:K);  
%[sortedCounts, idx] = sort(counts180,'descend');
% topCounts  = sortedCounts(1:K);
% topBins    = idx(1:K);           


% Print everything
fprintf('\n====== RESULTS for %s ======\n', InFileListShort{1});
fprintf('Total fractures: %d\n\n', numFractures);

fprintf('Length (%s):\n', scaleUnits);
fprintf('  Mean   = %.4f\n', meanLen);
fprintf('  Median = %.4f\n', medLen);
fprintf('  Std    = %.4f\n', stdLen);
fprintf('  Min    = %.4f\n', minLen);
fprintf('  Max    = %.4f\n\n', maxLen);

fprintf('\n====== WEIGHTED RESULTS for %s ======\n', InFileListShort{1});
fprintf('Top %d 5° bins by *summed fracture length* – mean orientation (° from North):\n',K);
for k = 1:K
  low        = (topBinsW(k)-1)*5;
  high       =  topBinsW(k)*5;
  inBin      = Orient_corr >= low & Orient_corr < high;
  topOrientsW(k) = mean( Orient_corr(inBin) );
  fprintf('  #%2d: %5.1f°  (%.0f %s)\n', k, topOrientsW(k), topCountsW(k), scaleUnits);
end
fprintf('=================================\n');

% fprintf('Top %d 5° bins by count – showing each bin’s mean orientation (° from North):\n', K);
% for i=1:K
%     % reconstruct the 0–180° bin edges from the bin index
%     lowEdge  = (topBins(i)-1)*5;
%     highEdge =  topBins(i)*5;
%     % pick the actual input angles in that bin
%     inBin     = Orient_corr >= lowEdge & Orient_corr < highEdge;
%     % compute their mean
%     avgAngle  = mean( Orient_corr(inBin) );
%     topOrients(i) = avgAngle;
%     fprintf('  #%d: %5.1f°  (%d fractures)\n', ...
%             i, avgAngle, topCounts(i));
% end
% fprintf('=================================\n');

% — Save summary stats to Excel —
% build a one‐row table of everything
% Tstats = table( ...
%     {InFileListShort{1}}, ...
%     numFractures, ...
%     meanLen, medLen, stdLen, minLen, maxLen, ...
%     topOrients(1), topCounts(1), ...
%     topOrients(2), topCounts(2), ...
%     topOrients(3), topCounts(3), ...
%     topOrients(4), topCounts(4), ...
%     topOrients(5), topCounts(5), ...
%     'VariableNames',{ ...
%       'ImageFilename', ...
%       'TotalFractures', ...
%       'MeanLength','MedianLength','StdLength','MinLength','MaxLength', ...
%       'DominantOrient1','Count1', ...
%       'DominantOrient2','Count2', ...
%       'DominantOrient3','Count3', ...
%       'DominantOrient4','Count4', ...
%       'DominantOrient5','Count5' } ...
% );

% — Save weighted summary stats to Excel —
Tstats = table( ...
    {InFileListShort{1}}, ...
    numFractures, ...
    meanLen, medLen, stdLen, minLen, maxLen, ...
    topOrientsW(1), topCountsW(1), ...
    topOrientsW(2), topCountsW(2), ...
    topOrientsW(3), topCountsW(3), ...
    topOrientsW(4), topCountsW(4), ...
    topOrientsW(5), topCountsW(5), ...
    topOrientsW(6), topCountsW(6), ...
    topOrientsW(7), topCountsW(7), ...
    topOrientsW(8), topCountsW(8), ...
    topOrientsW(9), topCountsW(9), ...
    topOrientsW(10), topCountsW(10), ...
    'VariableNames',{ ...
      'ImageFilename', ...
      'TotalFractures', ...
      'MeanLength','MedianLength','StdLength','MinLength','MaxLength', ...
      'DominantOrient1','Count1', ...
      'DominantOrient2','Count2', ...
      'DominantOrient3','Count3', ...
      'DominantOrient4','Count4', ...
      'DominantOrient5','Count5',...
      'DominantOrient6','Count6',...
      'DominantOrient7','Count7',...
      'DominantOrient8','Count8',...
      'DominantOrient9','Count9',...
      'DominantOrient10','Count10',...
      } ...
);
Tstats.WaveletEffSuppFraction = waveletEffSuppFraction_Combs;
Tstats.ScalesPerOctave        = scalesPerOctave_Combs;
Tstats.ShearLevel             = shearLevel_Combs;
Tstats.Alpha                  = alpha_Combs;
Tstats.MinContrast            = minContrast_Combs;
Tstats.Offset                 = offset_Combs;
Tstats.Threshold              = threshold;

% save into the same base folder as your roses & histograms
statsFile = fullfile(outfolder0, 'SummaryStatistics.xlsx');
writetable(Tstats, statsFile);

fprintf('Summary statistics saved to:\n  %s\n', statsFile);