%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: Ridge_Detection
% Modified by: Juan Manuel Losarcos (UIC), based on Prabhakaran et al. (2019)
%
% Description:
%   For each image in the input list, applies a set of precomputed shearlet 
%   systems and ridge extraction parameter combinations to compute ridge 
%   responses. All ridge realizations are summed and normalized into a 
%   probabilistic ridge map, resized to the original image dimensions.
%
% Inputs:
%   - InFileList: full paths of images to process
%   - InFileListShort: filenames without path for output naming
%   - no_Shearlets: number of precomputed shearlet systems
%   - RidgeExtract_Combs: matrix of [minContrast, offset] combinations
%   - onlyPositiveOrNegativeRidges: 1 or -1 for ridge polarity
%   - total_combs: total number of shearlet × ridge parameter combinations
%   - outfolder0: path where shearlet systems are stored
%   - output_folder: path to save the resulting ridge maps
%
% Output:
%   - C_Edges_norm: normalized summed ridge response for the last image
%
% Dependencies:
%   - CoSHREM toolbox functions: CSHRMsheardec, CSHRMgetRidges
%   - Requires shearlet systems generated at 1000x1000 resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C_Edges_norm] = Ridge_Detection(InFileList,InFileListShort,no_Shearlets,RidgeExtract_Combs,onlyPositiveOrNegativeRidges,total_combs,outfolder0,output_folder)


m = 1000;  % rows specified when the shearlet system was built
n = 1000;  % columns specified when the shearlet system was built
    % Note: If the dimensions of the image to be processed are not the same as
    % that of the rows and cols specified when shearlet system was built, the 
    % code doesn't work. The 70 shearlet systems that were used were set for
    % 1000 x 1000 pixel images. If the dimensions don't match, either resize 
    % the raw images (easier) or generate shearlets corresponding to raw image
    % size (takes more time if sizes of all raw images are different)

num_images = length(InFileList); % number of images that are selected



% Loop runs over each image
 for k=1:length(InFileList)
     
  counter = 0;
  tic   
  imageIN = imread(InFileList{k});
  
  % if images have geotiff tags then first 3 channels are extracted 
  % this statement can be commented if all images are grayscale
  %imageIN = rgb2gray(imageIN(:,:,1:3));  
  
  % converting image to grayscale
  imageIN = rgb2gray(imageIN);
  
  % storing the original size of the image     
  [imageIN_Xpixels,imageIN_Ypixels]=size(imageIN);
  
  outfilename = InFileList{k};
  outfilename = InFileListShort{k};
  
  % resizing the image to the rows and columns of shearlet systems
  imageIN = imresize(imageIN,[m n]);
  C_Edges = zeros(m,n,'double');
  
 for i=1:no_Shearlets % number of shearlet systems
     
    tic
     % loading shearlet systems from the shearlet folder
     load(strcat(outfolder0,'shearletSystem',num2str(i),'.mat'));
     disp(['Loaded Shearlet System Number: ',num2str(i)]);
    toc
    
    tic
     % creating coefficient matrix for each shearlet (CoSHREM function)
     coeffs = CSHRMsheardec(imageIN,shearletSystem);
     disp(['Created Coeffs Matrix for Shearlet ',num2str(i)]);
    toc
   
  numRidgeComb = size(RidgeExtract_Combs,1);

  for j=1:numRidgeComb % number of ridge realizations per shearlet
            
    tic
    % computing ridge image (CoSHREM function)
    [ridges,~] = CSHRMgetRidges(imageIN,shearletSystem,RidgeExtract_Combs(j,1),RidgeExtract_Combs(j,2),onlyPositiveOrNegativeRidges,coeffs);
    toc

    counter = counter + 1;
    
    disp(['Image: ',num2str(k),' ,Shearlet System:', num2str(i),'; Ridge Combination: ',  num2str(counter), ' of ',num2str(total_combs) ]);   

    % In order to check the effects of each ridge realization, you may want
    % to save each ridge to see the effects of each. Can comment these two
    % statements if each intermediate ridge is not required to be saved
    
        %ridges2 = imresize(ridges,[imageIN_Xpixels,imageIN_Ypixels]);
        %imwrite(ridges2, strcat([output_folder,'P_Ridges_',num2str(i),'_','MC',num2str(RidgeExtract_Combs(j,1)),'_Off_',num2str(RidgeExtract_Combs(j,2)),InFileListShort{k}]));


    % summation of computed ridge realization
    C_Edges = C_Edges + ridges;
     
    toc
     clearvars imageOUT ridges
     
  end  
 clearvars shearletSystem coeffs 
 end
 
 C_Edges_norm = mat2gray(C_Edges);
 
 % the summed up ridge is resized to the original image dimension
 C_Edges_norm = imresize(C_Edges_norm,[imageIN_Xpixels,imageIN_Ypixels]);
 
 % writing the ridge ensemble to folder
 imwrite(C_Edges_norm,strcat(output_folder,'P_Ridges_',outfilename));
 disp (['' newline ''])
 disp(['Writing Probabilistic Ridges Map for Image: ',num2str(k), ' of ', num2str(num_images)]);
%  clearvars C_Edges_norm C_Edges imageIN outfilename ridges2
 toc
 
 end
end