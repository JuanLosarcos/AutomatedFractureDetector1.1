%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: parameters
% Modified by: Juan Manuel Losarcos (UIC), based on Prabhakaran et al. (2019)
%
% Description:
%   Prompts the user to input shearlet and ridge extraction parameters.
%   Builds multiple shearlet systems and saves them to the output folder.
%   Returns all combinations for use in fracture detection workflows.
%
% Input:
%   - outfolder0: folder to save shearlet systems (.mat files)
%
% Output:
%   - Parameter arrays and counters for shearlet/ridge combinations
%
% Dependencies:
%   - CoSHREM toolbox (CSHRMgetContRidgeSystem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [no_Shearlets, RidgeExtract_Combs, onlyPositiveOrNegativeRidges, total_combs,...
    waveletEffSuppFraction_Combs, scalesPerOctave_Combs, shearLevel_Combs,...
    alpha_Combs, minContrast_Combs, offset_Combs] = parameters(outfolder0)
   
%-- Step 1: Ask user for shearlet-building parameters --
    disp(' ');
    disp('> Enter your desired values in square brackets [ ] for each parameter:');
    disp('  These parameters control the construction of the Shearlet System for fracture detection.');
    disp('  NOTE:');
    disp('    WavelengthEffSuppFraction (Divisor): higher values detect thinner fractures (min 2 values, >1 ).');
    disp('      e.g. [5 8 12]');
    disp('    ScalesPerOctave: higher values detect a wider range of fracture sizes (min 2 values, >0).');
    disp('      e.g. [1 2 3]');
    disp('    ShearLevel: smaller values improve straight-fracture resolution (single integer >= 2).');
    disp('      e.g. [3]');
    disp('    Alpha: lower values detect thinner fractures (one or more values).');
    disp('      e.g. [0.5]');
    disp(' ');

    waveletEffSuppFraction_Combs   = input('WavelengthEffSuppFraction (Divisor) = ');
    scalesPerOctave_Combs  = input('ScalesPerOctave = ');
    shearLevel_Combs       = input('ShearLevel = ');
    alpha_Combs            = input('Alpha = ');
    octaves = 3.5;
    tic

    % Generate all combinations of shearlet parameters
    [ca, cb, cc, cd] = ndgrid(waveletEffSuppFraction_Combs,scalesPerOctave_Combs,shearLevel_Combs,alpha_Combs);
    Shearlet_Combs =[ca(:),cb(:),cc(:),cd(:)];
    length_Shearlet_Combs = length(Shearlet_Combs);

    disp(['Total shearlet combinations: ', num2str(length_Shearlet_Combs)]);
    toc

   
    clearvars ca cb cc cd r
    disp('')

    %-- Step 2: Compute and save each shearlet system --
    % specify image dimensions (images are always regular 2D matrices) and the
    % resulting shearlet systems are also built for a specific image size
    % shearlet systems can be quite large for large image sizes and you may
    % run out of memory quite easily. Hence it is advisable to keep the image
    % sizes less than 1000 x 1000 pixels. 

    rows = 1000;  
    cols = 1000;

    tic
    for i=1:length(Shearlet_Combs)
        tic
        waveletEffSupp = ceil(rows/Shearlet_Combs(i,1));
        gaussianEffSupp = ceil(waveletEffSupp/2);
        scalesPerOctave = Shearlet_Combs(i,2);
        shearLevel = Shearlet_Combs(i,3);
        alpha = Shearlet_Combs(i,4);
        scales = scalesPerOctave*octaves;

        % build the continuous ridge/shearlet system
        %shearletSystem = CSHRMgetContRidgeSystem(rows,cols,ceil(rows/Shearlet_Combs(i,1)),ceil(ceil(rows/Shearlet_Combs(i,1))/2),Shearlet_Combs(i,2),Shearlet_Combs(i,3),Shearlet_Combs(i,4));
        shearletSystem = CSHRMgetContRidgeSystem(rows,cols,waveletEffSupp,gaussianEffSupp,scalesPerOctave,shearLevel,alpha,scales);
        save(strcat([outfolder0,'shearletSystem',num2str(i),'.mat']),'shearletSystem','-v7.3');
        disp(['> Calculating complex shearlet system based on the given inputs for the number of combinations: ', num2str(i)]);
        clearvars waveletEffSupp gaussianEffSupp scalesPerOctave shearLevel alpha scales shearletSystem
    end    
    toc

    %-- Step 3: Ask user for ridge-extraction parameters --
    disp(' ');
    disp('> Enter values in square brackets [ ] for ridge extraction:');
    disp('  NOTE:');
    disp('    Contrast: higher values detect high-contrast features and reduce noise.');
    disp('      e.g. [10 12 15]');
    disp('    Offset: smaller values reduce noise and improve connectivity of long fractures (<0.3 to avoid MEX crashes).');
    disp('      e.g. [0.001 0.01 0.05 0.1 0.2]');
    disp(' ');

    %% Setting the range for ridge extraction parameters
    % the ridge detection depends on both shearlet system parameter variation
    % (which we achieve by using multiple shearlet systems) and ridge parameter
    % variations. The total number of ridge realizations is the product of
    % 
    onlyPositiveOrNegativeRidges = -1;  % +1 -> positive ridges
                                    % -1 -> negative ridges

    minContrast_Combs = input('Contrast = ');
    offset_Combs      = input('Offset = ');
    [ce, cf] = ndgrid(minContrast_Combs, offset_Combs);
    RidgeExtract_Combs =[ce(:),cf(:)];

    % Shearlet Systems 
    no_Shearlets = length(Shearlet_Combs); 
    
    % total number or ridge realizations 
    total_combs = length(RidgeExtract_Combs) * no_Shearlets;
    disp(['total ridge realizations ', num2str(total_combs)]);
    
    % return 

end
