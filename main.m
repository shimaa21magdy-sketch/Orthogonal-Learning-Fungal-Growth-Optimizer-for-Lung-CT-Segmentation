clear all 
clc

% Load details of the selected benchmark function
% [lb, ub, dim, fobj] = Get_Functions_details(Function_name);
% Define the necessary parameters
fobj = @Otsu;  % Example objective function: sum of squares
lb = 1;  % Lower boundary for the search space
ub = 256;   % Upper boundary for the search space
level = 2;  % Parameter for the objective function (adjust as needed)
dim = level;  % Dimensionality of the problem
N = 50;    % Population size
MaxIter = 100;
I = imread('Dataset/IMG1.jpg');
I = rgb2gray(I);
MaxFEs = 1000 * dim;  % Maximum number of function evaluations;
[n_countR, x_valueR] = imhist(I(:,:,1));
Nt = size(I, 1) * size(I, 2); % Total number of pixels in the image (rows * cols)

Lmax = 256;   % Maximum number of levels in the image = 256 -> (0-255)
for i = 1:Lmax
    % Grayscale image
    probR(i) = n_countR(i) / Nt;
end
N_thresh = level;

bestOverallScore = -inf;  % Initialize the best score
bestThreshold = [];       % Initialize the best threshold

% Arrays to store performance metrics for each run
psnrValues = zeros(1, 30);
fsimValues = zeros(1, 30);
ssimValues = zeros(1, 30);
qualityRValues = zeros(1, 30);

% Repeat the optimization process 30 times
for run = 1:30
    fprintf('Running optimization iteration %d of 30...\n', run);
    
[Gb_Fit,Gb_Sol,cg_curve]=OLFGO(N,MaxIter,ub,lb,dim,fobj,I, level, Lmax, probR);
if Gb_Fit > bestOverallScore
        bestOverallScore = Gb_Fit;
        bestThreshold = Gb_Sol;
        bestCgCurve = cg_curve;
    end

    % Generate output image using the best threshold for this run
    Iout = imageGRAY(I, bestThreshold);

    % Calculate performance metrics for this run
    psnrValues(run) = PSNR(I, Iout);
    fsimValues(run) = FeatureSIM(I, Iout);
    ssimValues(run) = ssim(I, Iout);
    qualityRValues(run) = imageQualityIndex(I(:,:,1), Iout(:,:,1));
end

% Calculate the average of the metrics after 30 runs
avgPSNR = mean(psnrValues);
avgFSIM = mean(fsimValues);
avgSSIM = mean(ssimValues);
avgQualityR = mean(qualityRValues);

% Display the results after 30 runs
fprintf('Best Score after 30 runs: %.4f\n', bestOverallScore);
fprintf('Best Threshold found: ');
disp(bestThreshold);

% Display the average performance metrics
fprintf('Average PSNR: %.4f\n', avgPSNR);
fprintf('Average FSIM: %.4f\n', avgFSIM);
fprintf('Average SSIM: %.4f\n', avgSSIM);
fprintf('Average Image Quality Index: %.4f\n', avgQualityR);

% Generate the output image using the best threshold
Iout = imageGRAY(I, bestThreshold);

% Display the output image and original image
% Plot the original image on the left
% subplot(1, 2, 1); % 1 row, 2 columns, 1st position
% imshow(I);
% title('Original Image');
% 
% % Plot the optimized image on the right
% subplot(1, 2, 2); % 1 row, 2 columns, 2nd position
imshow(Iout);
title('Optimized Image');

% Plot convergence curve
figure;
plot(bestCgCurve);
title('Convergence Curve');
xlabel('Iterations');
ylabel('Best Fitness Value');
