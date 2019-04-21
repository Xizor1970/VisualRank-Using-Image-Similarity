%%
% Author: Hendrik A. Dreyer
% Due Date: 19 August 2018

%%
%%%%%%%%%%
% TASK 1 %
%%%%%%%%%%
% 1. Using the theory above, rank all 1400 MPEG7 shape images using VisualRank. 
% a. Load the ‘sim.mat’ data file into MATLAB. 
% This file contains the similarity matrix S for the entire image dataset,
% as well as a cell array of corresponding filenames. (0.5 marks) 

M = load('sim.mat');

%%
% b. From S, we can form the adjacency matrix A. Start by setting A=S,
% then remove elements from A to leave visual hyperlinks between only those
% images that are positively correlated. (1 mark) 

A = M.S;
A(A<0) = 0;

%%
% c. Finally, remove the elements from A that correspond to loop edges in
% the visual similarity graph. (1 mark)

A(A==1) = 0;

%%
% d. Use A to create the hyperlink matrix H. (0.5 marks)

% Normalise the adjacency matrix
H = A./sum(A);
% Then, transpose in order to have each row sum to 1
H = H';

%%
% e. Form the random jump matrix J. (0.5 marks) 

% 1400 images, therefore divide each entry by 1400
J = ones(1400)/1400;

%%
% f. Given a damping factor of ? = 0.85, create the modified visual
% hyperlink matrix Htilde. (0.5 marks) 
% ?H = ?H + (1 ? ?)J 

Htilde = 0.85*H + (1 - 0.85)*J;

%%
% g. Find the VisualRank vector r for all 1400 images. (3 marks) 

% Create a vector with initial eqaul probability of 1/1400 for each image
V = ones(1, 1400) / 1400;

% Create a matrix to store the final probabilities for each image
VFin = ones(1400) - 1;

% Iterate as per the VisualRank alg. 
% I've chosen a 100 iterations (bit of an overkill), but at least the 
% probabilities have settled by then

for i = 1:100
    VFin(i,:) = V;
    V = V*Htilde;
end

%Read the final probabilities
r = VFin(100,:);

%%
% h. Given this VisualRank vector, which image is ranked highest? 
% Display this image. (1 mark) 

% Find the maximim probability with its associated index.
[max_num, max_idx]=max(r(:));
[X,Y]=ind2sub(size(r),find(r==max_num));
% max_num = 0.0011
% Y = 713

most_viewed_file = M.filenames(Y);
% most_viewed_file = {'device9-20.png'}

[I,map] = imread(M.filenames{Y});
imshow(I,map);
title(['Highest Ranked Image - ', M.filenames{Y}]);


%%
%%%%%%%%%%
% TASK 2 %
%%%%%%%%%%
% 2. Rank the 20 heart-shaped images (indices 81 through 100) using
% VisualRank. 

% a. Make a smaller 20x20 adjacency matrix by indexing the necessary rows
% and columns of the full adjacency matrix from above. (1 mark) 
A20 = A(81:100,81:100);  % A is loaded from the sim.mat file above

%%
% b. Form the corresponding visual hyperlink matrix and find the 
% VisualRank vector for all 20 heart images. (2 marks) 

% Normalise the adjacency matrix
H20 = A20./sum(A20);

% Then, transpose in order to have each row sum to 1
H20 = H20';

% Form the random jump matrix J. 
J20 = ones(20)/20;

% Given a damping factor of ? = 0.85, create the modified visual
% hyperlink matrix H20tilde. 
H20tilde = 0.85*H20 + (1 - 0.85)*J20;

% Create a vector with initial eqaul probability of 1/20 for each image
V20 = ones(1, 20) / 20;

% Create a matrix to store the final probabilities for each image
V20Fin = ones(20) - 1;

% Iterate as per the VisualRank alg. 
% I've chosen a 100 (again) iterations 
for j = 1:100
    V20Fin(j,:) = V20;
    V20 = V20*H20tilde;
end

%Read the final probabilities
r20 = V20Fin(100,:);

% Visaul rank vector for the 20 heart images
%[0.0533427139113301,0.0479909572053565,0.0521449544317273,
% 0.0489477426126017,0.0513904692820485,0.0532249746650137,
% 0.0505492485963589,0.0472977583716867,0.0514912398048768,
% 0.0493367587774284,0.0510670925417854,0.0487639269372759,
% 0.0452888701951146,0.0524531067703493,0.0515672124054859,
% 0.0456319948596996,0.0521500344235125,0.0514203560984794,
% 0.0489371746018224,0.0470034135080401]

%%
% c. Given this ranking, which heart-shaped image is most representative 
% of the group of heart-shaped images? Display this image. (1 mark)

% Find the maximim probability with its associated index.
[max_num20, max_idx20]=max(r20(:));
[X20,Y20]=ind2sub(size(r20),find(r20==max_num20));
% max_num20 = 0.0533
% Y20 = 1
% Add offset 81 into A matric therefore,
Y20 = Y20 + 80;

most_viewed_file20 = M.filenames(Y20);
% most_viewed_file20 = {'Heart-1.png'}

%Display image
[I20,map20] = imread(M.filenames{Y20});
imshow(I20,map20);
title(['MOST representative heart shaped image - ', M.filenames{Y20}]);

%%
% d. Similarly, which heart-shaped image is the least heart-shaped
% (according to VisualRank)? Also display this image. (1 mark)
% Find the maximim probability with its associated index.

[min_num20, min_idx20]=min(r20(:));
[Xmin20,Ymin20]=ind2sub(size(r20),find(r20==min_num20));
% min_num20 = 0.0533
% Ymin20 = 1
% Add offset 81 into A matric therefore,
Ymin20 = Ymin20 + 80;

least_viewed_file20 = M.filenames(Ymin20);
% least_viewed_file20 = {'Heart-20.png'}

%Display image
[Imin20,mapmin20] = imread(M.filenames{Ymin20});
imshow(Imin20,mapmin20);
title(['LEAST representative heart shaped image - ', M.filenames{Ymin20}]);

%%
%%%%%%%%%%
% TASK 3 %
%%%%%%%%%%
% 3. Pretend you are a search engine that has been queried for an 
% image search of the following pentagon-looking shape that is present
% in the dataset as ‘device6-18.png’ (at index 650):

% a. We first need to create a graph from the similarity adjacency matrix A
% created in Task 1. However, we cannot use this matrix directly, as it 
% contains edge weights that are larger when images are more similar, 
% corresponding to a further distance. We would instead like images that 
% are similar to each other to be nearer in the graph. This can be achieved 
% simply forming a new adjacency matrix for which the elements are the 
% reciprocal of those in A. Do this, and form the graph G corresponding to 
% this new adjacency matrix. (3 marks)

% Calc reciprocal
Anew = 1./A; % Use A from Task 1-c.

% Remove all Inf and replace with 0
Anew(~isfinite(Anew)) = 0;

% Create graph G
%G = graph(Anew~=0);
G = graph(Anew);

%%
% b. Using MATLAB’s nearest function, find the 10 images nearest to
% ‘device618.png’ in the graph G. (2 marks) 

d = 10;
ind = 0;
%find index of 'device6-18.png'
for k = 1:1400
    if(strcmp(M.filenames{k}, 'device6-18.png') == 1)
        idx = k; % index = 650
        break
    end
end

% Calculate the nearest nodeIDs
[nodeIDs,dist] = nearest(G,idx,d);

% Copy 10 nearest node IDs
nearest10 = nodeIDs(1:10,:); %[645;654;657;647;653;649;644;643;655;656]

% Lookufilenames of nearest nodeIDs
nearestNames = strings(10,1);
for m = 1:10
    nearestNames(m,1) = M.filenames{nearest10(m)};
end

% The 10 nearest images to device6-18.png
% ["device6-13.png";"device6-3.png";"device6-6.png";"device6-15.png";
% "device6-20.png";"device6-17.png";"device6-12.png";"device6-11.png";
% "device6-4.png";"device6-5.png"]

%%
% c. Finally, rank these 10 nearest images using VisualRank and display
% them on a figure in their ranked order. This is the final search result. 
% (2 marks) 

% From the ranked vector r in Task 1-g.
%r = VFin(100,:);
% Find ranking number for each image in r
finalSearch = zeros(10,2);

for n = 1:10
    finalSearch(n,1) = r(nearest10(n));
    finalSearch(n,2) = nearest10(n)
end

% Sort the nearest 10
finalSearch = sortrows(finalSearch,1,'descend');

orderedFinalSearchNames = strings(1,10);

% Plot the final search
for p = 1:10
    [V,W] = imread(M.filenames{finalSearch(p,2)});
    orderedFinalSearchNames(1,p) = M.filenames{finalSearch(p,2)};
    subplot(10,1,p), imshow(V,W);
    title(['Rank # ',num2str(p),' - ', M.filenames{finalSearch(p,2)}]);
end

% Ordered nearest 10 images to ‘device618.png’
%["device6-20.png","device6-15.png","device6-6.png","device6-3.png",
% "device6-13.png","device6-5.png","device6-4.png","device6-12.png",
% "device6-17.png","device6-11.png"]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THE END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














