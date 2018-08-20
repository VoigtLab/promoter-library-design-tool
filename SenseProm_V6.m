% function SenseProm_V2

clear all; close all; clc;
tic
% This script will take input sequences of operators and promoters and
    % generate sequences for chip-synthesized oligos of a defined length.

% The program functions in modules, each of which independently perform a different design:
  %1. Assembles the promoters into background sequences, such that each background sequence is paired with each promoter. #total sequences = #promoters * #background seqs
  %2. Incorporates selected operators into each of the sequences from 1 above into every position.
    %#total sequences = #seqs from 1 * #operators to place * #positions each operators can be placed
  %3. Places selected operators in designated positions, taking into account "wobble" positions
  %4. Places  operators combinatorially into each of the designated positions, not accounting for wobble

% Before using, need to uncomment large sections of the code below (!!)

%Data input: Load an excel/csv file containing:
%   in Column 1: Sequence of each background sequence; must be 85bp if final length is 150bp
%   in Column 2: Sequence of the target promoter; must be 35bp in length
%   in Column 3: location of the -10 site in the promoter (counting from right (Transcription Start Site))
%   in Column 4: length of space between -10 and -35 (usually 17 bp for E coli)
%   in Column 5: Sequence of each operator
%   in Column 6: "1" or "0" input determines whether operator is used to saturate promoter to -101
%   in Column 7: "core" locations of operators rel to TSS, number string must end with a comma after the last number; MUST be from right(downstream-most) to left (upstream-most)
%   in Column 8: the number of operators  in column 7 to be combinatorially assembled; max of 4; ...
                 %must be in order from closest to TSS to farthest
                 %(e.g.'-41, -61, -71' NOT '-51, -41'); single numbers must
                 %be followed by a comma; NOT USED IN THIS ITERATION
%   in Column 9: the amount of "wobble" for every site in column 7; i.e.: input '3' will place single operator in all sites +/- 3bp with respect to
                 %each location in column 7. +/-3 covers both sides of the DNA strand (1 helical turn ~10bp), greater chance of functional "looping"
%   in Column 10: the sequence of the fw PCR primer recognition sequence
%   in Column 11: the sequence of the rev PCR primer recognition sequence
%   in Column 12: the desired location of the TSS relative to the end of the background sequence (e.g. -30 places +1 site -30 from the end of the backgound seq)
%   in Column 13: 1 or 0 input whether to integrate these operators combinatorially
%   in Column 14: 1 or 0 input whether to add "looping" feature to the  promoter saturation section

 %% Enter the inputs; reads excel file into 'input'

names = 'GluProm';  %enter names of output oligos

[a,b,input] = xlsread('SensePromInput_GluSensor_V1.xls');

bg = input(:,1); %makes column a separate cell array
bg(cellfun(@(x) any(isnan(x)),bg)) = [];  %removes NaN from column

prom = input(:,2); %makes column a separate cell array
prom(cellfun(@(x) any(isnan(x)),prom)) = [];  %removes NaN from column

neg10 = input(:,3); %makes column a separate cell array
neg10(cellfun(@(x) any(isnan(x)),neg10)) = []; %removes NaN from column

spacer = input(:,4); %makes column a separate cell array
spacer(cellfun(@(x) any(isnan(x)),spacer)) = []; %removes NaN from column

left = input(:,10); %makes column a separate cell array
left(cellfun(@(x) any(isnan(x)),left)) = []; %removes NaN from column

right = input(:,11); %makes column a separate cell array
right(cellfun(@(x) any(isnan(x)),right)) = []; %removes NaN from column

TSS = input(:,12); %makes column a separate cell array
TSS(cellfun(@(x) any(isnan(x)),TSS)) = []; %removes NaN from column
TSS = cell2mat(TSS);


%% First, concatenate the background sequences with each promoter.
%inputs are the background sequence, the promoter sequence, and
count = 1;
for i = 1:length(bg) %ie for each background sequence, ...
    for j = 1:length(prom) % and for each promoter sequence...
        bgseq{count,1} = bg{i}(1:end - length(prom{j}) + TSS); %pull out the left side of bg and put it in the first column of bgseq charmap; 10 is arbitrary lenght of 5'UTR
        bgseq{count,2} = prom{j}; %pull out the jth promoter seuquence and put it in the second column of bgseq charmp
        bgseq{count,3} = bg{i}(end +TSS +1:end); %pull out the right side of the background nad put it in third column of bgseq; 10 is arbitrary length chosen as part of 5'UTR to insulate promoter
        bgcat{count,1} = strcat(bgseq{count,1},bgseq{count,2},bgseq{count,3}); %concatenate the seqs of bgseq
        bgcat{count,2} = neg10{j}; %put corresponding -10 site info in second column of bgcat
        bgcat{count,3} = spacer{j}; %put corresponding spacer info in second column of bgcat
        count = count+1;
    end
end
clear i; clear j; clear count;

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
% % % Place an operator every position (from -35:-108ish) in a bg/promoter aggregate
% % % this assumes the operator is an activator; the sixth column in the excel file contains a logic input that designates which operators are used in this module
% satcat = {};
% for j = 1:length(bgcat) %for each concatenated background promoter sequence,...  assumes operator activates
%     for i = 1:length(input(:,5)) %for each operator in 'input'...
%         if input{i,6} == 1 %if that operator should be used to saturate the promoter; i.e. if input{i,6} == 1;
%             count = 1; %set count to 1; this resets EACH time the while loop is completed
%             index = 1 + length(satcat); %used for indexing rows made by the while loop
%             while length(bgcat{j}(1:end + TSS + bgcat{j,2} - bgcat{j,3} - 6 - round(length(input{i,5})/2) +1 -count))>0 %while the left part of bgcat - operator is > 0; inputs -10, spacer; assumes hexamer -10,-35;
%                 sat{index, 1} = bgcat{j}(1:end + TSS + bgcat{j,2} - bgcat{j,3} - 6 - round(length(input{i,5})/2) + 1 - count); %take left side of bgcat; removes seq up to -38; inputs -10, spacer; assumes hexamer -10,-35;
%                 sat{index, 2} = input{i,5}; %pull in operator seq into sat, column 2
%                 sat{index, 3} = bgcat{j}(end + TSS + bgcat{j,2} - bgcat{j,3} - 6 - round(length(input{i,5})/2) + 2 - count + length(input{i,5}):end); %pull out right side of bgcat
%                 satcat{index,1} = strcat(sat{index,1},sat{index,2},sat{index,3}); %concatenate the columns in sat and put them in their own row
%             count = count +1; %make new row of sat, satcat for each iteration of the while loop
%             index = index +1;
%             end
%         else
%             index = 1 + length(satcat); %used for indexing rows made by the while loop
%             satcat{index,1} = {}; %empty rows if 0's in column 6
%         end
%     end
% end
% % satcomp = satcat{:,index-1}(:)`; %put each row of satgather in a single column in satcomp, including the empty cells
% satcomp = cellfun(@(x) x(1:end),satcat(cellfun(@isempty,satcat) ~= 1),'un',0); %removes any empty cells from the satcat column and puts results in satcomp
% clear count; clear index;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
% Place an operator EVERY position  in a bg/promoter aggregate
% this assumes the operator is an activator; the sixth column in the excel file contains a logic input that designates which operators are used in this module
satcat = {};
loopcat = {};
for j = 1:length(bgcat) %for each concatenated background promoter sequence,...  assumes operator activates
    for i = 1:length(input(:,5)) %for each operator in 'input'...
        if input{i,6} == 1 %if that operator should be used to saturate the promoter; i.e. if input{i,6} == 1;
            count = 1; %set count to 1; this resets EACH time the while loop is completed
            index = 1 + length(satcat); %used for indexing rows made by the while loop
            while length(bgcat{j}(1:end -length(input{i,5}) +1 -count))>0 %while the left part of bgcat - operator is > 0; inputs -10, spacer; assumes hexamer -10,-35;
                sat{index, 1} = bgcat{j}(1:end -length(input{i,5}) +1 -count); %take left side of bgcat; removes seq up to -38; inputs -10, spacer; assumes hexamer -10,-35;
                sat{index, 2} = input{i,5}; %pull in operator seq into sat, column 2
                sat{index, 3} = bgcat{j}(end -length(input{i,5}) +2 -count+ length(input{i,5}):end); %pull out right side of bgcat
                satcat{index,1} = strcat(sat{index,1},sat{index,2},sat{index,3}); %concatenate the columns in sat and put them in their own row


                if input{i,14} == 1 %if want to incorporate "looping" into the saturation of this operator
                    n = 1; %n is the gap between the promoters
                    gap = 33; % gap is minimum size of gap; -1 has no gap
                    m = 1 + length(loopcat); %used for indexing rows made by the while loop
                    while length(satcat{index,1}(1:end -2*length(input{i,5}) +1 -count -n-gap))>0 %while the left part of bgcat - operator is > 0; inputs -10, spacer; assumes hexamer -10,-35;
                        loop{m, 1} = satcat{index,1}(1:end -2*length(input{i,5}) +1 -count -n-gap); %take left side of bgcat; removes seq up to -38; inputs -10, spacer; assumes hexamer -10,-35;
                        loop{m, 2} = input{i,5}; %pull in operator seq into sat, column 2
                        loop{m, 3} = satcat{index,1}(end -length(input{i,5}) +2 -count -n-gap:end); %pull out right side of bgcat
                        loopcat{m,1} = strcat(loop{m,1},loop{m,2},loop{m,3}); %concatenate the columns in sat and put them in their own row
                    n = n +1; % realize you can also change the spacing of the operators by changing +1 to +2, +3, etc
                    m = m +1;
                    end
                else
                    m = 1 + length(loopcat); %used for indexing rows made by the while loop
                    loopcat{m,1} = {}; %empty rows if 0's in column 6
                end

            count = count +1; %make new row of sat, satcat for each iteration of the while loop% realize you can also change the spacing of the operators by changing +1 to +2, +3, etc
            index = index +1;
            end
        else
            index = 1 + length(satcat); %used for indexing rows made by the while loop
            satcat{index,1} = {}; %empty rows if 0's in column 6
        end
    end
end
% satcomp = satcat{:,index-1}(:); %put each row of satgather in a single column in satcomp, including the empty cells
satcomp = cellfun(@(x) x(1:end),satcat(cellfun(@isempty,satcat) ~= 1),'un',0); %removes any empty cells from the satcat column and puts results in satcomp
loopcat = cellfun(@(x) x(1:end),loopcat(cellfun(@isempty,loopcat) ~= 1),'un',0);
satcomp = vertcat(satcomp,loopcat);
clear count; clear index;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Create operator location variable "loc"
%read the string in the "core location" column into new array; numbers must be text. use a comma after single numbers in excel to make them read as text
for i = 1:length(input(:,5))
    loc(i) = {strread(input{i,7},'%d','delimiter',',')};
end
clear i;

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Places operators at specified sites throughout the promoter background sequence.
% Operators listed in column 5 are placed in the locations listed in column 7 with "wobble" adjustment listed in column 9.
% e.g. input '3' in column 9 will place single operator in all sites +/- 3bp with respect to each location in column 7. +/-3 covers both sides of the DNA strand (1 helical turn ~10bp)

%First, check if the operator placement locations will go outside of
%the bounds of the promoter
for j = 1:length(bgcat)
        for i = 1:length(input(:,5));
            for k = 1:length(loc{i})
                if length(bgcat{j}) + TSS + loc{i}(k) - round(length(input{i,5})/2) <0 %if calculated location of operator position is negative, return error.
                    sprintf('ERROR: Operator %s cannot be placed in location %d because this location is outside (upstream) of the promoter.',char(input(i,5)),loc{i}(k))
                    return
                elseif length(bgcat{j}) + TSS + loc{i}(k) - round(length(input{i,5})/2) > length(bgcat{j}) %if calculated location of operator position is greater than length of the bgcat
                    sprintf('ERROR: Operator %s cannot be placed in location %d because this location is outside (downstream) of the promoter.',char(input(i,5)),loc{i}(k))
                    return
                else
                end
            end
        end
end
clear i; clear j; clear k;


%Insert a single operator at each designated location (loc(i)) from column 7, incorporates wobble position with variable m, pulled from column 9
%This mdoule is useful if the saturation term is unnecessary
count = 1;
% index = 1;
for j = 1:length(bgcat) %vary which background/promoter concatenation is used; assumes operator activates
    for i = 1:length(input(:,5)); %vary which operator is used
        if input{i,9}==0 %if statement to assess how many wobble bp to incorporate, if 0, place operator in only one position
            m = 0;
            for k = 1:length(loc{i})
                wob{count, 1} = bgcat{j}(1:end + TSS + loc{i}(k) - round(length(input{i,5})/2) + m); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                wob{count, 2} = input{i,5};
                wob{count, 3} = bgcat{j}(end + TSS + loc{i}(k) - round(length(input{i,5})/2) + length(input{i,5}) + m +1 :end);
                wobcat{count,1} = strcat(wob{count,1},wob{count,2},wob{count,3}); %concatenate the columns in sat and put them in their own row
                count = count +1;
            end
         elseif input{i,9}~=0  %if wobble is not 0, then incorporate the wobble positions for each designated operator location
            for k = 1:length(loc{i})
                for m = [-input{i,9}:1:input{i,9}]
                    wob{count, 1} = bgcat{j}(1:end + TSS + loc{i}(k) - round(length(input{i,5})/2) + m); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                    wob{count, 2} = input{i,5};
                    wob{count, 3} = bgcat{j}(end + TSS + loc{i}(k) - round(length(input{i,5})/2) + length(input{i,5}) + m +1 :end);
                    wobcat{count,1} = strcat(wob{count,1},wob{count,2},wob{count,3}); %concatenate the columns in wob and put them in their own row
                    count = count +1;
                end
            end
        else
        end
    end
end

wobcomp = vertcat(wobcat(:));
clear i; clear j; clear k; clear m; clear count;




%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Place individual operator in combinations of operator locations, both individually and in combination with other operators
%Inputs are the variables loc (operator location),  bgcat (background promoter sequence), 0 or 1 signals as to which operators to include, and other promoter sites lengths
%single operators are also included, so some redundancy with code above

%Essentially, these are multiple FOR loops that conditionally place operators in each designated site.
%WHILE loops are used to track of operator locations
%IF statements are embedded to determine which operators to combine

combcat2 = {};
count = 1; %to keep track of individual combinations
index = 1; %to keep track of multiple operator combinations
for j = 1:length(bgcat) %vary which background/promoter concatenation is used; assumes operator activates
    for i = 1:length(input(:,5)); %vary which operator is used
            for k = 1:length(loc{i})
                comb{count, 1} = bgcat{j}(1:end + TSS + loc{i}(k) - round(length(input{i,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                comb{count, 2} = input{i,5};
                comb{count, 3} = bgcat{j}(end + TSS + loc{i}(k) - round(length(input{i,5})/2) + length(input{i,5}) +1 :end);
                combcat{count,1} = strcat(comb{count,1},comb{count,2},comb{count,3}); %concatenate the columns in comb and put them in their own row

                if (i+1)<=length(loc) %condition to prevent i+1 from going out of bounds
                    if input{i,13}==1 && input{i+1,13}==1 %new if loop to integrate other operators;
                        for x = 1:length(loc{i+1})
                        comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                        comb2{index, 2} = input{i+1,5};
                        comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2) + length(input{i+1,5}) +1 :end);
                        combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3}); %concatenate the columns in comb2 and put them in their own row
                        index = index +1;
                        y = x+1;
                            while y <= length(loc{i+1})
                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                comb2{index, 2} = input{i+1,5};
                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                    end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                comb2{index, 4} = input{i+1,5};
                                comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1: end);
                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5}); %concatenate the columns in comb and put them in their own row
                                index = index +1;
                                z = y+1;
                                while z <= length(loc{i+1})
                                    comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                    comb2{index, 2} = input{i+1,5};
                                    comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                        end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2));
                                    comb2{index, 4} = input{i+1,5};
                                    comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                        end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                    comb2{index, 6} = input{i+1,5};
                                    comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1 :end);
                                    combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5},comb2{index,6},comb2{index,7}); %concatenate the columns in comb and put them in their own row
                                    index = index +1;
                                    xx = z+1;
                                    while xx <= length(loc{i+1})
                                        comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(xx) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                        comb2{index, 2} = input{i+1,5};
                                        comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(xx) - round(length(input{i+1,5})/2) + length(input{i+1,5}):...
                                            end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2));
                                        comb2{index, 4} = input{i+1,5};
                                        comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}):...
                                            end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2));
                                        comb2{index, 6} = input{i+1,5};
                                        comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}):...
                                            end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                        comb2{index, 8} = input{i+1,5};
                                        comb2{index, 9} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1 :end);
                                        combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},...
                                            comb2{count,5},comb2{index,6},comb2{index,7},comb2{index,8},comb2{index,9}); %concatenate the columns in comb and put them in their own row
                                        index = index +1;
                                        xx = xx+1;
                                    end
                                    z = z+1;
                                end
                                y = y+1;
                            end
                        end
                    end
                end

                if (i+2)<=length(loc) %condition to prevent i+1 from going out of bounds
                    if input{i,13}==1 && input{i+2,13}==1 %new if loop to integrate other operators;
                        for x = 1:length(loc{i+2})
                        comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                        comb2{index, 2} = input{i+2,5};
                        comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2) + length(input{i+2,5}) +1 :end);
                        combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3}); %concatenate the columns in comb2 and put them in their own row
                        index = index +1;
                        y = x+1;
                            while y <= length(loc{i+2})
                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                comb2{index, 2} = input{i+2,5};
                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+1,5})/2) + length(input{i+2,5})+1:...
                                    end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                comb2{index, 4} = input{i+2,5};
                                comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1: end);
                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5}); %concatenate the columns in comb and put them in their own row
                                index = index +1;
                                z = y+1;
                                while z <= length(loc{i+2})
                                    comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                    comb2{index, 2} = input{i,5};
                                    comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2) + length(input{i+2,5})+1:...
                                        end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2));
                                    comb2{index, 4} = input{i,5};
                                    comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                        end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                    comb2{index, 6} = input{i+2,5};
                                    comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1 :end);
                                    combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5},comb2{index,6},comb2{index,7}); %concatenate the columns in comb and put them in their own row
                                    index = index +1;
                                    xx = z+1;
                                    while xx <= length(loc{i+2})
                                        comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(xx) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                        comb2{index, 2} = input{i+2,5};
                                        comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(xx) - round(length(input{i+2,5})/2) + length(input{i+1,5})+1:...
                                            end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2));
                                        comb2{index, 4} = input{i+2,5};
                                        comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                            end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2));
                                        comb2{index, 6} = input{i+2,5};
                                        comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                            end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                        comb2{index, 8} = input{i+2,5};
                                        comb2{index, 9} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1 :end);
                                        combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},...
                                            comb2{count,5},comb2{index,6},comb2{index,7},comb2{index,8},comb2{index,9}); %concatenate the columns in comb and put them in their own row
                                        index = index +1;
                                        xx = xx+1;
                                    end
                                    z = z+1;
                                end
                                y = y+1;
                            end
                        end
                    end
                end

                count = count +1;
                l = k+1;
                while l <= length(loc{i})
                    comb{count, 1} = bgcat{j}(1:end + TSS + loc{i}(l) - round(length(input{i,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                    comb{count, 2} = input{i,5};
                    comb{count, 3} = bgcat{j}(end + TSS + loc{i}(l) - round(length(input{i,5})/2) + length(input{i,5})+1:...
                        end + TSS + loc{i}(k) - round(length(input{i,5})/2));
                    comb{count, 4} = input{i,5};
                    comb{count, 5} = bgcat{j}(end + TSS + loc{i}(k) - round(length(input{i,5})/2 )+ length(input{i,5}) +1: end);
                    combcat{count,1} = strcat(comb{count,1},comb{count,2},comb{count,3},comb{count,4},comb{count,5}); %concatenate the columns in comb and put them in their own row

                        if (i+1)<=length(loc) %condition to prevent i+1 from going out of bounds
                            if input{i,13}==1 && input{i+1,13}==1 %new if loop to integrate other operators;
                                for x = 1:length(loc{i+1})
                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                comb2{index, 2} = input{i+1,5};
                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2) + length(input{i+1,5}) +1 :end);
                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3}); %concatenate the columns in comb2 and put them in their own row
                                index = index +1;
                                y = x+1;
                                    while y <= length(loc{i+1})
                                        comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                        comb2{index, 2} = input{i+1,5};
                                        comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                            end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                        comb2{index, 4} = input{i+1,5};
                                        comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1: end);
                                        combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5}); %concatenate the columns in comb and put them in their own row
                                        index = index +1;
                                        z = y+1;
                                        while z <= length(loc{i+1})
                                            comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                            comb2{index, 2} = input{i+1,5};
                                            comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                                end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2));
                                            comb2{index, 4} = input{i+1,5};
                                            comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                                end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                            comb2{index, 6} = input{i+1,5};
                                            comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1 :end);
                                            combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5},comb2{index,6},comb2{index,7}); %concatenate the columns in comb and put them in their own row
                                            index = index +1;
                                            xx = z+1;
                                            while xx <= length(loc{i+1})
                                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(xx) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                                comb2{index, 2} = input{i+1,5};
                                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(xx) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2));
                                                comb2{index, 4} = input{i+1,5};
                                                comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2));
                                                comb2{index, 6} = input{i+1,5};
                                                comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                                comb2{index, 8} = input{i+1,5};
                                                comb2{index, 9} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1 :end);
                                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},...
                                                    comb2{count,5},comb2{index,6},comb2{index,7},comb2{index,8},comb2{index,9}); %concatenate the columns in comb and put them in their own row
                                                index = index +1;
                                                xx = xx+1;
                                            end
                                            z = z+1;
                                        end
                                        y = y+1;
                                    end
                                end
                            end
                        end

                        if (i+2)<=length(loc) %condition to prevent i+2 from going out of bounds
                            if input{i,13}==1 && input{i+2,13}==1 %new if loop to integrate other operators;
                                for x = 1:length(loc{i+2})
                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                comb2{index, 2} = input{i+2,5};
                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2) + length(input{i+2,5}) +1 :end);
                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3}); %concatenate the columns in comb2 and put them in their own row
                                index = index +1;
                                y = x+1;
                                    while y <= length(loc{i+2})
                                        comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                        comb2{index, 2} = input{i+2,5};
                                        comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+1,5})/2) + length(input{i+2,5})+1:...
                                            end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                        comb2{index, 4} = input{i+2,5};
                                        comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1: end);
                                        combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5}); %concatenate the columns in comb and put them in their own row
                                        index = index +1;
                                        z = y+1;
                                        while z <= length(loc{i+2})
                                            comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                            comb2{index, 2} = input{i,5};
                                            comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2) + length(input{i+2,5})+1:...
                                                end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2));
                                            comb2{index, 4} = input{i,5};
                                            comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                                end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                            comb2{index, 6} = input{i+2,5};
                                            comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1 :end);
                                            combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5},comb2{index,6},comb2{index,7}); %concatenate the columns in comb and put them in their own row
                                            index = index +1;
                                            xx = z+1;
                                            while xx <= length(loc{i+2})
                                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(xx) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                                comb2{index, 2} = input{i+2,5};
                                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(xx) - round(length(input{i+2,5})/2) + length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2));
                                                comb2{index, 4} = input{i+2,5};
                                                comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                                    end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2));
                                                comb2{index, 6} = input{i+2,5};
                                                comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                                    end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                                comb2{index, 8} = input{i+2,5};
                                                comb2{index, 9} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1 :end);
                                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},...
                                                    comb2{count,5},comb2{index,6},comb2{index,7},comb2{index,8},comb2{index,9}); %concatenate the columns in comb and put them in their own row
                                                index = index +1;
                                                xx = xx+1;
                                            end
                                            z = z+1;
                                        end
                                        y = y+1;
                                    end
                                end
                            end
                        end

                    count = count +1;
                    m = l+1;
                    while m <= length(loc{i})
                        comb{count, 1} = bgcat{j}(1:end + TSS + loc{i}(m) - round(length(input{i,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                        comb{count, 2} = input{i,5};
                        comb{count, 3} = bgcat{j}(end + TSS + loc{i}(m) - round(length(input{i,5})/2) + length(input{i,5})+1:...
                            end + TSS + loc{i}(l) - round(length(input{i,5})/2));
                        comb{count, 4} = input{i,5};
                        comb{count, 5} = bgcat{j}(end + TSS + loc{i}(l) - round(length(input{i,5})/2 )+ length(input{i,5})+1:...
                            end + TSS + loc{i}(k) - round(length(input{i,5})/2));
                        comb{count, 6} = input{i,5};
                        comb{count, 7} = bgcat{j}(end + TSS + loc{i}(k) - round(length(input{i,5})/2 )+ length(input{i,5}) +1 :end);
                        combcat{count,1} = strcat(comb{count,1},comb{count,2},comb{count,3},comb{count,4},comb{count,5},comb{count,6},comb{count,7}); %concatenate the columns in comb and put them in their own row

                        if (i+1)<=length(loc) %condition to prevent i+1 from going out of bounds
                            if input{i,13}==1 && input{i+1,13}==1 %new if loop to integrate other operators;
                                for x = 1:length(loc{i+1})
                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                comb2{index, 2} = input{i+1,5};
                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2) + length(input{i+1,5}) +1 :end);
                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3}); %concatenate the columns in comb2 and put them in their own row
                                index = index +1;
                                y = x+1;
                                    while y <= length(loc{i+1})
                                        comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                        comb2{index, 2} = input{i+1,5};
                                        comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                            end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                        comb2{index, 4} = input{i+1,5};
                                        comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1: end);
                                        combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5}); %concatenate the columns in comb and put them in their own row
                                        index = index +1;
                                        z = y+1;
                                        while z <= length(loc{i+1})
                                            comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                            comb2{index, 2} = input{i+1,5};
                                            comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                                end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2));
                                            comb2{index, 4} = input{i+1,5};
                                            comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                                end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                            comb2{index, 6} = input{i+1,5};
                                            comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1 :end);
                                            combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5},comb2{index,6},comb2{index,7}); %concatenate the columns in comb and put them in their own row
                                            index = index +1;
                                            xx = z+1;
                                            while xx <= length(loc{i+1})
                                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(xx) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                                comb2{index, 2} = input{i+1,5};
                                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(xx) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2));
                                                comb2{index, 4} = input{i+1,5};
                                                comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2));
                                                comb2{index, 6} = input{i+1,5};
                                                comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                                comb2{index, 8} = input{i+1,5};
                                                comb2{index, 9} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1 :end);
                                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},...
                                                    comb2{count,5},comb2{index,6},comb2{index,7},comb2{index,8},comb2{index,9}); %concatenate the columns in comb and put them in their own row
                                                index = index +1;
                                                xx = xx+1;
                                            end
                                            z = z+1;
                                        end
                                        y = y+1;
                                    end
                                end
                            end
                        end

                        if (i+2)<=length(loc) %condition to prevent i+2 from going out of bounds
                            if input{i,13}==1 && input{i+2,13}==1 %new if loop to integrate other operators;
                                for x = 1:length(loc{i+2})
                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                comb2{index, 2} = input{i+2,5};
                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2) + length(input{i+2,5}) +1 :end);
                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3}); %concatenate the columns in comb2 and put them in their own row
                                index = index +1;
                                y = x+1;
                                    while y <= length(loc{i+2})
                                        comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                        comb2{index, 2} = input{i+2,5};
                                        comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+1,5})/2) + length(input{i+2,5})+1:...
                                            end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                        comb2{index, 4} = input{i+2,5};
                                        comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1: end);
                                        combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5}); %concatenate the columns in comb and put them in their own row
                                        index = index +1;
                                        z = y+1;
                                        while z <= length(loc{i+2})
                                            comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                            comb2{index, 2} = input{i+2,5};
                                            comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2) + length(input{i+2,5})+1:...
                                                end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2));
                                            comb2{index, 4} = input{i+2,5};
                                            comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                                end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                            comb2{index, 6} = input{i+2,5};
                                            comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1 :end);
                                            combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5},comb2{index,6},comb2{index,7}); %concatenate the columns in comb and put them in their own row
                                            index = index +1;
                                            xx = z+1;
                                            while xx <= length(loc{i+2})
                                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(xx) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                                comb2{index, 2} = input{i+2,5};
                                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(xx) - round(length(input{i+2,5})/2) + length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2));
                                                comb2{index, 4} = input{i+2,5};
                                                comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                                    end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2));
                                                comb2{index, 6} = input{i+2,5};
                                                comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                                    end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                                comb2{index, 8} = input{i+2,5};
                                                comb2{index, 9} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1 :end);
                                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},...
                                                    comb2{count,5},comb2{index,6},comb2{index,7},comb2{index,8},comb2{index,9}); %concatenate the columns in comb and put them in their own row
                                                index = index +1;
                                                xx = xx+1;
                                            end
                                            z = z+1;
                                        end
                                        y = y+1;
                                    end
                                end
                            end
                        end

                        count = count +1;
                        n = m+1;
                        while n <= length(loc{i})
                            comb{count, 1} = bgcat{j}(1:end + TSS + loc{i}(n) - round(length(input{i,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                            comb{count, 2} = input{i,5};
                            comb{count, 3} = bgcat{j}(end + TSS + loc{i}(n) - round(length(input{i,5})/2) + length(input{i,5})+1:...
                                end + TSS + loc{i}(m) - round(length(input{i,5})/2));
                            comb{count, 4} = input{i,5};
                            comb{count, 5} = bgcat{j}(end + TSS + loc{i}(m) - round(length(input{i,5})/2 )+ length(input{i,5})+1:...
                                end + TSS + loc{i}(l) - round(length(input{i,5})/2));
                            comb{count, 6} = input{i,5};
                            comb{count, 7} = bgcat{j}(end + TSS + loc{i}(l) - round(length(input{i,5})/2 )+ length(input{i,5})+1:...
                                end + TSS + loc{i}(k) - round(length(input{i,5})/2));
                            comb{count, 8} = input{i,5};
                            comb{count, 9} = bgcat{j}(end + TSS + loc{i}(k) - round(length(input{i,5})/2 )+ length(input{i,5}) +1 :end);
                            combcat{count,1} = strcat(comb{count,1},comb{count,2},comb{count,3},comb{count,4},...
                                comb{count,5},comb{count,6},comb{count,7},comb{count,8},comb{count,9}); %concatenate the columns in comb and put them in their own row

                           if (i+1)<=length(loc) %condition to prevent i+1 from going out of bounds
                                if input{i,13}==1 && input{i+1,13}==1 %new if loop to integrate other operators;
                                    for x = 1:length(loc{i+1})
                                    comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                    comb2{index, 2} = input{i+1,5};
                                    comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2) + length(input{i+1,5}) +1 :end);
                                    combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3}); %concatenate the columns in comb2 and put them in their own row
                                    index = index +1;
                                    y = x+1;
                                        while y <= length(loc{i+1})
                                            comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                            comb2{index, 2} = input{i+1,5};
                                            comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                                end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                            comb2{index, 4} = input{i+1,5};
                                            comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1: end);
                                            combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5}); %concatenate the columns in comb and put them in their own row
                                            index = index +1;
                                            z = y+1;
                                            while z <= length(loc{i+1})
                                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                                comb2{index, 2} = input{i+1,5};
                                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2));
                                                comb2{index, 4} = input{i+1,5};
                                                comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                                    end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                                comb2{index, 6} = input{i+1,5};
                                                comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1 :end);
                                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5},comb2{index,6},comb2{index,7}); %concatenate the columns in comb and put them in their own row
                                                index = index +1;
                                                xx = z+1;
                                                while xx <= length(loc{i+1})
                                                    comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+1}(xx) - round(length(input{i+1,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                                    comb2{index, 2} = input{i+1,5};
                                                    comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+1}(xx) - round(length(input{i+1,5})/2) + length(input{i+1,5})+1:...
                                                        end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2));
                                                    comb2{index, 4} = input{i+1,5};
                                                    comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+1}(z) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                                        end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2));
                                                    comb2{index, 6} = input{i+1,5};
                                                    comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+1}(y) - round(length(input{i+1,5})/2 )+ length(input{i+1,5})+1:...
                                                        end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2));
                                                    comb2{index, 8} = input{i+1,5};
                                                    comb2{index, 9} = combcat{count,1}(end + TSS + loc{i+1}(x) - round(length(input{i+1,5})/2 )+ length(input{i+1,5}) +1 :end);
                                                    combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},...
                                                        comb2{count,5},comb2{index,6},comb2{index,7},comb2{index,8},comb2{index,9}); %concatenate the columns in comb and put them in their own row
                                                    index = index +1;
                                                    xx = xx+1;
                                                end
                                                z = z+1;
                                            end
                                            y = y+1;
                                        end
                                    end
                                end
                            end

                            if (i+2)<=length(loc) %condition to prevent i+2 from going out of bounds
                                if input{i,13}==1 && input{i+2,13}==1 %new if loop to integrate other operators;
                                    for x = 1:length(loc{i+2})
                                    comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                    comb2{index, 2} = input{i+2,5};
                                    comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2) + length(input{i+2,5}) +1 :end);
                                    combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3}); %concatenate the columns in comb2 and put them in their own row
                                    index = index +1;
                                    y = x+1;
                                        while y <= length(loc{i+2})
                                            comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                            comb2{index, 2} = input{i+2,5};
                                            comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+1,5})/2) + length(input{i+2,5})+1:...
                                                end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                            comb2{index, 4} = input{i+2,5};
                                            comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1: end);
                                            combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5}); %concatenate the columns in comb and put them in their own row
                                            index = index +1;
                                            z = y+1;
                                            while z <= length(loc{i+2})
                                                comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                                comb2{index, 2} = input{i,5};
                                                comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2) + length(input{i+2,5})+1:...
                                                    end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2));
                                                comb2{index, 4} = input{i,5};
                                                comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                                    end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                                comb2{index, 6} = input{i+2,5};
                                                comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1 :end);
                                                combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},comb2{index,5},comb2{index,6},comb2{index,7}); %concatenate the columns in comb and put them in their own row
                                                index = index +1;
                                                xx = z+1;
                                                while xx <= length(loc{i+2})
                                                    comb2{index, 1} = combcat{count,1}(1:end + TSS + loc{i+2}(xx) - round(length(input{i+2,5})/2)); %removes seq up to -38 (assumes hexamer for -10, -35), inputs are -10 loc and spacer size
                                                    comb2{index, 2} = input{i+2,5};
                                                    comb2{index, 3} = combcat{count,1}(end + TSS + loc{i+2}(xx) - round(length(input{i+2,5})/2) + length(input{i+1,5})+1:...
                                                        end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2));
                                                    comb2{index, 4} = input{i+2,5};
                                                    comb2{index, 5} = combcat{count,1}(end + TSS + loc{i+2}(z) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                                        end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2));
                                                    comb2{index, 6} = input{i+2,5};
                                                    comb2{index, 7} = combcat{count,1}(end + TSS + loc{i+2}(y) - round(length(input{i+2,5})/2 )+ length(input{i+2,5})+1:...
                                                        end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2));
                                                    comb2{index, 8} = input{i+2,5};
                                                    comb2{index, 9} = combcat{count,1}(end + TSS + loc{i+2}(x) - round(length(input{i+2,5})/2 )+ length(input{i+2,5}) +1 :end);
                                                    combcat2{index,1} = strcat(comb2{index,1},comb2{index,2},comb2{index,3},comb2{index,4},...
                                                        comb2{count,5},comb2{index,6},comb2{index,7},comb2{index,8},comb2{index,9}); %concatenate the columns in comb and put them in their own row
                                                    index = index +1;
                                                    xx = xx+1;
                                                end
                                                z = z+1;
                                            end
                                            y = y+1;
                                        end
                                    end
                                end
                            end

                            count = count +1;
                            n = n+1;
                        end
                        m = m+1;
                    end
                    l = l+1;
                end
            end
        end
    end


combcomp = vertcat(combcat(:),combcat2(:));
% combcomp = vertcat{:,index-1}(:); %put each row of satgather in a single column in satcomp, including the empty cells
combcomp = cellfun(@(x) x(1:end),combcomp(cellfun(@isempty,combcomp) ~= 1),'un',0); %removes any empty cells from the combcomp column

clear i; clear j; clear n; clear m; clear k; clear l; clear count;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Compile the controls (bgcat, bgw/O promoter, etc).



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
% Compile all the columns generated into a single, named set of columns;
% want first column to be names of oligos, second column to be oligo
% sequences; final results should be exported into excel

oligolist = vertcat(satcomp, wobcomp, combcomp, bg, bgcat(:,1)); %final list of oligos; incorporate the "negative controls", ie bg, bgcat

cuts = {};
for i = 1:length(oligolist); %create names of oligos and calculate lengths of each element of oligolist
    oligolist{i,1} = strcat(left{1}, oligolist{i,1},right{1}); %concatenate the left and right PCR primers onto the finished oligs
    f = int2str(i);
    f = cellstr(f);
    nameslist(i,:) = strcat(names,f); %make names list
%     [cuts{i,1},cuts{i,2}] = rebasecuts(oligolist{i,1}, {'BsaI'});
%     olen(i,:) = length(char(oligolist(i,:))); %get length of each oligo
%     prop(i,:) = oligoprop(char(oligolist(i,:))) %collect properties of each oligo
end

nameslist = cellstr(nameslist); %can only concatenate with oligolist as an array of strings
%     olen = int2str(olen);
%     olen = cellstr(olen);
%     cuts(:) = cellstr(cuts{:});

namesOligos = [nameslist, oligolist];
xlswrite('GluProm_OligoList.xls',namesOligos); %make excel file in SenseProm directory
toc
winopen('GluProm_OligoList.xls');
















%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%THE FOLLOWING IS LEFTOVER CODE FOR A GUI. MIGHT USE IT TO MAKE THIS MORE
%USER-FRIENDLY SOMEDAY.

% %create GUI that prompts input for: output names, which operators to
% %saturate, which operator cores to combine
% figure('units','pixels','position',[500,300,450,600],...
%              'toolbar','none','menu','none');
% neg10txt = uicontrol('style','text','position',[10,475,400,17],'string',...
%                 'Enter location of the E. coli -10 site (center of the -10 hexamer).');
% neg10lpos = uicontrol('style','edit','units','pixels',...
%                 'position',[20,450,400,17],'String','26',...
%                 'callback',@neg10_call);
% spacertxt = uicontrol('style','text','position',[10,425,400,17],'string',...
%                 'Enter length of spacer between -10 and -35 sites (center to center).');
% spacerlen = uicontrol('style','edit','units','pixels',...
%                 'position',[20,400,400,17],'String','17',...
%                 'callback',@spacer_call);
%
% % Create Repress/Activate checkboxes
% guitxt1 = uicontrol('style','text','position',[10,375,400,17],'string',...
%                 'Does the operator repress and/or activate?');
% guiin(1) = uicontrol('style','checkbox','units','pixels',...
%                 'position',[100,350,100,17],'string','Repress');
% guiin(2) = uicontrol('style','checkbox','units','pixels',...
%                 'position',[250,350,100,17],'string','Activate');
%
% % create edit text boxes for terminal ends
% cutsite5txt = uicontrol('style','text','position',[10,300,200,17],'string',...
%                 'Extra sequence @ 5-prime end?');
% cutsite5 = uicontrol('style','edit','units','pixels',...
%                 'position',[20,275,300,17],'string','',...
%                 'callback',@cutsite5_call);
% cutsite3txt = uicontrol('style','text','position',[10,250,200,17],'string',...
%                 'Extra sequence @ 3-prime end?');
% cutsite3 = uicontrol('style','edit','units','pixels',...
%                 'position',[20,225,300,17],'string','',...
%                 'callback',@cutsite3_call);
% oligonamestxt = uicontrol('style','text','position',[10,175,400,17],'string',...
%                 'Names of oligo outputs. Output will be ''string''1, ''string''2, etc.');
% oligonames = uicontrol('style','edit','units','pixels',...
%                 'position',[20,150,300,17],'string','name',...
%                 'callback',@oligonames_call);
%
% %create edit for max # constructs; can add later
% % constructstxt = uicontrol('style','text','position',[10,125,250,17],'string',...
% %                 'Max # of constructs to design?');
% % constructs = uicontrol('style','edit','units','pixels',...
% %                 'position',[20,100,300,17],'string','20',...
% %                 'callback',@constructs_call);
% % oligolentxt = uicontrol('style','text','position',[10,75,250,17],'string',...
% %                 'Max length of construction oligos?');
% % oligolen = uicontrol('style','edit','units','pixels',...
% %                 'position',[20,50,300,17],'string','60',...
% %                 'callback',@oligolen_call);
%
% % text inputs callbacks; collect global vars
% global prom
% global op
% global neg10
% global neg35
% global spacer
% global cutter5
% global cutter3
% global cons
% global oligolength
% global names
% global checkedsum
%
% %function callbacks for the GUI inputs
%     function promin_call(varargin)
%         prom = get(promin,'String');
%     end
%     function opin_call(varargin)
%         op = get(opin,'String');
%     end
%     function neg10_call(varargin)
%         neg10str = get(neg10lpos,'String');
%         neg10 = str2num(neg10str);
%     end
%     function spacer_call(varargin)
%         spacerstr = get(spacerlen,'String');
%         spacer = str2num(spacerstr);
%         neg35 = neg10 - spacer - 6;
%     end
%     function cutsite5_call(varargin)
%         cutter5 = get(cutsite5,'String');
%     end
%     function cutsite3_call(varargin)
%         cutter3= get(cutsite3,'String');
%     end
%     function constructs_call(varargin)
%         constr = get(constructs,'String');
%         cons = str2num(constr);
%     end
%     function oligolen_call(varargin)
%         oligoleng = get(oligolen,'String');
%         oligolength = str2num(oligoleng);
%     end
%     function oligonames_call(varargin)
%         names = get(oligonames,'string');
%         %oligos = strcat(name,int2str((1:cons*2).')); %chose to use
%         %variable names instead, but this line of code might be handy
%     end
%
%
%
%     % Create pushbutton; executes function design_call
%     guiok = uicontrol('style','pushbutton','units','pixels',...
%                 'position',[200,115,100,23],'string','Design!',...
%                 'callback',@design_call);
%
%
%     %pushbutton callback; first checks to see if necessary vars are filled,
%     %then calls @design for sequence design
%     function design_call(varargin)
%         clear nogo
%         vals = get(guiin,'Value');
%         checked = find([vals{:}]);
%         if isempty(checked)
%             checked = 'Please enter whether Repressor or Activator';
%             nogo = 1;
%         else nogo = 0;
%         end
%         if isempty(prom)
%             disp('Please enter a promoter sequence.')
%             nogo = 1;
%         else nogo = 0;
%         end
%         if isempty(op)
%             disp('Please enter an operator sequence.')
%             nogo = 1;
%         else nogo = 0;
%         end
%
%         if isempty(names)
%             disp('Please enter the names of the oligos.')
%             nogo = 1;
%         else nogo = 0;
%         end
% %         if isempty(cons)
% %             disp('Please enter the number of constructs you want to design')
% %             nogo = 1;
% %         else nogo = 0;
% %         end
% %         if isempty(oligolength)
% %             disp('Please enter the max length of the oligos.')
% %             nogo = 1;
% %         else nogo = 0;
% %         end
%         if nogo == 0;
%             checkedsum = sum(checked);
%             design (checkedsum)  %%call design function;
%             disp('Done.')
%         end
%     end
%
%
%     %function for sequence handeling; inputs are checkedsum and global vars
%     %are prom, neg10, neg35, neg;
%     %output is list of oligos
%     %General strategy is to collect "left" and "right" sides of target sequence and
%     %then concatenate them with the operator
%
%     function design(checkedsum)
%        if checkedsum == 1 %ie. if operator is repressor;
% %             disp('loop1'); %this line for checking
%
%             %loop to create the repressor designs; place full operator seq immediately downstream of -35 to
%             %end of promoter; index should be 1 -> #of times operator fits downstream
%             %of -35 +6 (length of operator at tail end of prom sequence)
%             for i = [1:(length(prom)-neg35 -3)];  % -3 here for half with of -35 half site
%                     repleft = cellstr(prom(1:(neg35+3+i-1))); %-1 corrects for inclusion of first bp in numbering (cause you can't use 0 to index)
%                     repright = cellstr(prom((neg35+3+i+length(op)):end));
%                     repleftop = strcat(repleft,op);
%                     repconstrarray = strcat(cutter5,repleftop,repright,cutter3);
%
%                     repconstr = char(repconstrarray); %must convert cell array constrarray to string again
%                     if (round(length(repconstr)/2)+20)<(length(repconstr));  %if loop to enable design when promoters are short (index becomes > than lenght of construct)
%                          reptop = repconstr(1:(round(length(repconstr)/2)+20));
%                          repbot = repconstr(round(length(repconstr)/2-20:end));
%                     else reptop = repconstr(1:end);
%                          repbot = repconstr(1:end);
%                     end
%
%                     reptopoligo(i,:) = cellstr(reptop); %design top oligo for annealing rx (20bp overlap)
%                     repbotoligo(i,:) = cellstr(seqrcomplement(repbot)); %design bottom oligo for annealing rx (20bp overlap); requires bioinfo toolbox
%
%                     acttopoligo = [];
%                     actbotoligo = [];
%              end
%
%         elseif checkedsum == 2 %ie. if operator is activator
% %              disp('loop2');
%              %loop to create activator designs; place full operator seq immediately
%              %upstream of -35, up to occupying half the -35 site w/ half the operator
%              for i = [1:(neg35 +1)];  % +1 here to account for construct where the last bp of the operator is -1 to left of input prom
%                  actright(i,:) = cellstr(prom(i:end)); %basically promoter sequence without the operator
%                  actopright(i,:) = strcat(op,actright(i,:));
%                  actoprightchar = char(actopright(i,:));
%                  if length(actoprightchar) < length(prom);  %if loop to add the 5' bp's back once operator is fully immersed in the promoter
%                     actleft(i,:) = cellstr(prom(1:i-1-length(op)));  %need to subtract length of operator - 1 (for numbering)
%                  else actleft(i,:) = cellstr('');
%                  end
%
%                 actconstrarray = strcat(cutter5,actleft(i,:),actopright(i,:),cutter3);
%                 actconstr = char(actconstrarray); %must convert cell array constrarray to string again
%
%                     if  length(actconstr) > (round(length(actconstr)/2)+20) %if loop in case actconstr is too short for long overlap
%                         acttopoligo(i,:) = cellstr(actconstr(1:(round(length(actconstr)/2)+20))); %design top oligo for annealing rx (20bp overlap)
%                     else
%                         acttopoligo(i,:) = cellstr(actconstr(1:end));
%                     end
%
%                     if  round(length(actconstr)/2)-20 > 0 %if loop in case actconstr is too short for long overlap
%                         actbotoligo(i,:) = cellstr(seqrcomplement(actconstr(round(length(actconstr)/2)-20:end))); %design bottom oligo for annealing rx (20bp overlap); requires bioinfo toolbox
%                     else
%                         actbotoligo(i,:) = cellstr(seqrcomplement(actconstr(1:end))); %design bottom oligo for annealing rx (20bp overlap); requires bioinfo toolbox
%                     end
%
%                 reptopoligo = [];
%                 repbotoligo = [];
%               end
%
%
%         elseif checkedsum == 3 %ie. if operator is repressor AND activator
% %             disp('loop1+2');
%             %loop to create the repressor designs; place full operator seq immediately downstream of -35 to
%             %end of promoter; index should be 1 -> #of times operator fits downstream
%             %of -35 +6 (length of operator at tail end of prom sequence)
%             for i = [1:(length(prom)-neg35 -3)];  % -3 here for half with of -35 half site
%                     repleft = cellstr(prom(1:(neg35+3+i-1))); %-1 corrects for inclusion of first bp in numbering (cause you can't use 0 to index)
%                     repright = cellstr(prom((neg35+3+i+length(op)):end));
%                     repleftop = strcat(repleft,op);
%                     repconstrarray = strcat(cutter5,repleftop,repright,cutter3);
%
%                     repconstr = char(repconstrarray); %must convert cell array constrarray to string again
%                     if (round(length(repconstr)/2)+20)<(length(repconstr));  %if loop to enable design when promoters are short (index becomes > than lenght of construct)
%                          reptop = repconstr(1:(round(length(repconstr)/2)+20));
%                          repbot = repconstr(round(length(repconstr)/2-20:end));
%                     else reptop = repconstr(1:end);
%                          repbot = repconstr(1:end);
%                     end
%
%                     reptopoligo(i,:) = cellstr(reptop); %design top oligo for annealing rx (20bp overlap)
%                     repbotoligo(i,:) = cellstr(seqrcomplement(repbot)); %design bottom oligo for annealing rx (20bp overlap); requires bioinfo toolbox
%
%             end
%
%             %loop to create activator designs; place full operator seq immediately
%             %upstream of -35, up to occupying half the -35 site w/ half the operator
%             for i = [1:(neg35 +1)];  % +1 here to account for construct where the last bp of the operator is -1 to left of input prom
%                  actright(i,:) = cellstr(prom(i:end)); %basically promoter sequence without the operator
%                  actopright(i,:) = strcat(op,actright(i,:));
%                  actoprightchar = char(actopright(i,:));
%                  if length(actoprightchar) < length(prom);  %if loop to add the 5' bp's back once operator is fully immersed in the promoter
%                     actleft(i,:) = cellstr(prom(1:i-1-length(op)));  %need to subtract length of operator - 1 (for numbering)
%                  else actleft(i,:) = cellstr('');
%                  end
%
%                 actconstrarray = strcat(cutter5,actleft(i,:),actopright(i,:),cutter3);
%                 actconstr = char(actconstrarray); %must convert cell array constrarray to string again
%
%                     if  length(actconstr) > (round(length(actconstr)/2)+20) %if loop in case actconstr is too short for long overlap
%                         acttopoligo(i,:) = cellstr(actconstr(1:(round(length(actconstr)/2)+20))); %design top oligo for annealing rx (20bp overlap)
%                     else
%                         acttopoligo(i,:) = cellstr(actconstr(1:end));
%                     end
%
%                     if  round(length(actconstr)/2)-20 > 0 %if loop in case actconstr is too short for long overlap
%                         actbotoligo(i,:) = cellstr(seqrcomplement(actconstr(round(length(actconstr)/2)-20:end))); %design bottom oligo for annealing rx (20bp overlap); requires bioinfo toolbox
%                     else
%                         actbotoligo(i,:) = cellstr(seqrcomplement(actconstr(1:end))); %design bottom oligo for annealing rx (20bp overlap); requires bioinfo toolbox
%                     end
%
%               end
%
%         end
%
%         oligolist = [reptopoligo; repbotoligo; acttopoligo; actbotoligo] %final list of oligos, in order: repTop,repBottom,actTop,actBottom
%
%         for i = 1:length(oligolist); %create names of oligos and calculate lengths of each element of oligolist
%             f = int2str(i);
%             f = cellstr(f);
%             nameslist(i,:) = strcat(names,f); %make names list
%             olen(i,:) = length(char(oligolist(i,:))); %get length of each oligo
% %             prop(i,:) = oligoprop(char(oligolist(i,:))) %collect properties of each oligo
%         end
%
%         nameslist = cellstr(nameslist); %can only concatenate with oligolist as an array of strings
%         olen = int2str(olen);
%         olen = cellstr(olen);
%
%         namesOligos = [nameslist, oligolist, olen];
%         xlswrite('SenseProm_OligoList.xls',namesOligos); %make excel file in SenseProm directory
%     end
%
% end
%
%
% % final = [oligos, oligolist]; %merge oligonames list and generated oligos
