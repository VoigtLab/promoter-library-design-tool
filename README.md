# promoter-library-design-tool

Voigt lab Promoter Library Design Tool

Used in Moser, et al. (201x) (unpublished as of 8/18)

Code is Matlab "SenseProm_V6".

Input is an excel or CSV file that detail different parameters of the promoter design task. 
Output is an excel file with a solid column of promoter sequences (150 bp) ready for chip synthesis. 

An example excel sheet for promoter synthesis is provided in the github folder.

Please note that large sections are commented (%) and therefore not functional until uncommented. 

Input parameters are as follows: 

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

