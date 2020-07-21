% [p, observeddifference, effectsize] = permutationTest(sample1, sample2, permutations [, varargin])
%
%       Permutation test (aka randomisation test), testing for a difference
%       in means between two samples. 
%
% In:
%       sample1 - vector of measurements representing one sample
%       sample2 - vector of measurements representing a second sample
%       permutations - the number of permutations
%
% Optional (name-value pairs):
%       sidedness - whether to test one- or two-sided:
%           'both' - test two-sided (default)
%           'smaller' - test one-sided, alternative hypothesis is that
%                       the mean of sample1 is smaller than the mean of
%                       sample2
%           'larger' - test one-sided, alternative hypothesis is that
%                      the mean of sample1 is larger than the mean of
%                      sample2
%       exact - whether or not to run an exact test, in which all possible
%               combinations are considered. this is only feasible for
%               relatively small sample sizes. the 'permutations' argument
%               will be ignored for an exact test. (1|0, default 0)
%       plotresult - whether or not to plot the distribution of randomised
%                    differences, along with the observed difference (1|0,
%                    default: 0)
%       showprogress - whether or not to show a progress bar. if 0, no bar
%                      is displayed; if showprogress > 0, the bar updates 
%                      every showprogress-th iteration.
%
% Out:  
%       p - the resulting p-value
%       observeddifference - the observed difference between the two
%                            samples, i.e. mean(sample1) - mean(sample2)
%       effectsize - the effect size
%
% Usage example:
%       >> permutationTest(rand(1,100), rand(1,100)-.25, 10000, ...
%          'plotresult', 1, 'showprogress', 250)
% 
%                    Copyright 2015-2018 Laurens R Krol
%                    Team PhyPA, Biological Psychology and Neuroergonomics,
%                    Berlin Institute of Technology

% 2020-07-14 lrk
%   - Added version-dependent call to hist/histogram
% 2019-02-01 lrk
%   - Added short description
%   - Increased the number of bins in the plot
% 2018-03-15 lrk
%   - Suppressed initial MATLAB:nchoosek:LargeCoefficient warning
% 2018-03-14 lrk
%   - Added exact test
% 2018-01-31 lrk
%   - Replaced calls to mean() with nanmean()
% 2017-06-15 lrk
%   - Updated waitbar message in first iteration
% 2017-04-04 lrk
%   - Added progress bar
% 2017-01-13 lrk
%   - Switched to inputParser to parse arguments
% 2016-09-13 lrk
%   - Caught potential issue when column vectors were used
%   - Improved plot
% 2016-02-17 toz
%   - Added plot functionality
% 2015-11-26 First version

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [p, observeddifference, effectsize] = permutationTest(sample1, sample2, permutations, varargin)

% parsing input
p = inputParser;

addRequired(p, 'sample1', @isnumeric);
addRequired(p, 'sample2', @isnumeric);
addRequired(p, 'permutations', @isnumeric);

addParamValue(p, 'sidedness', 'both', @(x) any(validatestring(x,{'both', 'smaller', 'larger'})));
addParamValue(p, 'exact' , 0, @isnumeric);
addParamValue(p, 'plotresult', 0, @isnumeric);
addParamValue(p, 'showprogress', 0, @isnumeric);

parse(p, sample1, sample2, permutations, varargin{:})

sample1 = rmmissing(p.Results.sample1);
sample2 = rmmissing(p.Results.sample2);
permutations = p.Results.permutations;
sidedness = p.Results.sidedness;
exact = p.Results.exact;
plotresult = p.Results.plotresult;
showprogress = p.Results.showprogress;

allobservations = [sample1, sample2];
num_variables = size(allobservations,1);
num_participants = size(allobservations,2);
observeddifference = nanmean(sample1,2) - nanmean(sample2,2); % avg diff b/w groups for all 90 regions
effectsize = observeddifference ./ nanmean([std(sample1,0,2), std(sample2,0,2)],2);

% running test
randomdifferences = zeros(1, permutations);
for n = 1:permutations
    permutation = cell2mat(arrayfun(@randperm, ones(size(allobservations,1),1)*(size(allobservations,2)), 'UniformOutput', false));
    linear_indices = sub2ind(size(permutation), repmat([1:num_participants]', num_variables, 1), permutation(:));
    randSample = reshape(allobservations(linear_indices), size(allobservations));
    
    % saving differences between the two samples
    randomdifferences(n) = nanmax(abs(nanmean(randSample(:,1:size(sample1,2)),2) - nanmean(randSample(:,1:size(sample2,2)),2)));
end

% getting probability of finding observed difference from random permutations
if strcmp(sidedness, 'both')
    p = (sum(abs(randomdifferences) > abs(observeddifference),2)+1) ./ (permutations+1);
elseif strcmp(sidedness, 'smaller')
    p = (sum(randomdifferences < observeddifference,2)+1) ./ (permutations+1);
elseif strcmp(sidedness, 'larger')
    p = (sum(randomdifferences > observeddifference,2)+1) / (permutations+1);
end

% plotting result
if plotresult
    figure;
    if verLessThan('matlab', '8.4')
        % MATLAB R2014a and earlier
        hist(randomdifferences, 20);
    else
        % MATLAB R2014b and later
        histogram(randomdifferences, 20);
    end
    hold on;
    xlabel('Random differences');
    ylabel('Count')
    od = plot(observeddifference, 0, '*r', 'DisplayName', sprintf('Observed difference.\nEffect size: %.2f,\np = %f', effectsize, p));
    legend(od);
end

end