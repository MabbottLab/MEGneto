function [p, observeddifference] = permutationTest(sample1, sample2, permutations, varargin)

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
[~, ~, ~, stats_observed] = ttest2(sample1', sample2');
observeddifference = abs(stats_observed.tstat);

% running test
randomdifferences = zeros(1, permutations);
for n = 1:permutations
    permutation = cell2mat(arrayfun(@randperm, ones(size(allobservations,1),1)*(size(allobservations,2)), 'UniformOutput', false));
    linear_indices = sub2ind(size(permutation), repmat([1:num_participants]', num_variables, 1), permutation(:));
    randSample = reshape(allobservations(linear_indices), size(allobservations));
    
    % saving differences between the two samples
    [~, ~, ~, stats] = ttest2(randSample(:,1:size(sample1,2))', randSample(:,1:size(sample2,2))');
    randomdifferences(n) = max(abs(stats.tstat));
end

p = arrayfun(@(x) sum(randomdifferences > x)/5000, observeddifference);

end