function p = MultiTestCorrection(p, nTests, Method, Inverse)
  % Correct p-values for multiple testing, by Bonferroni or Sidak method.
  % 
  % p = MultiTestCorrection(p, nTests, Method, Inverse)
  %
  % Given a certain number of tests or comparisons (nTests), correct
  % p-values (p) to obtain the family-wise error rate.  Uses the Bonferroni
  % correction by default (Method = 0), but can also use Sidak (Method =
  % 1).  Can go from FWER to single-test p-values by setting Inverse to
  % true (default false).
  %
  % 2014-05-12
  % Marc Lalancette, The Hospital for Sick Children, Toronto, Canada.

  if ~exist('Inverse', 'var') || isempty(Inverse)
    Inverse = false;
  end
  if ~exist('Method', 'var') || isempty(Method)
    Method = 0;
  end
  
  switch Method
    case 0 % Bonferroni correction.
      if ~Inverse
        p = p .* nTests;
      else
        p = p ./ nTests;
      end
      
    case 1 % Sidak correction.
      if ~Inverse
        p = 1 - (1 - p).^(nTests);
      else
        p = 1 - (1 - p).^(1/nTests);
      end
  
    otherwise
      error('Unknown multiple test correction method, should be 0 (Bonferroni) or 1 (Sidak).');
  end
  
end
