% ---------------------------------------
% (Self-explanatory)
function Outcome = isodd(Value)
  if any(mod(Value, 1))
    error('Value must be an integer.')
  end
  Outcome = mod(Value, 2);
end