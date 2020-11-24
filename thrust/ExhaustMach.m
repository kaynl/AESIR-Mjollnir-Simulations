function Me = ExhaustMach(opts,At)
%EXHAUSTMACH finds the 0 of the following function to find the exhaust mach
%number
Me_zero = fzero(@(Mach) ExhaustMach_fct(Mach,At,opts),2);
Me = max(0,Me_zero);


end

