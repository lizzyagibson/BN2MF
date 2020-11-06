funcs1 = {@loop1, @loop2} ;   % let fun1, fun2 be two functions 

% use parfor for parallel scripts
parfor (ii = 1:2)
      funcs1{ii}();
end
