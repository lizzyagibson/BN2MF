funcs1 = {@loop1, @loop2, @loop3, @loop4} ;   % let fun1, fun2 be two functions 

% use parfor for parallel scripts
parfor (ii = 1:4, 4)
      funcs1{ii}();
end
