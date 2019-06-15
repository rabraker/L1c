function l1c_assert_raises(fun, err_exp)
   err_thrown = '';
   try
       fun();
   catch ME
       err_thrown = ME.identifier;
   end
   
   status=true;
   if ~strcmp(err_thrown, err_exp);
     error('l1c:AssertionFail',...
           'Assertion Failed: expected: %s; Thrown: %s\n', err_exp, ...
           err_thrown)
   end
   
end
