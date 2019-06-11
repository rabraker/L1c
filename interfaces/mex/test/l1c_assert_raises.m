function [stat, err_thrown] = l1c_assert_raises(fun, err_exp)
   err_thrown = '';
   try
       fun();
   catch ME
       err_thrown = ME.identifier;
   end
   
   stat = strcmp(err_thrown, err_exp);

   st = dbstack();
   fprintf("%s: \n", st(2).name);
   if ~stat
    fprintf("\bFailed: expected: %s, Thrown: %s\n", err_exp, err_thrown)
  else
    fprintf("\bPassed\n");
  end
end