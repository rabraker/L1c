function stat = l1c_assert_error_free(fun)
    err_thrown = '';
   try 
       fun();
   catch ME
       err_thrown = ME.identifier;
   end
   stat = strcmp(err_thrown, '');
   
   st = dbstack();
   fprintf("%s: \n", st(2).name);
   if ~stat
    fprintf('\bFailed: expedected: n/a, Thrown: %s\n', err_thrown)
    % keyboard;
    fprintf('         Error Message: %s\n', ME.message);
  else
    fprintf('\bPassed\n');
  end
end