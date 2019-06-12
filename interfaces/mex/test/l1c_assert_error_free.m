function stat = l1c_assert_error_free(fun)
  err_thrown = '';
   try 
       fun();
   catch ME
       err_thrown = ME.identifier;
       message = ME.message;
   end
   stat = strcmp(err_thrown, '');
   st = dbstack();
   fprintf("%s: ", st(2).name);
   if ~stat
    fprintf([clrs.fail_str('Failed: '), 'expected: n/a\n     Thrown: %s\n'], err_thrown)
    % keyboard;
    fprintf('         Error Message: %s\n', message);
  else
    fprintf(clrs.pass_str('Passed\n')); 
  end
end