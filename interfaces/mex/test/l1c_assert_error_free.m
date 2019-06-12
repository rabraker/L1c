function l1c_assert_error_free(fun)
  err_thrown = '';
   try 
       fun();
   catch ME
       err_thrown = ME.identifier;
       message = ME.message;
       error('l1c:AssertionFail',...
             'Assertion Failed: function %s raised:\n    %s', func2str(fun), message);
   end
   
end