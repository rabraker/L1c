classdef L1cMexTesting
  
  
methods (Static)
  function [stat, msg] = assert_eq(x, y)
  % try 
    stat = x == y;
    msg = '';
    if ~stat 
      error('l1c:AssertionFail', 'Assertion failed: x !< y: x=%f, y=%f\n', x, y);
    end
  end


  function assert_lt(x, y)
  % try 
    stat = x < y;
    msg = '';
    if ~stat 
      error('l1c:AssertionFail', 'Assertion failed: x !< y: x=%f, y=%f\n', x, y);
    end
  end

  function assert_array_eq_tol(x, y, tol)
    x = x(:);
    y = y(:);
    
    [nx, mx] = size(x);
    [ny, my] = size(x);
    if (nx ~= ny)
      error('l1c:AssertionFail',...
            'Array size mismatch: nx=%d, ny=%d\n', nx, ny);
    end
    
    z = abs(x - y);
    idx_fail = find(z > tol);
    
    if ~isempty(idx_fail)
      idx_end = min(10, length(idx_fail));
      str = '';
      for k=1:idx_end
        str = sprintf('%s|x[%d] - y[%d]| !<TOL, (|%f - %f| !< %f)\n', str, k, ...
                      k, x(k), y(k), tol);
      end
      error('l1c:AssertionFail', 'Assertion failed: \n%s', str);
    end
    
  end  
  function assert_error_free(fun)
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
  function assert_raises(fun, err_exp)
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
  
  function pass_fail = run_suite(case_funcs)
    MEX_RUN_CASE = getenv('MEX_RUN_CASE');
    pass_fail = true;
    for case_k = case_funcs
      if ~isempty(MEX_RUN_CASE)
        name = func2str(case_k{1});
        name = regexprep(name, '@|\([a-zA-Z0-9]*\)', '');
        
        if ~strcmp(name, MEX_RUN_CASE)
          fprintf('%s: ', name);
          status_str = sprintf('SKIP\n');
          status_str = clrs.skip_str(status_str);
          fprintf('%s\n', status_str);
          continue;
        end
      end
      
      pass_fail = pass_fail & L1cMexTesting.run_case(case_k{1});
    end

    st1 = '---------------------------------------------------------';
    st2 = sprintf('Check in %s ', mfilename);
    if pass_fail
      status_str = sprintf('\n%s\n%s PASSED\n', st1, st2);
      status_str = clrs.pass_str(status_str);
    else
      status_str = sprintf('%s\n%s FAILED\n', st1, st2);
      status_str = clrs.fail_str(status_str);
    end
    fprintf('%s', status_str);    
    
    
  end
  
  function status = run_case(f)
  % status = TMI_runner(f)
  % (TMI=Test Mex InterFace). Runs the function f and reports 
  % pass or fail status. WHen f fails, it should throw an error.
  % 
    name = func2str(f);
    fprintf('%s: ', name);
    status = true;

    try
      f();
      fprintf(clrs.pass_str('Passed\n'));
    catch ME
      status = false;
      fprintf([clrs.fail_str('FAILED: '), '%s\n'], ME.message);
      if ~strcmp(ME.identifier, 'l1c:AssertionFail')
        for k=1:length(ME.stack)
          fprintf('    In function %s, line %d of file %s\n', ME.stack(k).name, ...
                  ME.stack(k).line, ME.stack(k).file);
        end
      end
      status = false;
    end
  end

  
  
end

  
end
