function status = TMI_runner(f)
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
