if ~exist('data_dir', 'var')
  data_dir = get_data_dir();
end

MEX_RUN_CASE = getenv('MEX_RUN_SUITE');

suites = {@()test_mex_interface(data_dir),...
          @test_mex_utils
         };

pf = true;
total = length(suites);
n_skip = 0;
n_fail = 0;
n_pass = 0;
for suite = suites
  if ~isempty(MEX_RUN_CASE)
    name = func2str(suites{1});
    name = regexprep(name, '@|\([a-zA-Z0-9]*\)', '');
      
    if ~strcmp(name, MEX_RUN_CASE)
      fprintf('%s: ', name);
      status_str = sprintf('SKIP\n');
      status_str = clrs.skip_str(status_str);
      fprintf('%s\n', status_str);
      n_skip = n_skip + 1;
      continue;
    end
  end
    status = suite{1}();
    if status
      n_pass = n_pass + 1;
    else
      n_fail = n_fail + 1;
    end
    
    pf = pf & status;
end

st1 = '---------------------------------------------------------';
st2 = sprintf('Check in %s ', mfilename);
if pf
    status_str = sprintf('\n%s\n%s PASSED\n', st1, st2);
    status_str = clrs.pass_str(status_str);
    status = 0;
else
    status_str = sprintf('%s\n%s FAILED\n', st1, st2);
    status_str = clrs.fail_str(status_str);
    status = 1;
end
fprintf('%s', status_str);


exit(status);



function data_dir = get_data_dir()
  data_dir = getenv('TEST_DATA_DIR');
end
