if ~exist('data_dir', 'var')
  data_dir = get_data_dir();
end

MEX_RUN_SUITE = getenv('MEX_RUN_SUITE');

suites = {@()test_mex_interface(data_dir),...
          @test_mex_utils
         };

pf = true;
total = length(suites);
n_skip = 0;
n_fail = 0;
n_pass = 0;
for suite = suites
  if ~isempty(MEX_RUN_SUITE)
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

st1 = '----------------------------------------';
st2 = sprintf('Summary of %s :', mfilename);

sn_pass = sprintf('# Passed:  %d / %d', n_pass, total);
sn_fail = sprintf('# Failed:  %d / %d', n_fail, total);
sn_skip = sprintf('# Skipped: %d / %d', n_skip, total);


if n_pass > 0
  sn_pass = clrs.pass_str(sn_pass);
end
if n_skip > 0
  sn_skip = clrs.skip_str(sn_skip);
end
if n_fail > 0
  sn_fail = clrs.fail_str(sn_fail);
end

if n_fail > 0
    st1 = clrs.fail_str(st1);
    st2 = sprintf('%s FAILED', st2);
    st2 = clrs.fail_str(st2);
    status = 1;
else
    st1 = clrs.pass_str(st1);
    st2 = sprintf('%s PASSED', st2);
    st2 = clrs.pass_str(st2);
    status = 0;
end

fprintf('%s\n%s\n%s\n', st1, st2, st1);
fprintf('%s\n', sn_pass);
fprintf('%s\n', sn_skip);
fprintf('%s\n', sn_fail);
fprintf('%s\n\n\n', st1);


exit(status);



function data_dir = get_data_dir()
  data_dir = getenv('TEST_DATA_DIR');
end
