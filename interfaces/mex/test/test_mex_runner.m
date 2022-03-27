if ~exist('data_dir', 'var')
  data_dir = get_data_dir();
end

MEX_RUN_SUITE = getenv('MEX_RUN_SUITE');

suites = {@(skip_all)test_mex_interface(skip_all, data_dir),...
          @(skip_all)test_mex_utils(skip_all)
         };

pf = true;
n_skip = 0;
n_fail = 0;
n_pass = 0;
for suite = suites
  if ~isempty(MEX_RUN_SUITE)
    name = func2str(suite{1});
    name = regexprep(name, '@\([a-zA-Z0-9_,]*\)|\([a-zA-Z0-9_,]*\)', '');

    if ~strcmp(name, MEX_RUN_SUITE)
      fprintf('%s: ', name);
      status_str = sprintf('SKIP\n');
      status_str = clrs.skip_str(status_str);
      fprintf('%s\n', status_str);
      [~,~,nskip_k] = suite{1}(true);
      n_skip = n_skip + nskip_k;
      continue;
    end
  end

  [npass_k, nfail_k, nskip_k] = suite{1}(false);
  n_pass = n_pass + npass_k;
  n_fail = n_fail + nfail_k;
  n_skip = n_skip + nskip_k;  
    
    pf = pf & (nfail_k ==0);
end

total = n_pass + n_fail + n_skip;
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
