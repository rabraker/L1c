
if ~exist('fpath', 'var')
  fpath = fullfile(get_data_dir(), 'example_img_data.json');
end

MEX_RUN_CASE = getenv('MEX_RUN_CASE');

checks = {@check_mex_assert_num_inputs,...
    @check_mex_get_double_scalar_or_fail,...
    @check_get_double_array_or_fail,...
    @check_get_int_array_or_fail,...
    @check_assert_2Darray_with_size,...
    @()check_breg_TV_mex(fpath),...
    @()test_l1qc_dct_mex(fpath),...
    @()test_nesta_dctTV(fpath),...
         };

pf = true;
for check = checks
  if ~isempty(MEX_RUN_CASE)
    name = func2str(check{1});
    name = regexprep(name, '@|\([a-zA-Z0-9]*\)', '');
      
    if ~strcmp(name, MEX_RUN_CASE)
      fprintf('%s: ', name);
      status_str = sprintf('SKIP\n');
      status_str = clrs.skip_str(status_str);
      fprintf('%s\n', status_str);
      continue;
    end
  end
  
    pf = pf & runner(check{1});
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




function status = runner(f)
    name = func2str(f);
    fprintf('%s: ', name);
    status = true;

    try
        f();
        fprintf(clrs.pass_str('Passed\n'));
    catch ME
        status = false;
        fprintf([clrs.pass_str('FAILED: '), '%s\n'], ME.message);
        if ~strcmp(ME.identifier, 'l1c:AssertionFail')
          for k=1:length(ME.stack)
            fprintf('    In function %s, line %d of file %s\n', ME.stack(k).name, ...
                    ME.stack(k).line, ME.stack(k).file);
          end
        end
        status = false;
    end
end


function data_dir = get_data_dir()
  data_dir = getenv('TEST_DATA_DIR');
end


function test_nesta_dctTV(fpath)
  verbose = 0;
  fid = fopen(fpath, 'r');
  dat_json = fscanf(fid, '%s');
  fclose(fid);
  
  dat = jsondecode(dat_json);
  
  opts = nesta_opts('alpha_v', 1, 'alpha_h', 1);

  if dat.one_based_index == 1
    pix_idx = dat.pix_idx;
  else
    pix_idx = dat.pix_idx+1;
  end
  fprintf('mrow:%f, mcol:%f\n', dat.mrow, dat.mcol);
  [x_est, status] = nesta_dctTV(dat.mrow, dat.mcol, dat.b, pix_idx, opts);

  l1c_assert_eq(status, 0);
  
  % TODO: implement the property check we use in the c-testsuite.
end


function test_l1qc_dct_mex(fpath)
  verbose = 0;
  fid = fopen(fpath, 'r');
  dat_json = fscanf(fid, '%s');
  fclose(fid);
  
  dat = jsondecode(dat_json);
  
  opts = l1qc_dct_opts('verbose', verbose, 'l1_tol', 1e-5,...
                       'epsilon', 0.1, 'mu', 10);

  if dat.one_based_index == 1
    pix_idx = dat.pix_idx;
  else
    pix_idx = dat.pix_idx+1;
  end
  [x_est, LBRes]= l1qc_dct(dat.mtot, 1, dat.b, pix_idx, opts);

  l1c_assert_eq(LBRes.status, 0);
  
  % TODO: implement the property check we use in the c-testsuite.
end

function check_breg_TV_mex(fpath)
  fid = fopen(fpath, 'r');
  dat_json = fscanf(fid, '%s');
  fclose(fid);
  
  dat = jsondecode(dat_json);
  
  n = sqrt(dat.mtot);
  m = n;
  x = dat.x_orig(:) + rand(dat.mtot, 1);
  img = reshape(x, n, n);
  mu = 5;
  tol = 0.001;
  
  Xopt = breg_anistropic_TV(img, mu, tol, 1000);
  
  l1c_assert_eq(size(Xopt, 1), n);
  l1c_assert_eq(size(Xopt, 2), m);
  
  err = norm(Xopt(:) - x(:))/norm(x(:));
  l1c_assert_lt(err, 1.0);
end

function [stat, msg] = l1c_assert_eq(x, y)
  % try 
    stat = x == y;
    msg = '';
    if ~stat 
      error('l1c:AssertionFail', 'Assertion failed: x !< y: x=%f, y=%f\n', x, y);
    end
end


function l1c_assert_lt(x, y)
  % try 
    stat = x < y;
    msg = '';
    if ~stat 
      error('l1c:AssertionFail', 'Assertion failed: x !< y: x=%f, y=%f\n', x, y);
    end
end

function l1c_assert_array_eq_tol(x, y, tol)
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

function check_assert_2Darray_with_size()
    err_exp = 'l1c:not2DArray';
    f = @()TMU_mex_assert_2Darray_with_size(5);
    l1c_assert_raises(f, err_exp);

    err_exp = 'l1c:ArraySize';
    f = @()TMU_mex_assert_2Darray_with_size([5, 6], 1);
    l1c_assert_raises(f, err_exp);

    f = @()TMU_mex_assert_2Darray_with_size([5, 6; 4, 5], 1);
    l1c_assert_raises(f, err_exp);

    A = [1,2,3,4; 1,2,3,4; 5,6,7,8];
    f = @()TMU_mex_assert_2Darray_with_size(A, 1);
    l1c_assert_error_free(f);
    
end

function  check_get_int_array_or_fail()
    err_exp = 'l1c:not1DArray';
    f = @()TMU_mex_get_int_array_or_fail(int32(5));
    l1c_assert_raises(f, err_exp);

    f = @()TMU_mex_get_int_array_or_fail(int32([5, 6]), 1);
     l1c_assert_error_free(f);
    
    err_exp = 'l1c:notInt32';
    f = @()TMU_mex_get_int_array_or_fail([5, 6], 1);
    l1c_assert_raises(f, err_exp);
    
end

function check_get_double_array_or_fail()

    err_exp = 'l1c:not1DArray';
    f = @()TMU_mex_get_double_array_or_fail(5);
    l1c_assert_raises(f, err_exp);

    f = @()TMU_mex_get_double_array_or_fail([5, 6], 1);
    l1c_assert_error_free(f);
    
    err_exp = 'l1c:notDouble';
    f = @()TMU_mex_get_double_array_or_fail(int32([5, 6]), 1);
    l1c_assert_raises(f, err_exp);    
end

function check_mex_get_double_scalar_or_fail()
    
    f = @()TMU_mex_get_double_scalar_or_fail(5, [1, 2]);
    
    err_exp = 'l1c:notDoubleScaler';
    l1c_assert_raises(f, err_exp);
    
    f = @()TMU_mex_get_double_scalar_or_fail(5, 1);
    
    l1c_assert_error_free(f);
end


function check_mex_assert_num_inputs()
  
  f = @()TMU_mex_assert_num_inputs(5);
  
  err_exp = 'l1c:nlhs';
  l1c_assert_raises(f, err_exp);
  
  
  f = @()TMU_mex_assert_num_inputs(5, 2);
  l1c_assert_error_free(f);
end









