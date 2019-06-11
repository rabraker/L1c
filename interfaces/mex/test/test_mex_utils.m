
checks = {@check_mex_assert_num_inputs,...
    @check_mex_get_double_scalar_or_fail,...
    @check_get_double_array_or_fail,...
    @check_get_int_array_or_fail,...
    @check_assert_2Darray_with_size};

pf = true;
for check = checks
    pf = pf & check{1}();
end

fprintf("---------------------------------------------------------\n");
if pf
    status_str = 'PASSED';
    status = 0;
else
    status = 1;
    status_str = 'FAILED';
end
fprintf("Checks in file %s %s\n", mfilename, status_str);


function pass_fail = check_assert_2Darray_with_size()
    err_exp = 'l1c:not2DArray';
    f = @()TMU_mex_assert_2Darray_with_size(5);
    pass_fail = l1c_assert_raises(f, err_exp);

    err_exp = 'l1c:ArraySize';
    f = @()TMU_mex_assert_2Darray_with_size([5, 6], 1);
    pass_fail = pass_fail & l1c_assert_raises(f, err_exp);

    f = @()TMU_mex_assert_2Darray_with_size([5, 6; 4, 5], 1);
    pass_fail = pass_fail & l1c_assert_raises(f, err_exp);

    A = [1,2,3,4; 1,2,3,4; 5,6,7,8];
    f = @()TMU_mex_assert_2Darray_with_size(A, 1);
    pass_fail = pass_fail & l1c_assert_error_free(f);
    
end

function pass_fail = check_get_int_array_or_fail()
    err_exp = 'l1c:not1DArray';
    f = @()TMU_mex_get_int_array_or_fail(int32(5));
    pass_fail = l1c_assert_raises(f, err_exp);

    f = @()TMU_mex_get_int_array_or_fail(int32([5, 6]), 1);
    pass_fail = pass_fail & l1c_assert_error_free(f);
    
    err_exp = 'l1c:notInt32';
    f = @()TMU_mex_get_int_array_or_fail([5, 6], 1);
    pass_fail = pass_fail & l1c_assert_raises(f, err_exp);
    
end

function pass_fail = check_get_double_array_or_fail()

    err_exp = 'l1c:not1DArray';
    f = @()TMU_mex_get_double_array_or_fail(5);
    pass_fail = l1c_assert_raises(f, err_exp);

    f = @()TMU_mex_get_double_array_or_fail([5, 6], 1);
    pass_fail = pass_fail & l1c_assert_error_free(f);
    
    err_exp = 'l1c:notDouble';
    f = @()TMU_mex_get_double_array_or_fail(int32([5, 6]), 1);
    pass_fail = pass_fail & l1c_assert_raises(f, err_exp);    
end

function pass_fail = check_mex_get_double_scalar_or_fail()
    
    f = @()TMU_mex_get_double_scalar_or_fail(5, [1, 2]);
    
    err_exp = 'l1c:notDoubleScaler';
    pass_fail = l1c_assert_raises(f, err_exp);
    
    f = @()TMU_mex_get_double_scalar_or_fail(5, 1);
    
    pass_fail = pass_fail & l1c_assert_error_free(f);
end


function pass_fail = check_mex_assert_num_inputs()
  
  f = @()TMU_mex_assert_num_inputs(5);
  
  err_exp = 'l1c:nlhs';
  pass_fail = l1c_assert_raises(f, err_exp);
  
  
  f = @()TMU_mex_assert_num_inputs(5, 2);
  pass_fail = pass_fail & l1c_assert_error_free(f);
end









