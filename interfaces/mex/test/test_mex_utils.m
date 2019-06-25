function [npass, nfail, nskip] = test_mex_utils(skip_all)
  if nargin < 1
    skip_all = false;
  end
  
  cases = {@check_mex_assert_num_inputs,...
            @check_mex_get_double_scalar_or_fail,...
            @check_get_double_array_or_fail,...
            @check_get_int_array_or_fail,...
            @check_assert_2Darray_with_size,...
            @check_assert_scalar_struct,...
            @check_get_double_from_struct_or_fail,...
           };

  if ~skip_all
    [npass, nfail, nskip] = L1cMexTesting.run_suite(cases);
  else
    npass = 0;
    nfail = 0;
    nskip = length(cases);
  end

end


function check_assert_2Darray_with_size()
    err_exp = 'l1c:not2DArray';
    f = @()TMU_mex_assert_2Darray_with_size(5);
    L1cMexTesting.assert_raises(f, err_exp);

    err_exp = 'l1c:ArraySize';
    f = @()TMU_mex_assert_2Darray_with_size([5, 6], 1);
    L1cMexTesting.assert_raises(f, err_exp);

    f = @()TMU_mex_assert_2Darray_with_size([5, 6; 4, 5], 1);
    L1cMexTesting.assert_raises(f, err_exp);

    A = [1,2,3,4; 1,2,3,4; 5,6,7,8];
    f = @()TMU_mex_assert_2Darray_with_size(A, 1);
    L1cMexTesting.assert_error_free(f);
    
end

function  check_get_int_array_or_fail()
    err_exp = 'l1c:not1DArray';
    f = @()TMU_mex_get_int_array_or_fail(int32(5));
    L1cMexTesting.assert_raises(f, err_exp);

    f = @()TMU_mex_get_int_array_or_fail(int32([5, 6]), 1);
     L1cMexTesting.assert_error_free(f);
    
    err_exp = 'l1c:notInt32';
    f = @()TMU_mex_get_int_array_or_fail([5, 6], 1);
    L1cMexTesting.assert_raises(f, err_exp);
    
end

function  check_assert_scalar_struct()
    err_exp = 'l1c:notScalar';
    ST = struct('field1', 55, 'field2', pi);
    ST_array = [ST; ST];
    f = @()TMU_mex_assert_scalar_struct(ST_array);
    L1cMexTesting.assert_raises(f, err_exp);

    f = @()TMU_mex_assert_scalar_struct(ST);
     L1cMexTesting.assert_error_free(f);
    
    err_exp = 'l1c:notStruct';
    f = @()TMU_mex_assert_scalar_struct(5);
    L1cMexTesting.assert_raises(f, err_exp);
    
end

function  check_get_double_from_struct_or_fail()
    err_exp = 'l1c:notAField';
    ST = struct('field0', [55, 25], 'field2', pi);
    
    f = @()TMU_mex_get_double_from_struct_or_fail(ST);
    L1cMexTesting.assert_raises(f, err_exp);

    err_exp = 'l1c:notDoubleScalar';
    ST = struct('field1', [55, 25], 'field2', pi);
    f = @()TMU_mex_get_double_from_struct_or_fail(ST);
    L1cMexTesting.assert_raises(f, err_exp);
    
    
    ST = struct('field1', [55], 'field2', pi);
    f = @()TMU_mex_get_double_from_struct_or_fail(ST);
    L1cMexTesting.assert_error_free(f);
    
end

function check_get_double_array_or_fail()

    err_exp = 'l1c:not1DArray';
    f = @()TMU_mex_get_double_array_or_fail(5);
    L1cMexTesting.assert_raises(f, err_exp);

    f = @()TMU_mex_get_double_array_or_fail([5, 6], 1);
    L1cMexTesting.assert_error_free(f);
    
    err_exp = 'l1c:notDouble';
    f = @()TMU_mex_get_double_array_or_fail(int32([5, 6]), 1);
    L1cMexTesting.assert_raises(f, err_exp);    
end

function check_mex_get_double_scalar_or_fail()
    
    f = @()TMU_mex_get_double_scalar_or_fail(5, [1, 2]);
    
    err_exp = 'l1c:notDoubleScalar';
    L1cMexTesting.assert_raises(f, err_exp);
    
    f = @()TMU_mex_get_double_scalar_or_fail(5, 1);
    
    L1cMexTesting.assert_error_free(f);
end


function check_mex_assert_num_inputs()
  
  f = @()TMU_mex_assert_num_inputs(5);
  
  err_exp = 'l1c:nlhs';
  L1cMexTesting.assert_raises(f, err_exp);
  
  
  f = @()TMU_mex_assert_num_inputs(5, 2);
  L1cMexTesting.assert_error_free(f);
end









