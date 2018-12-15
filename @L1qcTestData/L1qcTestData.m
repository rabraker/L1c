classdef L1qcTestData
% The static class provides an interface to build the test data needed by the
% c-code unit tests. To build all the data, call, e.g.,
%
% L1qcTestData.build_all_test_data(test_data_root);
  properties
  end
  
  methods (Static)
    function build_all_test_data(test_data_root)
      L1qcTestData.build_newton_test_data(test_data_root);
      L1qcTestData.cgsolve_newton_test_data(test_data_root);
      L1qcTestData.dct_test_data(test_data_root);
    end
    build_l1qc_newton_test_data(test_data_root);
    build_cgsolve_test_data(test_data_root);
    build_dct_test_data(test_data_root);
  end
  
end