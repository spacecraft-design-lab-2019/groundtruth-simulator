import pytest
import os, sys, inspect

# add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

import sun_sensor_math as LA

# usual convention is to name the function "test_<function_to_test>"
def test_matTimesVec():
    x = [1, 1, 1]
    M = [x, x, x]

    assert LA.matTimesVec(M, x) == [3, 3, 3]

    M = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    assert LA.matTimesVec(M, x) == x

    # ... more and more tests go here.

    # For function that have a known and desired error (i.e. incorrect bounds in multiplication)
    # you can test for that error with
    # with pytest.raises(Exception):
        # something_wrong(...)
