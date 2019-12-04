
import pytest
import os, sys, inspect
import math

# add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

import sun_sensor_math as LA

def pytest_addoption(parser):
    parser.addoption(
        "--cmdopt", action="store", default="type1", help="my option: type1 or type2"
    )

@pytest.fixture
def cmdopt(request):
    return request.config.getoption("--cmdopt")

# usual convention is to name the function "test_<function_to_test>"
def test_test_vecTimesMat():
    x = [1, 1, 1]
    M = [[1,2,3],[4,5,6],[7,8,9]]

    assert LA.vecTimesMat(x, M) == [12, 15, 18]

    M = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    assert LA.vecTimesMat(x, M) == x

    x = [1, 2, 3]
    M = [[1,1,1],[2,2,2],[0,0,0]]
    
    assert LA.vecTimesMat(x, M) == [5,5,5]
    
def test_matTimesVec():
    x = [1, 1, 1]
    M = [[1,0,0],[0,1,0],[0,0,1]]

    assert LA.matTimesVec(M, x) == x

    x = [1, 2, 3]
    M = [[1,0,0],[0,1,0],[0,0,1]]

    assert LA.matTimesVec(M, x) == x
    
    x = [1, 2, 3]
    M = [[8,1,6],[3,5,7],[4,9,2]]
    
    assert LA.matTimesVec(M, x) == [28,34,28]
    
def test_normalize():
    vec = [1,0,0]
    assert LA.normalize(vec) == vec
    
    vec = [0,1,0]
    assert LA.normalize(vec) == vec
    
    vec = [0,0,-1]
    assert LA.normalize(vec) == vec
    
    vec = [0,0,0]
    assert LA.normalize(vec) == vec
    
    vec = [.1,0,0]
    assert LA.normalize(vec) == [1,0,0]
    
    vec = [1,2,3]
    assert LA.normalize(vec) == [1/math.sqrt(14), 2/math.sqrt(14), 3/math.sqrt(14)]
    
    vec = [1,-2,3]
    assert LA.normalize(vec) == [1/math.sqrt(14), -2/math.sqrt(14), 3/math.sqrt(14)]
    
    vec = [3,0,4]
    assert LA.normalize(vec) == [3/5, 0, 4/5]
    
    vec = [-3,0,-4]
    assert LA.normalize(vec) == [-3/5, 0, -4/5]

def test_norm(): 
    vec = [1,0,0]
    assert LA.norm(vec) == 1
    
    vec = [3,0,4]
    assert LA.norm(vec) == 5
    
    vec = [0,0,1]
    assert LA.norm(vec) == 1
    
    vec = [1,2,3]
    assert LA.norm(vec) == math.sqrt(14)

def test_transpose():
    M = [[1,2,3],[4,5,6],[7,8,9]]
    assert LA.transpose(M) == [[1,4,7],[2,5,8],[3,6,9]]
    
    M = [[1,0,0],[0,1,0],[0,0,1]]
    assert LA.transpose(M) == M
    
def test_dot():
    v1 = [1, 2, 3]
    v2 = [1, 0, 0]
    assert LA.dot(v1, v2) == 1
    
    v1 = [1, 2, 3]
    v2 = [0, 1, 0]
    assert LA.dot(v1, v2) == 2
    
    v1 = [1, 2, 3]
    v2 = [0, 0, 1]
    assert LA.dot(v1, v2) == 3
    
    v1 = [1, 2, 3]
    v2 = [1, 1, 1]
    assert LA.dot(v1, v2) == 6
    
    v1 = [2, 2, 4]
    v2 = [1, 6, 9]
    assert LA.dot(v1, v2) == 50
    
    v1 = [-2, -2, 4]
    v2 = [1, -6, 9]
    assert LA.dot(v1, v2) == 46 
    
    v1 = [-2, 2, 4]
    v2 = [1, 6, -9]
    assert LA.dot(v1, v2) == -26
    
def test_scale():
    vec = [1, 2.02, 3.04]
    scalar = 1
    assert LA.scale(vec, scalar) == vec
    
    vec = [1.09, -2.02, 3.04]
    scalar = 0
    assert LA.scale(vec, scalar) == [0, 0, 0]
    
    vec = [1, 2.02, -3.04]
    scalar = 100
    assert LA.scale(vec, scalar) == [100, 202, -304]
    
def test_scale_error():
    with pytest.raises(TypeError):
        vec = [1, 2.02, 3.04]
        scalar = [.1,1]
        assert LA.scale(vec, scalar)
        
        scalar = [.1,1,4]
        assert LA.scale(vec, scalar)
        
        scalar = [1]
        assert LA.scale(vec, scalar)
    
def test_sub():
    v1 = [3,4,5]
    v2 = [1,1,1]
    assert LA.sub(v1,v2) == [2,3,4]
    
    v1 = [-3,-4,-5]
    v2 = [1,1,1]
    assert LA.sub(v1,v2) == [-4,-5,-6]
    
    v1 = [3,2,9]
    v2 = [0,0,0]
    assert LA.sub(v1,v2) == v1
    
    v1 = [0,0,0]
    v2 = [0,4,3]
    assert LA.sub(v1,v2) == [0,-4,-3]
    
    v1 = [0,0,0]
    v2 = [0,-4,-3]
    assert LA.sub(v1,v2) == [0,4,3]
    
def test_add():
    v1 = [3,4,5]
    v2 = [1,1,1]
    assert LA.add(v1,v2) == [4,5,6]
    
    v1 = [-3,-4,-5]
    v2 = [1,1,1]
    assert LA.add(v1,v2) == [-2,-3,-4]
    
    v1 = [3,2,9]
    v2 = [0,0,0]
    assert LA.add(v1,v2) == v1
    
    v1 = [0,0,0]
    v2 = [0,4,3]
    assert LA.add(v1,v2) == v2
    
    v1 = [0,0,0]
    v2 = [0,-4,-3]
    assert LA.add(v1,v2) == v2    
    
    # ... more and more tests go here.

    # For function that have a known and desired error (i.e. incorrect bounds in multiplication)
    # you can test for that error with
    # with pytest.raises(Exception):
        # something_wrong(...)
