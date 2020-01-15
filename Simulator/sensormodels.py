import numpy as np

class SpacecraftSensors():
    """
    A class to initialize and store spacecraft sensors
    """
    def __init__(self, mag_params, gyro_params, sun_params):
        self.magnetometer = withParams(mag_params)
        self.gyroscope = withParams(gyro_params)
        self.sunsensor = withParams(sun_params)

    def withParams(self, params):
        sensor = Sensor()
        sensor.errormodel.T = getTmatrix(params["scalefactor"], params["crossaxis_sensitivity"])
        sensor.errormodel.b = params["b"]
        sensor.errormodel.cov = params["cov"]
        return sensor
        

class Sensor():
    """
    Generic sensor. Implements measure(x), which returns the
    measured value associated with ground-truth value x.
    Defaults to returning ground truth measurements unless an error model is specified.
    To specify the error model, see the parameter `errormodel` and the LinearErrorModel class.
    NOTE, we may have to implement other error models to account for temp-dependence and whatever else.

    - dim = dimension of values the sensor returns
    - errormodel - the actual model used to "jitter" the measurement. Defaults to identity via the default LinearErrorModel.
    - name - optional name for the sensor

    Example Usage:
        S = Sensor(errormodel = LinearErrorModel.withDim(3, b = np.ones(3), cov = 0.0005))
        S.measure(np.zeros(3)) -> returns e.g.: array([0.99649388, 0.9849549 , 1.04930531])
    """
    def __init__(self, dim = 3, errormodel = None, name = None):
        if errormodel == None:
            self.errormodel = LinearErrorModel.withDim(N = dim)
        else:
            self.errormodel = errormodel

        self.name = name

    def measure(self, x):
        return self.errormodel.measure(x)


class LinearErrorModel():
    """
    Class that handles a linear error model of the form: (I+T)x + b + W

    T - the "T" matrix in the error model (i.e. the combined misalignment + scale factor matrices)
    b - bias vector
    cov - covariance of the white noise `W`. May be scalar or matrix valued.

    In addition to the regular constructor, may also be initialized with:
    LinearErrorModel.withDim(dim = N, [optional_args])
    """
    def __init__(self, T=np.zeros((3,3)), b=0, cov=0, dim=3):
        self.T = T
        self.b = b
        self.cov = cov
        self.dim = dim

    @classmethod
    def withDim(cls, N = 3, T = None, b = None, cov = None):

        T = np.zeros((N, N)) if T is None else T
        b = np.zeros(N) if b is None else b
        cov = 0 if cov is None else cov

        return cls(T, b, cov)

    def measure(self, x):
        """
        Add measurement error to the value `x`. x must be vector valued.
        TODO consider whether we should also have a scalar valued version
        """
        assert len(x.shape) == 1, "x is not vector valued. Got {}".format(x)
        assert x.shape[0] == self.T.shape[1], "x is not compatible with the T matrix. Got {}".format(x)

        I = np.eye(self.T.shape[0])
        n = x.shape[0]
        return (I + self.T) @ x + self.b + whitenoise(self.cov, dims = n)



def whitenoise(cov = 0, dims = None):
    """
    Basically just a wrapper for np.random.multivariate_normal that returns zero mean noise.
    """
    if isinstance(cov, int) or isinstance(cov, float):
        if dims == None:
            raise(Exception("Bad inputs to whitenoise, cov = {}, dims = {}".format(cov, dims)))
        cov = cov * np.eye(dims)
    else:
        dims = cov.shape[0]

    return np.random.multivariate_normal(np.zeros(dims), cov)

def whitenoise2(mean = 0, cov = 0, dims = None):
    """
    Basically just a wrapper for np.random.multivariate_normal that returns nonzero mean noise.
    """
    if isinstance(cov, int) or isinstance(cov, float):
        if dims == None:
            raise(Exception("Bad inputs to whitenoise, cov = {}, dims = {}".format(cov, dims)))
        cov = cov * np.eye(dims)
    else:
        dims = cov.shape[0]

    return np.random.multivariate_normal(mean*np.ones(dims), cov)
    

def getTmatrix(scalefactor, crossaxis_sensitivity):
    """
    Inputs:
        scalefactor: scale factor for a particular sensor
        crossaxis_sensitivity: cross-axis sensitivity for particular sensor
    Outputs:
        T-matrix, combined misalignment and scaling matrix for linear error model
    """
    scaleFmat = np.eye(3) + np.diag(whitenoise(scalefactor,3))
    misalign = np.reshape(np.random.multivariate_normal(np.zeros(9),crossaxis_sensitivity*np.eye(9)),(3,3))
    np.fill_diagonal(misalign, 0)
    T = np.dot(scaleFmat,misalign)
    return T