import numpy as np

class SpacecraftSensors():
    """
    A class to initialize and store spacecraft sensors from config parameter files.
    """
    def __init__(self, mag_params, gyro_params, sun_params):

        self.magnetometer = LinearErrorModel(T = getTmatrix(mag_params["scalefactor"], mag_params["crossaxis_sensitivity"]),
                                             b = mag_params["b"],
                                             cov = mag_params["cov"])

        self.sunsensor = LinearErrorModel(T = getTmatrix(sun_params["scalefactor"], sun_params["crossaxis_sensitivity"]),
                                          b = sun_params["b"],
                                          cov = sun_params["cov"])

        self.gyroscope = LinearErrorModel(T = getTmatrix(gyro_params["scalefactor"], gyro_params["crossaxis_sensitivity"]),
                                          b = gyro_params["b"],
                                          cov = gyro_params["cov"],
                                          random_walk_cov = gyro_params["random_walk_cov"])


class LinearErrorModel():
    """
    Class that handles a sensor with linear error model of the form: (I+T)x + b + W

    T - the "T" matrix in the error model (i.e. the combined misalignment + scale factor matrices)
    b - bias vector
    cov - covariance of the white noise `W`. May be scalar or matrix valued.

    In addition to the regular constructor, may also be initialized with:
    LinearErrorModel.withDim(dim = N, [optional_args])
    """
    def __init__(self, T=np.zeros((3, 3)), b=0, cov=0, random_walk_cov = 0):
        self.T = T
        self.b = b
        self.cov = cov

        self.random_walk_cov = random_walk_cov
        self.b_original = b

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

    def update(self):
        self.b += whitenoise(self.random_walk_cov, dims=self.T.shape[0])


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