import numpy as np

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
        S = Sensor(errormodel = LinearErrorModel.withDim(3, cov = 0.0005))
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
    Class that handles a linear error model of the form: Tx + b + W

    T - the "T" matrix in the error model (i.e. the combined misalignment + scale factor matrices)
    b - bias vector
    cov - covariance of the white noise `W`. May be scalar or matrix valued.

    In addition to the regular constructor, may also be initialized with:
    LinearErrorModel.withDim(dim = N, [optional_args])
    """
    def __init__(self, T=np.zeros((3,3)), b=np.zeros(3), cov=0):
        """
        T - the "T" matrix in the error model (i.e. the combined misalignment + scale factor matrices)
        b - bias vector
        cov - covariance of the noise. May be scalar or matrix valued.
        """
        assert T.shape[0] == T.shape[1], "T is not square."
        assert b.shape[0] == T.shape[0], "b is not compatible with T"

        self.T = T
        self.b = b
        self.cov = cov

    @classmethod
    def withDim(cls, N = 3, T = None, b = None, cov = None):
        M = cls(T = np.zeros((N, N)), b = np.zeros(N))
        if T != None:
            M.T = T
        if b != None:
            M.b = b
        if cov != None:
            M.cov = cov
        return M

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
