import numpy as np

class Sensor():
    """
    Generic sensor. Implements measure(x), which returns the
    measured value associated with ground-truth value x.
    Defaults to returning ground truth measurements unless an error model is specified.
    To specify the error model, see the parameter `errormodel` and the LinearErrorModel class.
    NOTE, we may have to implement other error models to account for temp-dependence and whatever else.

    Example Usage:

    S = Sensor(errormodel = LinearErrorModel(T = np.eye(3), b = np.ones(3), noisecov = 0.0005))
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
    class that handles a linear error model of the form:
        Tx + b + N
    """
    def __init__(self, T=np.zeros((3,3)), b=np.zeros(3), noisecov=0):
        """
        T - the "T" matrix in the error model
        """
        assert T.shape[0] == T.shape[1], "T is not square."
        assert b.shape[0] == T.shape[0], "b is not compatible with T"

        self.T = T
        self.b = b
        self.cov = noisecov # TODO do this properly

    @classmethod
    def withDim(cls, N = 3):
        return cls(T = np.zeros((N, N)), b = np.zeros(N))

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
