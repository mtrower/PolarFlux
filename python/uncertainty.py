import numpy as np

class Measurement:
    #Set array priority to override some of ndarray's ufunc binary relations
    __array_priority__ = 100

    def __init__(self, value, uncertainty):
        self.v = value
        self.u = uncertainty

    def __getitem__(self, index):
        row = index[0]
        col = index[1]
        print ("{}+/-{}".format(self.v[row, col], self.u[row, col]))
        return

    def __repr__(self):
        return 'Measurement({}, {})'.format(self.v, self.u)

    def __str__(self):
        return '{}+/-{}'.format(self.v, self.u)

    def __add__(self, right):
        if isinstance(right, (float, np.ndarray, int)):
            return Measurement(self.v + right, self.u)
        unc = np.sqrt(self.u**2 + right.u**2)
        return Measurement(self.v + right.v, unc)

    def __radd__(self, left):
        self + left

    def __sub__(self, right):
        if isinstance(right, (float, np.ndarray, int)):
            return Measurement(self.v - right, self.u)
        unc = np.sqrt(self.u**2 + right.u**2)
        return Measurement(self.v - right.v, unc)

    def __mul__(self, right):
        if isinstance(right, (float, np.ndarray, int)):
            val = self.v*right
            unc = self.u*right
            return Measurement(val, unc)
        unc = np.sqrt((self.u/self.v)**2 + (right.u/right.v)**2)
        val = self.v*right.v
        return Measurement(val, unc*val)

    def __rmul__(self, left):
        return Measurement(left*self.v, left*self.u)

    def __truediv__(self, right):
        if isinstance(right, (float, np.ndarray, int)):
            return Measurement(self.v/right, self.u/right)
        unc = np.sqrt((self.u/self.v)**2 + (right.u/right.v)**2)
        val = self.v/right.v
        return Measurement(val, unc*val)


    def __rtruediv__(self, left):
        #TODO: think about this one, it might be a power function for floats
        if isinstance(left, (float, np.ndarray, int)):
            return Measurement(self.v/left, self.u/left)
        unc = np.sqrt((self.u/self.v)**2 + (left.u/left.v)**2)
        val = left.v/left.v
        return Measurement(val, unc*val)

    def cos(x):
        return Measurement(np.cos(x.v), np.abs(np.sin(x.v)*x.u))

    def sin(x):
        return Measurement(np.sin(x.v), np.abs(np.cos(x.v)*x.u))

    def sqrt(x):
        return Measurement(np.sqrt(x.v), 0.5*np.abs(x.u/np.sqrt(x.v)))