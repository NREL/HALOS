"""
HALOS model

Classes for Points and Vectors
and their operations


"""
import numpy as np

class Point(object):
    def __init__(self, nparray):
        self._point = nparray
        self._x = nparray[0]
        self._y = nparray[1]
        self._z = nparray[2]
    
    @property
    def point(self): return self._point
    @property
    def x(self): return self._x
    @property
    def y(self): return self._y
    @property
    def z(self): return self._z

    # use setPoint to update Point values
    def setPoint(self, nparray):
        self._point = nparray
        self._x = self._point[0]
        self._y = self._point[1]
        self._z = self._point[2]

    def Add(self, nparray):
        self._point += nparray
        self._x = self._point[0]
        self._y = self._point[1]
        self._z = self._point[2]


class Vector(object):
    def __init__(self, nparray):
        self._vec = nparray
        self._i = nparray[0]
        self._j = nparray[1]
        self._k = nparray[2]
    
    @property
    def vec(self): return self._vec
    @property
    def i(self): return self._i
    @property
    def j(self): return self._j
    @property
    def k(self): return self._k

    # use setVector to update vector values
    def setVector(self, nparray):
        self._vec = nparray
        self._i = nparray[0]
        self._j = nparray[1]
        self._k = nparray[2]

class PointVect(object):
    def __init__(self):
        self.x = 0.
        self.y = 0.
        self.z = 0.

        self.i = 0.
        self.j = 0.
        self.k = 0.

        # TODO: add methods as required



def rotation(theta, axis, PointVector):
    """
	This method takes a point, a rotation angle, and the axis of rotation and 
	rotates the point about the origin in the specified direction. 

	The inputs are:
	theta	| [rad]	| Angle of rotation
	axis	| none	| X=0, Y=1, Z=2 : Axis to rotate about

	The method returns a modified point "P" that has been rotated according to the inputs.

	Rotation is clockwise about each axis (left hand rule)(I think this is wrong). In other words, 
	positive rotation for each axis is defined by the apparent motion when positive end 
	of the axis points toward the viewer. 
    """

    # TODO: check through code to see if this should be left or right hand rule
    # These Rotation metrices follow right hand rules
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    
    R = np.zeros((3,3))     # Rotation matrix
    """
    The rotation vectors are entered in as the transverse for convenience of multiplication.
    The correct matrix form for each are:
    """
    if axis == 0:
        # X axis
        """
        [1,					0,				0,
            0,				cos(theta),		-sin(theta),
            0,				sin(theta),		cos(theta)]
        """
        R[0,0] = 1
        R[1,1] = costheta
        R[1,2] = -sintheta
        R[2,1] = sintheta
        R[2,2] = costheta
    elif axis == 1:
        # Y axis
        """
        [cos(theta),		0,			sin(theta),
            0,				1,				0,
        -sin(theta),		0,			cos(theta)]
        """
        R[0,0] = costheta
        R[0,2] = sintheta
        R[1,1] = 1
        R[2,0] = -sintheta
        R[2,2] = costheta
    elif axis == 2:
        # Z axis
        """
        [cos(theta),	-sin(theta),		0,
           sin(theta),	cos(theta),			0,
            0,				0,				1	]
        """
        R[0,0] = costheta
        R[0,1] = -sintheta
        R[1,0] = sintheta
        R[1,1] = costheta
        R[2,2] = 1
    else:
        print("Internal error: invalid axis number specified in rotation() method.")
    
    if isinstance(PointVector, Point):
        PointVector.setPoint(np.dot(R,PointVector.point))
    elif isinstance(PointVector, Vector):
        PointVector.setVector(np.dot(R,PointVector.vec))

    return



if __name__ == "__main__":

    # testing point, vector, and rotation
    test = np.array([0,1,0])
    x = Point(test)
    y = Vector(test)
    theta = np.pi

    rotation(theta,0,x)

    pass
