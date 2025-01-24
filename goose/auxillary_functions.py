from abc import ABC, abstractmethod
from typing import Any, Tuple, Dict, Callable, get_type_hints
from scipy.interpolate import RegularGridInterpolator, CubicSpline, Akima1DInterpolator, PchipInterpolator
import sparrow
import numpy as np
import finches
import inspect
from goose.goose_exceptions import GooseError, GooseInputError, GooseException



'''
-----------------
Constants
-----------------
'''

#These are the allowable methods for 1D interpolation
#These will need updating if you add a new method
INTERP_METHODS_1D = ['linear', 'cubic spline', 'akima', 'makima', 'pchip']





'''
----------------
Functions that are useful for defining comparison functions
----------------
'''
def regular_grid_normalized_interpolator(n_dim_data : np.ndarray, interpolation_shape : tuple[int], method : str = "linear") -> np.ndarray:
    '''Interpolates a new N-dimensional matrix on the input matrixes values with a given shape

    Note: the data gets squeezed to ensure it can avoid problems
    
    Parameters
    ----------
    n_dim_data : numpy.ndarray
        This is the data that you are wishing to interpolate. This is an N dimensional matrix.
    interpolation_shape : tuple[int]
        This is the shape of the interpolated matrix that you wish to have. If you have a 
        matrix you are trying to match, matched_mat.shape is all you need in this parameter.
    method : str
        This is the method you want to use for your interpolation. Any method that can be used in
        scipy.interpolate's RegularGridInterpolator can be used on this. If an invalid method is passed
        the error will be the same as if that method was passed to RegularGridInterpolator.

    Returns
    -------
    np.ndarray
        The new N dimentional matrix will have the shape specified in interpolation_shape
        except it will be squeezed.
    '''
    #check that the proper types were passed to this function
    if not isinstance(n_dim_data, np.ndarray):
        raise GooseInputError(f"The n_dim_data is not a numpy.ndarray ({type(n_dim_data)}).")

    #remove dimensions that have a value of 1 (aka redundant dimensions)
    squeezed_interp_shape = tuple([val for val in interpolation_shape if val != 1])
    #get the interpolation number of dimensions
    interp_num_dimensions = len(squeezed_interp_shape)

    #squeeze the dim_data
    squeezed_data = np.squeeze(n_dim_data)

    #determine the number of dimensions
    original_dimensions = squeezed_data.shape
    original_num_dimensions = len(original_dimensions)

    #check that the new interpolation shape has the same number of dimensions
    if original_num_dimensions != interp_num_dimensions:
        raise GooseInputError(f"The number of dimensions in n_dim_data ({original_num_dimensions}) is not the same as interpolation_shape ({interp_num_dimensions}).")

    #generate the x1,x2,x3,x4,... grid numbers for the original data
    original_grid_tup = tuple([np.linspace(0,1,val) for val in original_dimensions]) #need a tuple for the interpolator

    #generate the interpolation function
    #grid_interp_function = RegularGridInterpolator(original_grid_tup, squeezed_data, method=method)
    grid_interp_function = RegularGridInterpolator(original_grid_tup, squeezed_data)

    #All of this manipulation is just to get the new data positions in a form they can be worked with by the interpolator
    #get grid tuple for the interpolated values
    interp_grid_tup = tuple([np.linspace(0,1,val) for val in squeezed_interp_shape])
    #meshgrid the interpolated values so that it can be fed into the interpolation function
    interp_meshgrid = np.meshgrid(*interp_grid_tup, indexing='ij')
    #ravel each dimension
    raveled_positions = np.array([dim.ravel() for dim in interp_meshgrid])
    #take the transpose swap the inner and outter arrays
    raveled_transpose = raveled_positions.T


    #compute the new points in the matrix that are interpolated
    new_matrix = grid_interp_function(raveled_transpose)

    #return the new matrix
    return new_matrix.reshape(squeezed_interp_shape)


def vector_interpolator(data_vec : np.ndarray, num_interp_pts : int, method : str = 'linear') -> np.ndarray:
    '''Interpolates the vector passed based on the number of points you want the vector to be


    This is effectively a wrapper on a whole host of interpolation methods offered by numpy and
    scipy. This allows the user to find a method that will be stable for them.


    Allowed methods:
    "linear", "cubic spline", "akima", "makima", "pchip"

    
    Parameters
    ----------
    data_vec : numpy.ndarray
        This is the data that you want to interpolate. It should only have 1 dimension
    num_interp_pts : int
        This is the total number of points that you wish to get back in the returned data
    method : str
        This is the interpolation method that you wish to use
    '''
    
    #check thtat the number of points is greater than 2
    if num_interp_pts < 3:
        raise GooseInputError(f"The number of interpolated points must be >= 3.")
    
    #check that the data_vec is a numpy.ndarray
    if not isinstance(data_vec, np.ndarray):
        raise GooseInputError(f"Must pass a numpy array to data_vec not a {type(data_vec)}")

    #check that the array is n dimensional
    if data_vec.ndim != 1:
        raise GooseException(f"The data_vec numpy.array must be 1 dimension. It is currently {data_vec.ndim} dimensions.")
    
    #check that the method passed is in the allowed list
    if method not in INTERP_METHODS_1D:
        raise GooseInputError(f"A proper method must be specified. The method passed was: {method}.\nPlease choose from the allowable list: {INTERP_METHODS_1D}")

    #get the set of old and new x values to pull from
    old_x = np.linspace(0,1,len(data_vec))
    new_x = np.linspace(0,1,num_interp_pts)

    #perform the interpolation based on the method specified
    new_vals = None #initialize a value
    if method == 'linear':
        #just directly compute from the numpy function for linear interpolation
        new_vals = np.interp(new_x, old_x, data_vec)
    elif method == 'cubic spline':
        #get the cubic spline interpolation function
        cubic_spline_func = CubicSpline(old_x, data_vec)
        #perform the interpolation using the function
        new_vals = cubic_spline_func(new_x)
    elif method == 'akima' or method == 'makima':
        #get the either the akima or modified akima interpolation method based on the specified method
        ak_func = Akima1DInterpolator(old_x, data_vec, method=method)
        #perform the akima interpolation
        new_vals = ak_func(new_x)
    elif method == 'pchip':
        #get the pchip interpolator
        pchip_func = PchipInterpolator(old_x, data_vec)
        #calculate the interpolated value based on the pchip function
        new_vals = pchip_func(new_x)

    #return the interpolated data
    return new_vals