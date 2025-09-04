'''
Various functions for the optimizer properties.
'''
import numpy as np

class MatrixManipulation:
    '''
    class to hold matrix manipulation functions.
    '''
    @staticmethod
    def multiply_values_of_matrix(matrix: np.ndarray, multiplier: float, only_positive=False, only_negative=False) -> np.ndarray:
        """
        Multiply all values in the matrix by a specified multiplier.

        Parameters
        ----------
        matrix : np.ndarray
            The input matrix with values to be modified.
        multiplier : float
            The value to multiply all positive elements.
        only_positive : bool
            If only_positive is True, only positive elements will be multiplied.
        only_negative : bool    
            If only_negative is True, only negative elements will be multiplied.

        Returns
        -------
        np.ndarray
            A new matrix with positive values increased.
        """
        # Create a copy of the matrix to avoid modifying the original
        modified_matrix = np.copy(matrix)
        # Increase positive values
        if only_positive:
            modified_matrix[modified_matrix > 0] *= multiplier
        elif only_negative:
            modified_matrix[modified_matrix < 0] *= multiplier
        else:
            modified_matrix *= multiplier
        return modified_matrix

    @staticmethod
    def increment_values_of_matrix(matrix: np.ndarray, increment: float, only_positive=False, only_negative=False) -> np.ndarray:
        """
        Increment all values in the matrix by a specified amount.

        Parameters
        ----------
        matrix : np.ndarray
            The input matrix with values to be modified.
        increment : float
            The value to add to all elements.
        only_positive : bool
            If only_positive is True, only positive elements will be multiplied.
        only_negative : bool    
            If only_negative is True, only negative elements will be multiplied.    

        Returns
        -------
        np.ndarray
            A new matrix with all values incremented.
        """
        modified_matrix = np.copy(matrix)
        if only_positive:
            modified_matrix[modified_matrix > 0] += increment
        elif only_negative:
            modified_matrix[modified_matrix < 0] += increment
        else:
            modified_matrix += increment
        return modified_matrix
    
    @staticmethod
    def scale_matrix_to_size(original_matrix: np.ndarray, target_size: tuple) -> np.ndarray:
        """
        Scale the original matrix to the target size using interpolation for maximum information preservation.
        Follows GOOSE's vectorized operation patterns for performance.

        Parameters
        ----------
        original_matrix : np.ndarray
            The input matrix to be scaled.
        target_size : tuple
            The desired size for the output matrix (rows, cols).

        Returns
        -------
        np.ndarray
            A new matrix scaled to the target size with preserved information.
            
        Raises
        ------
        ValueError
            If target_size is not a 2-tuple or contains non-positive values.
        """
        # Input validation following GOOSE patterns
        if not isinstance(target_size, tuple) or len(target_size) != 2:
            raise ValueError("target_size must be a tuple of (rows, cols)")
        if target_size[0] <= 0 or target_size[1] <= 0:
            raise ValueError("target_size dimensions must be positive")
        
        # Early return if no scaling needed
        if original_matrix.shape == target_size:
            return np.copy(original_matrix)
        
        # Use scipy's zoom for high-quality scaling if available (GOOSE uses scipy elsewhere)
        try:
            from scipy.ndimage import zoom
            
            # Calculate zoom factors
            zoom_factors = (target_size[0] / original_matrix.shape[0], 
                        target_size[1] / original_matrix.shape[1])
            
            # Use order=1 (bilinear) for good balance of speed and quality
            # order=3 (cubic) would be highest quality but slower
            scaled_matrix = zoom(original_matrix, zoom_factors, order=1, mode='nearest')
            
            # Ensure exact target size (zoom might have slight rounding differences)
            if scaled_matrix.shape != target_size:
                # Crop or pad as needed
                scaled_matrix = scaled_matrix[:target_size[0], :target_size[1]]
                
            return scaled_matrix
            
        except ImportError:
            # Fallback to numpy-based interpolation if scipy not available
            return MatrixManipulation._scale_matrix_numpy_fallback(original_matrix, target_size)

    @staticmethod
    def _scale_matrix_numpy_fallback(original_matrix: np.ndarray, target_size: tuple) -> np.ndarray:
        """
        Fallback matrix scaling using pure NumPy with bilinear interpolation.
        Follows GOOSE's vectorized operation patterns.
        
        Parameters
        ----------
        original_matrix : np.ndarray
            The input matrix to be scaled.
        target_size : tuple
            The desired size for the output matrix.
            
        Returns
        -------
        np.ndarray
            Scaled matrix using bilinear interpolation.
        """
        old_rows, old_cols = original_matrix.shape
        new_rows, new_cols = target_size
        
        # Create coordinate grids for interpolation (vectorized approach)
        row_indices = np.linspace(0, old_rows - 1, new_rows)
        col_indices = np.linspace(0, old_cols - 1, new_cols)
        
        # Create meshgrid for all target positions
        col_grid, row_grid = np.meshgrid(col_indices, row_indices)
        
        # Get integer parts and fractional parts for bilinear interpolation
        row_floor = np.floor(row_grid).astype(int)
        col_floor = np.floor(col_grid).astype(int)
        row_ceil = np.minimum(row_floor + 1, old_rows - 1)
        col_ceil = np.minimum(col_floor + 1, old_cols - 1)
        
        # Fractional parts
        row_frac = row_grid - row_floor
        col_frac = col_grid - col_floor
        
        # Bilinear interpolation (vectorized)
        # Get the four corner values for each target pixel
        top_left = original_matrix[row_floor, col_floor]
        top_right = original_matrix[row_floor, col_ceil]
        bottom_left = original_matrix[row_ceil, col_floor]
        bottom_right = original_matrix[row_ceil, col_ceil]
        
        # Interpolate horizontally first, then vertically
        top = top_left * (1 - col_frac) + top_right * col_frac
        bottom = bottom_left * (1 - col_frac) + bottom_right * col_frac
        
        # Final vertical interpolation
        scaled_matrix = top * (1 - row_frac) + bottom * row_frac
        
        return scaled_matrix

    @staticmethod
    def scale_matrix_block_average(original_matrix: np.ndarray, target_size: tuple) -> np.ndarray:
        """
        Alternative scaling method using block averaging for downscaling.
        Best for preserving overall matrix statistics when making matrices smaller.
        Follows GOOSE's information preservation patterns.
        
        Parameters
        ----------
        original_matrix : np.ndarray
            The input matrix to be scaled down.
        target_size : tuple
            The desired smaller size for the output matrix.
            
        Returns
        -------
        np.ndarray
            Matrix scaled using block averaging.
            
        Notes
        -----
        This method is particularly good for epsilon matrices where you want to preserve
        the overall interaction strength patterns rather than individual pixel values.
        """
        if not isinstance(target_size, tuple) or len(target_size) != 2:
            raise ValueError("target_size must be a tuple of (rows, cols)")
        
        old_rows, old_cols = original_matrix.shape
        new_rows, new_cols = target_size
        
        # Only use block averaging for downscaling
        if new_rows >= old_rows or new_cols >= old_cols:
            return MatrixManipulation.scale_matrix_to_size(original_matrix, target_size)
        
        # Calculate block sizes
        row_block_size = old_rows / new_rows
        col_block_size = old_cols / new_cols
        
        scaled_matrix = np.zeros(target_size)
        
        for i in range(new_rows):
            for j in range(new_cols):
                # Define the block boundaries
                row_start = int(i * row_block_size)
                row_end = int((i + 1) * row_block_size)
                col_start = int(j * col_block_size)
                col_end = int((j + 1) * col_block_size)
                
                # Average the block
                scaled_matrix[i, j] = np.mean(original_matrix[row_start:row_end, col_start:col_end])
        
        return scaled_matrix

class VectorManipulation:
    """
    Class for various vector manipulation techniques.
    """
    @staticmethod
    def resize_vector(original_vector: list[float], new_size: int) -> list[float]:
        """
        Resizes a 1D vector to a new size using linear interpolation.

        Parameters
        ---------
        original_vector: The input list of floats to be resized.
        new_size: The desired size (number of elements) of the new vector.
                    Must be a positive integer.

        Returns
        -------
        A new list of floats representing the resized vector.

        Raises:
            ValueError: If the original_vector is empty or if new_size is not a
                        positive integer.

        """
        # --- Input Validation ---
        if not isinstance(new_size, int) or new_size <= 0:
            raise ValueError("new_size must be a positive integer.")
        if not original_vector:
            raise ValueError("original_vector cannot be empty.")

        original_len = len(original_vector)

        # --- Handle Edge Cases ---
        if original_len == new_size:
            return original_vector[:]  # Return a copy of the original vector

        # If the target is a single value, the average is most representative.
        if new_size == 1:
            return [sum(original_vector) / original_len]

        # If the original has only one value, we can only repeat it.
        if original_len == 1:
            return [original_vector[0]] * new_size

        # --- Main Interpolation Logic ---
        new_vector = [0.0] * new_size

        # The scaling factor maps an index in the new vector to its corresponding
        # floating-point position in the old vector. We subtract 1 from the lengths
        # to ensure the first and last elements of the old and new vectors align.
        scale_factor = (original_len - 1) / (new_size - 1)

        for i in range(new_size):
            # Calculate the floating-point position in the original vector.
            original_pos = i * scale_factor

            # Find the two surrounding indices in the original vector.
            floor_pos = int(original_pos)

            # Protect against floating-point inaccuracies that might push the
            # final index slightly out of bounds.
            if floor_pos >= original_len - 1:
                new_vector[i] = original_vector[-1]
                continue

            # Get the fractional part, which represents how far we are between
            # the 'floor_pos' and the next point.
            fraction = original_pos - floor_pos

            # Get the values of the two surrounding points.
            val1 = original_vector[floor_pos]
            val2 = original_vector[floor_pos + 1]

            # Perform the linear interpolation.
            # The new value is the first point's value plus a fraction of the
            # difference between the two points.
            interpolated_value = val1 + (val2 - val1) * fraction
            new_vector[i] = interpolated_value

        return new_vector

    def multiply_attractive_force(vector: list[float], factor: float) -> list[float]:
        """
        Multiply the attractive force of a vector by a given factor.

        Parameters
        ----------
        vector : list[float]
            The input vector to be modified.
        factor : float
            The factor by which to increase the attractive force.

        Returns
        -------
        list[float]
            The modified vector with increased attractive force.
        """
        if not isinstance(factor, (int, float)):
            raise TypeError("Factor must be a numerical value.")

        # only change attractive values, which are negative.
        return [x * factor if x < 0 else x for x in vector]
    
    def multiply_repulsive_force(vector: list[float], factor: float) -> list[float]:
        """
        Multiply the repulsive  force of a vector by a given factor.

        Parameters
        ----------
        vector : list[float]
            The input vector to be modified.
        factor : float
            The factor by which to increase the repulsive force.

        Returns
        -------
        list[float]
            The modified vector with increased repulsive force.
        """
        if not isinstance(factor, (int, float)):
            raise TypeError("Factor must be a numerical value.")

        # only change repulsive values, which are positive.
        return [x * factor if x > 0 else x for x in vector]