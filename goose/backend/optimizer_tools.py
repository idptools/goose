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