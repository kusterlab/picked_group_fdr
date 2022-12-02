Converting numpy array to vector: https://stackoverflow.com/questions/17855032/passing-and-returning-numpy-arrays-to-c-methods-via-cython

In the `.pyx` file:

```
cdef divide_wrapper(np.ndarray[np.double_t, ndim=1] col_i, np.ndarray[np.double_t, ndim=1] col_j, vector[double]& ratios):
    divide(&col_i[0], &col_j[0], len(col_i), ratios)
```

Speed ups for vector functions:

```
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
def vector_function():
```
