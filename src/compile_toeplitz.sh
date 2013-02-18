echo '*******************************'
echo 'compile toeplitz.f90'
echo '*******************************'
f2py -c -m _toeplitz toeplitz.f90
mv _toeplitz.so ../