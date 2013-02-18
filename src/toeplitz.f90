!http://people.sc.fsu.edu/~jburkardt/f_src/toeplitz/toeplitz.html
!http://www.osti.gov/energycitations/product.biblio.jsp?osti_id=5116519
function c4_abs1 ( x )

!*****************************************************************************80
!
!! C4_ABS1 computes the L1 absolute value of a complex number.
!
!  Modified:
!
!    29 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex X, the number whose L1 absolute value is desired.
!
!    Output, real C4_ABS1, the L1 absolute value of X.
!
  implicit none

  real c4_abs1
  complex x

  c4_abs1 = abs ( real ( x ) ) + abs ( imag ( x ) )

  return
end
subroutine c4_swap ( x, y )

!*****************************************************************************80
!
!! C4_SWAP swaps two complex values.
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, complex X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  complex x
  complex y
  complex z

  z = x
  x = y
  y = z

  return
end
subroutine caxpy ( n, ca, cx, incx, cy, incy )

!*****************************************************************************80
!
!! CAXPY adds a constant times one vector to another.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, complex CA, the multiplier.
!
!    Input, complex CX(*), the vector to be scaled and added to Y.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries of X.
!
!    Input/output, complex Y(*), the vector to which a multiple of
!    X is to be added.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive entries of Y.
!
  implicit none

  complex ca
  complex cx(*)
  complex cy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n

  if ( n <= 0 ) then

  else if ( ca == cmplx ( 0.0E+00, 0.0E+00 ) ) then

  else if ( incx == 1 .and. incy == 1 ) then

    cy(1:n) = cy(1:n) + ca * cx(1:n)

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      cy(iy) = cy(iy) + ca * cx(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine cbto_sl ( m, l, a1, a2, b, x )

!*****************************************************************************80
!
!! CBTO_SL solves the complex block Toeplitz linear system A * X = B.
!
!  Modified:
!
!    29 August 2002
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, complex A1(M**2,L), the first row of blocks of the matrix.
!    Each block is represented by columns.
!
!    Input, complex A2(M**2,L-1), the first column of blocks of the matrix,
!    beginning with the second block.  Each block is represented by columns.
!
!    Input, complex B(M*L), the right hand side vector.
!
!    Output, complex X(M*L), the solution vector.  X may coincide with B.
!
!f2py intent(in) m
!f2py intent(in) l
!f2py intent(in) a1
!f2py intent(in) a2
!f2py intent(in) b
!f2py intent(out) x

  implicit none

  integer ( kind = 4 ) l, m
  complex a1(m*m,l), a2(m*m,l-1), b(m,l), c1(m,m,l-1), c2(m,m,l-1)
  integer ( kind = 4 ) i, i1, i2, i3, ii, j, n
  integer ( kind = 4 ) pivot(m)
  complex r1(m,m), r2(m,m), r3(m,m), r5(m,m), r6(m,m), x(m,l)
!
!  Solve the system with the principal minor of order M.
!
  i3 = 1
  do j = 1, m
    do i = 1, m
      c1(i,j,1) = a1(i3,1)
      r1(i,j) = a1(i3,1)
      r3(i,j) = r1(i,j)
      i3 = i3 + 1
    end do
    x(j,1) = b(j,1)
  end do

  call cgefa ( r3, m, m, pivot, ii )

  call cgesl ( r3, m, m, pivot, x(1,1), 0 )

  if ( l == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the block Toeplitz matrix
!  for N = 2, L.
!
  do n = 2, l
!
!  Compute multiples of the first and last block columns of
!  the inverse of the principal minor of order M*N.
!
    i3 = 1
    do j = 1, m
      do i = 1, m
        r5(i,j) = a2(i3,n-1)
        r6(i,j) = a1(i3,n)
        i3 = i3 + 1
      end do
    end do

    if ( 2 < n ) then

      c1(1:m,1:m,n-1) = r2(1:m,1:m)

      do i1 = 1, n-2
        i2 = n - i1
        do j = 1, m
          i3 = 1
          do i = 1, m
            call caxpy ( m, c1(i,j,i2), a2(i3,i1), 1, r5(1,j), 1 )
            call caxpy ( m, c2(i,j,i1), a1(i3,i1+1), 1, r6(1,j), 1 )
            i3 = i3 + m
          end do
        end do
      end do

    end if

    do j = 1, m
      r2(1:m,j) = -r5(1:m,j)
      call cgesl ( r3, m, m, pivot, r2(1,j), 0 )
    end do

    r3(1:m,1:m) = r6(1:m,1:m)
    r6(1:m,1:m) = -c1(1:m,1:m,1)

    do j = 1, m
      do i = 1, m
        call caxpy ( m, r2(i,j), r3(1,i), 1, c1(1,j,1), 1 )
      end do
    end do

    call cgefa ( r6, m, m, pivot, ii )

    do j = 1, m
      call cgesl ( r6, m, m, pivot, r3(1,j), 0 )
      do i = 1, m
        call caxpy ( m, r3(i,j), r5(1,i), 1, r1(1,j), 1 )
      end do
    end do

    if ( 2 < n ) then

      r6(1:m,1:m) = c2(1:m,1:m,1)

      do i1 = 2, n-1

        if ( i1 /= n-1 ) then
          r5(1:m,1:m) = c2(1:m,1:m,i1)
        end if

        do j = 1, m
          c2(1:m,j,i1) = r6(1:m,j)
          do i = 1, m
            call caxpy ( m, r3(i,j), c1(1,i,i1), 1, c2(1,j,i1), 1 )
          end do
        end do

        do j = 1, m
          do i = 1, m
            call caxpy ( m, r2(i,j), r6(1,i), 1, c1(1,j,i1), 1 )
          end do
        end do

        r6(1:m,1:m) = r5(1:m,1:m)

      end do

    end if

    c2(1:m,1:m,1) = r3(1:m,1:m)
!
!  Compute the solution of the system with the principal minor of order M*N.
!
    r3(1:m,1:m) = r1(1:m,1:m)
    x(1:m,n) = b(1:m,n)

    do i1 = 1, n-1
      i2 = n - i1
      i3 = 1
      do i = 1, m
        call caxpy ( m, -x(i,i2), a2(i3,i1), 1, x(1,n), 1 )
        i3 = i3 + m
      end do
    end do

    call cgefa ( r3, m, m, pivot, ii )

    call cgesl ( r3, m, m, pivot, x(1,n), 0 )

    do i1 = 1, n-1
      do i = 1, m
        call caxpy ( m, x(i,n), c2(1,i,i1), 1, x(1,i1), 1 )
      end do
    end do

  end do

  return
end
subroutine cccc_sl ( a, x, r, m, l, k, lda )

!*****************************************************************************80
!
!! CCCC_SL solves the complex double column circulant system A * X = B.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input/output, complex A(M*L,K), the first row of outer blocks of the
!    CCC matrix.  Each outer block is represented by its first row
!    of inner blocks.  Each inner block is represented by its first row.
!    On return, A has been destroyed.
!
!    Input/output, complex X(M*L*K)
!    On input, the right hand side vector.
!    On output, the solution.
!
!    Workspace, complex R(max(M,2*L,2*K)).
!
!    Input, integer ( kind = 4 ) M, the order of the inner blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of inner blocks in a row or column
!    of an outer block of A.
!
!    Input, integer ( kind = 4 ) K, the number of outer blocks in a row or column
!    of the matrix A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  complex a(lda,k)
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ml
  complex r(1)
  real rk
  complex x(m,l,k)

  rk = real ( k )
  ml = m * l
!
!  Reduce the CCC matrix to a block-diagonal matrix
!  by the inverse discrete Fourier transformation.
!
  call csalw ( a, ml, k, lda, -1 )
!
!  Compute the discrete Fourier transformation of the right hand side vector.
!
  call csalw ( x, ml, k, ml, 1 )
!
!  Solve the block-diagonal system, blocks of which are CC matrices.
!
  do i3 = 1, k
    call ccc_sl ( a(1,i3), x(1,1,i3), r, m, l, m )
  end do
!
!  Solve the system by the inverse discrete Fourier transformation.
!
  call csalw ( x, ml, k, ml, -1 )

  x(1:m,1:l,1:k) = x(1:m,1:l,1:k) / rk

  return
end
subroutine cccg_sl ( a, x, r, m, l, k, lda )

!*****************************************************************************80
!
!! CCCG_SL solves the complex CCG linear system A * X = B.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input/output, complex A(M**2*L,K).
!    On input, the first row of outer blocks of the CCG matrix.
!    Each outer block is represented by its first row
!    of inner blocks.  Each inner block is represented by columns.
!    On return, A has been destroyed.
!
!    Input/output, complex X(M*L*K).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Workspace, complex R(max(M,2*L,2*K)).
!
!    Input, integer ( kind = 4 ) M, the order of the inner blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of inner blocks in a row or column
!    of an outer block of A.
!
!    Input, integer ( kind = 4 ) K, the number of outer blocks in a row or column
!    of the matrix A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  complex a(lda,k)
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mm
  complex r(1)
  real rk
  complex x(m,l,k)

  rk = real ( k )
  mm = m**2
  ml = m * l
!
!  Reduce the CCG matrix to a block-diagonal matrix
!  by the inverse discrete Fourier transformation.
!
  call csalw ( a, mm*l, k, lda, -1 )
!
!  Compute the discrete Fourier transformation of
!  the right hand side vector.
!
  call csalw ( x, ml, k, ml, 1 )
!
!  Solve the block-diagonal system, blocks of which are CG matrices.
!
  do i3 = 1, k
    call ccg_sl ( a(1,i3), x(1,1,i3), r, m, l, mm )
  end do
!
!  Compute the solution of the given system by
!  the inverse discrete Fourier transformation.
!
  call csalw ( x, ml, k, ml, -1 )

  x(1:m,1:l,1:k) = x(1:m,1:l,1:k) / rk

  return
end
subroutine ccc_sl ( a, x, r, m, l, lda )

!*****************************************************************************80
!
!! CCC_SL solves the complex column circulant system A * X = B.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input/output, complex A(M,L)
!    On input, the first row of blocks of the CC matrix.
!    Each block is represented by its first row.
!    On output, A has been destroyed.
!
!    Input/output, complex X(M*L)
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Workspace, complex R(max(M,2*L)).
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  complex a(lda,l)
  integer ( kind = 4 ) i2
  complex r(*)
  complex x(m,l)
!
!  Reduce the CC matrix to a block-diagonal matrix
!  by the inverse discrete Fourier transformation.
!
  call csalw ( a, m, l, lda, -1 )
!
!  Compute the discrete Fourier transformation of
!  the right hand side vector.
!
  call csalw ( x, m, l, m, 1 )
!
!  Solve the block-diagonal system, blocks of which are circulant matrices.
!
  do i2 = 1, l
    call cci_sl ( m, a(1,i2), x(1,i2) )
  end do
!
!  Compute the solution of the given system by
!  the inverse discrete Fourier transformation.
!
  call csalw ( x, m, l, m, -1 )

  x(1:m,1:l) = x(1:m,1:l) / real ( l )

  return
end
subroutine ccct_sl ( a, x, r, m, l, k, lda )

!*****************************************************************************80
!
!! CCCT_SL solves the complex CCT linear system A * X = B..
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input/output, complex A((2*M-1)*L,K).
!    On input, the first row of outer blocks of the CCT matrix.
!    Each outer block is represented by its first row of inner blocks.
!    Each inner block is represented by its first row followed by its
!    first column beginning with the second element.
!    On output, A has been destroyed.
!
!    Input/output, complex X(M*L*K).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Workspace, complex R(max(2*M - 2,2*L,2*K)).
!
!    Input, integer ( kind = 4 ) M, the order of the inner blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of inner blocks in a row or column
!    of an outer block of A.
!
!    Input, integer ( kind = 4 ) K, the number of outer blocks in a row or column
!    of the matrix A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  complex a(lda,k)
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) ml
  complex r(*)
  real rk
  complex x(m,l,k)

  rk = real ( k )
  m2 = 2 * m - 1
  ml = m * l
!
!  Reduce the CCT matrix to a block-diagonal matrix
!  by the inverse discrete Fourier transformation.
!
  call csalw ( a, m2*l, k, lda, -1 )
!
!  Compute the discrete Fourier transformation of
!  the right hand side vector.
!
  call csalw ( x, ml, k, ml, 1 )
!
!  Solve the block-diagonal system, blocks of which are CT matrices.
!
  do i3 = 1, k
    call cct_sl ( a(1,i3), x(1,1,i3), r, m, l, m2 )
  end do
!
!  Compute the solution of the given system by
!  the inverse discrete Fourier transformation.
!
  call csalw ( x, ml, k, ml, -1 )

  x(1:m,1:l,1:k) = x(1:m,1:l,1:k) / rk

  return
end
subroutine ccg_sl ( a, x, r, m, l, lda )

!*****************************************************************************80
!
!! CCG_SL solves the complex CG linear system A * X = B.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input/output, complex A(M**2,L).
!    On input, the first row of blocks of the CG matrix.
!    Each block is represented by columns.
!    On output, A has been destroyed.
!
!    Input/output, complex X(M*L).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Workspace, complex R(max(M,2*L)).
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  complex a(lda,l)
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ii
  complex r(*)
  complex x(m,l)
!
!  Reduce the CG matrix to a block-diagonal matrix
!  by the inverse discrete Fourier transformation.
!
  call csalw ( a, m**2, l, lda, -1 )
!
!  Compute the discrete Fourier transformation of
!  the right hand side vector.
!
  call csalw ( x, m, l, m, 1 )
!
!  Solve the block-diagonal system, blocks of which are G matrices.
!
  do i2 = 1, l
    call cgefa ( a(1,i2), m, m, r, ii )
    call cgesl ( a(1,i2), m, m, r, x(1,i2), 0 )
  end do
!
!  Compute the solution of the given system by
!  the inverse discrete Fourier transformation.
!
  call csalw ( x, m, l, m, -1 )

  x(1:m,1:l) = x(1:m,1:l) / real ( l )

  return
end
subroutine cci_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! CCI_MXV multiplies a complex circulant matrix times a vector.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex A(N), the entries of the first row of the circulant matrix.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(n)
  complex b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex x(n)

  do i = 1, n

    b(i) = cmplx ( 0.0E+00, 0.0E+00 )

    do j = 1, i-1
      b(i) = b(i) + a(n+j+1-i) * x(j)
    end do

    do j = i, n
      b(i) = b(i) + a(j+1-i) * x(j)
    end do

  end do

  return
end
subroutine cci_print ( n, a, title )

!*****************************************************************************80
!
!! CCI_PRINT prints a complex circulant matrix.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A(N), the N by N circulant matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(n)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call cci_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine cci_print_some ( n, a, ilo, jlo, ihi, jhi )

!*****************************************************************************80
!
!! CCI_PRINT_SOME prints some of a complex circulant matrix.
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A(N), the N by N circulant matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) n

  complex a(n)
  complex aij
  character ( len = 20 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) 'Columns', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+j+1-i)
        end if

        if ( aij == cmplx ( 0.0E+00, 0.0E+00 ) ) then
          ctemp(j2) = '     0.0            '
        else if ( aimag ( aij ) == 0.0E+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( aij )
        else
          write ( ctemp(j2), '(2g10.3)' ) aij
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine cci_random ( n, a )

!*****************************************************************************80
!
!! CCI_RANDOM randomizes a complex circulant matrix.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Output, complex A(N), the randomized matrix, with entries between
!    0 and 1.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(n)
  real ai(n)
  real ar(n)

  call random_number ( harvest = ar(1:n) )
  call random_number ( harvest = ai(1:n) )

  a(1:n) = cmplx ( ar(1:n), ai(1:n) )

  return
end
subroutine cci_sl ( m, a, x )

!*****************************************************************************80
!
!! CCI_SL solves the complex circulant system A * X = B.
!
!  Modified:
!
!    11 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the matrix A.
!
!    Input, complex A(M), the first row of the circulant matrix.
!
!    Input/output, complex X(M)
!    On input, the right hand side vector.
!    On output, the solution vector.
!
  implicit none

  integer ( kind = 4 ) m

  complex a(m)
  complex e
  complex e1
  complex f
  complex f1
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real p
  complex r(m)
  real ri
  real rm
  complex t
  complex t1
  real v1
  real v2
  complex x(m)

  t1 = x(1)
  x(1) = t1 / a(1)

  if ( m == 1 ) then
    return
  end if

  rm = real ( m )
!
!  Compute the inverse discrete Fourier transformation
!  of the first row of the matrix and the discrete
!  Fourier transformation of the right hand side vector.
!
  t = cmplx ( 0.0E+00, 0.0E+00 )

  do i1 = 1, m

    ri = real ( i1 - 1 )
!
!  Minimize error in forming multiples of 2 PI.
!
    p = ( ( 201.0E+00 / 32.0E+00 ) * ri &
      + 1.93530717958647692528E-03 * ri ) / real ( m )

    v1 = cos ( p )
    v2 = sin ( p )
    e = cmplx ( v1, -v2 )
    e1 = cmplx ( v1, v2 )
    f = a(1)
    f1 = t1
    do i2 = 2, m
      f = e * f + a(i2)
      f1 = e1 * f1 + x(i2)
    end do
    r(i1) = ( e1 * f1 ) / ( e * f )
    t = t + r(i1)

  end do
!
!  Compute the solution of the given system by
!  the inverse discrete Fourier transformation.
!
  x(1) = t / rm

  do i1 = 2, m

    ri = real ( i1 - 1 )
!
!  Minimize error in forming multiples of 2 PI.
!
    p = ( ( 201.0E+00 / 32.0E+00 ) * ri &
      + 1.93530717958647692528E-03 * ri ) / real ( m )

    v1 = cos ( p )
    v2 = sin ( p )
    e = cmplx ( v1, -v2 )

    f = r(1)
    do i2 = 2, m
      f = e * f + r(i2)
    end do

    x(i1) = e * f / real ( m )

  end do

  return
end
subroutine cctg_sl ( a, x, r, m, l, k, lda )

!*****************************************************************************80
!
!! CCTG_SL solves the complex CTG linear system A * X = B.
!
!  Modified:
!
!    29 August 2002
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input/output, complex A(M**2*(2*L - 1),K).
!    On input, the first row of outer blocks of the CTG matrix.
!    Each outer block is represented by its first row of inner blocks
!    followed by its first column of inner blocks beginning with the
!    second block.  Each inner block is represented by columns.
!    On output, A has been destroyed.
!
!    Input/output, complex X(M*L*K).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Workspace, complex R(max(M**2*(2*L + 3) + M,2*K)).
!
!    Input, integer ( kind = 4 ) M, the order of the inner blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of inner blocks in a row or column
!    of an outer block of A.
!
!    Input, integer ( kind = 4 ) K, the number of outer blocks in a row or column
!    of the matrix A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  complex a(lda,k)
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mm
  complex r(*)
  real rk
  complex x(m,l,k)

  rk = real ( k )
  mm = m**2
  ml = m * l
!
!  Reduce the CTG matrix to a block-diagonal matrix
!  by the inverse discrete Fourier transformation.
!
  call csalw ( a, mm*(2*l - 1), k, lda, -1 )
!
!  Compute the discrete Fourier transformation of
!  the right hand side vector.
!
  call csalw ( x, ml, k, ml, 1 )
!
!  Solve the block-diagonal system, blocks of which are block Toeplitz matrices.
!
  do i3 = 1, k
    call ctg_sl ( a(1,i3), x(1,1,i3), r, m, l, lda )
  end do
!
!  Compute the solution of the given system by
!  the inverse discrete Fourier transformation.
!
  call csalw ( x, ml, k, ml, -1 )

  x(1:m,1:l,1:k) = x(1:m,1:l,1:k) / rk

  return
end
subroutine cct_sl ( a, x, r, m, l, lda )

!*****************************************************************************80
!
!! CCT_SL solves the complex CT linear system A * X = B.
!
!  Discussion:
!
!    This routine can handle linear systems in which the matrix is
!    a complex CT matrix.  The entries are complex numbers.  The matrix
!    has a block structure.
!
!    The matrix has order L*M by L*M.
!    As a block matrix, the matrix has order L by L.
!    Each block has order M by M.
!
!    Each block is a Toeplitz matrix.
!    The blocks appear in the matrix in circulant form.
!
!    In other words, the block matrix is a circulant.  The individual
!    "entries" of the block matrix are Toeplitz matrices.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input/output, complex A(2*M-1,L).
!    On input, the first row of blocks of the CT matrix.
!    Each block is represented by its first row followed by its first
!    column beginning with the second element.
!    On output, A has been destroyed.
!
!    Input/output, complex X(M*L).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Workspace, complex R(max(2*M - 2,2*L)).
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  complex a(lda,l)
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) job
  complex r(*)
  complex x(m,l)
!
!  Reduce the CT matrix to a block-diagonal matrix by the inverse discrete
!  Fourier transform.
!
  call csalw ( a, 2*m - 1, l, lda, -1 )
!
!  Compute the discrete Fourier transform of the right hand side.
!
  call csalw ( x, m, l, m, 1 )
!
!  Solve the block-diagonal system, blocks of which are Toeplitz matrices.
!
  job = 0
  do i2 = 1, l
    call cto_sl ( m, a(1,i2), x(1,i2), r, job )
  end do
!
!  Compute the solution of the given system by the inverse discrete
!  Fourier transform.
!
  call csalw ( x, m, l, m, -1 )

  x(1:m,1:l) = x(1:m,1:l) / real ( l )

  return
end
function cdotc ( n, cx, incx, cy, incy )

!*****************************************************************************80
!
!! CDOTC forms the dot product of two vectors, conjugating the first.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, complex CX(*), one of the vectors to be multiplied.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of CX.
!
!    Input, complex CY(*), one of the vectors to be multiplied.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    elements of CY.
!
!    Output, real CDOTC, the (conjugated) dot product of CX and CY.
!
  implicit none

  complex cdotc
  complex ctemp
  complex cx(*)
  complex cy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n

  ctemp = cmplx ( 0.0E+00, 0.0E+00 )

  cdotc = cmplx ( 0.0E+00, 0.0E+00 )

  if ( n <= 0 ) then
    return
  end if

  if ( incx == 1 .and. incy == 1 ) then

    ctemp = dot_product ( conjg ( cx(1:n) ), cy(1:n) )

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      ctemp = ctemp + conjg ( cx(ix) ) * cy(iy)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  cdotc = ctemp

  return
end
subroutine cgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! CGEFA factors a complex matrix by gaussian elimination.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN 0-89871-172-X.
!
!  Parameters:
!
!    Input/output, complex A(LDA,N).
!    On input, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers used to
!    obtain it.
!    The factorization can be written A = L * U where L is a product of
!    permutation and unit lower triangular matrices and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity indicator.
!    0, normal value.
!    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
!    but it does indicate that CGESL or CGEDI will divide by zero
!    if called.  Use RCOND in CGECO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex a(lda,n)
  real c4_abs1
  integer ( kind = 4 ) icamax
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  complex t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find L = pivot index.
!
    l = icamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l
!
!  A zero pivot implies this column already triangularized.
!
    if ( c4_abs1 ( a(l,k) ) == 0.0E+00 ) then
      info = k
      cycle
    end if
!
!  Interchange if necessary.
!
    if ( l /= k ) then
      call c4_swap ( a(l,k), a(k,k) )
    end if
!
!  Compute multipliers.
!
    t = - cmplx ( 1.0E+00, 0.0E+00 ) / a(k,k)
    call cscal ( n-k, t, a(k+1,k), 1 )
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        call c4_swap ( a(l,j), a(k,j) )
      end if

      t = a(k,j)
      call caxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )

    end do

  end do

  ipvt(n) = n

  if ( c4_abs1 ( a(n,n) ) == 0.0E+00 ) then
    info = n
  end if

  return
end
subroutine cgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! CGESL solves the complex general system A * X = B.
!
!  Discussion:
!
!    The system matrix must have been factored by CGECO or CGEFA.
!
!    The routine can also solve the system (A*) * X = B.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if CGECO has set 0.0 < RCOND
!    or CGEFA has set INFO == 0.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN 0-89871-172-X.
!
!  Parameters:
!
!    Input, complex A(LDA,N), the output from CGECO or CGEFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CGECO or CGEFA.
!
!    Input/output, complex B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB, specifies the task.
!    0, solve A * X = B,
!    nonzero, solve (A*) * X = B, where (A*) is the conjugate transpose.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex a(lda,n)
  complex b(n)
  complex cdotc
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  complex t
!
!  Solve A * X = B.
!
  if ( job ==  0 ) then

    do k = 1, n - 1

      l = ipvt(k)

      if ( l /= k ) then
        call c4_swap ( b(l), b(k) )
      end if

      t = b(k)
      call caxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      call caxpy ( k-1, t, a(1,k), 1, b(1), 1 )
    end do
!
!  Solve (A*) * X = B.
!
  else

    do k = 1, n
      t = cdotc ( k-1, a(1,k), 1, b(1), 1 )
      b(k) = ( b(k) - t ) / conjg ( a(k,k) )
    end do

    do k = n-1, 1, -1

      b(k) = b(k) + cdotc ( n-k, a(k+1,k), 1, b(k+1), 1 )
      l = ipvt(k)

      if ( l /=  k ) then
        call c4_swap ( b(l), b(k) )
      end if

    end do

  end if

  return
end
subroutine csalw ( a, m, l, lda, job )

!*****************************************************************************80
!
!! CSALW Fourier transforms the rows of a complex rectangular matrix.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input/output, complex A(M,L)
!    On input, the matrix to be transformed.
!    On output, the transformed matrix.
!
!    Input, integer ( kind = 4 ) M, L, the number of rows and columns of A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) JOB.
!    +1, for direct Fourier transform.
!    -1, for inverse Fourier transform.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda

  complex a(lda,l)
  complex e
  complex f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) job
  integer ( kind = 4 ) m
  real p
  complex r1(l)
  complex r2(l)
  real ri
  real v1
  real v2

  if ( l <= 1 ) then
    return
  end if

  r1(1) = cmplx ( 1.0E+00, 0.0E+00 )

  do i1 = 2, l

    ri = real ( i1 - 1 )
!
!  Minimize error in forming multiples of 2 PI.
!
    p = ( ( 201.0E+00 / 32.0E+00 ) * ri &
      + 1.93530717958647692528E-03 * ri ) / real ( l )

    v1 = cos ( p )
    v2 = sin ( p )
    if ( job == -1 ) then
      v2 = -v2
    end if

    r1(i1) = cmplx ( v1, v2 )

  end do

  do i = 1, m

    do i1 = 1, l
      e = r1(i1)
      f = a(i,1)
      do i2 = 2, l
        f = e * f + a(i,i2)
      end do
      r2(i1) = e * f
    end do

    a(i,1:l) = r2(1:l)

  end do

  return
end
subroutine cscal ( n, ca, cx, incx )

!*****************************************************************************80
!
!! CSCAL scales a complex vector by a constant.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, complex CA, the multiplier.
!
!    Input/output, complex CX(*), the vector to be multiplied.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of CX.
!
  implicit none

  complex ca
  complex cx(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    return
  end if

  if ( incx == 1 ) then

    cx(1:n) = ca * cx(1:n)

  else

    do i = 1, n * incx, incx
      cx(i) = ca * cx(i)
    end do

  end if

  return
end
subroutine csscal ( n, sa, cx, incx )

!*****************************************************************************80
!
!! CSSCAL scales a complex vector by a real constant.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real SA, the multiplier.
!
!    Input/output, complex CX(N), the vector to be multiplied.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of CX.
!
  implicit none

  complex cx(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nincx
  real sa

  if ( n <= 0 ) then
    return
  end if

  if ( incx == 1 ) then

    do i = 1, n
      cx(i) = cmplx ( sa * real ( cx(i) ), sa * aimag ( cx(i) ) )
    end do

  else

    nincx = n * incx
    do i = 1, nincx, incx
      cx(i) = cmplx ( sa * real ( cx(i) ), sa * aimag ( cx(i) ) )
    end do

  end if

  return
end
subroutine ctg_sl ( a, x, r, m, l, lda )

!*****************************************************************************80
!
!! CTG_SL solves a linear system involving a complex TG matrix.
!
!  Modified:
!
!    29 August 2002
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input, complex A(LDA,2*L-1), an M**2 by 2*L-2 array, containing
!    the first row of blocks of the CTG matrix, followed by its first
!    column of blocks beginning with the second block.  Each block is
!    represented by columns.
!
!    Input/output, complex A(M*L).  On input, the right hand side vector,
!    and on output the solution.
!
!    Workspace, complex R(M**2*(2*L+3)+M).
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be at least M**2.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  complex a(lda,2*l-1)
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mml
  integer ( kind = 4 ) mml1
  integer ( kind = 4 ) mml2
  integer ( kind = 4 ) mml3
  integer ( kind = 4 ) mml4
  integer ( kind = 4 ) mml5
  complex r(*)
  complex x(m,l)

  mm = m * m
  mml = mm*(l - 1) + 1
  mml1 = 2*mml - 1
  mml2 = mml1 + mm
  mml3 = mml2 + mm
  mml4 = mml3 + mm
  mml5 = mml4 + mm

  call ctg_sl1 ( a, a(1,l+1), x, x, r, r(mml), r(mml1), r(mml2), &
    r(mml3), r(mml4), r(mml5), m, l, lda )

  return
end
subroutine ctg_sl1 ( a1, a2, b, x, c1, c2, r1, r2, r3, r5, r6, m, &
  l, lda )

!*****************************************************************************80
!
!! CTG_SL1 solves a CTG linear system.
!
!  Modified:
!
!    29 August 2002
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input, complex A1(M**2,L), the first row of blocks of the TG matrix.
!    Each block is represented by columns.
!
!    Input, complex A2(M**2,L-1), the first column of blocks of the TG
!    matrix, beginning with the second block.  Each block is
!    represented by columns.
!
!    Input, complex B(M*L), the right hand side vector.
!
!    Output, complex X(M*L), the solution vector, which may coincide
!    with B.
!
!    Workspace, complex C1(M,M,L-1), C2(M,M,L-1), R1(M,M), R2(M,M),
!    R3(M,M), R5(M,M), R6(M,M).
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column
!    of A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  complex a1(lda,l)
  complex a2(lda,l-1)
  complex b(m,l)
  complex c1(m,m,l-1)
  complex c2(m,m,l-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) pivot(m)
  complex r1(m,m)
  complex r2(m,m)
  complex r3(m,m)
  complex r5(m,m)
  complex r6(m,m)
  complex x(m,l)
!
!  Solve the system with the principal minor of order M.
!
  i3 = 1

  do j = 1, m
    do i = 1, m
      c1(i,j,1) = a1(i3,1)
      r1(i,j) = a1(i3,1)
      r3(i,j) = r1(i,j)
      i3 = i3 + 1
    end do
    x(j,1) = b(j,1)
  end do

  call cgefa ( r3, m, m, pivot, ii )

  call cgesl ( r3, m, m, pivot, x(1,1), 0 )

  if ( l == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system
!  with the TG matrix for N = 2 through L.
!
  do n = 2, l
!
!  Compute multiples of the first and last block columns of
!  the inverse of the principal minor of order M*N.
!
    n1 = n - 1
    n2 = n - 2
    i3 = 1

    do j = 1, m
      do i = 1, m
        r5(i,j) = a2(i3,n1)
        r6(i,j) = a1(i3,n)
        i3 = i3 + 1
      end do
    end do

    if ( 2 < n ) then

      do j = 1, m
        do i = 1, m
          c1(i,j,n1) = r2(i,j)
        end do
      end do

      do i1 = 1, n2
        i2 = n - i1
        do j = 1, m
          i3 = 1
          do i = 1, m
            do ii = 1, m
              r5(ii,j) = r5(ii,j) + c1(i,j,i2)*a2(i3,i1)
              r6(ii,j) = r6(ii,j) + c2(i,j,i1)*a1(i3,i1+1)
              i3 = i3 + 1
            end do
          end do
        end do
      end do

    end if

    do j = 1, m
      do i = 1, m
        r2(i,j) = -r5(i,j)
      end do
      call cgesl ( r3, m, m, pivot, r2(1,j), 0 )
    end do

    do j = 1, m
      do i = 1, m
        r3(i,j) = r6(i,j)
        r6(i,j) = -c1(i,j,1)
      end do
    end do

    do j = 1, m
      do i = 1, m
        do ii = 1, m
          c1(ii,j,1) = c1(ii,j,1) + r2(i,j)*r3(ii,i)
        end do
      end do
    end do

    call cgefa ( r6, m, m, pivot, ii )

    do j = 1, m
      call cgesl ( r6, m, m, pivot, r3(1,j), 0 )
      do i = 1, m
        do ii = 1, m
          r1(ii,j) = r1(ii,j) + r3(i,j)*r5(ii,i)
        end do
      end do
    end do

    if ( 2 < n ) then

      do j = 1, m
        do i = 1, m
          r6(i,j) = c2(i,j,1)
        end do
      end do

      do i1 = 2, n1

        if ( i1 /= n1 ) then
          do j = 1, m
            do i = 1, m
              r5(i,j) = c2(i,j,i1)
            end do
          end do
        end if

        do j = 1, m
          do i = 1, m
            c2(i,j,i1) = r6(i,j)
          end do
          do i = 1, m
            do ii = 1, m
              c2(ii,j,i1) = c2(ii,j,i1) + r3(i,j)*c1(ii,i,i1)
            end do
          end do
        end do

        do j = 1, m
          do i = 1, m
            do ii = 1, m
              c1(ii,j,i1) = c1(ii,j,i1) + r2(i,j)*r6(ii,i)
            end do
          end do
        end do

        do j = 1, m
          do i = 1, m
            r6(i,j) = r5(i,j)
          end do
        end do

      end do

    end if

    do j = 1, m
      do i = 1, m
        c2(i,j,1) = r3(i,j)
      end do
    end do
!
!  Compute the solution of the system with the
!  principal minor of order M*N.
!
    do j = 1, m
      do i = 1, m
        r3(i,j) = r1(i,j)
      end do
      x(j,n) = b(j,n)
    end do

    do i1 = 1, n1
      i2 = n - i1
      i3 = 1
      do i = 1, m
        do ii = 1, m
          x(ii,n) = x(ii,n) - x(i,i2)*a2(i3,i1)
          i3 = i3 + 1
        end do
      end do
    end do

    call cgefa ( r3, m, m, pivot, ii )
    call cgesl ( r3, m, m, pivot, x(1,n), 0 )

    do i1 = 1, n1
      do i = 1, m
        do ii = 1, m
          x(ii,i1) = x(ii,i1) + x(i,n)*c2(ii,i,i1)
        end do
      end do
    end do

  end do

  return
end
subroutine cto_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! CTO_MXV multiplies a complex Toeplitz matrix times a vector.
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex A(2*N-1), the entries of the first row of the Toeplitz
!    matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(2*n-1)
  complex b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex x(n)

  do i = 1, n

    b(i) = cmplx ( 0.0E+00, 0.0E+00 )

    do j = 1, i-1
      b(i) = b(i) + a(n+i-j) * x(j)
    end do

    do j = i, n
      b(i) = b(i) + a(j+1-i) * x(j)
    end do

  end do

  return
end
subroutine cto_print ( n, a, title )

!*****************************************************************************80
!
!! CTO_PRINT prints a complex Toeplitz matrix.
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(2*n-1)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call cto_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine cto_print_some ( n, a, ilo, jlo, ihi, jhi )

!*****************************************************************************80
!
!! CTO_PRINT_SOME prints some of a complex Toeplitz matrix.
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) n

  complex a(2*n-1)
  complex aij
  character ( len = 20 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) 'Columns', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+i-j)
        end if

        if ( aij == cmplx ( 0.0E+00, 0.0E+00 ) ) then
          ctemp(j2) = '    0.0'
        else if ( aimag ( aij ) == 0.0E+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( aij )
        else
          write ( ctemp(j2), '(2g10.3)' ) aij
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine cto_random ( n, a )

!*****************************************************************************80
!
!! CTO_RANDOM randomizes a complex Toeplitz matrix.
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Output, complex A(2*N-1), the randomized matrix, with entries between
!    0 and 1.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(2*n-1)
  real ai(2*n-1)
  real ar(2*n-1)

  call random_number ( harvest = ar(1:2*n-1) )
  call random_number ( harvest = ai(1:2*n-1) )

  a(1:2*n-1) = cmplx ( ar(1:2*n-1), ai(1:2*n-1) )

  return
end
subroutine cto_sl ( n, a, b, x, job )

!*****************************************************************************80
!
!! CTO_SL solves the complex Toeplitz system A * X = B.
!
!  Modified:
!
!    11 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex A(2*N-1), the first row of the Toeplitz matrix, followed
!    by the first column of the Toeplitz matrix, beginning with the second
!    element.
!
!    Input, complex B(N) the right hand side vector.
!
!    Output, complex X(N), the solution vector.  X and B may share the
!    same storage.
!
!    Input, integer ( kind = 4 ) JOB,
!    0 to solve A*X=B,
!    nonzero to solve A'*X=B.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(2*n-1)
  complex b(n)
  complex c1(n-1)
  complex c2(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nsub
  complex r1
  complex r2
  complex r3
  complex r5
  complex r6
  complex x(n)

  if ( n < 1 ) then
    return
  end if
!
!  Solve the system with the principal minor of order 1.
!
  r1 = a(1)
  x(1) = b(1) / r1

  if ( n == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the Toeplitz matrix.
!
  do nsub = 2, n
!
!  Compute multiples of the first and last columns of the inverse of
!  the principal minor of order NSUB.
!
    if ( job == 0 ) then
      r5 = a(n+nsub-1)
      r6 = a(nsub)
    else
      r5 = a(nsub)
      r6 = a(n+nsub-1)
    end if

    if ( 2 < nsub ) then

      c1(nsub-1) = r2

      do i = 1, nsub-2
        if ( job == 0 ) then
          r5 = r5 + a(n+i) * c1(nsub-i)
          r6 = r6 + a(i+1) * c2(i)
        else
          r5 = r5 + a(i+1) * c1(nsub-i)
          r6 = r6 + a(n+i) * c2(i)
        end if
      end do

    end if

    r2 = -r5 / r1
    r3 = -r6 / r1
    r1 = r1 + r5 * r3

    if ( 2 < nsub ) then

      r6 = c2(1)
      c2(nsub-1) = cmplx ( 0.0E+00, 0.0E+00 )
      do i = 2, nsub-1
        r5 = c2(i)
        c2(i) = c1(i) * r3 + r6
        c1(i) = c1(i) + r6 * r2
        r6 = r5
      end do

    end if

    c2(1) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
    r5 = cmplx ( 0.0E+00, 0.0E+00 )

    do i = 1, nsub-1
      if ( job == 0 ) then
        r5 = r5 + a(n+i) * x(nsub-i)
      else
        r5 = r5 + a(i+1) * x(nsub-i)
      end if
    end do

    r6 = ( b(nsub) - r5 ) / r1

    x(1:nsub-1) = x(1:nsub-1) + c2(1:nsub-1) * r6
    x(nsub) = r6

  end do

  return
end
subroutine cto_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! CTO_VXM multiplies a vector times a complex Toeplitz matrix.
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex A(2*N-1), the entries of the first row of the Toeplitz
!    matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(2*n-1)
  complex b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex x(n)

  do i = 1, n

    b(i) = cmplx ( 0.0E+00, 0.0E+00 )

    do j = 1, i
      b(i) = b(i) + a(i+1-j) * x(j)
    end do

    do j = i+1, n
      b(i) = b(i) + a(n+j-i) * x(j)
    end do

  end do

  return
end
subroutine ctrdi ( t, ldt, n, det, job, info )

!*****************************************************************************80
!
!! CTRDI computes the determinant and inverse of a complex triangular matrix.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN 0-89871-172-X.
!
!  Parameters:
!
!    Input/output, complex T(LDT,N).
!    On input, T contains the triangular matrix.  The zero elements of the
!    matrix are not referenced, and the corresponding elements of the array
!    can be used to store other information.
!    On output, T contains the inverse of the matrix, if requested.
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of T.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix T.
!
!    Output, complex DET(2), the determinant if requested.
!    The determinant = DET(1) * 10.0**DET(2)
!    with 1.0 <= c4_abs1(DET(1)) < 10.0 or DET(1) == 0.0.
!
!    Input, integer ( kind = 4 ) JOB, indicates the shape of the matrix,
!    and the task.
!    010, inverse of lower triangular matrix only.
!    011, inverse of upper triangular matrix only.
!    100, determinant only.
!    110, determinant and inverse of lower triangular matrix.
!    111, determinant and inverse of upper triangular matrix.
!
!    Output, integer ( kind = 4 ) INFO, inverse information.
!    If the inverse was requested, then INFO is:
!    0, if the system is nonsingular,
!    nonzero, if the system is singular.
!
  implicit none

  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) n

  real c4_abs1
  complex det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  complex t(ldt,n)
  complex temp
  real, parameter :: ten = 10.0E+00
!
!  Determinant
!
  if ( job / 100 /= 0 ) then

    det(1) = cmplx ( 1.0E+00, 0.0E+00 )
    det(2) = cmplx ( 0.0E+00, 0.0E+00 )

    do i = 1, n

      det(1) = t(i,i) * det(1)

      if ( c4_abs1 ( det(1) ) == 0.0E+00 ) then
        exit
      end if

      do while ( c4_abs1 ( det(1) ) < 1.0E+00 )
        det(1) = cmplx ( ten, 0.0E+00 ) * det(1)
        det(2) = det(2) - cmplx ( 1.0E+00, 0.0E+00 )
      end do

      do while ( ten <= c4_abs1 ( det(1) ) )
        det(1) = det(1) / cmplx ( ten, 0.0E+00 )
        det(2) = det(2) + cmplx ( 1.0E+00, 0.0E+00 )
      end do

    end do

  end if

  if ( mod ( job / 10, 10 ) == 0 ) then
    return
  end if
!
!  Inverse of upper triangular matrix.
!
  if ( mod ( job, 10 ) /= 0 ) then

    info = 0

    do k = 1, n

      if ( c4_abs1 ( t(k,k) ) == 0.0E+00 ) then
        info = k
        exit
      end if

      t(k,k) = cmplx ( 1.0E+00, 0.0E+00 ) / t(k,k)
      temp = -t(k,k)
      call cscal ( k-1, temp, t(1,k), 1 )

      do j = k + 1, n
        temp = t(k,j)
        t(k,j) = cmplx ( 0.0E+00, 0.0E+00 )
        call caxpy ( k, temp, t(1,k), 1, t(1,j), 1 )
      end do

    end do

  else
!
!  Inverse of lower triangular matrix.
!
    info = 0

    do k = n, 1, -1

      if ( c4_abs1 ( t(k,k) ) == 0.0E+00 ) then
        info = k
        exit
      end if

      t(k,k) = cmplx ( 1.0E+00, 0.0E+00 ) / t(k,k)
      temp = -t(k,k)

      if ( k /= n ) then
        call cscal ( n-k, temp, t(k+1,k), 1 )
      end if

      do j = 1, k-1
        temp = t(k,j)
        t(k,j) = cmplx ( 0.0E+00, 0.0E+00 )
        call caxpy ( n-k+1, temp, t(k,k), 1, t(k,j), 1 )
      end do

    end do

  end if

  return
end
subroutine c4vec_indicator ( n, a )

!*****************************************************************************80
!
!! C4VEC_INDICATOR sets a C4VEC to the indicator vector.
!
!  Discussion:
!
!    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
!
!  Modified:
!
!    04 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = cmplx ( i, -i )
  end do

  return
end
subroutine c4vec_print ( n, a, title )

!*****************************************************************************80
!
!! C4VEC_PRINT prints a C4VEC, with an optional title.
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, complex A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2g14.6)' ) i, a(i)
  end do

  return
end
subroutine c4vec_print_some ( n, x, max_print )

!*****************************************************************************80
!
!! C4VEC_PRINT_SOME prints some of a C4VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    14 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex X(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines to print.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  complex x(n)

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i8,2x,2g14.6)' ) i, x(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(i8,2x,2g14.6)' ) i, x(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i8,2x,2g14.6)' ) i, x(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,2x,2g14.6)' ) i, x(i)
    end do
    i = max_print
    write ( *, '(i8,2x,2g14.6,2x,a)' ) i, x(i), '...more entries...'

  end if

  return
end
subroutine c4vec_random ( alo, ahi, n, a )

!*****************************************************************************80
!
!! C4VEC_RANDOM returns a random C4VEC in a given range.
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the range allowed for the entries.
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Output, complex A(N), the vector of randomly chosen values.
!
  implicit none

  integer ( kind = 4 ) n

  complex a(n)
  real ahi
  real ai(n)
  real alo
  real ar(n)

  call random_number ( harvest = ai(1:n) )
  ai(1:n) = alo + ai(1:n) * ( ahi - alo )

  call random_number ( harvest = ar(1:n) )
  ar(1:n) = alo + ar(1:n) * ( ahi - alo )

  a(1:n) = cmplx ( ar(1:n), ai(1:n) )

  return
end
function icamax ( n, cx, incx )

!*****************************************************************************80
!
!! ICAMAX finds the index of element having maximum absolute value.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, complex C(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of CX.
!
!    Output, integer ( kind = 4 ) ICAMAX, the index of the element of CX of
!    maximum absolute value.
!
  implicit none

  real c4_abs1
  complex cx(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icamax
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real smax

  if ( n < 1 ) then
    icamax = 0
    return
  end if

  icamax = 1

  if ( n == 1 ) then
    return
  end if

  if ( incx == 1 ) then

    smax = c4_abs1 ( cx(1) )

    do i = 2, n
      if ( smax < c4_abs1 ( cx(i) ) ) then
        icamax = i
        smax = c4_abs1 ( cx(i) )
      end if
    end do

  else

    ix = 1
    smax = c4_abs1 ( cx(1) )
    ix = ix + incx
    do i = 2, n
      if ( smax < c4_abs1 ( cx(ix) ) ) then
        icamax = i
        smax = c4_abs1 ( cx(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function isamax ( n, x, incx )

!*****************************************************************************80
!
!! ISAMAX finds the index of the vector element of maximum absolute value.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of SX.
!
!    Output, integer ( kind = 4 ) ISAMAX, the index of the element of SX of
!    maximum absolute value.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) isamax
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real samax
  real x(*)

  if ( n <= 0 ) then

    isamax = 0

  else if ( n == 1 ) then

    isamax = 1

  else if ( incx == 1 ) then

    isamax = 1
    samax = abs ( x(1) )

    do i = 2, n

      if ( samax < abs ( x(i) ) ) then
        isamax = i
        samax = abs ( x(i) )
      end if

    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    isamax = 1
    samax = abs ( x(ix) )

    ix = ix + incx

    do i = 2, n
      if ( samax < abs ( x(ix) ) ) then
        isamax = i
        samax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
subroutine r4_random ( rlo, rhi, r )

!*****************************************************************************80
!
!! R4_RANDOM returns a random real in a given range.
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Output, real R, the randomly chosen value.
!
  implicit none

  logical, save :: seed = .false.
  real r
  real rhi
  real rlo
  real t

  if ( .not. seed ) then
    call random_seed ( )
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R.
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r4vec_indicator ( n, a )

!*****************************************************************************80
!
!! R4VEC_INDICATOR sets a real vector to the indicator vector.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i )
  end do

  return
end
subroutine r4vec_print ( n, a, title )

!*****************************************************************************80
!
!! R4VEC_PRINT prints a real vector.
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
  end do

  return
end
subroutine r4vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R4VEC_PRINT_SOME prints "some" of a real vector.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    10 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i8,2x,g14.6)' ) i, a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(i8,2x,g14.6)' ) i, a(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i8,2x,g14.6)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,2x,g14.6)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(i8,2x,g14.6,2x,a)' ) i, a(i), '...more entries...'

  end if

  return
end
subroutine r4vec_random ( alo, ahi, n, a )

!*****************************************************************************80
!
!! R4VEC_RANDOM returns a random real vector in a given range.
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the range allowed for the entries.
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Output, real A(N), the vector of randomly chosen values.
!
  implicit none

  integer ( kind = 4 ) n

  real a(n)
  real ahi
  real alo

  call random_number ( harvest = a(1:n) )

  a(1:n) = alo + a(1:n) * ( ahi - alo )

  return
end
function samax ( n, x, incx )

!*****************************************************************************80
!
!! SAMAX returns the maximum absolute value of the entries in a vector.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Output, real SAMAX, the maximum absolute value of an element of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real samax
  real x(*)

  if ( n <= 0 ) then

    samax = 0.0E+00

  else if ( n == 1 ) then

    samax = abs ( x(1) )

  else if ( incx == 1 ) then

    samax = abs ( x(1) )

    do i = 2, n
      if ( samax < abs ( x(i) ) ) then
        samax = abs ( x(i) )
      end if
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    samax = abs ( x(ix) )
    ix = ix + incx

    do i = 2, n
      if ( samax < abs ( x(ix) ) ) then
        samax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
subroutine saxpy ( n, sa, x, incx, y, incy )

!*****************************************************************************80
!
!! SAXPY adds a constant times one vector to another.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real SA, the multiplier.
!
!    Input, real X(*), the vector to be scaled and added to Y.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Input/output, real Y(*), the vector to which a multiple of X is to
!    be added.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    entries of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n
  real sa
  real x(*)
  real y(*)

  if ( n <= 0 ) then

  else if ( sa == 0.0E+00 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    y(1:n) = y(1:n) + sa * x(1:n)

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      y(iy) = y(iy) + sa * x(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine sbto_mxv ( m, l, a1, a2, x, b )

!*****************************************************************************80
!
!! SBTO_MXV computes the real block Toeplitz matrix product A * X = B.
!
!  Discussion:
!
!    The full matrix has order M * L, and can be regarded
!    as comprising L by L blocks.  Each block is of order M.
!
!    Example:
!
!      M = 2, L = 3
!
!      1 2 | 3 4 | 5 6
!      5 5 | 6 6 | 7 7
!      ----+-----+-----
!      7 8 | 1 2 | 3 4
!      8 8 | 5 5 | 6 6
!      ----+-----+-----
!      9 0 | 7 8 | 1 2
!      9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ 91, 134, 73, 125, 97, 129 /)
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!
!    Input, real X(M*L), the vector to be multiplied.
!
!    Output, real B(M*L), the product vector, A * X.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real a1(m,m,l)
  real a2(m,m,l-1)
  real b(m,l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real x(m,l)
!
!  Construct the right hand side by blocks.
!
  do i = 1, l

    b(1:m,i) = 0.0E+00

    do j = 1, i-1
      b(1:m,i) = b(1:m,i) + matmul ( a2(1:m,1:m,i-j), x(1:m,j) )
    end do

    do j = i, l
      b(1:m,i) = b(1:m,i) + matmul ( a1(1:m,1:m,j+1-i), x(1:m,j) )
    end do

  end do

  return
end
subroutine sbto_print ( m, l, a1, a2, title )

!*****************************************************************************80
!
!! SBTO_PRINT prints a block Toeplitz matrix.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real a1(m,m,l)
  real a2(m,m,l-1)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sbto_print_some ( m, l, a1, a2, 1, 1, m*l, m*l )

  return
end
subroutine sbto_print_some ( m, l, a1, a2, ilo, jlo, ihi, jhi )

!*****************************************************************************80
!
!! SBTO_PRINT_SOME prints some of a block Toeplitz matrix.
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real a1(m,m,l)
  real a2(m,m,l-1)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3hi
  integer ( kind = 4 ) i3lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) j3hi
  integer ( kind = 4 ) j3lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) n

  n = m * l
!
!  Print the columns of the matrix, in strips of 5.
!
  do j3lo = jlo, jhi, incx

    j3hi = j3lo + incx - 1
    j3hi = min ( j3hi, n )
    j3hi = min ( j3hi, jhi )

    inc = j3hi + 1 - j3lo

    write ( *, '(a)' ) ' '

    do j = j3lo, j3hi
      j3 = j + 1 - j3lo
      write ( ctemp(j3), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j3), j3 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i3lo = max ( ilo, 1 )
    i3hi = min ( ihi, n )

    do i = i3lo, i3hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j3 = 1, inc

        j = j3lo - 1 + j3
!
!  i = M * ( i1 - 1 ) + i2
!  j = M * ( j1 - 1 ) + j2
!
        i1 = ( i - 1 ) / m + 1
        i2 = i - m * ( i1 - 1 )
        j1 = ( j - 1 ) / m + 1
        j2 = j - m * ( j1 - 1 )

        if ( i1 <= j1 ) then
          aij = a1(i2,j2,j1+1-i1)
        else
          aij = a2(i2,j2,i1-j1)
        end if

        write ( ctemp(j3), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j3), j3 = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine sbto_sl ( m, l, a1, a2, b, x )

!*****************************************************************************80
!
!! SBTO_SL solves the real block Toeplitz linear system A * X = B.
!
!  Discussion:
!
!    The full matrix has order M * L, and can be regarded
!    as comprising L by L blocks.  Each block is of order M.
!
!    Example:
!
!      M = 2, L = 3
!
!      1 2 | 3 4 | 5 6
!      5 5 | 6 6 | 7 7
!      ----+-----+-----
!      7 8 | 1 2 | 3 4
!      8 8 | 5 5 | 6 6
!      ----+-----+-----
!      9 0 | 7 8 | 1 2
!      9 9 | 8 8 | 5 5
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, real A1(M*M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  Each block is represented by columns.
!
!    Input, real A2(M*M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!    Each block is represented by columns.
!
!    Input, real B(M*L), the right hand side vector.
!
!    Output, real X(M*L), the solution vector.  X and B may share storage.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real a1(m*m,l)
  real a2(m*m,l-1)
  real b(m,l)
  real c1(m,m,l-1)
  real c2(m,m,l-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) pivot(m)
  real r1(m,m)
  real r2(m,m)
  real r3(m,m)
  real r5(m,m)
  real r6(m,m)
  real x(m,l)
!
!  Solve the system with the principal minor of order M.
!
  i3 = 1
  do j = 1, m
    do i = 1, m
      c1(i,j,1) = a1(i3,1)
      r1(i,j) = a1(i3,1)
      i3 = i3 + 1
    end do
  end do

  r3(1:m,1:m) = r1(1:m,1:m)
  x(1:m,1) = b(1:m,1)

  call sgefa ( r3, m, m, pivot, info )

  call sgesl ( r3, m, m, pivot, x(1,1), 0 )

  if ( l == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system
!  with the block Toeplitz matrix for N = 2 through L.
!
  do n = 2, l
!
!  Compute multiples of the first and last block columns of
!  the inverse of the principal minor of order M*N.
!
    i3 = 1
    do j = 1, m
      do i = 1, m
        r5(i,j) = a2(i3,n-1)
        r6(i,j) = a1(i3,n)
        i3 = i3 + 1
      end do
    end do

    if ( 2 < n ) then

      c1(1:m,1:m,n-1) = r2(1:m,1:m)

      do i1 = 1, n-2
        i2 = n - i1
        do j = 1, m
          i3 = 1
          do i = 1, m
            call saxpy ( m, c1(i,j,i2), a2(i3,i1), 1, r5(1,j), 1 )
            call saxpy ( m, c2(i,j,i1), a1(i3,i1+1), 1, r6(1,j), 1 )
            i3 = i3 + m
          end do
        end do
      end do

    end if

    do j = 1, m
      r2(1:m,j) = -r5(1:m,j)
      call sgesl ( r3, m, m, pivot, r2(1,j), 0 )
    end do

    r3(1:m,1:m) = r6(1:m,1:m)
    r6(1:m,1:m) = -c1(1:m,1:m,1)

    do j = 1, m
      do i = 1, m
        call saxpy ( m, r2(i,j), r3(1,i), 1, c1(1,j,1), 1 )
      end do
    end do

    call sgefa ( r6, m, m, pivot, info )

    do j = 1, m
      call sgesl ( r6, m, m, pivot, r3(1,j), 0 )
      do i = 1, m
        call saxpy ( m, r3(i,j), r5(1,i), 1, r1(1,j), 1 )
      end do
    end do

    if ( 2 < n ) then

      r6(1:m,1:m) = c2(1:m,1:m,1)

      do i1 = 2, n-1

        if ( i1 /= n-1 ) then
          r5(1:m,1:m) = c2(1:m,1:m,i1)
        end if

        do j = 1, m
          c2(1:m,j,i1) = r6(1:m,j)
          do i = 1, m
            call saxpy ( m, r3(i,j), c1(1,i,i1), 1, c2(1,j,i1), 1 )
          end do
        end do

        do j = 1, m
          do i = 1, m
            call saxpy ( m, r2(i,j), r6(1,i), 1, c1(1,j,i1), 1 )
          end do
        end do

        r6(1:m,1:m) = r5(1:m,1:m)

      end do

    end if

    c2(1:m,1:m,1) = r3(1:m,1:m)
!
!  Compute the solution of the system with the principal minor of order M*N.
!
    r3(1:m,1:m) = r1(1:m,1:m)
    x(1:m,n) = b(1:m,n)

    do i1 = 1, n-1
      i2 = n - i1
      i3 = 1
      do i = 1, m
        call saxpy ( m, -x(i,i2), a2(i3,i1), 1, x(1,n), 1 )
        i3 = i3 + m
      end do
    end do

    call sgefa ( r3, m, m, pivot, info )

    call sgesl ( r3, m, m, pivot, x(1,n), 0 )

    do i1 = 1, n-1
      do i = 1, m
        call saxpy ( m, x(i,n), c2(1,i,i1), 1, x(1,i1), 1 )
      end do
    end do

  end do

  return
end
subroutine sbto_to_sge ( m, l, a1, a2, lda, n, a )

!*****************************************************************************80
!
!! SBTO_TO_SGE converts a block Toeplitz matrix to a Linpack General matrix.
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the SBTO matrix.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of the
!    SBTO matrix.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the SBTO matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the SBTO matrix, beginning with the second block.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the GE matrix.
!
!    Output, integer ( kind = 4 ) N, the order of the GE matrix.
!
!    Output, real A(LDA,N), the N by N GE matrix.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m

  real a(lda,m*l)
  real a1(m,m,l)
  real a2(m,m,l-1)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n

  n = m * l

  do i = 1, n

    i1 = ( i - 1 ) / m + 1
    i2 = i - m * ( i1 - 1 )

    if ( debug ) then
      write ( *, '(a,3i8)' ) 'I:', i, i1, i2
    end if

    do j = 1, n

      j1 = ( j - 1 ) / m + 1
      j2 = j - m * ( j1 - 1 )

      if ( debug ) then
        write ( *, '(a,3i8)' ) 'J:', j, j1, j2
      end if

      if ( i1 <= j1 ) then
        a(i,j) = a1(i2,j2,j1+1-i1)
      else
        a(i,j) = a2(i2,j2,i1-j1)
      end if

    end do

  end do

  return
end
subroutine sbto_vxm ( m, l, a1, a2, x, b )

!*****************************************************************************80
!
!! SBTO_VXM computes the real block Toeplitz matrix product X * A = B.
!
!  Discussion:
!
!    The full matrix has order M * L, and can be regarded
!    as comprising L by L blocks.  Each block is of order M.
!
!    Example:
!
!      M = 2, L = 3
!
!      1 2 | 3 4 | 5 6
!      5 5 | 6 6 | 7 7
!      ----+-----+-----
!      7 8 | 1 2 | 3 4
!      8 8 | 5 5 | 6 6
!      ----+-----+-----
!      9 0 | 7 8 | 1 2
!      9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ ? /)
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column of A.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!
!    Input, real X(M*L), the vector to be multiplied.
!
!    Output, real B(M*L), the product vector, X * A.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real a1(m,m,l)
  real a2(m,m,l-1)
  real b(m,l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real x(m,l)
!
!  Construct the right hand side by blocks.
!
  do i = 1, l

    b(1:m,i) = 0.0E+00

    do j = 1, i
      b(1:m,i) = b(1:m,i) + matmul ( transpose ( a1(1:m,1:m,i+1-j) ), x(1:m,j) )
    end do

    do j = i+1, l
      b(1:m,i) = b(1:m,i) + matmul ( transpose ( a2(1:m,1:m,j-i) ), x(1:m,j) )
    end do

  end do

  return
end
subroutine scc_qr ( a, q, s, m, l, ldq, lds )

!*****************************************************************************80
!
!! SCC_QR computes the QR factorization of a real M by L column circulant matrix.
!
!  Discussion:
!
!    The factorization has the form A * inverse(R) = Q.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input, real A(M), the first column of the column-circulant matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrices A and Q.
!    M must be at least as large as L.
!
!    Input, integer ( kind = 4 ) L, the number of columns of the matrices A and Q
!    and the order of the upper triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDQ, the leading dimension of Q.
!
!    Input, integer ( kind = 4 ) LDS, the leading dimension of S.
!
!    Output, real Q(LDQ,L), the M by L matrix Q of the factorization.
!    The columns of Q are orthonormal.
!
!    Output, real S(LDS,L), the L by L inverse of the R matrix of the
!    factorization.  Elements below the main diagonal are not accessed.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) ldq
  integer ( kind = 4 ) lds
  integer ( kind = 4 ) m

  real a(m)
  real c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) ji
  real q(ldq,l)
  real s(lds,l)
  real scale
  real sdot
  real snrm2
!
!  Initialization.
!  The last column of Q is used as a work vector.
!
  q(1:m,1) = a(1:m)
  q(1:m,l) = a(1:m)
!
!  Recurrent process for the lattice algorithm with normalization.
!
  do j1 = 1, l

    j = j1 + 1
    scale = 1.0E+00 / snrm2 ( m, q(1,j1), 1 )

    if ( j1 /= l ) then

      c = - scale * ( q(m,j1) * q(1,l) + &
        sdot ( m-1, q(1,j1), 1, q(2,l), 1 ) ) / snrm2 ( m, q(1,l), 1 )

      q(1,j) = q(m,j1) + c * q(1,l)
      do i = 2, m
        q(i,j) = q(i-1,j1) + c * q(i,l)
      end do

      if ( j /= l ) then
        q(1,l) = q(1,l) + c * q(m,j1)
        call saxpy ( m-1, c, q(1,j1), 1, q(2,l), 1 )
      end if

      s(1,j) = c

      if ( 2 < j ) then
        do i = 2, j1
          ji = j - i
          s(i,j) = s(i-1,j1) + c * s(ji,j1)
        end do
      end if

    end if

    call sscal ( m, scale, q(1,j1), 1 )
    s(j1,j1) = 1.0E+00
    call sscal ( j1, scale, s(1,j1), 1 )

  end do

  return
end
function scnrm2 ( n, x, incx )

!*****************************************************************************80
!
!! SCNRM2 computes the unitary norm of a complex vector.
!
!  Discussion:
!
!    The original SCNRM2 algorithm is accurate but written in a bizarre,
!    unreadable and obsolete format.  This version goes for clarity.
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, complex X(*), the vector whose norm is to be computed.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Output, real SCNRM2, the unitary norm of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real scnrm2
  real stemp
  complex x(*)

  if ( n <= 0 ) then
    scnrm2 = 0.0E+00
    return
  end if

  if ( 0 <= incx ) then
    ix = 1
  else
    ix = ( - n + 1 ) * incx + 1
  end if

  stemp = 0.0E+00

  do i = 1, n
    stemp = stemp + conjg ( x(ix) ) * x(ix)
    ix = ix + incx
  end do

  scnrm2 = sqrt ( stemp )

  return
end
function sdot ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! SDOT forms the dot product of two vectors.
!
!  Modified:
!
!    02 June 2000
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real X(*), one of the vectors to be multiplied.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Input, real Y(*), one of the vectors to be multiplied.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    elements of Y.
!
!    Output, real SDOT, the dot product of X and Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n
  real sdot
  real stemp
  real x(*)
  real y(*)

  if ( n <= 0 ) then

    sdot = 0.0E+00

  else if ( incx == 1 .and. incy == 1 ) then

    sdot = dot_product ( x(1:n), y(1:n) )

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    stemp = 0.0E+00
    do i = 1, n
      stemp = stemp + x(ix) * y(iy)
      ix = ix + incx
      iy = iy + incy
    end do

    sdot = stemp

  end if

  return
end
subroutine sgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! SGEFA factors a real matrix.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN 0-89871-172-X.
!
!  Parameters:
!
!    Input/output, real A(LDA,N).
!    On intput, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers used to
!    obtain it.
!    The factorization can be written A=L*U, where L is a product of
!    permutation and unit lower triangular matrices, and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity indicator.
!    0, normal value.
!    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
!    but it does indicate that SGESL or SGEDI will divide by zero if called.
!    Use RCOND in SGECO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real a(lda,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) isamax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find L = pivot index.
!
    l = isamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0E+00 ) then
      info = k
      cycle
    end if
!
!  Interchange if necessary.
!
    if ( l /= k ) then
      t = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Compute multipliers.
!
    t = -1.0E+00 / a(k,k)
    call sscal ( n-k, t, a(k+1,k), 1 )
!
!  Row elimination with column indexing.
!
    do j = k+1, n
      t = a(l,j)
      if ( l /= k ) then
        a(l,j) = a(k,j)
        a(k,j) = t
      end if
      call saxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
    end do

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0E+00 ) then
    info = n
  end if

  return
end
subroutine sgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! SGESL solves a real general linear system A * X = B.
!
!  Discussion:
!
!    SGESL can solve either of the systems A * X = B or transpose ( A ) * X = B.
!
!    The system matrix must have been factored by SGECO or SGEFA.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if SGECO has set 0.0E+00 < RCOND
!    or SGEFA has set INFO == 0.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN 0-89871-172-X.
!
!  Parameters:
!
!    Input, real A(LDA,N), the output from SGECO or SGEFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from SGECO or SGEFA.
!
!    Input/output, real B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, solve A * X = B;
!    nonzero, solve transpose ( A ) * X = B.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real a(lda,n)
  real b(n)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real sdot
  real t
!
!  Solve A * X = B.
!
  if ( job == 0 ) then

    do k = 1, n-1

      l = ipvt(k)
      t = b(l)

      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if

      call saxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      call saxpy ( k-1, t, a(1,k), 1, b(1), 1 )
    end do

  else
!
!  Solve transpose ( A ) * X = B.
!
    do k = 1, n
      t = sdot ( k-1, a(1,k), 1, b(1), 1 )
      b(k) = ( b(k) - t ) / a(k,k)
    end do

    do k = n-1, 1, -1

      b(k) = b(k) + sdot ( n-k, a(k+1,k), 1, b(k+1), 1 )
      l = ipvt(k)

      if ( l /= k ) then
        t = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
function snrm2 ( n, x, incx )

!*****************************************************************************80
!
!! SNRM2 computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    The original SNRM2 algorithm is accurate but written in a bizarre,
!    unreadable and obsolete format.  This version goes for clarity.
!
!  Modified:
!
!    01 June 2000
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real X(*), the vector whose norm is to be computed.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Output, real SNRM2, the Euclidean norm of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real samax
  real snrm2
  real stemp
  real x(*)
  real xmax

  if ( n <= 0 ) then

    snrm2 = 0.0E+00

  else

    xmax = samax ( n, x, incx )

    if ( xmax == 0.0E+00 ) then

      snrm2 = 0.0E+00

    else

      if ( 0 <= incx ) then
        ix = 1
      else
        ix = ( - n + 1 ) * incx + 1
      end if

      stemp = 0.0E+00
      do i = 1, n
        stemp = stemp + ( x(ix) / xmax )**2
        ix = ix + incx
      end do

      snrm2 = xmax * sqrt ( stemp )

    end if

  end if

  return
end
subroutine sscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! SSCAL scales a vector by a constant.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real SA, the multiplier.
!
!    Input/output, real X(*), the vector to be scaled.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real sa
  real x(*)

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    x(1:m) = sa * x(1:m)

    do i = m+1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
subroutine sto_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! STO_MXV multiplies a Toeplitz matrix times a vector.
!
!  Modified:
!
!    06 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real A(2*N-1), the entries of the first row of the Toeplitz
!    matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real a(2*n-1)
  real b(n)
  integer ( kind = 4 ) i
  real x(n)

  do i = 1, n
    b(i) = dot_product ( a(n+i-1:n+1:-1), x(1:i-1) ) &
         + dot_product ( a(1:n+1-i), x(i:n) )
  end do

  return
end
subroutine sto_print ( n, a, title )

!*****************************************************************************80
!
!! STO_PRINT prints a Toeplitz matrix.
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) n

  real a(2*n-1)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sto_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine sto_print_some ( n, a, ilo, jlo, ihi, jhi )

!*****************************************************************************80
!
!! STO_PRINT_SOME prints some of a Toeplitz matrix.
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real a(2*n-1)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) 'Columns', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+i-j)
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine sto_random ( n, a )

!*****************************************************************************80
!
!! STO_RANDOM randomizes a Toeplitz matrix.
!
!  Modified:
!
!    06 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Output, real A(2*N-1), the randomized matrix, with entries between
!    0 and 1.
!
  implicit none

  integer ( kind = 4 ) n

  real a(2*n-1)

  call random_number ( harvest = a(1:2*n-1) )

  return
end
subroutine sto_sl ( n, a, b, x, job )

!*****************************************************************************80
!
!! STO_SL solves the real Toeplitz system A * X = B.
!
!  Modified:
!
!    11 March 2001
!
!  Author:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Oleg Arushanian, MK Samarin, Valentin Voevodin, Evgeny Tyrtyshnikov,
!    Burton Garbow, James Boyle, Wayne Cowell, Kenneth Dritz,
!    The TOEPLITZ Package User's Guide,
!    Argonne National Laboratory,
!    ANL-83-16, 1983.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real A(2*N-1), the first row of the Toeplitz matrix, followed by
!    the first column of the Toeplitz matrix beginning with the second element.
!
!    Input, real B(N) the right hand side vector.
!
!    Output, real X(N), the solution vector.  X and B may share the
!    same storage.
!
!    Input, integer ( kind = 4 ) JOB,
!    0 to solve A*X=B,
!    nonzero to solve Transpose(A)*X=B.
!
  implicit none

  integer ( kind = 4 ) n

  real a(2*n-1)
  real b(n)
  real c1(n-1)
  real c2(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nsub
  real r1
  real r2
  real r3
  real r5
  real r6
  real x(n)

  if ( n < 1 ) then
    return
  end if
!
!  Solve the system with the principal minor of order 1.
!
  r1 = a(1)
  x(1) = b(1) / r1

  if ( n == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the Toeplitz matrix.
!
  do nsub = 2, n
!
!  Compute multiples of the first and last columns of the inverse of
!  the principal minor of order NSUB.
!
    if ( job == 0 ) then
      r5 = a(n+nsub-1)
      r6 = a(nsub)
    else
      r5 = a(nsub)
      r6 = a(n+nsub-1)
    end if

    if ( 2 < nsub ) then

      c1(nsub-1) = r2

      do i = 1, nsub-2
        if ( job == 0 ) then
          r5 = r5 + a(n+i) * c1(nsub-i)
          r6 = r6 + a(i+1) * c2(i)
        else
          r5 = r5 + a(i+1) * c1(nsub-i)
          r6 = r6 + a(n+i) * c2(i)
        end if
      end do

    end if

    r2 = - r5 / r1
    r3 = - r6 / r1
    r1 = r1 + r5 * r3

    if ( 2 < nsub ) then

      r6 = c2(1)
      c2(nsub-1) = 0.0E+00

      do i = 2, nsub-1
        r5 = c2(i)
        c2(i) = c1(i) * r3 + r6
        c1(i) = c1(i) + r6 * r2
        r6 = r5
      end do

    end if

    c2(1) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
    r5 = 0.0E+00
    do i = 1, nsub-1
      if ( job == 0 ) then
        r5 = r5 + a(n+i) * x(nsub-i)
      else
        r5 = r5 + a(i+1) * x(nsub-i)
      end if
    end do

    r6 = ( b(nsub) - r5 ) / r1

    x(1:nsub-1) = x(1:nsub-1) + c2(1:nsub-1) * r6
    x(nsub) = r6

  end do

  return
end
subroutine sto_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! STO_VXM multiplies a vector times a Toeplitz matrix.
!
!  Modified:
!
!    06 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real A(2*N-1), the entries of the first row of the Toeplitz
!    matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  real a(2*n-1)
  real b(n)
  integer ( kind = 4 ) i
  real x(n)

  do i = 1, n

    b(i) = dot_product ( a(i:1:-1), x(1:i) ) + &
           dot_product ( a(n+1:2*n-i), x(i+1:n) )

  end do

  return
end
subroutine strdi ( t, ldt, n, det, job, info )

!*****************************************************************************80
!
!! STRDI computes the determinant and inverse of a real triangular matrix.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN 0-89871-172-X.
!
!  Parameters:
!
!    Input/output, real T(LDT,N).
!    On input, T contains the triangular matrix. The zero elements of the
!    matrix are not referenced, and the corresponding elements of the array
!    can be used to store other information.
!    On output, T contains the inverse matrix, if it was requested.
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of T.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix T.
!
!    Output, real DET(2), the determinant of the matrix, if requested.
!    The determinant = DET(1) * 10.0**DET(2), with 1.0 <= abs ( DET(1) ) < 10.0,
!    or DET(1) == 0.
!
!    Input, integer ( kind = 4 ) JOB, specifies the shape of T, and the task.
!    010, inverse of lower triangular matrix.
!    011, inverse of upper triangular matrix.
!    100, determinant only.
!    110, determinant and inverse of lower triangular.
!    111, determinant and inverse of upper triangular.
!
!    Output, integer ( kind = 4 ) INFO.
!    If the inverse was requested, then
!    0, if the system was nonsingular;
!    nonzero, if the system was singular.
!
  implicit none

  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) n

  real det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real t(ldt,n)
  real temp
  real, parameter :: ten = 10.0E+00
!
!  Determinant.
!
  if ( job / 100 /= 0 ) then

    det(1) = 1.0E+00
    det(2) = 0.0E+00

    do i = 1, n

      det(1) = t(i,i) * det(1)

      if ( det(1) == 0.0E+00 ) then
        exit
      end if

      do while ( abs ( det(1) ) < 1.0E+00 )
        det(1) = ten * det(1)
        det(2) = det(2) - 1.0E+00
      end do

      do while ( ten <= abs ( det(1) ) )
        det(1) = det(1) / ten
        det(2) = det(2) + 1.0E+00
      end do

    end do

  end if

  if ( mod ( job / 10, 10 ) == 0 ) then
    return
  end if
!
!  Inverse of an upper triangular matrix.
!
  if ( mod ( job, 10 ) /= 0 ) then

    info = 0

    do k = 1, n

      if ( t(k,k) == 0.0E+00 ) then
        info = k
        exit
      end if

      t(k,k) = 1.0E+00 / t(k,k)
      temp = -t(k,k)
      call sscal ( k-1, temp, t(1,k), 1 )

      do j = k + 1, n
        temp = t(k,j)
        t(k,j) = 0.0E+00
        call saxpy ( k, temp, t(1,k), 1, t(1,j), 1 )
      end do

    end do
!
!  Inverse of a lower triangular matrix.
!
  else

    info = 0

    do k = n, 1, -1

      if ( t(k,k) == 0.0E+00 ) then
        info = k
        exit
      end if

      t(k,k) = 1.0E+00 / t(k,k)
      temp = -t(k,k)

      if ( k /= n ) then
        call sscal ( n-k, temp, t(k+1,k), 1 )
      end if

      do j = 1, k-1
        temp = t(k,j)
        t(k,j) = 0.0E+00
        call saxpy ( n-k+1, temp, t(k,k), 1, t(k,j), 1 )
      end do

    end do

  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
