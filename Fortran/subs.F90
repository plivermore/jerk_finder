MODULE SUBS
USE sort
REAL( KIND = 8), PARAMETER :: PI = 3.14159265358979_8


CONTAINS

SUBROUTINE Find_linear_interpolated_values( k, x_min, x_max, pt, endpt, nd, grid, interpolated_signal)
IMPLICIT none
INTEGER :: K, nd
REAL( KIND = 8) ::linear_description_time(1:k+2),linear_description_intensity(1:k+2),interpolated_signal(1:nd)
REAL( KIND = 8) :: x_max, x_min, pt(:,:), endpt(:), grid(:)
INTENT(OUT) :: interpolated_signal

linear_description_time(1) = x_min
linear_description_time(2:k+1) = pt(1:k,1)
linear_description_time(k+2) = x_max

linear_description_intensity(1) = endpt(1)
linear_description_intensity(2:k+1) = pt(1:k,2)
linear_description_intensity(k+2) = endpt(2)




call interp_linear( 1, k+2, linear_description_time, linear_description_intensity, nd, &
grid(1:nd), interpolated_signal )

return
END SUBROUTINE Find_linear_interpolated_values

REAL( KIND = 8) FUNCTION randn()
IMPLICIT none
REAL( KIND = 8):: RANDOM_NUMBERS(2)

CALL RANDOM_NUMBER( RANDOM_NUMBERS(1:2) )
! Use Box_Muller transform
RANDN = SQRT( -2.0_8 * LOG( RANDOM_NUMBERS(1) ) ) * COS( 2.0_8 * Pi * RANDOM_NUMBERS(2) )
!  ANOTHER IS = SQRT( -2.0_LONG_REAL * LOG( RANDOM_NUMBERS(1) )) * SIN( 2.0_LONG_REAL * Pi * RANDOM_NUMBERS(2) )
RETURN
END FUNCTION

subroutine interp_linear ( m, data_num, t_data, p_data, interp_num, &
t_interp, p_interp )

!*****************************************************************************80
!
!! INTERP_LINEAR: piecewise linear interpolation to a curve in M dimensions.
!
!  Discussion:
!
!    From a space of M dimensions, we are given a sequence of
!    DATA_NUM points, which are presumed to be successive samples
!    from a curve of points P.
!
!    We are also given a parameterization of this data, that is,
!    an associated sequence of DATA_NUM values of a variable T.
!    The values of T are assumed to be strictly increasing.
!
!    Thus, we have a sequence of values P(T), where T is a scalar,
!    and each value of P is of dimension M.
!
!    We are then given INTERP_NUM values of T, for which values P
!    are to be produced, by linear interpolation of the data we are given.
!
!    Note that the user may request extrapolation.  This occurs whenever
!    a T_INTERP value is less than the minimum T_DATA or greater than the
!    maximum T_DATA.  In that case, linear extrapolation is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the value of the
!    independent variable at the sample points.  The values of T_DATA
!    must be strictly increasing.
!
!    Input, real ( kind = 8 ) P_DATA(M,DATA_NUM), the value of the
!    dependent variables at the sample points.
!
!    Input, integer ( kind = 4 ) INTERP_NUM, the number of points
!    at which interpolation is to be done.
!
!    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the value of the
!    independent variable at the interpolation points.
!
!    Output, real ( kind = 8 ) P_INTERP(M,DATA_NUM), the interpolated
!    values of the dependent variables at the interpolation points.
!
implicit none

integer ( kind = 4 ) data_num
integer ( kind = 4 ) m
integer ( kind = 4 ) interp_num

integer ( kind = 4 ) interp
integer ( kind = 4 ) left
real ( kind = 8 ) p_data(m,data_num)
real ( kind = 8 ) p_interp(m,interp_num)
integer ( kind = 4 ) right
real ( kind = 8 ) t
real ( kind = 8 ) t_data(data_num)
real ( kind = 8 ) t_interp(interp_num)

if ( .not. r8vec_ascends_strictly ( data_num, t_data ) ) then
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'INTERP_LINEAR - Fatal error!'
write ( *, '(a)' ) &
'  Independent variable array T_DATA is not strictly increasing. T_DATA WRITTEN TO FORT.99'
WRITE(99,*) t_data
stop 1
end if

do interp = 1, interp_num

t = t_interp(interp)
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
call r8vec_bracket ( data_num, t_data, t, left, right )

p_interp(1:m,interp) = &
( ( t_data(right) - t                ) * p_data(1:m,left)   &
+ (                 t - t_data(left) ) * p_data(1:m,right) ) &
/ ( t_data(right)     - t_data(left) )

end do

return
end subroutine interp_linear

function r8vec_ascends_strictly ( n, x )

!*****************************************************************************80
!
!! R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!    Notice the effect of entry number 6 in the following results:
!
!      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
!      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
!      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )
!
!      R8VEC_ASCENDS_STRICTLY ( X ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Y ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Z ) = TRUE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, real ( kind = 8 ) X(N), the array to be examined.
!
!    Output, logical R8VEC_ASCENDS_STRICTLY, is TRUE if the
!    entries of X strictly ascend.
!
implicit none

integer ( kind = 4 ) n

integer ( kind = 4 ) i
logical r8vec_ascends_strictly
real ( kind = 8 ) x(n)

do i = 1, n - 1
if ( x(i+1) <= x(i) ) then
r8vec_ascends_strictly = .false.
return
end if
end do

r8vec_ascends_strictly = .true.

return
end function r8vec_ascends_strictly


subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
implicit none

integer ( kind = 4 ) n

integer ( kind = 4 ) i
integer ( kind = 4 ) left
integer ( kind = 4 ) right
real ( kind = 8 ) x(n)
real ( kind = 8 ) xval

do i = 2, n - 1

if ( xval < x(i) ) then
left = i - 1
right = i
return
end if

end do

left = n - 1
right = n

return
end subroutine r8vec_bracket

FUNCTION CHECK_DIFFERENT(x1, x2, vector )
! returns 0 if any two elements are the same, 1 if they are all different.
IMPLICIT NONE
REAL( KIND = 8) :: vector(:), vector2(1:size(vector)+2),x1,x2, vector3(1:size(vector)+2)
INTEGER :: CHECK_DIFFERENT, ORDER(1: SIZE(VECTOR)+2), I

vector2(1:size(vector)) = vector(:)
vector2(size(vector)+1) = x1
vector2(size(vector)+2) = x2

order = rargsort(vector2)
do i = 1, SIZE(vector2)
vector3(i) = vector2( order(i))
enddo
vector2 = vector3

CHECK_DIFFERENT = 0
DO I = 1, SIZE(VECTOR2)-1
IF( VECTOR2(I) .eq. VECTOR2(I+1)) THEN
CHECK_DIFFERENT = 1
EXIT
ENDIF
ENDDO

END FUNCTION CHECK_DIFFERENT

END MODULE SUBS

