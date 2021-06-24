MODULE Jerks
!
!  This program defines a Python-callable function that returns Bayesian statistics on fitting a piecewise linear segment though a dataset.
! The algorithm is based on a transdimensional random walk through model space (reverse-jump Markov chain Monte Carlo)

! **********
! Inputs:
! **********
!
! sigmas: the perturbations defining the random walk used in the RJMCMC method
! in the order sigma_change_value, sigma_move, sigma_birth
!
! burn_in: the length of burn in period which is ignored for collecting statistics
!
! nsample: the length of the Markov chain including the burn_in period
!
! num_data: the size of the dataset
!
! times: an array of times that define the dataset
!
! Y: an array defining the data at the given time (see above)
!
! delta_Y: the error in Y, assumed to be the standard deviation of a normal distribution
!
! Y_min, Y_max: the bounds of the prior distribution (assumed uniform) for the range of the internal vertices
!
! times_min, times_max: the bounds on the model period, which needs to include all the data.
!
! k_min, k_max: the bounds on the prior for the number of internal vertices
!
! discretise_size: the number of grid points that defines the resolution of some of the outputs (e.g. ensemble mean, median timeseries).
!
! cp_bins: the number of bins defining the histogram of internal vertices
!
! thin: the amount of chain thinning. E.g. a value of n means that statistics are only evaluated only at every nth model
!
! nbins: the number of bins defining the model density histogram
!
! credible: credible interval bounds as percentage e.g. 95. A value of 0 means that credible intervals are not calculated.
!
! RUNNING_MODE: 1 indicates that the posterior is returned; otherwise the prior distribution is returned.
!
! **********
! OUTPUTS
! **********
!
! Acceptance_rates:  acceptance rates for change value, move (in time), birth, death
!
! CREDIBLE_SUP: the upper bound for the credible interval (array size: discretize_size), or an array of zeros if credible interval is not calculated.
!
! CREDIBLE_INF: the lower bound for the credible interval (array size: discretize_size), or an array of zeros if credible interval is not calculated.
!
! ENSEMBLE_AV: the ensemble average model (array size: discretize_size)
!
! ENSEMBLE_MEDIAN: the ensemble median model (array size: discretize_size)
!
! ENSEMBLE_MODE: the ensemble modal model (array size: discretize_size)
!
! CHANGE_POINTS: a histogram of the timing of internal vertices (array, size CP_BINS)
!
! MARGINAL_DENSITY: a 2D marginal density for the linear ensemble or zeros (if not calculated); (array, size: discretise_size,NBINS)
!
! N_changepoint_hist: a histogram of the number of internal vertices (1D array, dimensions K_MIN:K_MAX)


USE SORT
USE SUBS

CONTAINS
SUBROUTINE RJMCMC(sigmas, burn_in, NSAMPLE, NUM_DATA, TIMES, &
           Y, delta_Y, Y_MIN, Y_MAX, TIMES_MIN, &
           TIMES_MAX, K_MIN, K_MAX, discretise_size, CP_NBINS, &
           THIN, NBINS, credible, RUNNING_MODE, &
           Acceptance_rates, CREDIBLE_SUP, CREDIBLE_INF, ENSEMBLE_AV, ENSEMBLE_MEDIAN, &
           ENSEMBLE_MODE, CHANGE_POINTS, MARGINAL_DENSITY, &
           N_changepoint_hist )


IMPLICIT NONE
INTEGER, PARAMETER :: dp = 8!

! INPUT ARGUMENTS:
integer, intent(in)   :: BURN_IN, RUNNING_MODE
integer, intent(in)   :: NUM_DATA
real(dp), intent(in)  :: TIMES(NUM_DATA), Y(NUM_DATA), delta_Y(NUM_DATA)
real(dp), intent(in) :: Y_MIN, Y_MAX, TIMES_MIN, TIMES_MAX, credible
integer, intent(in)   :: thin, NBINS, K_MIN, K_MAX, discretise_size
integer, intent(in) :: NSAMPLE
real(dp), intent(in) :: sigmas(3)
INTEGER, intent(in) :: CP_NBINS

! OUTPUTS
real(dp), intent(out) :: Acceptance_rates(4)
REAL( KIND = dp), DIMENSION(1:discretise_size), intent(out) :: ENSEMBLE_AV, CREDIBLE_SUP, &
        CREDIBLE_INF, ENSEMBLE_MEDIAN, ENSEMBLE_MODE
INTEGER, intent(out) :: CHANGE_POINTS(CP_NBINS)

REAL( KIND = dp), intent(out) :: MARGINAL_DENSITY(discretise_size,NBINS)
INTEGER, intent(out) :: N_changepoint_hist(K_MIN:K_MAX)


! Internal variables

INTEGER :: b, AB, AD, PD, PB, ACV, PCV, AP, PP, PA, AA, P_sd, A_sd, k_max_array_bound, out, NUM
INTEGER :: s, birth, death, move, ind, k_prop, accept,BIN_INDEX, CHANGE_value, I, K, j, k_init
REAL(dp) :: ENDPT(2), ENDPT_PROP(2), like, RAND(2), alpha, like_prop, prob, &
            like_init, pt_death(2), u, sigma_move, sigma_change_value, sigma_birth
REAL(dp), ALLOCATABLE :: VAL_MIN(:), VAL_MAX(:),   MINI(:,:), MAXI(:,:), PT(:,:), PT_PROP(:,:)
REAL(dp), ALLOCATABLE :: interpolated_signal(:), TIME(:), PTS_NEW(:,:)
INTEGER, ALLOCATABLE, DIMENSION(:) :: IND_MIN, IND_MAX, ORDER
INTEGER, ALLOCATABLE :: discrete_history(:,:)
LOGICAL :: CALC_CREDIBLE
INTEGER :: input_random_seed
INTEGER, ALLOCATABLE :: SEED(:)

IF( credible > 0.0_dp) THEN
CALC_CREDIBLE = .TRUE.
ELSE
CALC_CREDIBLE = .FALSE.
CREDIBLE_SUP(:) = 0.0_dp
CREDIBLE_INF(:) = 0.0_dp
ENDIF

input_random_seed = 1
CALL RANDOM_SEED(SIZE = i)
ALLOCATE( SEED(1:i) )
SEED(:) = input_random_seed
CALL RANDOM_SEED(PUT = SEED)

sigma_change_value = sigmas(1)
sigma_move = sigmas(2)
sigma_birth = sigmas(3)

! Other parameters are fixed here
k_max_array_bound = k_max + 1;

NUM = ceiling((nsample-burn_in)*(100.0_8-credible)/200.0_8/thin) ! number of collected samples for credible intervals

ALLOCATE( TIME(discretise_size) )
DO I=1, discretise_size
TIME(I) = TIMES_MIN + REAL(I-1, KIND = 8)/REAL(discretise_size-1, KIND = 8) * (TIMES_MAX - TIMES_MIN)
ENDDO

ALLOCATE( val_min(discretise_size), val_max(discretise_size), &
        ind_min(discretise_size), ind_max(discretise_size) )
ALLOCATE( MINI(discretise_size, NUM), MAXI(discretise_size, NUM) )
ALLOCATE( discrete_history(discretise_size,NBINS) )
ALLOCATE( pt(k_max_array_bound,2), pt_prop(k_max_array_bound,2)  )

discrete_history(:,:) = 0
CHANGE_POINTS(:) = 0
val_min(:) = 0.0_dp
val_max(:) = 0.0_dp
ind_min(:) = 0
ind_max(:) = 0
MINI(:,:) = 0.0_dp
MAXI(:,:) = 0.0_dp
pt(:,:) = 0.0_dp
b = 0
AB=0
AD=0
PB=0
PD=0
ACV=0
PCV=0
AP=0
PP=0
PA = 0
P_sd = 0
A_sd = 0
AA = 0

! initial condition
CALL RANDOM_NUMBER( RAND(1))
k_init = floor(RAND(1) * (K_max - K_min+1)) + k_min
k = k_init

DO i=1,k_init
CALL RANDOM_NUMBER( RAND(1:2))
pt(i,1)=TIMES_min+rand(1) * (TIMES_MAX-TIMES_Min)  ! position of internal vertex
pt(i,2)=Y_min+rand(2) * (Y_MAX-Y_MIN)  ! magnitude of vertices
enddo

! Now randomly choose the end points:
CALL RANDOM_NUMBER (RAND(1:2))
endpt(1) = Y_min+RAND(1) * (Y_MAX-Y_MIN)
endpt(2) = Y_min+RAND(2) * (Y_MAX-Y_MIN)


! make sure the positions are sorted in ascending order.

ALLOCATE( ORDER(k_init), pts_new(k_init,1:2) )
! Find the order that sorts the first row of the array pt:
order = rargsort(pt(1:k_init,1))
! Sort both the values and times based on this ordering:
do i = 1, k_init
pts_new(i,1:2) = pt( order(i), 1:2)
enddo
pt(1:k_init,1:2) = pts_new(:,:)
DEALLOCATE (ORDER, pts_new)

! COMPUTE INITIAL MISFIT
k = k_init
like=0;
! interpolate. First, assemble the complete linear description
ALLOCATE( interpolated_signal(max(NUM_DATA,discretise_size)) ) !generic output space for interpolation.

CALL Find_linear_interpolated_values( k, TIMES_min, TIMES_max, pt, endpt, NUM_DATA, TIMES, interpolated_signal)

! Running mode
IF( RUNNING_MODE .eq. 1) THEN

!compute likelihood based on normal distribution
do i=1,NUM_DATA
like=like+(Y(i) - interpolated_signal(i))**2/(2.0_8 * delta_Y(i)**2)
enddo

ELSE !set likelihood as 1. The prior distribution is then returned.
like = 1.0_dp
endif

like_init=like


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%% START RJ-MCMC SAMPLING %%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do s=1,nsample



birth = 0
move  = 0
death = 0
change_value = 0


pt_prop = pt
endpt_prop = endpt
like_prop = like
k_prop = k
prob = 1.0_8
out = 1


!----------------------------------------------------------------------
! Every other iteration, propose a new value
!======================================================================
if ( mod(s,2) == 0 ) then ! Change Value
!======================================================================
if (s > burn_in)  PCV = PCV + 1
change_value = 1
k_prop = k
CALL RANDOM_NUMBER( RAND(1))
ind=ceiling(RAND(1)*(k+2))

! Check bounds to see if outside prior

if(ind == k+1)  then! change left end point
endpt_prop(1) = endpt(1) + randn()*sigma_change_value
if( endpt_prop(1) < Y_min .or. endpt_prop(1) > Y_max ) out = 0
elseif( ind == k+2) then ! change right end point
endpt_prop(2) = endpt(2) + randn()*sigma_change_value
if( endpt_prop(2) < Y_min .or. endpt_prop(2) > Y_max ) out = 0
else ! change interior point
pt_prop(ind,2) = pt(ind,2) + randn()*sigma_change_value
if( pt_prop(ind,2)>Y_max .or. pt_prop(ind,2)<Y_min) out = 0
endif

!-----------------------------------------------------------------------
! Every other iteration iteration change the vertex positions
!======================================================================
else ! Change position
!======================================================================
CALL RANDOM_NUMBER( RAND(1))
u = RAND(1) !Chose randomly between 3 different types of moves
if ( u < 0.333 ) then ! BIRTH ++++++++++++++++++++++++++++++++++++++
birth = 1
if (s>burn_in) PB=PB+1
k_prop = k+1
CALL RANDOM_NUMBER( RAND(1))
pt_prop(k+1,1)=TIMES_min+RAND(1)*(TIMES_max-TIMES_min)

! check that all the time-points are different:
DO WHILE ( CHECK_DIFFERENT(TIMES_MIN, TIMES_MAX, pt_prop(1:k+1,1)) .EQ. 1)
CALL RANDOM_NUMBER( RAND(1))
pt_prop(k+1,1)=TIMES_min+RAND(1)*(TIMES_max-TIMES_min)
END DO

!interpolate to find magnitude as inferred by current state
CALL Find_linear_interpolated_values(k, TIMES_min, TIMES_max, pt, endpt, 1, &
pt_prop(k+1:k+1,1), interpolated_signal)
pt_prop(k+1,2)=interpolated_signal(1)+randn()*sigma_birth

!GET prob
prob = (1.0_dp/(sigma_birth*sqrt(2.0_dp*pi))) * &
exp(-(interpolated_signal(1)-pt_prop(k+1,2))**2/(2.0_dp*sigma_birth**2))

!Check BOUNDS to see if outside prior
out = 1
if ((pt_prop(k+1,2)>Y_max) .OR. (pt_prop(k+1,2)<Y_min))  out = 0
if ((pt_prop(k+1,1)>TIMES_max) .OR. (pt_prop(k+1,1)<TIMES_min))  out = 0
if (k_prop>k_max) out=0

! make sure the positions are sorted in ascending order.
ALLOCATE( ORDER(k_prop), pts_new(k_prop,1:2) )
! Find the order that sorts the first row of the array pt:
order = rargsort(pt_prop(1:k_prop,1))
! Sort both the values and ages based on this ordering:
do i = 1, k_prop
pts_new(i,1:2) = pt_prop( order(i), 1:2)
enddo
pt_prop(1:k_prop,1:2) = pts_new(:,:)
DEALLOCATE (ORDER, pts_new)


elseif ( u < 0.666 ) then !  DEATH +++++++++++++++++++++++++++++++++++++++++

death = 1
if ( s > burn_in) PD = PD + 1
out = 1
k_prop = k-1

if (k_prop<k_min) out=0

if (out == 1) then
CALL RANDOM_NUMBER( RAND(1))
ind=ceiling(RAND(1)*k)
pt_death(1:2) = pt(ind,1:2)
! remove point to be deleted, shifting everything to the left.
pt_prop(1:ind-1,1:2)=pt(1:ind-1,1:2)
pt_prop(ind:k-1,1:2) = pt(ind+1:k,1:2)

!GET prob
!interpolate to find magnitude of current state as needed by birth

CALL Find_linear_interpolated_values(k_prop, TIMES_min, TIMES_max, pt_prop, endpt_prop, &
1, pt_death(1:1), interpolated_signal)

prob = 1.0_dp/(sigma_birth*sqrt(2.0_dp*pi))  * &
exp(-(interpolated_signal(1)-pt_death(2))**2/(2.0_dp*sigma_birth**2))

endif !if (out==1)

else ! MOVE +++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (s > burn_in ) PP = PP + 1
move = 1
k_prop = k
CALL RANDOM_NUMBER( RAND(1:2))
ind=ceiling(RAND(1)*k)
if(k > 0) then ! If there are no points to move, then we can't move any
pt_prop(ind,1) = pt(ind,1)+randn()*sigma_move         !Normal distribution of move destination
else
pt_prop = pt
endif  !

!Check BOUNDS
out = 1
if ( (pt_prop(ind,1) > TIMES_max) .OR. (pt_prop(ind,1) < TIMES_min) )  out = 0

!! check that all the time-points are different:
DO WHILE ( CHECK_DIFFERENT( TIMES_MIN, TIMES_MAX, pt_prop(1:k,1)) .EQ. 1)
pt_prop(ind,1) = pt(ind,1)+randn()*sigma_move
END DO


! make sure the positions are sorted in ascending order.
ALLOCATE( ORDER(k_prop), pts_new(k_prop,1:2) )

! Find the order that sorts the first row of the array pt:
order = rargsort(pt_prop(1:k_prop,1))

! Sort both the values and ages based on this ordering:
do i = 1, k_prop
pts_new(i,1:2) = pt_prop( order(i), 1:2)
enddo
pt_prop(1:k_prop,1:2) = pts_new(:,:)

DEALLOCATE (ORDER, pts_new)

ENDIF !change vertex structure
ENDIF ! decide on what proposal to make
!----------------------------------------------------------------------


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! COMPUTE MISFIT OF THE PROPOSED MODEL
! If the proposed model is not outside the bounds of the uniform prior,
! compute its misfit : "like_prop"
if (out == 1)  then
like_prop = 0.0_dp
CALL Find_linear_interpolated_values( k_prop, TIMES_min, TIMES_max, pt_prop, endpt_prop, NUM_DATA, TIMES, interpolated_signal)

IF( RUNNING_MODE .eq. 1) THEN
do i=1,NUM_DATA
like_prop=like_prop+(Y(i) - interpolated_signal(i))**2/(2.0_dp * delta_Y(i)**2)
enddo

else
like_prop = 1.0_dp
endif

endif !if (out==1)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %% SEE WHETHER MODEL IS ACCEPTED

accept = 0
alpha = 0
! The acceptance term takes different
! values according the the proposal that has been made.

if (birth == 1) then
if(out ==1) alpha = ((1.0_dp/((Y_max-Y_min)*prob))*exp(-like_prop+like))
CALL RANDOM_NUMBER( RAND(1))
if (RAND(1) < alpha) THEN
accept = 1
if ( s > burn_in) AB = AB + 1
endif

elseif (death==1) then
if(out == 1) alpha = ((Y_max-Y_min)*prob)*exp(-like_prop+like)
CALL RANDOM_NUMBER( RAND(1))
if (RAND(1)<alpha) then
accept=1
if (s>burn_in) AD = AD + 1
endif

else ! NO JUMP, i.e no change in dimension
if(out ==1) alpha = exp(-like_prop+like)
CALL RANDOM_NUMBER( RAND(1))
if (RAND(1)<alpha) then
accept=1
if (s > burn_in) then
if ( change_value .eq. 1) then
ACV = ACV + 1
elseif( move .eq. 1) then
AP = AP + 1

endif
endif !if (s>burn_in)
endif !accept?
endif !birth/death/no jump?

! If accept, update the values
if (accept == 1 ) then
k = k_prop
pt = pt_prop
like = like_prop
endpt = endpt_prop
endif

! Update acceptance rates
Acceptance_rates(1) = 100.0*ACV/PCV
Acceptance_rates(2) = 100.0*AP/PP
Acceptance_rates(3) = 100.0*AB/PB
Acceptance_rates(4) = 100.0*AD/PD


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Collect samples for the ensemble solution

if ( s > burn_in .AND. mod( s-burn_in, thin) == 0 ) THEN

b = b + 1

CALL Find_linear_interpolated_values( k, TIMES_min, TIMES_max, pt, endpt, discretise_size, &
TIME(1:discretise_size), interpolated_signal)


! build the sum of the models to calculate the average
ENSEMBLE_AV(:)=ENSEMBLE_AV(:)+interpolated_signal(1:discretise_size)


! build marginal intensity density
DO i=1,discretise_size
BIN_INDEX = FLOOR( (interpolated_signal(i)-Y_MIN)/REAL(Y_MAX-Y_MIN, KIND = dp) * NBINS ) + 1
IF( BIN_INDEX < 0 .OR. BIN_INDEX > NBINS) THEN
PRINT*, 'FATAL ERROR, BIN_INDEX IS OUT OF RANGE'
PRINT*, ' MODEL POINT ', I, ' VALUE ',interpolated_signal(i)
PRINT*, 'INTENSITY MIN/MAX ', Y_MIN, Y_MAX
STOP
ENDIF
discrete_history(i,BIN_INDEX) = discrete_history(i,BIN_INDEX) + 1
enddo


IF ( CALC_CREDIBLE ) THEN
do i=1 , discretise_size

! Do the (e.g.) 95% credible interval by keeping the lowest and greatest 2.5% of
! all models at each sample point. We could either keep ALL the data and at
! the end determine these regions, or (better) keep a running list of the
! num data points we need. At the end of the algorithm, simply take the
! maximum of the smallest points, and the min of the largest, to get the
! bounds on the credible intervals.
! Method:
! num is the number of data points corresponding to 2.5% of the total number
! of samples (after thinning).
! Collect num datapoints from the first num samples.
! For each subsequent sample, see if the value should actually be inside
! the 2.5% tail. If it is, replace an existing value by the current value.
! Repeat.

! Interestingly, this is by far the speed-bottleneck in the algorithm. As the algorithm progresses, the credible intervals converge and the
! tails need not be changed very much. This leads to an increase in speed of the algorithm as the iteration count increases.

if (b<=num) then
MINI(i,b) = interpolated_signal(i)
MAXI(i,b) = interpolated_signal(i)
if (b == num) then
val_min(i) = MINVAL(MAXI(i,:) ); ind_min(i:i) = MINLOC( MAXI(i,:) )
val_max(i) = MAXVAL(MINI(i,:) ); ind_max(i:i) = MAXLOC( MINI(i,:) )
endif
else
if (interpolated_signal(i)>val_min(i)) then
MAXI(i,ind_min(i))=interpolated_signal(i);
val_min(i) = MINVAL(MAXI(i,:) ); ind_min(i:i) = MINLOC( MAXI(i,:) )
endif
if (interpolated_signal(i)<val_max(i)) then
MINI(i,ind_max(i))=interpolated_signal(i)
val_max(i) = MAXVAL(MINI(i,:) ); ind_max(i:i) = MAXLOC( MINI(i,:) )
endif
endif

enddo !i

ENDIF !CALC_CREDIBLE


! Build histogram of number of changepoints: k
N_changepoint_hist(k)=N_changepoint_hist(k) + 1


!calculate the histogram on change points

DO i=1,k
BIN_INDEX = FLOOR( (pt(i,1)-TIMES_MIN)/REAL(TIMES_MAX-TIMES_MIN, KIND = dp) * CP_NBINS ) + 1
IF( BIN_INDEX < 0 .OR. BIN_INDEX > CP_NBINS) THEN
PRINT*, 'FATAL ERROR, BIN_INDEX IS OUT OF RANGE'
PRINT*, ' MODEL POINT ', I, ' VALUE ',pt(i,1)
PRINT*, 'TIMES MIN/MAX ', TIMES_MIN, TIMES_MAX
STOP
ENDIF
CHANGE_POINTS(BIN_INDEX) = CHANGE_POINTS(BIN_INDEX) + 1
enddo


ENDIF !collect stats

enddo ! the Sampling of the mcmc

! Compute the average
ENSEMBLE_AV(:)=ENSEMBLE_AV(:)/b

do i=1, discretise_size

IF ( CALC_CREDIBLE ) THEN
! Compute the credible intervals
CREDIBLE_SUP(i) = MINVAL(MAXI(i,:) )
CREDIBLE_INF(i) = MAXVAL(MINI(i,:) )
ENDIF

! Compute the mode
ENSEMBLE_MODE(i) = (0.5_dp + REAL(MAXLOC( discrete_history(i,:),DIM = 1)-1, KIND = dp))/NBINS * (Y_MAX-Y_MIN) + Y_MIN

! Compute the median. Get the first instance of the count from the left being greater than half the total:
do j=1, NBINS
if( sum( discrete_history(i,1:j)) .GE. sum( discrete_history(i,1:NBINS))/2) then
ENSEMBLE_MEDIAN(i) = (REAL(j-1, KIND = dp)+0.5_dp)/NBINS * (Y_MAX-Y_MIN) + Y_MIN
exit
endif
enddo

enddo

! normalise marginal distributions
MARGINAL_DENSITY(:,:) = REAL(discrete_history(:,:), KIND = dp)/ sum( discrete_history(1,:) )



RETURN
END SUBROUTINE RJMCMC



END MODULE Jerks
