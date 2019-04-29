C
C      This is the LINES()/POWDER subroutines
C      these again do full octants, but use the powder algorithm of ALderman
C
C      SUBROUTINES AND FUNCTIONS INCLUDED
C
C          SUBROUTINE POWDERSACE(SPEC,POINTS,FSTART,FWIDTH)
C          SUBROUTINE POWDERS(SPEC,POINTS,FSTART,FWIDTH)
C          SUBROUTINE POWDERSLU(SPEC,POINTS,FSTART,FWIDTH,PARAM)
C          SUBROUTINE PREPOWDERS(PARAM)
C
C
        SUBROUTINE POWDERS(SPEC,POINTS,FSTART,FWIDTH)
C
C            DONALD W. ALDERMAN
C            Department of Chemistry
C            University of Utah
C
C       VERSION:
C            January 7, 1986
C
C       CALLING ARGUMENTS:
C            SPEC:    REAL  Array to which computed spectrum is added.
C            POINTS:  INTEGER  Number of pooints in SPEC array,
C                                must be greater than three,
C                                typically 1024 or greater.
C            FSTART:  REAL  Frequency of start of spectrum.
C            FWIDTH:  REAL  Width of spectrum, must be positive.
C
C       SYSTEM:
C            This program has been compiled and tested under Digital
C            VAX/VMS Fortran Version 4.2.
C
C       DESCRIPTION:
C            POWDER computes the powder sample NMR spectrum which results
C            from a single line whose position and amplitude is a function
C            of the orientation of each small single crystal.  The user
C            specifies the nature of this dependence by supplying a
C            subroutine LINE(COSX,COSY,COSZ,FREQ,AMP).  The arguments
C            COSX, COSY, and COSZ, supplied by POWDER as inputs to LINE,
C            are the direction cosines specifying the orientation of the
C            magnetic field vector relative to the single crystal.
C            LINE must then calculate the frequency and amplitude of the
C            NMR line when the magnetic field is so oriented and return
C            the result to POWDER in the arguments FREQ and AMP.
C            The arguments COSX, COSY, COXZ, FREQ, and AMP, are all REAL.
C            If LINE returns a frequency which is beyond the limits of the
C            spectrum specified by the arguments FSTART and FWIDTH,
C            POWDER types a message and aborts.
C            POWDER does not zero the array SPEC before computing the 
C            spectrum.  Thus multiple calls to POWDER can be used to
C            accumulate multiline powder spectra.
C            POWDER works in frequency, rather than field, coordinates.
C            Thus the initial points of the spectrum corespond to the
C            lowest frequency and the final points to the highest frequency.
C            If the user wants to display the spectrum in the usual 
C            old fashioned, backwards manner with field increasing to the 
C            right, he must reverse the spectrum.
C            POWDER sums contributions only over a half-sphere
C            on the assumption that the sign of the direction of the
C            magnetic field does not alter the line position or amplitude.
C            The vertical scale of the spectrum produced by POWDER is
C            arbitrary and depends both on the number of points in the
C            spectrum, specified by the calling argument POINTS, and on the
C            parameter NT described below.  It is up to the user to scale
C            the resultant spectrum to meet his requirements.  Two spectra 
C            produced by two calls to POWDER with the same values of POINTS 
C            and NT bear the correct relation to one another.
C            POWDER achieves its efficiency by using a two dimensional
C            interpolation scheme and by dividing the sphere using the
C            symmetry of an octahedron.
       IMPLICIT NONE
C
C            NT is the number of intervals into which the edge of the 
C            octahedron is divided to produce triangular grid on each face.
C            The total number of fequency-amplitude calculations required
C            is 2*NT*NT+1.   A value of 32 is typical for this parameter.
C            The calcuation may be done faster at the expense of accuracy
C            by decreasing NT.  In fact, if significant broadening is to be
C            done, NT can be made surprisingly small and the spectrum will
C            still be acceptable.  In any case, if NT is dereased it is 
C            advised that the spectrum be carefully examined to insure 
C            that it meets requirements for accuracy.
C            The calculation may be done more accurately at the expense 
C            of slower execution by increasing NT.
C
C Declaration of argument variables.
        INTEGER NT,POINTS
        PARAMETER (NT=32)
	REAL SPEC(0:POINTS-1),FSTART,FWIDTH
C
C Declaration of internal variables.
        REAL FREQ(0:NT,0:2*NT),AMP(0:NT,0:2*NT)
        INTEGER I,J
        REAL X,Y,Z,R,FINC,R2
C
C FIRST EXECUTABLE STATEMENT
C CHECK ARGUMENT VALUES.                          
C       IF(POINTS.GE.3.AND.FWIDTH.GT.0.) GO TO 10
C       TYPE*,'POINTS.LT.3 or FWIDTH.LE.0, POWDER aborting'
C       CALL EXIT
C COMPUTE FREQUENCIES AND AMPLITUDES AT TRIANGULAR GRID INTERSECTIONS
C ON FACES OF OCTAHEDRON
        DO 40 J=0,NT-1
         DO 20 I=0,NT-J
C        THIS SET OF I & J LOOP THROUGH A FACE IN THE FIRST OCTANT
C        THAT IS  X,Y & Z ALL >0
          X=NT-I-J
          Y=I
          Z=J
          R=SQRT(X*X+Y*Y+Z*Z)
	  R2=SQRT(X*X+Y*Y)
          CALL LINES(X/R2,Y/R2,Z/R,R2/R,FREQ(I,J),AMP(I,J))
          AMP(I,J)=AMP(I,J)/R/R/R
 20      CONTINUE
         DO 30 I=NT-J+1,NT
          X=NT-I-J
          Y=NT-J           
          Z=NT-I
          R=SQRT(X*X+Y*Y+Z*Z)
	  R2=SQRT(X*X+Y*Y)
          CALL LINES(X/R2,Y/R2,Z/R,R2/R,FREQ(I,J),AMP(I,J))
          AMP(I,J)=AMP(I,J)/R/R/R
 30      CONTINUE
 40     CONTINUE
        DO 70 J=NT,2*NT-1
         DO 50 I=J-NT+1,NT-1
          X=-NT-I+J
          Y=NT-J
          Z=NT-I
          R=SQRT(X*X+Y*Y+Z*Z)
	  R2=SQRT(X*X+Y*Y)
          CALL LINES(X/R2,Y/R2,Z/R,R2/R,FREQ(I,J),AMP(I,J))
          AMP(I,J)=AMP(I,J)/R/R/R
 50      CONTINUE
         DO 60 I=1,J-NT
          X=-NT-I+J
          Y=-I     
          Z=2*NT-J
          R=SQRT(X*X+Y*Y+Z*Z)
	  R2=SQRT(X*X+Y*Y)
          CALL LINES(X/R2,Y/R2,Z/R,R2/R,FREQ(I,J),AMP(I,J))
          AMP(I,J)=AMP(I,J)/R/R/R
 60      CONTINUE
 70     CONTINUE
        CALL LINES(0.,0.,1.,0.,FREQ(0,NT),AMP(0,NT))
        AMP(0,NT)=AMP(0,NT)/NT/NT/NT
        DO 80 J=0,NT-1
        FREQ(0,2*NT-J)=FREQ(0,J)
        AMP(0,2*NT-J)=AMP(0,J)
 80    CONTINUE
       DO 90 I=0,NT
        FREQ(NT,NT+I)=FREQ(I,0)
        AMP(NT,NT+I)=AMP(I,0)
 90    CONTINUE
       DO 100 J=1,NT-1
        FREQ(NT-J,2*NT)=FREQ(NT,J)
        AMP(NT-J,2*NT)=AMP(NT,J)
 100   CONTINUE
C
C      FORM SPECTRUM FROM FREQUENCIES AND AMPLITUDES AT TRIANGULAR GRID 
C      INTERSECTIONS ON FACES OF OCTAHEDRON BY ADDING "TENTS" TO SPECTRUM
       FINC=FWIDTH/REAL(POINTS-1)
       DO 130 I=0,NT-1
        DO 110 J=0,NT-1
         CALL TENT(FREQ(I+1,J),FREQ(I,J+1),FREQ(I,J),
     1    AMP(I+1,J)+AMP(I,J+1)+AMP(I,J),SPEC,
     1    POINTS,FSTART,FWIDTH,FINC)
         CALL TENT(FREQ(I+1,J),FREQ(I,J+1),FREQ(I+1,J+1),
     1    AMP(I+1,J)+AMP(I,J+1)+AMP(I+1,J+1),SPEC,
     1    POINTS,FSTART,FWIDTH,FINC)
 110    CONTINUE
        DO 120 J=NT,2*NT-1
         CALL TENT(FREQ(I,J),FREQ(I+1,J+1),FREQ(I+1,J),
     1    AMP(I,J)+AMP(I+1,J+1)+AMP(I+1,J),SPEC,
     1    POINTS,FSTART,FWIDTH,FINC)
         CALL TENT(FREQ(I,J),FREQ(I+1,J+1),FREQ(I,J+1),
     1    AMP(I,J)+AMP(I+1,J+1)+AMP(I,J+1),SPEC,
     1    POINTS,FSTART,FWIDTH,FINC)
 120    CONTINUE
 130   CONTINUE
       RETURN
       END
