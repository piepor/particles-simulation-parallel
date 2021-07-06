!   
!                          Program Particles 
!  mimics the behaviour of a system of particles affected by mutual forces
!
!  Final application of the Course "Parallel Computing using MPI and OpenMP"
!
!  This program is meant to be used by course participants for demonstrating
!  the abilities acquired in optimizing and parallelising programs.
!
!  The techniques learnt at the course should be extensively used in order to
!  improve the program response times as much as possible, while gaining
!  the same or very closed results, i.e. the series of final produced images
!  and statistical results.
!
!  The code implemented herewith has been written for course exercise only,
!  therefore source code, algorithms and produced results must not be 
!  trusted nor used for anything different from their original purpose. 
!  
!  Description of the program:
!  a squared grid is hit by a field whose result is the distribution of particles
!  with different properties.
!  After having been created the particles move under the effect of mutual
!  forces.
!
!  Would you please send comments to m.cremonesi@cineca.it 
!
!  Program outline:
!  1 - the program starts reading the initial values (InitGrid)
!  2 - the generating field is computed (GeneratingField)
!  3 - the set of created particles is computed (ParticleGeneration)
!  4 - the evolution of the system of particles is computed (SystemEvolution)
!

MODULE Utilities
   IMPLICIT NONE
   
   
CONTAINS

   FUNCTION char2int(txt) RESULT(ival)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: txt
      INTEGER :: ival, l
      
      l = LEN_TRIM(txt)
      IF (l < 2) THEN
         READ(txt,"(I1)") ival
      ELSE IF (l < 3) THEN
         READ(txt,"(I2)") ival
      ELSE IF (l < 4) THEN
         READ(txt,"(I3)") ival
      ELSE IF (l < 5) THEN
         READ(txt,"(I4)") ival
      ELSE IF (l < 6) THEN
         READ(txt,"(I5)") ival
      ELSE IF (l < 7) THEN
         READ(txt,"(I6)") ival
      ELSE IF (l < 8) THEN
         READ(txt,"(I7)") ival
      ELSE IF (l < 9) THEN
         READ(txt,"(I8)") ival
      ELSE IF (l < 10) THEN
         READ(txt,"(I9)") ival
      ELSE 
         READ(txt,"(I10)") ival
      END IF
   
      RETURN
   END FUNCTION char2int
   
   function char2float(s) result(f)
      implicit none
	  integer, parameter :: mn=ICHAR('0'),mx=ICHAR('9')
	  character*(*) :: s
	  real(8) :: f
	  integer :: l, c, cc, e
	  
	  l = LEN_TRIM(s)
	  e = 0
	  f = 0
	  do c = l, 1, -1
	     cc = ICHAR(s(c:c))
		 cc = cc - mn
		 if ( cc >= 0 .and. cc <= 9 ) then
		    f = f+cc*10**e
			e = e + 1
		 else if ( s(c:c) == '.' ) then
		    f = f / (10**e); e = 0
		 else if ( s(c:c) == '+' ) then
		    continue
		 else if ( s(c:c) == '-' ) then
		    f = -1.0 * f
		 else
		    write(*,*) "Error: invalid character in numeric value: ",TRIM(s)
		 endif
	  enddo	  
	  
	  return
   end function char2float

   
   SUBROUTINE get_date(time_date,time)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN OUT) :: time_date
      REAL(8), INTENT(OUT) :: time
      INTEGER, DIMENSION(8) :: time_data  ! year, month, day, UTC diff, hh, mm, ss, msec
      
      CALL DATE_AND_TIME(VALUES=time_data)
      WRITE(time_date,"(2(I2.2,':'),I2.2,' - ',2(I2.2,'/'),I4.4)") time_data(5), &
   				& time_data(6),time_data(7), &
   				& time_data(3), time_data(2),time_data(1)
      time = time_data(8)/1000.0D0 + time_data(7) + time_data(6)*60.0D0 + time_data(5)*3600.0D0
   END SUBROUTINE get_date
   

   SUBROUTINE IntVal2ppm(s1, s2, idata, vmin, vmax, name) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: s1, s2
      INTEGER, DIMENSION(s1,s2), INTENT(IN) :: idata
      INTEGER, INTENT(IN OUT) :: vmin, vmax
      CHARACTER(*), INTENT(IN) :: name
      !
      ! Simple subroutine to dump double data with fixed min & max values 
      !   in a PPM format
      !
   INTEGER ::  i, j, st
   INTEGER, DIMENSION(3,256) :: cm  ! R,G,B, Colour Map 
   INTEGER, PARAMETER :: ouni=33, ColMap=44
   INTEGER ::   vp, vs
   INTEGER ::   rmin, rmax, value
   CHARACTER(80) ::  fname, jname, command
   !
   !  Define color map: 256 colours
   !
   OPEN(ColMap, FILE="ColorMap.txt", STATUS="OLD",IOSTAT=st)
   if ( st /= 0 ) THEN
   	   WRITE(*,'("Error read opening file ColorMap.txt")')
   	   STOP
   ENDIF
   cm = 0
   DO i = 1, 256
      !READ(ColMap,'(A)', IOSTAT=st) command
      !print*,TRIM(command)
      READ(ColMap,'(3(1X,I3.3))', IOSTAT=st) cm(1,i) ,cm(2,i) ,cm(3,i)
      IF ( st /= 0 ) THEN
    	  	  WRITE(*,'("Error reading colour map at line ",I3,": r, g, b =")',ADVANCE="NO") i
   	  	  WRITE(*,'(3(1X,I3))') cm(:,i)
   	  	  STOP
   	  ENDIF
   ENDDO
   CLOSE(ColMap)

   !
   !   Write on unit ouni with  PPM format
   !
   fname = TRIM(name)//".ppm"
   OPEN(ouni, FILE=TRIM(fname), STATUS="UNKNOWN",IOSTAT=st)
   if ( st /= 0 ) THEN
   	   WRITE(*,'("!!!! Error write access to file ",A)') TRIM(fname)
   	   STOP
   ENDIF

   !   Magic code 
   WRITE(ouni,'("P3")')
   !   Dimensions 
   WRITE(ouni,'(2(I6,1X))') s1, s2
   !   Maximum value 
   WRITE(ouni,'("255")')
   !   Values from 0 to 255 
   rmin = MINVAL(idata); rmax = MAXVAL(idata)
   if ( (vmin == vmax) .AND. (vmin == 0) ) THEN
      vmin = rmin; vmax = rmax
   ELSE
  	  rmin = vmin; rmax = vmax
   ENDIF
   vs = 0
     !! print*,"cm=",cm
   DO j = 1, s2
      DO i = 1, s1
      	 value = idata(i,j)
      	 if ( value < rmin ) value = rmin
      	 if ( value > rmax ) value = rmax
         vp = 1 + (DBLE(value - rmin) * 255.0D0 / DBLE(rmax - rmin))
         vs = vs + 1
         WRITE(ouni,'(3(1X,I4))', ADVANCE="NO") cm(:,vp)
         if ( vs >= 10 ) THEN
            WRITE(ouni,'(1X)')
            vs = 0
         ENDIF
      ENDDO
      WRITE(ouni,'(1X)')
      vs = 0
   ENDDO
   CLOSE(ouni)
   ! the following instructions require ImageMagick tool: comment out if not available
   jname = TRIM(name)//".jpg"
   command = "convert "//TRIM(fname)//" "//TRIM(jname)
   CALL system(command)

   END SUBROUTINE IntVal2ppm

END MODULE Utilities

MODULE CompCode
   USE Utilities
   IMPLICIT NONE
   
   ! Structures and global objects
   
   TYPE i2dGrid
      INTEGER :: EX, EY ! extensions in X and Y directions 
	  REAL(8) :: Xs, Xe, Ys, Ye ! initial and final value for X and Y directions
	  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Values ! 2D matrix of values
   END TYPE i2dGrid
   
   TYPE particle
      REAL(8) :: weight, x, y, vx, vy, fx, fy
   END TYPE particle

   TYPE Population
      INTEGER :: np
      REAL(8), DIMENSION(:), ALLOCATABLE :: weight, x, y, vx, vy; ! particles have a position and few other properties
   END TYPE Population



   TYPE(Population) :: Particles
	
   TYPE(i2dGrid) :: GenFieldGrid, ParticleGrid
   
   
   INTEGER :: MaxIters, MaxSteps
   REAL(8) :: TimeBit  ! Evolution time steps


CONTAINS


   SUBROUTINE newparticle(p, weight, x, y, vx, vy)
      IMPLICIT NONE
      TYPE(particle), INTENT(IN OUT) :: p
      REAL(8), INTENT(IN) :: weight, x, y, vx, vy
      !
	  ! define a new object with passed parameters
	  !
    	p%weight = weight;  p%x = x;  p%y = y;  p%vx = vx;  p%vy = vy;

   END SUBROUTINE newparticle


   SUBROUTINE print_i2dGrid(g) 
      IMPLICIT NONE
      TYPE(i2dGrid), INTENT(IN) :: g

      WRITE(*,'("i2dGrid: EX, EY = ",I6,", ",I6)') g%EX, g%EY
	  WRITE(*,'("         Xs, Xe = ",F12.6,", ",F12.6,"; Ys, Ye = ",F12.6,", ",F12.6)') &
	     & g%Xs,g%Xe,g%Ys,g%Ye
   END SUBROUTINE print_i2dGrid

   SUBROUTINE print_particle(p) 
      IMPLICIT NONE
      TYPE(particle), INTENT(IN) :: p

      WRITE(*,'("particle: weight=",F10.6,", x,y=(",F10.6,",",F10.6,"), vx,vy=(",F10.6,",",F10.6, &
         & "), fx,fy=(",F10.6,",",F10.6,")")')	p%weight, p%x, p%y, p%vx, p%vy, p%fx, p%fy
	    
   END SUBROUTINE print_particle


   SUBROUTINE print_Population(p) 
      IMPLICIT NONE
      TYPE(Population), INTENT(IN) :: p
      
      WRITE(*,'("Population: np = ",I6)') p%np
      
   END SUBROUTINE print_Population
   
   
   
   SUBROUTINE DumpPopulation(p, t)
      IMPLICIT NONE
      TYPE(Population), INTENT(IN) :: p
      INTEGER, INTENT(IN) :: t
      !
      ! save population values on file 
      !
      CHARACTER(80) :: fname
      INTEGER, PARAMETER :: dump=11
      INTEGER :: st
   
      WRITE(fname,'("Population",I4.4,".dmp")') t
      OPEN(dump,FILE=TRIM(fname),FORM="UNFORMATTED",IOSTAT=st)
      IF ( st /= 0 ) THEN
         WRITE(*,'("Error write open file ",A)') TRIM(fname)
         STOP
      ENDIF
      WRITE(dump) p%np, p%weight, p%x, p%y
      CLOSE(dump)
      
   END SUBROUTINE DumpPopulation



   SUBROUTINE ParticleStats(p, t)
      IMPLICIT NONE
      TYPE(Population), INTENT(IN) :: p
      INTEGER, INTENT(IN) :: t
      !
	  ! write a file with statistics on population 
	  !
      INTEGER, PARAMETER :: stats=22
      REAL(8) :: w, xg, yg, wmin, wmax
      INTEGER :: i, st
      
      IF ( t <= 0 ) THEN
         OPEN(stats,FILE="Population.sta",IOSTAT=st) 
      ELSE
         OPEN(stats,FILE="Population.sta",POSITION="APPEND",IOSTAT=st) 
      ENDIF
      IF ( st /= 0 ) THEN
	     WRITE(*,'("Error append/open file Population.sta")')
	     STOP
	  ENDIF

      w = 0.0; xg = 0.0; yg = 0.0
      wmin = p%weight(1)
      wmax = p%weight(1)
      DO i = 1, p%np
	     IF ( wmin > p%weight(i) ) wmin = p%weight(i)
	     IF ( wmax < p%weight(i) ) wmax = p%weight(i)
         w = w + p%weight(i)
	     xg = xg + (p%weight(i) * p%x(i));
	     yg = yg + (p%weight(i) * p%y(i));
      ENDDO
      xg = xg / w;  yg = yg / w;
      WRITE(stats,'("At iteration ",I4," particles: ",I6,"; wmin, wmax = ",F12.2,", ",F12.2,";")') &
             &  t, p%np, wmin, wmax
      WRITE(stats,'("   total weight = ",F12.2,"; CM = (",F10.4,",",F10.4,")")') w,xg,yg
      CLOSE(stats)

   END SUBROUTINE ParticleStats



   SUBROUTINE ForceCompt(f, p1, p2)
      IMPLICIT NONE
      REAL(8), DIMENSION(2), INTENT(INOUT) :: f
      TYPE(particle), INTENT(IN) :: p1, p2     
      !
      ! Compute force acting on p1 by p1-p2 interactions 
	  !
	   REAL(8) :: force, d, d2, dx, dy
       REAL(8), PARAMETER :: k=0.001, tiny=1.0D0/1000000.0D0
	
	   dx = p2%x - p1%x; dy = p2%y - p1%y;
	   d2 = dx*dx + dy*dy  ! what if particles get in touch? Simply avoid the case
	   IF ( d2 < tiny ) d2 = tiny
	   force = (k * p1%weight * p2%weight) / d2;
	   f(1) = force * dx / sqrt(d2); f(2) = force * dy / sqrt(d2);
	
   END SUBROUTINE ForceCompt



   SUBROUTINE ComptPopulation(p, forces)
      IMPLICIT NONE
      TYPE(Population), INTENT(IN OUT) :: p
      REAL(8), DIMENSION(:,:), INTENT(IN) :: forces
      !
	  ! compute effects of forces on particles in a interval time
	  !
      INTEGER :: i
      REAL(8) :: x0, x1, y0, y1
      
      DO i = 1, p%np

		x0 = p%x(i); y0 = p%y(i); 
				
		p%x(i) = p%x(i) + (p%vx(i)*TimeBit) + &
		     (0.5*forces(1,i)*TimeBit*TimeBit/p%weight(i));
		p%vx(i) = p%vx(i) + forces(1,i)*TimeBit/p%weight(i);

		p%y(i) = p%y(i) + (p%vy(i)*TimeBit) + &
		     (0.5*forces(2,i)*TimeBit*TimeBit/p%weight(i));		
		p%vy(i) = p%vy(i) + forces(2,i)*TimeBit/p%weight(i);

	  ENDDO
	  
   END SUBROUTINE ComptPopulation



   SUBROUTINE InitGrid(InputFile)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: InputFile
      !
      ! Output:
      ! GenFieldGrid, ParticleGrid initialization
      ! Maxiters, TimeBit
      !
      INTEGER :: nv, iv, st
      REAL(8) :: dv
      CHARACTER(80) :: filerow
      INTEGER :: inpunit=44
      
      WRITE(*,'("Initializing grids ...")')
      
      OPEN(inpunit,FILE=TRIM(InputFile), IOSTAT=st)
      IF ( st /= 0 ) THEN
         WRITE(*,'("!!!! Error read access to file",A)') TRIM(InputFile)
         STOP
      ENDIF

      ! Now read measured values; they are read in the followin order:
      ! GenFieldGrid.EX, GenFieldGrid.EY, 
      ! GenFieldGrid.Xs, GenFieldGrid.Xe, GenFieldGrid.Ys, GenFieldGrid.Ye
      ! ParticleGrid.Xs, ParticleGrid.Xe, ParticleGrid.Ys, ParticleGrid.Ye
   
      nv = 0; iv = 0; dv = 0.0
      DO
         READ(inpunit,'(A)',IOSTAT=st) filerow
         IF ( st /= 0 ) THEN
      	  	  WRITE(*,'("Error reading input file")')
      	  	  STOP
      	 ENDIF
      	 IF ( filerow(1:1) == '#' ) CYCLE
      	 IF ( nv <= 0 ) THEN
            iv = char2int(filerow)
			GenFieldGrid%EX = iv; nv = 1;
			CYCLE
		 ENDIF		 
         IF ( nv == 1 ) THEN
            iv = char2int(filerow)
			GenFieldGrid%EY = iv; nv = nv+1;
			CYCLE
		 ENDIF
         IF ( nv == 2 ) THEN
            dv = char2float(filerow)
			GenFieldGrid%Xs = dv; nv = nv + 1

			CYCLE
		 ENDIF
         IF ( nv == 3 ) THEN
            dv = char2float(filerow)
			GenFieldGrid%Xe = dv; nv = nv + 1

			CYCLE
		 ENDIF
         IF ( nv == 4 ) THEN
            dv = char2float(filerow)
			GenFieldGrid%Ys = dv; nv = nv + 1

			CYCLE
		 ENDIF
         IF ( nv == 5 ) THEN
            dv = char2float(filerow)
			GenFieldGrid%Ye = dv; nv = nv + 1

			CYCLE
		 ENDIF
         IF ( nv <= 6 ) THEN
            iv = char2int(filerow)
			ParticleGrid%EX = iv; nv = nv + 1

			CYCLE
		 ENDIF
         IF ( nv == 7 ) THEN
            iv = char2int(filerow)
			ParticleGrid%EY = iv; nv = nv + 1

			CYCLE
		 ENDIF
         IF ( nv == 8 ) THEN
            dv = char2float(filerow)
			ParticleGrid%Xs = dv; nv = nv + 1

			CYCLE
		 ENDIF		 
         IF ( nv == 9 ) THEN
            dv = char2float(filerow)
			ParticleGrid%Xe = dv; nv = nv + 1

			CYCLE
		 ENDIF		 
         IF ( nv == 10 ) THEN
            dv = char2float(filerow)
			ParticleGrid%Ys = dv; nv = nv + 1

			CYCLE
		 ENDIF
         IF ( nv == 11 ) THEN
            dv = char2float(filerow)
			ParticleGrid%Ye = dv; nv = nv + 1

			EXIT
		 ENDIF
      ENDDO
	  
      !
      !  Now read MaxIters
      !
      MaxIters = 0;
      DO
         READ(inpunit,'(A)',IOSTAT=st) filerow
         IF ( st /= 0 ) THEN
      	  	  WRITE(*,'("Error reading MaxIters from input file")')
      	  	  STOP
      	 ENDIF
      	 IF ( filerow(1:1) == '#' .OR. LEN_TRIM(filerow) < 1 ) CYCLE
      	 MaxIters = char2int(filerow)
         
      	  WRITE(*,'("MaxIters = ",I6)') MaxIters
         EXIT    	  
      ENDDO    

      !
      !  Now read MaxSteps
      !
      MaxSteps = 0;
      DO
         READ(inpunit,'(A)',IOSTAT=st) filerow
         IF ( st /= 0 ) THEN
      	  	  WRITE(*,'("Error reading MaxIters from input file")')
      	  	  STOP
      	 ENDIF
      	 IF ( filerow(1:1) == '#'  .OR. LEN_TRIM(filerow) < 1 ) CYCLE
      	 MaxSteps = char2int(filerow)
         
      	  WRITE(*,'("MaxSteps = ",I6)') MaxSteps
         EXIT    	  
      ENDDO

      !
      !  Now read TimeBit
      !
      TimeBit = 0.0;
      DO
         READ(inpunit,'(A)',IOSTAT=st) filerow
         IF ( st /= 0 ) THEN
      	  	  WRITE(*,'("Error reading MaxIters from input file")')
      	  	  STOP
      	 ENDIF
      	 IF ( filerow(1:1) == '#'  .OR. LEN_TRIM(filerow) < 1 ) CYCLE
      	 TimeBit = char2float(filerow)
         
      	  WRITE(*,'("TimeBit = ",F10.6)') TimeBit
         EXIT    	  
      ENDDO

      
      CLOSE(inpunit)

      ! Grid allocations
      ALLOCATE(GenFieldGrid%Values(GenFieldGrid%EX,GenFieldGrid%EY),STAT=st)
	  if ( st /= 0 ) THEN
              WRITE(*,'("Error allocating GenFieldGrid%Values(:,:)")')
      	  	  STOP
  	  ENDIF
	  iv = ParticleGrid%EX * ParticleGrid%EY;
	  
      ALLOCATE(ParticleGrid%Values(ParticleGrid%EX,ParticleGrid%EY),STAT=st)
	  if ( st /= 0 ) THEN
              WRITE(*,'("Error allocating ParticleGrid%Values(:,:)")')
      	  	  STOP
  	  ENDIF

      WRITE(*,'("GenFieldGrid ")',ADVANCE="NO")
	  CALL print_i2dGrid(GenFieldGrid)
      WRITE(*,'("ParticleGrid ")',ADVANCE="NO")
	  CALL print_i2dGrid(ParticleGrid)

END SUBROUTINE InitGrid



   SUBROUTINE GeneratingField( grid, MaxIt)
      IMPLICIT NONE
      TYPE(i2dGrid), INTENT(IN OUT) :: grid
      INTEGER, INTENT(IN) :: MaxIt
      !
      !  Compute "generating" points 
      !  Output:
      !    *grid.Values 
      !

   INTEGER :: ix, iy, iz
   REAL(8) :: ca, cb, za, zb    
   REAL(8) ::  rad, zan, zbn
   REAL(8) ::  Xinc, Yinc, Sr, Si, Ir, Ii    
   INTEGER ::  izmn, izmx
   INTEGER ::  Xdots, Ydots
   
   WRITE(*,'("Computing generating field ...")')
   
   Xdots = grid%EX; Ydots = grid%EY
   Sr = grid%Xe - grid%Xs 
   Si = grid%Ye - grid%Ys
   Ir = grid%Xs 
   Ii = grid%Ys
   Xinc = Sr / DBLE(Xdots)
   Yinc = Si / DBLE(Ydots)

   izmn=9999; izmx=-9
   DO iy = 1, Ydots
      DO ix = 1, Xdots
         ca = Xinc * ix + Ir 
         cb = Yinc * iy + Ii 
         rad = SQRT( ca * ca * ( 1.0D0 + (cb/ca)*(cb/ca) ) )
         zan = 0.0D0
         zbn = 0.0D0
         DO iz = 1, MaxIt
            if ( rad > 2.0D0 ) EXIT  
            za = zan
            zb = zbn
            zan = ca + (za-zb)*(za+zb)    
            zbn = 2.0 * ( za*zb + cb/2.0D0 )   
            rad = SQRT( zan * zan * ( 1.0D0 + (zbn/zan)*(zbn/zan) ) )
         ENDDO
         if (izmn > iz) izmn=iz
         if (izmx < iz) izmx=iz
		 if ( iz >= MaxIt ) iz = 0
       	 grid%Values(ix,iy) = iz
      ENDDO
   ENDDO

   END SUBROUTINE GeneratingField


   SUBROUTINE ParticleGeneration( grid, pgrid, pp)
      IMPLICIT NONE
      TYPE(i2dGrid), INTENT(IN) :: grid, pgrid
      TYPE(Population), INTENT(IN OUT) :: pp
 
	  !! A system of particles is generated according to the value distribution of 
	  !! grid.Values
	  
	INTEGER :: vmax, vmin, v
	INTEGER :: Xdots, Ydots
	INTEGER :: ix, iy, np, n, st
	REAL(8) :: p
	
	Xdots = grid%EX; Ydots = grid%EY
	vmax = MAXVAL(grid%Values)
    vmin = MINVAL(grid%Values)
    
    !! Just count number of particles to be generated
    vmin = DBLE(1*vmax + 29*vmin) / 30.0D0
    np = 0
    DO ix = 1, Xdots
       DO iy = 1, Ydots
			v = grid%Values(ix,iy)
			IF ( v <= vmax .AND. v >= vmin ) np = np + 1
		ENDDO
	ENDDO

	!! allocate memory space for particles
	pp%np = np
	ALLOCATE(pp%weight(np), STAT=st)
	ALLOCATE(pp%x(np), STAT=st)
	ALLOCATE(pp%y(np), STAT=st)
	ALLOCATE(pp%vx(np), STAT=st)
	ALLOCATE(pp%vy(np), STAT=st)

	!! Population initialization
	n = 1
	DO ix = 1, grid%EX
	   DO iy = 1, grid%EY
			v = grid%Values(ix,iy)
			IF ( v <= vmax .AND. v >= vmin ) THEN
				pp%weight(n) = v*10.0D0 

				p = (pgrid%Xe-pgrid%Xs) * ix / (grid%EX * 2.0D0)
				pp%x(n) = pgrid%Xs + ((pgrid%Xe-pgrid%Xs)/4.0D0) + p
				
				p = (pgrid%Ye-pgrid%Ys) * iy / (grid%EY * 2.0D0)
				pp%y(n) = pgrid%Ys + ((pgrid%Ye-pgrid%Ys)/4.0D0) + p
				
				pp%vx(n) = 0.0D0 ! at start particles are still 
				pp%vy(n) = 0.0D0 

				n = n + 1
				IF ( n > np ) EXIT
			ENDIF
		ENDDO
		if ( n > np ) EXIT
	ENDDO
	
      CALL print_Population(pp)
   END SUBROUTINE ParticleGeneration


   SUBROUTINE SystemEvolution( pgrid, pp, mxiter)
      IMPLICIT NONE
      TYPE(i2dGrid), INTENT(IN OUT) :: pgrid
      TYPE(Population), INTENT(IN OUT) :: pp
      INTEGER, INTENT(IN) :: mxiter
 
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: forces
    REAL(8) :: vmin, vmax
	TYPE(particle) :: p1, p2
	REAL(8), DIMENSION(2) :: f
	INTEGER :: i, j, t, st

     !! temporary array of forces 
     ALLOCATE(forces(2,pp%np), STAT=st)
     IF ( st /= 0 ) THEN
		 WRITE(*,'("Error mem alloc of forces!")');
		 STOP
	 ENDIF

	 !! compute forces acting on each particle step by step
	 DO t = 0, mxiter-1
		 WRITE(*,'("Step ",I6," of ",I6)') t,mxiter
		 CALL ParticleScreen(pgrid,pp, t)
		 ! ParticleDump call frequency may be changed
		 if ( mod(t,4) == 0 ) CALL DumpPopulation(pp, t)
		 CALL ParticleStats(pp,t)
         forces = 0.0D0
         DO i = 1, pp%np
            CALL newparticle(p1,pp%weight(i),pp%x(i),pp%y(i),pp%vx(i),pp%vy(i))
            DO j = 1, pp%np
				IF ( j /= i ) THEN
                 CALL newparticle(p2,pp%weight(j),pp%x(j),pp%y(j),pp%vx(j),pp%vy(j))
				 CALL ForceCompt(f,p1,p2)
				 forces(1,i) = forces(1,i) + f(1)
				 forces(2,i) = forces(2,i) + f(2)
				ENDIF
			ENDDO
		 ENDDO
		 CALL ComptPopulation(pp,forces)
	 ENDDO
	 DEALLOCATE(forces)    
    
   END SUBROUTINE SystemEvolution



   SUBROUTINE ParticleScreen( pgrid,  pp, step)
      IMPLICIT NONE
      TYPE(i2dGrid), INTENT(IN OUT) :: pgrid
      TYPE(Population), INTENT(IN) :: pp
      INTEGER, INTENT(IN) :: step

	!! Distribute a particle population in a grid for visualization purposes

	INTEGER :: ix, iy, Xdots, Ydots
	INTEGER :: np, n, wp
	REAL(8) :: rmin, rmax
	INTEGER, SAVE :: vmin, vmax
	REAL(8) :: Dx, Dy, wint, wv
	CHARACTER(40) :: name
	
	Xdots = pgrid%EX; Ydots = pgrid%EY
	pgrid%Values = 0

    rmin = MINVAL(pp%weight)
    rmax = MAXVAL(pp%weight)
    wint = rmax - rmin
	Dx = pgrid%Xe - pgrid%Xs
	Dy = pgrid%Ye - pgrid%Ys
	DO n = 1, pp%np
		!! keep a tiny border free anyway
		ix = Xdots * pp%x(n) / Dx; IF ( ix >= Xdots .OR. ix <= 1 ) CYCLE
		iy = Ydots * pp%y(n) / Dy; IF ( iy >= Ydots .OR. iy <= 1 ) CYCLE
		wv = pp%weight(n) - rmin; wp = 10.0*wv/wint
		pgrid%Values(ix,iy) = wp
 	    pgrid%Values(ix-1,iy) = wp
        pgrid%Values(ix+1,iy) = wp
        pgrid%Values(ix,iy-1) = wp
        pgrid%Values(ix,iy+1) = wp
	    
	ENDDO
	WRITE(name,'("stage",I3.3)') step
    IF ( step <= 0 ) THEN
       vmin = 0; vmax = 0
    ENDIF
    CALL IntVal2ppm(pgrid%EX, pgrid%EY, pgrid%Values, vmin, vmax, name)
    
   END SUBROUTINE ParticleScreen

END MODULE CompCode

PROGRAM P_PARTICLES
   USE CompCode
   IMPLICIT NONE

   CHARACTER(80) :: asctime
   REAL(8) :: time0, time1
   
   CALL get_date(asctime,time0)
   WRITE(*,'("Starting at: ",a)') TRIM(asctime)
   
   CALL InitGrid("Particles.inp")

   ! GenFieldGrid initialization
   WRITE(*,'("GeneratingField...")')
   CALL GeneratingField(GenFieldGrid,MaxIters)

   ! Particle population initialization
   WRITE(*,'("ParticleGeneration...")')
   CALL ParticleGeneration(GenFieldGrid, ParticleGrid, Particles)
   
   ! Compute evolution of the particle population
   WRITE(*,'("SystemEvolution...")')
   CALL SystemEvolution(ParticleGrid, Particles, MaxSteps)
   
   CALL get_date(asctime,time1)
   WRITE(*,'("Ending   at: ",a)') TRIM(asctime)
   WRITE(*,'("Computations ended in ",F12.2," seconds")') time1-time0

   WRITE(*,'("End of program!")')

END PROGRAM P_PARTICLES
