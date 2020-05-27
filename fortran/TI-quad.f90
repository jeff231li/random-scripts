PROGRAM TI_quadrature
    IMPLICIT NONE

    ! Initialize Variables
    INTEGER, PARAMETER                   :: mx_line=2000000
    CHARACTER(LEN=80)                    :: filename_inp, arg
    CHARACTER(LEN=10)                    :: col1
    INTEGER                              :: ndx, step, length, data_nr, data, sample_point
    INTEGER                              :: i, j, k, narg, Nlen
    INTEGER                              :: skip, av_pnt, Ind
    REAL                                 :: e(8), l1, l2, dl
    REAL                                 :: av_fe, av_fee, av_few, av_dg
    REAL                                 :: w5(5), w3(3), w7(7), w9(9), w12(12), W15(15)
    REAL                                 :: a5(5), a3(3), a7(7), a9(9), a12(12), a15(15)

    ! Initialize Allocatable Array Variables
    REAL,    ALLOCATABLE, DIMENSION(:,:) :: dv, adv, dva, dve, dvw, fe, sigma_fe
    REAL,    ALLOCATABLE, DIMENSION(:,:) :: my_dvdl, my_dvdel, my_dvdwl, my_dgdl, dvea, dvwa
    REAL,    ALLOCATABLE, DIMENSION(:)   :: sigma_dvdl, sigma_dvdel, sigma_dvdwl, sigma_dgdl
    REAL,    ALLOCATABLE, DIMENSION(:)   :: w, lambda
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: time_step

    ! Gaussian Quadrature Coefficients
    DATA w3  /0.27777, 0.44444, 0.27777/
    DATA w5  /0.11846, 0.23931, 0.28444, 0.23931, 0.11846/
    DATA w7  /0.06474, 0.13985, 0.19091, 0.20897, 0.19091, 0.13985, 0.06474/
    DATA w9  /0.04064, 0.09032, 0.13031, 0.15617, 0.16512, 0.15617, 0.13031, 0.09032, 0.04064/
    DATA w12 /0.02359, 0.05347, 0.08004, 0.10158, 0.11675, 0.12457, 0.12457, 0.11675, 0.10158, &
            & 0.08004, 0.05347, 0.02359/
    DATA w15 /0.01538, 0.03518, 0.05358, 0.06979, 0.08313, 0.09308, 0.09922, 0.10129, 0.09922, &
            & 0.09308, 0.08313, 0.06979, 0.05358, 0.03518, 0.01538/

    DATA a3  /0.11270, 0.50000, 0.88730/
    DATA a5  /0.04691, 0.23077, 0.50000, 0.76923, 0.95309/
    DATA a7  /0.02545, 0.12923, 0.29708, 0.50000, 0.70292, 0.87077, 0.97455/
    DATA a9  /0.01592, 0.08198, 0.19331, 0.33787, 0.50000, 0.66213, 0.80669, 0.91802, 0.98408/
    DATA a12 /0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738, 0.56262, 0.68392, 0.79366, &
            & 0.88495, 0.95206, 0.99078/
    DATA a15 /0.00600, 0.03136, 0.07590, 0.13779, 0.21451, 0.30292, 0.39940, 0.50000, 0.60060, &
            & 0.69708, 0.78549, 0.86221, 0.92410, 0.96864, 0.99400/

    ! Help information
    narg = COMMAND_ARGUMENT_COUNT()
    IF (narg == 0) THEN
        WRITE(*,*) '==============='
        WRITE(*,*) ' TI QUADRATURE    -    Author: Jeffry Setiadi, USyd (2016)'
        WRITE(*,*) '==============='
        WRITE(*,*) ''
        WRITE(*,*) ' * Fortran program to calculate free energy profile from Thermodynamic'
        WRITE(*,*) '   Integration (TI) simulation using Gaussian quadrature from NAMD output.'
        WRITE(*,*) ''
        WRITE(*,*) ' * Requires 3 command line argument:'
        WRITE(*,*) '       -> input file with list of *.fepout filenames'
        WRITE(*,*) '       -> Integer value for block averaging'
        WRITE(*,*) '       -> Integer value number of points to skip'
        WRITE(*,*) ''
        WRITE(*,*) ' * Currently only Gaussian Quadrature with 3, 5, 7, 9, 12 and 15 nodes'
        WRITE(*,*) '   are supported.'
        WRITE(*,*) ''
        WRITE(*,*) ' * Example Usage:'
        WRITE(*,*) '        TI_quad files.inp 100 0'
        WRITE(*,*) ''
        WRITE(*,*) '        files.inp:'
        WRITE(*,*) '             out_t1l1_f.fepout'
        WRITE(*,*) '             out_t1l2_f.fepout'
        WRITE(*,*) '             out_t1l3_f.fepout'
        WRITE(*,*) '             ...'
        WRITE(*,*) ''
        WRITE(*,*) '        Comments are ignored in input file (i.e. lines starting with #)'
        WRITE(*,*) ''
        WRITE(*,*) ' * Output files:'
        WRITE(*,*) '      -> "Ensemble-Average.dat"'
        WRITE(*,*) '         ->> Integrated free energy profile running average'
        WRITE(*,*) '      -> "Ensemble-Block.dat"'
        WRITE(*,*) '         ->> Integrated free energy profile block average'
        WRITE(*,*) '      -> "Block-Average.dat"'
        WRITE(*,*) '         ->> Block data analysis for each lambda values'
        WRITE(*,*) '      -> "Statistics-Block.dat"'
        WRITE(*,*) '         ->> Average values of block data series'
        WRITE(*,*) '      -> "Window-Derivative.dat"'
        WRITE(*,*) '         ->> Integrated free energy profile running & block average'
        WRITE(*,*) '      -> "Free-Energy.dat"'
        WRITE(*,*) '         ->> Free energy values with error'
        WRITE(*,*) ''
        STOP
    END IF

    ! Command line argument
    CALL getarg(1, arg)
    READ(arg,*) filename_inp
    CALL getarg(2, arg)
    READ(arg,*) av_pnt
    CALL getarg(3, arg)
    READ(arg,*) skip

    WRITE(*,*) '==============='
    WRITE(*,*) ' TI QUADRATURE    -    Author: Jeffry Setiadi, USyd (2016)'
    WRITE(*,*) '==============='
    WRITE(*,*) ''

    ! Open input file and determine number of files
    OPEN(UNIT=7, FILE=filename_inp, STATUS='old')
    data_nr = 0
    Nlen = 0
    DO DATA=1,mx_line
        READ(7,*,END=1) arg
        Nlen = Nlen + 1
        IF (INDEX(arg, "#") /= 0) CYCLE
        data_nr = data_nr + 1
    END DO
    1 WRITE(*,'(A,I3)'), ' * Number of Quadrature nodes:', data_nr
    WRITE(*,*) ''
    REWIND(7)

    ! Determine Quadrature Length
    WRITE(*,'(2A10,A11)'), 'Window', 'Lambda', 'Weights'
    WRITE(*,'(2A10,A11)'), '======', '======', '======='
    SELECT CASE (data_nr)
        CASE(3)
            ALLOCATE(w(3))
            DO i=1,3
                w(i) = w3(i)
                WRITE(*,'(I10,F10.4,F11.4)'), i, a3(i), w(i)
            END DO
        CASE(5)
            ALLOCATE(w(5))
            DO i=1,5
                w(i) = w5(i)
                WRITE(*,'(I10,F10.4,F11.4)'), i, a5(i), w(i)
            END DO
        CASE(7)
            ALLOCATE(w(7))
            DO i=1,7
                w(i) = w7(i)
                WRITE(*,'(I10,F10.4,F11.4)'), i, a7(i), w(i)
            END DO
        CASE(9)
            ALLOCATE(w(9))
            DO i=1,9
                w(i) = w9(i)
                WRITE(*,'(I10,F10.4,F11.4)'), i, a9(i), w(i)
            END DO
        CASE(12)
            ALLOCATE(w(12))
            DO i=1,12
                w(i) = w12(i)
                WRITE(*,'(I10,F10.4,F11.4)'), i, a12(i), w(i)
            END DO
        CASE(15)
            ALLOCATE(w(15))
            DO i=1,15
                w(i) = w15(i)
                WRITE(*,'(I10,F10.4,F11.4)'), i, a15(i), w(i)
            END DO
        CASE DEFAULT
            WRITE(*,*) &
            & 'Number of files (nodes) not currently supported [only 3, 5, 7, 9, 12 and 15].'
            STOP
    END SELECT
    WRITE(*,*) ''

    ! Read individual files and determine number of sample points
    WRITE(*,*) '* Reading Input file'
    sample_point = 0
    DO data=1,data_nr
       READ(7,*) arg
       IF (INDEX(arg, "#") /= 0) CYCLE
       ndx = 0
       OPEN(UNIT=8, FILE=TRIM(arg), STATUS='old')
       WRITE(*,'(3A)'), '   -> fep file opened "', TRIM(arg), '"'
       DO j=1, mx_line
          READ(8,*,END=2) col1
          IF (col1 == 'FepEnergy:') THEN
             ndx = ndx + 1
          END IF
       END DO
       2 sample_point = MAX(sample_point, ndx)
       WRITE(*,'(A,2I7)'), '     * Number of sampling sample_points ', sample_point, ndx
    CLOSE(8)
    END DO
    REWIND(7)
    WRITE(*,*) ''

    ! Query the user for number of points for Averaging and Equilibration
    length = (sample_point - skip) / av_pnt
    WRITE(*,*) '* Statistics'
    WRITE(*,'(A15,I7)'), '   -> Average  :', length
    WRITE(*,'(A15,I7)'), '   -> Skip     :', skip
    WRITE(*,'(A15,I7)'), '   -> Sampling :', sample_point
    WRITE(*,*)

    ! Allocate Arrays
    ALLOCATE(my_dvdl(data_nr,av_pnt), my_dvdel(data_nr,av_pnt), my_dvdwl(data_nr,av_pnt))
    ALLOCATE(my_dgdl(data_nr,av_pnt))
    ALLOCATE(sigma_dvdl(data_nr), sigma_dvdel(data_nr), sigma_dvdwl(data_nr))
    ALLOCATE(sigma_fe(5, av_pnt), sigma_dgdl(data_nr))
    ALLOCATE(dv(data_nr,sample_point), adv(data_nr,sample_point))
    ALLOCATE(dva(data_nr,sample_point), dvea(data_nr,sample_point), dvwa(data_nr,sample_point))
    ALLOCATE(fe(5,sample_point), lambda(data_nr))
    ALLOCATE(time_step(sample_point), dve(data_nr,sample_point), dvw(data_nr,sample_point))

    ! Open Files
    OPEN(UNIT=11, FILE="Ensemble-Average.dat",   STATUS='unknown')
    OPEN(UNIT=12, FILE="Ensemble-Block.dat",     STATUS='unknown')
    OPEN(UNIT=13, FILE="Block-Average.dat",      STATUS='unknown')
    OPEN(UNIT=14, FILE="Statistics-Block.dat",   STATUS='unknown')
    OPEN(UNIT=15, FILE="Windows-Derivative.dat", STATUS='unknown')
    OPEN(UNIT=16, FILE="Free-Energy.dat",        STATUS='unknown')

    ! Read FEP data from fepout files and Average
    WRITE(*,*) '* Quadrature Integration'
    WRITE(14,'(A,A9,10A15)') '#', 'Lambda', 'dG/dl', 'dE/dl', 'Ave dE/dl', 'sig dE/dl', &
        & 'dElec/dl', 'Ave dElec/dl', 'sig dElec/dl', 'dvdW/dl', 'Ave dvdW/dl', 'sig dvdW/dl'
    Ind=1
    DO data=1, Nlen
        READ(7,*) arg
        IF (INDEX(arg, "#") /= 0) CYCLE
        OPEN(UNIT=8, FILE=TRIM(arg), STATUS='old')
        WRITE(*,'(3A)'), '   -> fep file opened "', TRIM(arg), '"'
        ndx = 0

        ! READ through to the end of the file
        DO j=1, mx_line
            READ(8,*,END=3) col1

            ! Obtain the L1, L2 and determine dL
            IF (col1 == '#NEW') THEN
                BACKSPACE(8)
                READ(8,*) (col1, k=1,6), l1, col1, l2
                dl = l2 - l1
                lambda(Ind) = l1
                WRITE(*,'(A,3F8.4)'),    '     * L1, L2 and DL ', l1, l2, dl
                WRITE(15,'(A,I3)')       '# ... Window ', Ind
                WRITE(15,'(A1,A4,7A15)') '#', 'i', 'dE/dl', 'Ave dE/dl', 'dElec/dl', &
                    & 'Ave dElec/dl','dvdW/dl', 'Ave dvdW/dl', 'dG/dl'
            END IF

            ! Obtain Alchemical Data
            IF (col1 == 'FepEnergy:') THEN
                ! Read line
                BACKSPACE(8)
                READ(8,*) col1, step, (e(k), k=1,8)
                ndx = ndx + 1

                ! Store data
                IF (ndx <= skip) CYCLE
                dve (Ind,ndx)  = (e(2) - e(1)) / dl                    ! dElec/dl
                dvw (Ind,ndx)  = (e(4) - e(3)) / dl                    ! dvdW/dl
                dv  (Ind,ndx)  = e(5) / dl                             ! dE/dl
                adv (Ind,ndx)  = e(8) / dl                             ! dG/dl
                dvea(Ind,ndx)  = SUM(dve(Ind,skip+1:ndx)) / (ndx-skip) ! Mean(dElec/dl)
                dvwa(Ind,ndx)  = SUM(dvw(Ind,skip+1:ndx)) / (ndx-skip) ! Mean(dvdW/dl)
                dva (Ind,ndx)  = SUM(dv (Ind,skip+1:ndx)) / (ndx-skip) ! Mean(dE/dl)
                time_step(ndx) = step                                  ! dt(i)

                ! Print -> i | dE/dl    | Ave dE/dl    |
                !            | dElec/dl | Ave dElec/dl |
                !            | dvdW/dl  | Ave dvdw/dl  |
                !            | dG/dl    |
                WRITE(15,'(I5,7F15.5)') ndx, dv(Ind,ndx), dve(Ind,ndx), dvw(Ind,ndx), &
                    & adv(Ind,ndx), dva(Ind,ndx), dvea(Ind,ndx), dvwa(Ind,ndx)
            END IF
        END DO
        3 CONTINUE
        WRITE(15,*)

        ! Sum over average points (Block data)
        j=0
        WRITE(13,'(A,I1)')       '# Window - ', Ind
        WRITE(13,'(A1,A5,4A15)') '#', 'i', 'Blk dE/dl', 'Blk dElec/dk', 'Blk dvdW/dl', 'Blk dG/dl'
        DO i=skip,sample_point,length
            j = j+1
            IF (j > av_pnt) EXIT
            my_dvdl (Ind,j) = SUM(dv (Ind,i+1:i+length)) / length
            my_dvdel(Ind,j) = SUM(dve(Ind,i+1:i+length)) / length
            my_dvdwl(Ind,j) = SUM(dvw(Ind,i+1:i+length)) / length
            my_dgdl (Ind,j) = SUM(adv(Ind,i+1:i+length)) / length

            ! Print -> i | Average dE/dl | dG/dl
            WRITE(13,'(I6,4F15.4)') skip+(j-1)*length+1, my_dvdl(Ind,j), &
                & my_dvdel(Ind,j), my_dvdwl(Ind,j), adv(Ind,j)
        END DO
        WRITE(13,*)

        ! Calculate standard deviation
        sigma_dgdl (Ind) = SQRT(SUM((my_dgdl (Ind,:) - adv (Ind,ndx))**2) / av_pnt)
        sigma_dvdl (Ind) = SQRT(SUM((my_dvdl (Ind,:) - dva (Ind,ndx))**2) / av_pnt)
        sigma_dvdel(Ind) = SQRT(SUM((my_dvdel(Ind,:) - dvea(Ind,ndx))**2) / av_pnt)
        sigma_dvdwl(Ind) = SQRT(SUM((my_dvdwl(Ind,:) - dvwa(Ind,ndx))**2) / av_pnt)

        ! Print to screen dE/dl for each window
        WRITE(*,'(A,I5,10F12.4)'), '     * sigma_dvdl(i) ', Ind, sigma_dvdl(Ind), &
            & SQRT(SUM((my_dvdl(Ind,:) - dva(Ind,ndx))**2) / av_pnt)

        ! Print -> i | Average dE/dl    | Average Block dEdl     | sigma dE/dl    |
        !            | Average dElec/dl | Average Block dElec/dl | sigma dElec/dl |
        !            | Average dvdW/dl  | Average Block dvdW/dl  | sigma dvdW/dl  |
        WRITE(14,'(F10.5,10F15.5)') lambda(Ind), adv(Ind,ndx), &
            & dva (Ind,ndx), SUM(my_dvdl (Ind,:)) / av_pnt, sigma_dvdl (Ind), &
            & dvea(Ind,ndx), SUM(my_dvdel(Ind,:)) / av_pnt, sigma_dvdel(Ind), &
            & dvwa(Ind,ndx), SUM(my_dvdwl(Ind,:)) / av_pnt, sigma_dvdwl(Ind)

        Ind = Ind + 1
        CLOSE(8)
    END DO
    WRITE(14,*)

    ! Perform Quadrature Calculation
    WRITE(11,'(A1,A5,4A12)') '#', 'i','dE/dl', 'dElec/dl', 'dvdW/dl', 'dG/dl'
    DO i=skip+1,sample_point
        fe(2,i) = SUM(w(:) * adv (:,i)) ! dG/dl
        fe(3,i) = SUM(w(:) * dva (:,i)) ! dE/dl
        fe(4,i) = SUM(w(:) * dvea(:,i)) ! dElec/dl
        fe(5,i) = SUM(w(:) * dvwa(:,i)) ! dvdW/dl

        ! Print -> i | dE/dl | dElec/dl | dvdW/dl | dG/dl |
        WRITE(11,'(I6,4F12.4)') i, fe(3,i), fe(4,i), fe(5,i), fe(2,i)
    END DO
    WRITE(11,*)

    ! Perform Quadrature on Block-Data
    WRITE(12,'(A1,A5,4A12)') '#', 'i', 'dE/dl', 'dElec/dl', 'dvdW/dl', 'dG/dl'
    DO i=1,av_pnt
        ! Block Data Quadrature
        fe(2,i) = SUM(w(:) * my_dgdl (:,i)) ! dG/dl
        fe(3,i) = SUM(w(:) * my_dvdl (:,i)) ! dE/dl
        fe(4,i) = SUM(w(:) * my_dvdel(:,i)) ! dElec/dl
        fe(5,i) = SUM(w(:) * my_dvdwl(:,i)) ! dvdW/dl

        ! Standard deviation
        sigma_fe(2,i) = SQRT(SUM((w(:) * sigma_dgdl (:))**2)) ! dG/dl
        sigma_fe(3,i) = SQRT(SUM((w(:) * sigma_dvdl (:))**2)) ! dE/dl
        sigma_fe(4,i) = SQRT(SUM((w(:) * sigma_dvdel(:))**2)) ! dElec/dl
        sigma_fe(5,i) = SQRT(SUM((w(:) * sigma_dvdwl(:))**2)) ! dvdW/dl

        ! Print -> i | Average Block dE/dl    |
        !            | Average Block dElec/dl |
        !            | Average Block dvdW/dl  |
        WRITE(12,'(I6,4F12.4)') i*length+skip, fe(3,i), fe(4,i), fe(5,i), fe(2,i)
    END DO
    WRITE(12,*)

    ! Average of Block-Data Analysis
    av_dg  = SUM(fe(2,1:av_pnt)) / av_pnt  ! dG/dl
    av_fe  = SUM(fe(3,1:av_pnt)) / av_pnt  ! dE/dl
    av_fee = SUM(fe(4,1:av_pnt)) / av_pnt  ! dElec/dl
    av_few = SUM(fe(5,1:av_pnt)) / av_pnt  ! dvdW/dl

    ! Print -> | Delta E    | sigma Delta E    |
    !          | Delta Elec | sigma Delta Elec |
    !          | Delta vdW  | sigma Delta vdW  |
    !          | Delta G    | sigma Delta G    |
    WRITE(*,'(/,A)')            ' * Results'
    WRITE(*,'(A,F8.3,A3,F6.3)') '   -> Delta E    (kcal/mol):', av_fe,  '+-', sigma_fe(3,1)
    WRITE(*,'(A,F8.3,A3,F6.3)') '   -> Delta Elec (kcal/mol):', av_fee, '+-', sigma_fe(4,1)
    WRITE(*,'(A,F8.3,A3,F6.3)') '   -> Delta vdW  (kcal/mol):', av_few, '+-', sigma_fe(5,1)
    WRITE(*,'(A,F8.3,A3,F6.3)') '   -> Delta G    (kcal/mol):', av_dg,  '+-', sigma_fe(2,1)
    WRITE(*,*)

    ! Print -> length | dE/dl    | sigma dE/dl    |
    !                 | dElec/dl | sigma dElec/dl |
    !                 | dvdW/dl  | sigma dvdw/dl  |
    WRITE(16,'(A1,A5,2A12)') '#', 'N', 'Average', 'Sigma'
    WRITE(16,'(I6,3F12.4)') (i-1)*length+skip, av_fe,  sigma_fe(3,1)
    WRITE(16,'(I6,3F12.4)') (i-1)*length+skip, av_fee, sigma_fe(4,1)
    WRITE(16,'(I6,3F12.4)') (i-1)*length+skip, av_few, sigma_fe(5,1)
    WRITE(16,'(I6,3F12.4)') (i-1)*length+skip, av_dg,  sigma_fe(2,1)
    WRITE(16,*)

    ! Close Files
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
    CLOSE(14)
    CLOSE(15)
    CLOSE(16)

END PROGRAM TI_quadrature
