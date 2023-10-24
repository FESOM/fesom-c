!==============================================================================|
!   Input parameters and flags, which control the simulation                               |
!==============================================================================|
!VF, updated 02.09
SUBROUTINE READ_DATA_RUN()
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE o_UTILIT

use g_parsup

   IMPLICIT NONE
   real(kind=WP) :: REALVEC(150)
   INTEGER ::  INTVEC,ISCAN,KTEMP, NOUT_VARS
   real(kind=WP) :: KSwTMP,ZKLTMP, DAYS
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=3) :: STRINGVEC(16)
   CHARACTER(LEN=80) :: type_vd
   CHARACTER(LEN=80) :: restartc
   INTEGER :: i, j

!==============================================================================|
!   READ IN VARIABLES AND SET VALUES                                           |
!==============================================================================|

   FNAME = "./fv_run.dat"


!==============================================================================|
!   CASENAME                                         |
!==============================================================================|

   ISCAN = SCAN_FILE(TRIM(FNAME),"TITLE",CVAL = TITLE)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING TITLE: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!   Barotropic time step (dt_2d)
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"dt_2D",FSCAL = dt_2D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING dt_2D: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!    Baroclinic time step (dt)
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"dt",FSCAL = dt)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING dt: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!   Julian day (beginning of simulation period)
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"JD0",FSCAL = Time_jd0)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING JD0: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!  Type of coordinates
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"cartesian",LVAL = cartesian)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING cartesian: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "cartesian" TO INPUT FILE'
     END IF
     STOP
   END IF
!------------------------------------------------------------------------------|
!   Total number of barotropic steps
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"nsteps",ISCAL = nsteps)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING nsteps: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "INDIR"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"INDIR",CVAL = INDIR)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING INDIR: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "MESHPATH"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"MESHPATH",CVAL = MESHPATH)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING MESHPATH: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "OUTDIR"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"OUTDIR",CVAL = OUTDIR)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING OUTDIR: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "IREP"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"IREP",ISCAL = IREP)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING IREP: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "IREC"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"IREC",ISCAL = IREC)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING IREC: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "IRESTART"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"IRESTART",ISCAL = IRESTART)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING IRESTART: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "RESTART"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"RESTART",CVAL = RESTARTC)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING RESTART: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "VERBOSITY"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"VERBOSITY",ISCAL = iverbosity)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'NO VERBOSITY VALUE GIVEN, SETTING VERBOSITY=1 '
     iverbosity = 1
   ELSEIF (iverbosity < 0 ) THEN
     WRITE(*,*)'VERBOSITY TO SMALL, ADJUSTING TO VERBOSITY=0 '
     iverbosity = 0
   ENDIF
!------------------------------------------------------------------------------|
!     "mom_adv_2D"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"mom_adv_2D",ISCAL = mom_adv_2D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING mom_adv_2D: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "mom_adv_3D"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"mom_adv_3D",ISCAL = mom_adv_3D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING mom_adv_3D: ',ISCAN
     STOP
   END IF


!------------------------------------------------------------------------------|
!     "tracer_adv"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"tracer_adv",ISCAL = tracer_adv)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING tracer_adv: ',ISCAN
     STOP
   END IF
if (tracer_adv==3) then
!------------------------------------------------------------------------------|
!     "mu_w" - weight for muscl scheme
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"mu_w",FSCAL = mu_w)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING mu_w: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "up_w" - weight for upwind scheme
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"up_w",FSCAL = up_w)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING up_w: ',ISCAN
     STOP
   END IF
if (.not.((up_w+mu_w)==1.0_WP)) WRITE(*,*)'up_w+mu_w should be equal 1. ', up_w+mu_w
END IF
!------------------------------------------------------------------------------|
!     "vert_adv"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"vert_adv",ISCAL = vert_adv)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING vert_adv: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "ex_vert_visc"
!------------------------------------------------------------------------------|
    ISCAN = SCAN_FILE(TRIM(FNAME),"ex_vert_visc",LVAL = ex_vert_visc)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING ex_vert_visc: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "ex_vert_visc" TO INPUT FILE'
     END IF
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "im_vert_visc"
!------------------------------------------------------------------------------|
    ISCAN = SCAN_FILE(TRIM(FNAME),"im_vert_visc",LVAL = im_vert_visc)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING im_vert_visc: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "im_vert_visc" TO INPUT FILE'
     END IF
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "im_vert_visc_fic"
!------------------------------------------------------------------------------|
    ISCAN = SCAN_FILE(TRIM(FNAME),"im_vert_visc_fic",LVAL = im_vert_visc_fic)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING im_vert_visc_fic: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "im_vert_visc_fic" TO INPUT FILE'
     END IF
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "Cd"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"Cd",FSCAL = C_d)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING Cd: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "BFC"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BFC",ISCAL = BFC)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING BFC: ',ISCAN
     STOP
   END IF
   if ((BFC/=1).and.(BFC/=2).and.(BFC/=3)) then
     WRITE(*,*)'ERROR READING BFC, it should be 1 or 2 or 3: ',BFC
     STOP
   END IF


!------------------------------------------------------------------------------|
!     "fic_point"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"fic_point",ISCAL = fic_point)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING fic_point: ',ISCAN
     STOP
   END IF
   if ((fic_point<0).or.(fic_point>2)) then
     WRITE(*,*)'ERROR READING fic_point, it should be 0,1 or 2: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "z0b_min"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"z0b_min",FSCAL = z0b_min)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING z0b_min: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "z0s_min"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"z0s_min",FSCAL = z0s_min)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING z0s_min: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "za"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"za",FSCAL = za)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING za: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "tau_inv"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"tau_inv_filt",FSCAL = tau_inv_filt)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING tau_inv: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "Abh0"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"Abh0",FSCAL = Abh0)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING Abh0: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "C_Smag"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"C_Smag",FSCAL = C_Smag)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING C_Smag: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "tau_c"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"tau_c",FSCAL = tau_c)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING tau_c: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "snul"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"snul",FSCAL = snul)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING snul: ',ISCAN
     STOP
   END IF
     ISCAN = SCAN_FILE(TRIM(FNAME),"filt_2D",LVAL = filt_2D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING filt_2D: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "filt_2D" TO INPUT FILE'
     END IF
     STOP
   END IF
        ISCAN = SCAN_FILE(TRIM(FNAME),"filt_bi_2D",LVAL = filt_bi_2D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING filt_bi_2D: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "filt_bi_2D" TO INPUT FILE'
     END IF
     STOP
   END IF
           ISCAN = SCAN_FILE(TRIM(FNAME),"bih_2D",LVAL = bih_2D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING bih_2D: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "bih_2D" TO INPUT FILE'
     END IF
     STOP
   END IF
   ISCAN = SCAN_FILE(TRIM(FNAME),"filt_3D",LVAL = filt_3D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING filt_3D: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "filt_3D" TO INPUT FILE'
     END IF
     STOP
   END IF
        ISCAN = SCAN_FILE(TRIM(FNAME),"filt_bi_3D",LVAL = filt_bi_3D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING filt_bi_3D: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "filt_bi_3D" TO INPUT FILE'
     END IF
     STOP
   END IF
           ISCAN = SCAN_FILE(TRIM(FNAME),"bih_3D",LVAL = bih_3D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING bih_3D: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "bih_3D" TO INPUT FILE'
     END IF
     STOP
   END IF


!------------------------------------------------------------------------------|
!  BAROTROPIC FLAG
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"BAROTROPIC",LVAL = BAROTROPIC)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING BAROTROPIC: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "BAROTROPIC" TO INPUT FILE'
     END IF
     STOP
   END IF

!------------------------------------------------------------------------------|
!  BAROTROPIC_3D FLAG
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"BAROTROPIC_3D",LVAL = BAROTROPIC_3D)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING BAROTROPIC_3D: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "BAROTROPIC_3D" TO INPUT FILE'
     END IF
     STOP
   END IF

!------------------------------------------------------------------------------|
!  TEMP_ON
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"TEMP_ON",LVAL = TEMP_ON)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING TEMP_ON: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "TEMP_ON" TO INPUT FILE'
     END IF
     STOP
   END IF

!------------------------------------------------------------------------------|
!  SALINITY_ON
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"SALINITY_ON",LVAL = SALINITY_ON)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING SALINITY_ON: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "SALINITY_ON" TO INPUT FILE'
     END IF
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "S_const"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"S_const",FSCAL = S_const)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING S_const: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "T_const"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"T_const",FSCAL = T_const)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING T_const: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!  Coriolis_TF
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"Coriolis_TF",LVAL = Coriolis_TF)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING Coriolis_TF: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(*,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
       WRITE(*,*)'PLEASE ADD LOGICAL (T/F) VARIABLE "Coriolis_TF" TO INPUT FILE'
     END IF
     STOP
   END IF


if (cartesian) then
!------------------------------------------------------------------------------|
!     "lat_cartesian"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"lat_cartesian",FSCAL = lat_cartesian)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING lat_cartesian: ',ISCAN
     STOP
   END IF

lat_cartesian=lat_cartesian*pi/180.0_WP
endif
!------------------------------------------------------------------------------|
!     "type_vd"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"type_vd",CVAL = type_vd)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING type_vert_dis: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "nsigma"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"nsigma",ISCAL = nsigma)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING nsigma: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "ind_sigma"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"ind_sigma",ISCAL = ind_sigma)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING ind_sigma: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "pow_sig"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"pow_sig",FSCAL = pow_sig)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING pow_sig: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "lev_bot"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"lev_bot",FSCAL = lev_bot)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING lev_bot: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "lev_sur"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"lev_sur",FSCAL = lev_sur)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING lev_sur: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "KSw"
!------------------------------------------------------------------------------|

     ISCAN = SCAN_FILE(FNAME,"KSw",FVEC = REALVEC,NSZE = KTEMP)
     IF(ISCAN /= 0)THEN
       WRITE(*,*)'ERROR READING KSw: ',ISCAN
       STOP
     END IF

     if (REALVEC(1).eq.0.0) then
	  KTEMP = 0
	  ALLOCATE(KSw(1)); KSw=0.0_WP
	 else
	 ALLOCATE(KSw(KTEMP)); KSw=0.0_WP
     KSw(1:KTEMP)= REALVEC(1:KTEMP)
	 endif


!------------------------------------------------------------------------------|
!      "Kbw"
!------------------------------------------------------------------------------|
     ISCAN = SCAN_FILE(FNAME,"KBw",FVEC = REALVEC,NSZE = KTEMP)
     IF(ISCAN /= 0)THEN
       WRITE(*,*)'ERROR READING KBw: ',ISCAN
       STOP
     END IF

     if (REALVEC(1).eq.0.0) then
	  KTEMP = 0
	  ALLOCATE(KBw(1)); KBw=0.0_WP
	 else
	 ALLOCATE(KBw(KTEMP))
     KBw(1:KTEMP)= REALVEC(1:KTEMP)
	 endif

!------------------------------------------------------------------------------|
!     "D_ref_max"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"D_ref_max",FSCAL = D_ref_max)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING D_ref_max: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "Vert_turb"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"Vert_turb",ISCAL = ver_mix)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING Vert_turb: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!  Rivers
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"key_rivers",LVAL = key_rivers)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING key_rivers: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "cyclic_length"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"cyclic_length",FSCAL = cyclic_length)
   IF((ISCAN /= 0).or. (cyclic_length>360).or. (cyclic_length<=0)) THEN
     WRITE(*,*)'ERROR READING cyclic_length: ',ISCAN
     STOP
   END IF
   cyclic_length=cyclic_length*rad
!------------------------------------------------------------------------------|
!     "Dcr"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"Dcr",FSCAL = Dcr)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING Dcr: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "Dmin"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"Dmin",FSCAL = Dmin)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING Dmin: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "TF_presence" - turn on (T), turn off (F) tidal forcing
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"TF_presence",LVAL = TF_presence)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING TF: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "T_potential"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"T_potential",LVAL = T_potential)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING T_potential: ',ISCAN
     STOP
   END IF

 If ((TF_presence).or.(T_potential)) THEN
!------------------------------------------------------------------------------|
!     "Harmonics_tf"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"Harmonics_tf",CVEC = STRINGVEC,NSZE = NOUT_VARS)
   if(iscan /= 0)then
     WRITE(*,*)'error reading Harmonics_tf: ',iscan
     STOP
   end if
   if(nout_vars <= 0)then
     WRITE(*,*)'incorrect number of netcdf Harmonics_tf variables specified'
     WRITE(*,*)'in input file',nout_vars
     STOP
   end if
   allocate(Harmonics_tf(NOUT_VARS))
   Harmonics_tf(1:NOUT_VARS)= STRINGVEC(1:NOUT_VARS)
n_th=0
do i=1,NOUT_VARS
do j=1,12
if (Harmonics_tf(i).eq.Harmonics_tf_base(j)) then
a_th(j)=1
n_th=n_th+1
exit
endif
end do
end do

!write(*,*) a_th
endif
!------------------------------------------------------------------------------|
!     "T_out"   !!CONTROLS TIDAL OUTPUT
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"T_out",LVAL = T_out)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING T_out: ',ISCAN
     STOP
   END IF
   if (T_out) then
!------------------------------------------------------------------------------|
!     "T_out_el_va"   !!CONTROLS TIDAL OUTPUT
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"T_out_el_va",LVAL = T_out_el_va)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING T_out_el_va: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "T_out_3Dvel"   !!CONTROLS TIDAL OUTPUT
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"T_out_3Dvel",LVAL = T_out_3Dvel)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING T_out_3Dvel: ',ISCAN
     STOP
   END IF
   If ((T_out).and.(.not.(T_out_el_va)).and.(.not.(T_out_3Dvel))) then
        WRITE(*,*)'T_out is T, but there is no value for output'
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "T_period"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"T_period",FSCAL = T_period)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING T_period: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "T_step"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"T_step",FSCAL = T_step)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING T_step: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "T_aval"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"T_aval",ISCAL = T_aval)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING T_aval: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "T_num"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"T_num",ISCAL = T_num)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING T_num: ',ISCAN
     STOP
   END IF

   If (mod(real(T_period*3600.0_WP/(T_aval*T_step)),1.0_WP)/= 0.0_WP) then
        WRITE(*,*)'T_step/T_period*3600/(T_aval*i)', T_period*3600.0_WP/(real(T_aval)*T_step)
     STOP
   END IF
endif
!------------------------------------------------------------------------------|
!     "key_nc_output"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"key_nc_output",LVAL = key_nc_output)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING key_nc_output: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "key_atm"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"key_atm",LVAL = key_atm)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING key_atm: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "key_obc"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"key_obc",LVAL = key_obc)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING key_obc: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "key_ic"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"key_ic",LVAL = key_ic)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING key_ic: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "type_swr_body"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"type_swr_body",ISCAL = type_swr_body)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING type_swr_body: ',ISCAN
     STOP
   END IF

   If ((.not.(type_swr_body==1)).and.(.not.(type_swr_body==2))) then
        WRITE(*,*)'type_swr_body should be 1 or 2'
     STOP
   END IF

!------------------------------------------------------------------------------|
!     "swr_bot_refl_part"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"swr_bot_refl_part",FSCAL = swr_bot_refl_part)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING swr_bot_refl_part: ',ISCAN
     STOP
   END IF
   If ((swr_bot_refl_part>1.0_WP).or.(swr_bot_refl_part<0.0_WP)) then
        WRITE(*,*)'swr_bot_refl_part should be between 0 and 1'
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "swr_bot_refl_part"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"swr_bot_min",FSCAL = swr_bot_min)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING swr_bot_min: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "jer_const"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"jer_const",LVAL = jer_const)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING jer_const: ',ISCAN
     STOP
   END IF
if (jer_const) then
!------------------------------------------------------------------------------|
!     "R"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"R",FSCAL = Rjc)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING R: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "a"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"a",FSCAL = ajc)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING a: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "b"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"b",FSCAL = bjc)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING a: ',ISCAN
     STOP
   END IF

endif
!------------------------------------------------------------------------------|
!     "cl_solid_boder_relax"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"cl_solid_boder_relax",LVAL = cl_solid_boder_relax)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING cl_solid_boder_relax: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "cl_relax"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"cl_relax",LVAL = cl_relax)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING cl_relax: ',ISCAN
     STOP
   END IF
if (cl_relax) then
!------------------------------------------------------------------------------|
!     "clim_relax"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"clim_relax",FSCAL = clim_relax)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING clim_relax: ',ISCAN
     STOP
   END IF
 clim_relax=1.0_WP/(clim_relax*3600.0_WP*24.0_WP)
!write(*,*) 'clim_relax',clim_relax
!------------------------------------------------------------------------------|
!     "relax2clim_ac"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"relax2clim_ac",FSCAL = relax2clim_ac)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING relax2clim_ac: ',ISCAN
     STOP
   END IF
 relax2clim_ac=1.0_WP/(relax2clim_ac*3600.0_WP*24.0_WP)

endif

 if (mype==0) write(*,*) 'relax2clim_ac',relax2clim_ac

!------------------------------------------------------------------------------|
!     "WET_DRY_ON"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"WET_DRY_ON",LVAL = WET_DRY_ON)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING WET_DRY_ON: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "OB_nudg" -Temp/Salt Series Nudging
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"OB_nudg",LVAL = OB_nudg)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING OB_nudg: ',ISCAN
     STOP
   END IF

!------------------------------------------------------------------------------|
! obn stability regulation
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"mixed_BP_grad_OB",LVAL = mixed_BP_grad_OB)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING mixed_BP_grad_OB: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "Cd_bz"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"Cd_bz",FSCAL = Cd_bz)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING Cd_bz: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "Bp_bz"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"Bp_bz",FSCAL = Bp_bz)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING Bp_bz: ',ISCAN
     STOP
   END IF

   if (OB_nudg) then
   !------------------------------------------------------------------------------|
!     "ALPHA_nudg"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"ALPHA_nudg",FSCAL = ALPHA_nudg)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING ALPHA_nudg: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "RADIUS"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"RADIUS",FSCAL = RADIUS)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING RADIUS: ',ISCAN
     STOP
   END IF

   endif

!------------------------------------------------------------------------------|
!     "Sponge layer"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"SL_obn",LVAL = SL_obn)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING SL_obn: ',ISCAN
     STOP
   END IF
!------------------------------------------------------------------------------|
!     "Sponge layer, user defined"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"SL_obn_user_def",LVAL = SL_obn_user_def)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING SP_obn_user_def: ',ISCAN
     STOP
   END IF
   if (SL_obn_user_def) then
   if (.not.(SL_obn)) then
   WRITE(*,*)'ERROR, SL_obn_user_def - true, but SL_obn - false: '
   endif
   endif
   if (SL_obn) then

!------------------------------------------------------------------------------|
!     "Sponge layer radius"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"SL_radius",FSCAL = SL_radius)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING SL_radius: ',ISCAN
     STOP
   END IF

   endif
!------------------------------------------------------------------------------|
!     "Sediments"
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"comp_sediment",LVAL = comp_sediment)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING comp_sediment: ',ISCAN
     STOP
   END IF
   if (comp_sediment) then

   ISCAN = SCAN_FILE(TRIM(FNAME),"sed_boun_flux",LVAL = sed_boun_flux)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING sed_boun_flux: ',ISCAN
     STOP
   END IF

   ISCAN = SCAN_FILE(TRIM(FNAME),"cl_relax_sed",LVAL = cl_relax_sed)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING cl_relax_sed: ',ISCAN
     STOP
   END IF

  if (sed_boun_flux) then

   ISCAN = SCAN_FILE(FNAME,"sed_relax",FSCAL = sed_relax)
   IF(ISCAN /= 0) THEN
   WRITE(*,*)'ERROR READING sed_relax: ',ISCAN
   STOP
   END IF
   sed_relax=1/(sed_relax*3600.0_WP*24.0_WP)
   ISCAN = SCAN_FILE(FNAME,"relax2sed",FSCAL = relax2sed)
   IF(ISCAN /= 0) THEN
   WRITE(*,*)'ERROR READING relax2sed: ',ISCAN
   STOP
   END IF
   relax2sed=1/(relax2sed*3600.0_WP*24.0_WP)

  endif

   ISCAN = SCAN_FILE(FNAME,"d_m",FSCAL = d_m)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING d_m: ',ISCAN
     STOP
   END IF
   ISCAN = SCAN_FILE(FNAME,"plop_s",FSCAL = plop_s)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING plop_s: ',ISCAN
     STOP
   END IF
   ISCAN = SCAN_FILE(FNAME,"e_p",FSCAL = e_p)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING e_p: ',ISCAN
     STOP
   END IF
   ISCAN = SCAN_FILE(FNAME,"Er0",FSCAL = Er0)
   IF(ISCAN /= 0) THEN
     WRITE(*,*)'ERROR READING Er0: ',ISCAN
     STOP
   END IF
   endif

!==============================================================================|
!            SET PHYSICAL PARAMETERS                                           !
!==============================================================================|

   Mt=int(dt/dt_2D)
   ini_time=1.0
   DAYS=NSTEPS*dt/24.0_WP/3600.0_WP


!==============================================================================|
!            ERROR CHECKING                                                    !
!==============================================================================|


   IF(trim(RESTARTC) /= 'cold_start' .AND. trim(RESTARTC) /= 'hot_start')THEN
     WRITE(*,*) 'RESTART NOT CORRECT --->',RESTARTC
     WRITE(*,*) 'SHOULD BE "cold_start" or "hot_start"'
     STOP
   END IF

   IF (trim(RESTARTC)=='cold_start') then
   RESTART=.false.
   else
   RESTART=.true.
   endif

   IF(trim(type_vd) /= 'sigma') THEN
     WRITE(*,*) 'type_vert_dis is not correct --->',type_vd
     WRITE(*,*) 'SHOULD BE "sigma"'
     STOP
   END IF

 !  IF(WINDTYPE /= 'stress' .AND. WINDTYPE /= 'speed') THEN
 !    WRITE(*,*) 'WINDTYPE NOT CORRECT,--->',WINDTYPE
 !    WRITE(*,*) 'SHOULD BE "stress" or "speed"'
 !    STOP
 !  END IF
 !  IF(M_TYPE /= 'uniform' .AND. M_TYPE /= 'non-uniform') THEN
 !    WRITE(*,*) 'M_TYPE NOT CORRECT,--->',M_TYPE
 !    WRITE(*,*) 'SHOULD BE "uniform" or "non-uniform"'
 !    STOP
 !  END IF
   IF(ver_mix<1 .AND. ver_mix<3) THEN
     WRITE(*,*) 'Vert_turb NOT CORRECT,--->',ver_mix
     WRITE(*,*) 'SHOULD BE "1", "2" or "3" '
     STOP
   END IF
 !  IF(H_TYPE /= 'body_h' .AND. H_TYPE /= 'flux_h') THEN
 !    WRITE(*,*) 'H_TYPE NOT CORRECT,--->',H_TYPE
 !    WRITE(*,*) 'SHOULD BE "body_h" or "flux_h"'
 !    STOP
 !  END IF

  ! =========================================================================
  ! Inizialize type of task (barotropic 2D, barotropic 3D, baroclinic......)
  !==========================================================================
  type_task=0
  if (Barotropic) then
  type_task=1
  endif
  if (Barotropic_3D) then
  type_task=2
  endif
  if ((SALINITY_ON).or.(TEMP_ON)) then
  type_task=3
  else
  if (type_task==0)then
  write(*,*) 'No tracer, but task is set to baroclinic!'
  STOP
  endif
  endif
  if (mype==0) write(*,*) 'type_task = ', type_task
  if (mype==0) write(*,*) 'T_const, S_const ', T_const, S_const
  !==============================================================================|
  !            SCREEN REPORT OF SET VARIABlES                                    !
  !==============================================================================|

  if (mype==0) then
     if (cartesian) then
        WRITE(*,*)'!  # Type of the coordinates             : cartesian '
     else
        WRITE(*,*)'!  # Type of the coordinates             : spherical '
     endif
     WRITE(*,*)'!  # Nsteps                              : ',NSTEPS
     WRITE(*,*)'!  # Julian day (start of sim.)          : ',Time_jd0
     WRITE(*,*)'!  # Restart                             : ',TRIM(RESTARTC)
     WRITE(*,*)'!  # A_hor                               : ',A_hor
     WRITE(*,*)'!  # Min_depth                           : ',Dmin
     if (BFC/=2) then
        WRITE(*,*)'!  # Cd                                  : ',C_d
     else
        WRITE(*,*)'!  # Z0                                  : ', z0b_min +  za
     endif
     WRITE(*,*)'!  # surface roughness height(min value) : ',z0s_min
     WRITE(*,*)'!  # bottom roughness height(min value)  : ',z0b_min
     WRITE(*,*)'!  # roughness caused by suspended sed.  : ', za
     if (comp_sediment) then
        WRITE(*,*)'!  # Sediment module                     : Activated'
        WRITE(*,*)'!  # grain size                          : ',d_m
        WRITE(*,*)'!  # density of individual grain         : ',plop_s
        WRITE(*,*)'!  # porous of the grain                 : ',e_p
     else
        WRITE(*,*)'!  # comp_sediment                        : ', comp_sediment
     endif
     if (TF_presence) then
        WRITE(*,*)'!  # Harmonics_tf                        : ', Harmonics_tf
     endif
     if (T_potential) then
        WRITE(*,*)'!  # Tidal potential                     : Activated'
     endif
     WRITE(*,*)'!  # IRESTART                            : ',IRESTART
     WRITE(*,*)'!  # Ireport                             : ',IREP
     WRITE(*,*)'!  # Irecord                             : ',IREC
     WRITE(*,*)'!  # CASE TITLE                          : ',trim(TITLE)
     WRITE(*,*)'!  # OUTDIR                              : ',trim(OUTDIR)
     WRITE(*,*)'!  # MESHPATH                            : ',trim(MESHPATH)

     IF(WET_DRY_ON)THEN
        WRITE(*,*)'!  # Wetting/Drying                      : ACTIVE'
     ELSE
        WRITE(*,*)'!  # Wetting/Drying                      : INACTIVE'
     END IF
     IF(BAROTROPIC)THEN
        WRITE(*,*)'!  # BAROTROPIC RUN                      : ACTIVE'
     END IF
     IF(BAROTROPIC_3D)THEN
        WRITE(*,*)'!  # BAROTROPIC 3D RUN                   : ACTIVE'
     END IF
     IF(TEMP_ON .AND. .NOT. BAROTROPIC)THEN
        WRITE(*,*)'!  # TEMPERATURE EQUATION                : ACTIVE'
     END IF
     IF(SALINITY_ON .AND. .NOT. BAROTROPIC)THEN
        WRITE(*,*)'!  # SALINITY EQUATION                   : ACTIVE'
     END IF
     if (type_task>1) then
        WRITE(*,*)'!  # Number sigma levels                 : ',nsigma
        if (ver_mix==1) then
           WRITE(*,*)'!  # Vert_turb                           : Pacanowski and Philander, 1981 (this is schema for deep ocean)'
        else
           if (ver_mix==2) then
              WRITE(*,*)'!  # Vert_turb                           : GOTM'
           else
              if (ver_mix==3) then
                 WRITE(*,*)'!  # Vert_turb                           : b-l'
              else
                 WRITE(*,*)'!  # Vert_turb                           : no vertical mixing'
              endif
           endif
	endif
        IF(OB_nudg)THEN
           WRITE(*,*)'!  # Open boundary T/S nudging           :  ACTIVE'
           WRITE(*,*)'!  # Nudging  coefficient                : ',ALPHA_nudg
           WRITE(*,*)'!  # Nudging  coeff                      : ',RADIUS
        END IF
        WRITE(*,*)'!  # Baroclinic time steps               : ',dt
     endif
     WRITE(*,*)'!  # Barotropic time steps               : ',dt_2d
     WRITE(*,*)'!  # Amount of simulated days            : ',DAYS
     IF (Coriolis_TF) then
        WRITE(*,*) '!  # Coriolis                            : Turn on'
     else
        WRITE(*,*) '!  # Coriolis                            : Turn off'
     endif
  end if
  RETURN

END SUBROUTINE READ_DATA_RUN

