!=============================================================================
! MODULE: Module Name
!
!> @author
!> Felix Dietzsch
!
! DESCRIPTION:
!> Defines some user specified input parameters
!>
! REVISION HISTORY
! 25 07 2011 - Initial Version
!=============================================================================
MODULE INIT
    USE NRTYPE
    INTEGER(SP) :: DIM
    REAL(DP) :: DX
    NAMELIST /NML_INIT/ DIM,DX
END MODULE INIT
