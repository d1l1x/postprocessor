!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> @details
!>Defines some basic kinds for variable definitions.
!> It is strongly recommended to use this module in every
!> module to be defined.
!> @brief
!>Defines some basic kinds for variable definitions.
!
! REVISION HISTORY:
! 25 07 2011 Initial Version
!
!=============================================================================
MODULE NRTYPE
    INTEGER,PARAMETER :: SPC =KIND((1.0,1.0))
    INTEGER,PARAMETER :: DPC =KIND((1.0D0,1.0D0))
    INTEGER,PARAMETER :: SP = KIND(1.0)
    INTEGER,PARAMETER :: DP = KIND(1.0D0)
    INTEGER,PARAMETER :: I4B = SELECTED_INT_KIND(R=9)
    INTEGER, PARAMETER :: LGT = KIND(.true.)
    REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
    REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
END MODULE NRTYPE
!
