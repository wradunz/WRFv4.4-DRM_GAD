!WRF:MEDIATION_LAYER:PHYSICS

!==============================================================================
!
! Copyright 2009. Lawrence Livermore National Security, LLC. All rights reserved.
! This work was produced at the Lawrence Livermore National Laboratory (LLNL) under
! contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE)
! and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL. Copyright
! is reserved to Lawrence Livermore National Security, LLC for purposes of controlled
! dissemination, commercialization through formal licensing, or other disposition under
! terms of Contract 44; DOE policies, regulations and orders; and U.S. statutes. The rights
! of the Federal Government are reserved under Contract 44.
!
! DISCLAIMER
! This work was prepared as an account of work sponsored by an agency of the United States
! Government. Neither the United States Government nor Lawrence Livermore National
! Security, LLC nor any of their employees, makes any warranty, express or implied, or
! assumes any liability or responsibility for the accuracy, completeness, or usefulness of
! any information, apparatus, product, or process disclosed, or represents that its use
! would not infringe privately-owned rights. Reference herein to any specific commercial
! products, process, or service by trade name, trademark, manufacturer or otherwise does
! not necessarily constitute or imply its endorsement, recommendation, or favoring by the
! United States Government or Lawrence Livermore National Security, LLC. The views and
! opinions of authors expressed herein do not necessarily state or reflect those of the
! United States Government or Lawrence Livermore National Security, LLC, and shall not be
! used for advertising or product endorsement purposes.
!
! LICENSING REQUIREMENTS
! Any use, reproduction, modification, or distribution of this software or documentation
! for commercial purposes requires a license from Lawrence Livermore National Security,
! LLC. Contact: Lawrence Livermore National Laboratory, Industrial Partnerships Office,
! P.O. Box 808, L-795, Livermore, CA 94551
!
!=============================================================================
!
! Modification History: 
!
! Implemented 12/2009 by Jeff Mirocha, jmirocha@llnl.gov
!
! 01/22: modified by Xiaoming Shi (https://github.com/shixm-cloud/WRF-RNA)
!        - ASW summarizing XS modifications
!        - added code for DRM to work with moisture and microphysics
!
! 03/24: modified by Adam Wise (adamwise@berkeley.edu)
!        - added logic to decouple DWL and DRM (DWL-SGS + RSFS)
!
!=============================================================================

MODULE module_sfs_driver

CONTAINS

!=============================================================================

SUBROUTINE sfs_driver( grid, config_flags, &
                       nba_mij, n_nba_mij, & 
                       nba_rij, n_nba_rij, & 
          			   rsfs_rec, n_rsfs_rec,     &   ! XS
          			   rsfs_filt, n_rsfs_filt,   &   ! XS
          			   rsfs_prec, n_rsfs_prec,   &   ! XS
          			   rsfs_rmij, n_rsfs_rmij,   &   ! XS
          			   rsfs_qmij, n_rsfs_qmij,   &   ! XS
          			   rsfs_qsij, n_rsfs_qsij,   &   ! XS
          			   tfu_drm, n_tfu_drm,       &   ! XS
          			   uc_drm, n_uc_drm,         &   ! XS
          			   tfuc_drm, n_tfuc_drm,     &   ! XS
          			   sijc_drm, n_sijc_drm,     &   ! XS
          			   tfsij_drm, n_tfsij_drm,   &   ! XS
          			   tfsijc_drm, n_tfsijc_drm, &   ! XS
          			   hij_drm, n_hij_drm,       &   ! XS
          			   hijc_drm, n_hijc_drm,     &   ! XS
          			   tfrmij_drm, n_tfrmij_drm, &   ! XS
          			   tgrmij_drm, n_tgrmij_drm, &   ! XS
          			   evc_drm, n_evc_drm,       &   ! XS
          			   moist, n_moist,           &   ! XS
          			   scalar, n_scalar, bn2)        ! XS

!-----------------------------------------------------------------------------
!
! PURPOSE: Calls turbulence subfilter stress model subroutines and handles
!          all MPI and OMP operations
!
!-----------------------------------------------------------------------------


! Driver layer modules
  USE module_domain
  USE module_model_constants
  USE module_configure
  USE module_tiles
  USE module_machine
  USE module_state_description
! Model layer modules
  USE module_bc

!! *** add new modules of schemes here

  USE module_sfs_nba
  USE module_sfs_drm    ! XS 
  USE module_sfs_rsfs   ! XS   
#ifdef DM_PARALLEL
   USE module_dm
   USE module_comm_dm, ONLY : &
                           HALO_EM_NBA_RIJ_sub   &
                          ,PERIOD_EM_NBA_RIJ_sub   &
                          ,HALO_EM_NBA_MIJ_sub   &
                          ,PERIOD_EM_NBA_MIJ_sub &
                          ,HALO_EM_RSFS_REC_sub  &    !---------- XS --------- 
                          ,PERIOD_EM_RSFS_REC_sub  &                           
                          ,HALO_EM_RSFS_FILT_sub  &                            
                          ,PERIOD_EM_RSFS_FILT_sub  &                          
                          ,HALO_EM_RSFS_PREC_sub  &                            
                          ,PERIOD_EM_RSFS_PREC_sub  &                          
                          ,HALO_EM_RSFS_RMIJ_sub   &                           
                          ,PERIOD_EM_RSFS_RMIJ_sub  &                          
                          ,HALO_EM_RSFS_QMIJ_sub   &                           
                          ,PERIOD_EM_RSFS_QMIJ_sub  &                          
                          ,HALO_EM_RSFS_QSIJ_sub   &                           
                          ,PERIOD_EM_RSFS_QSIJ_sub  &                          
                          ,HALO_EM_RSFS_UVWS_sub  &                            
                          ,PERIOD_EM_RSFS_UVWS_sub  &                          
                          ,HALO_EM_RSFS_PUVWS_sub  &                           
                          ,PERIOD_EM_RSFS_PUVWS_sub  &                         
                          ,HALO_EM_TFU_DRM_sub  &                              
                          ,PERIOD_EM_TFU_DRM_sub  &                            
                          ,HALO_EM_UC_DRM_sub  &                               
                          ,PERIOD_EM_UC_DRM_sub  &                             
                          ,HALO_EM_TFUC_DRM_sub  &                             
                          ,PERIOD_EM_TFUC_DRM_sub  &                           
                          ,HALO_EM_TFRMIJ_DRM_sub  &                           
                          ,PERIOD_EM_TFRMIJ_DRM_sub  &                         
                          ,HALO_EM_SIJC_DRM_sub   &                            
                          ,PERIOD_EM_SIJC_DRM_sub   &                          
                          ,HALO_EM_HIJ_DRM_sub   &                             
                          ,PERIOD_EM_HIJ_DRM_sub  &                            
                          ,HALO_EM_HIJC_DRM_sub  &                             
                          ,PERIOD_EM_HIJC_DRM_sub  &
                          !,HALO_EM_EVC_DRM_H_sub & ! Julia
                          !,PERIOD_EM_EVC_DRM_H_sub & ! Julia
                          ,period_bdy_em_phy_bc_sub & ! Julia
                          ,HALO_EM_TFSIJ_DRM_sub  &                            
                          ,PERIOD_EM_TFSIJ_DRM_sub  &                          
                          ,HALO_EM_TFSIJC_DRM_sub  &                           
                          ,PERIOD_EM_TFSIJC_DRM_sub   &                        
                          ,HALO_EM_EVC_DRM_sub   &                             
                          ,PERIOD_EM_EVC_DRM_sub   &                           
                          ,HALO_EM_EV_sub   &                                  
                          ,PERIOD_EM_EV_sub  &                                 
                          ,HALO_EM_NWTAU_sub   &  !ASW                         
                          ,PERIOD_EM_NWTAU_sub    !ASW                         

#endif

  IMPLICIT NONE

! Input data.

  TYPE(domain) , TARGET          :: grid

  TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags

  INTEGER, INTENT(  IN ) :: n_nba_mij, n_nba_rij

  INTEGER, INTENT(  IN ) :: n_rsfs_rec, n_rsfs_filt, n_rsfs_prec, n_rsfs_rmij, n_rsfs_qmij, n_rsfs_qsij, &  
	  n_tfu_drm, n_uc_drm, n_tfuc_drm, &                                                                  
	  n_sijc_drm, n_tfsij_drm, n_tfsijc_drm, n_hij_drm,  &                                                
	  n_hijc_drm, n_tfrmij_drm, n_tgrmij_drm, n_evc_drm, &                                                
	  n_moist, n_scalar  ! XS                                                                             
								                                                                                                              

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_nba_mij) &
  :: nba_mij

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_nba_rij) &
  :: nba_rij

  !--- XS ---
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_rsfs_rec) &
  :: rsfs_rec

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_rsfs_filt) &
  :: rsfs_filt

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_rsfs_prec) &
  :: rsfs_prec

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_rsfs_rmij) &
  :: rsfs_rmij

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_rsfs_qmij) &
  :: rsfs_qmij
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_rsfs_qsij) &
  :: rsfs_qsij

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_tfu_drm) &
  :: tfu_drm

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_uc_drm) &
  :: uc_drm

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_tfuc_drm) &
  :: tfuc_drm

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_sijc_drm) &
  :: sijc_drm
  
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_tfsij_drm) &
  :: tfsij_drm
  
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_tfsijc_drm) &
  :: tfsijc_drm
  
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_hij_drm) &
  :: hij_drm
  
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_hijc_drm) &
  :: hijc_drm
  
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_tfrmij_drm) & 
  :: tfrmij_drm
  
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_tgrmij_drm) &
  :: tgrmij_drm
  
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_evc_drm) &
  :: evc_drm
  
  ! XS added moist variables below. 9/4/2019
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist), &
	  INTENT(IN) ::  moist

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar), &
      INTENT(IN) ::  scalar
  
  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33), &
      INTENT(IN) ::  bn2
!-----------------------------------------


! Local data

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33) &     
  :: temp, qvapor                                                                    
  ! qvapor below is a dummy for all moist variables. It is just for historical reason
  ! we call it qvapor. XS                                                           


  INTEGER :: k_start , k_end, its, ite, jts, jte
  INTEGER :: ids , ide , jds , jde , kds , kde , &
             ims , ime , jms , jme , kms , kme , &
             ips , ipe , jps , jpe , kps , kpe

  INTEGER :: imsx, imex, jmsx, jmex, kmsx, kmex, &
             ipsx, ipex, jpsx, jpex, kpsx, kpex, &
             imsy, imey, jmsy, jmey, kmsy, kmey, &
             ipsy, ipey, jpsy, jpey, kpsy, kpey
 
  INTEGER :: ij, i, j, k, nq
                                                       
  ! RSFS variables-----------------------------------  
  INTEGER :: recn                                      
  REAL :: tmp !jss                                     
  REAL, DIMENSION(0:10) :: rcoef                       
  ! DRM variables------------------------------------- 
  INTEGER :: itst                                      
  INTEGER :: iexp                                      
  INTEGER :: switch                                    

  CALL get_ijk_from_grid ( grid ,                              &
                           ids, ide, jds, jde, kds, kde,       &
                           ims, ime, jms, jme, kms, kme,       &
                           ips, ipe, jps, jpe, kps, kpe,       &
                           imsx, imex, jmsx, jmex, kmsx, kmex, &
                           ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                           imsy, imey, jmsy, jmey, kmsy, kmey, &
                           ipsy, ipey, jpsy, jpey, kpsy, kpey  )

  k_start         = kps
  k_end           = kpe

  iexp = 1   ! XS
  itst = 3   ! XS
                 
! Compute these starting and stopping locations for each tile and number of tiles.
! See: http://www.mmm.ucar.edu/wrf/WG2/topics/settiles

! Solve_em has already called this, so should not be necessary to reset tiles here
  CALL set_tiles ( ZONE_SFS, grid , ids , ide , jds , jde , ips , ipe , jps , jpe )

  IF ( (config_flags%sfs_opt .EQ. 1) .OR. (config_flags%sfs_opt .EQ. 2) ) THEN

!=======================================================================
!
!                                BEGIN NBA
!
!=======================================================================
      
!  IF ( grid%itimestep .EQ. 1 ) THEN
!         
!    IF ( (config_flags%sfs_opt .EQ. 2) .AND. (config_flags%km_opt .NE. 2)) THEN
!  
!    CALL wrf_error_fatal( 'Must use km_opt=2 with sfs_opt=2' )
!    
!    ENDIF
!
!  ENDIF

!_______________________________________________________________________
!
! Compute NBA model constants
!_______________________________________________________________________


    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !---------------------------------------- 

        CALL calc_mij_constants( )

    ENDDO !-------------------------------------------------------------
    !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Compute Smn*Smn
!_______________________________________________________________________

    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !---------------------------------------- 

        CALL calc_smnsmn( nba_rij(ims,kms,jms,P_smnsmn),    &
                          grid%defor11, grid%defor22,       &
                          grid%defor33, grid%defor12,       &
                          grid%defor13, grid%defor23,       &
                          config_flags,                     &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start    , k_end                )

    ENDDO !-------------------------------------------------------------
    !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for R12, R13, R23 and smnsmn
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_NBA_RIJ.inc"
#      include "PERIOD_EM_NBA_RIJ.inc"
#endif

    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !----------------------------------------

        CALL set_physical_bc3d( nba_rij(ims,kms,jms,P_r12), 'd',  &
                                config_flags,                     &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start    , k_end                )


        CALL set_physical_bc3d( nba_rij(ims,kms,jms,P_r13), 'e',  &
                                config_flags,                     &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start    , k_end                )

        CALL set_physical_bc3d( nba_rij(ims,kms,jms,P_r23), 'f',  &
                                config_flags,                     &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start    , k_end                )

        CALL set_physical_bc3d( nba_rij(ims,kms,jms,P_smnsmn), 'c', &
                                config_flags,                       &
                                ids, ide, jds, jde, kds, kde,       &
                                ims, ime, jms, jme, kms, kme,       &
                                ips, ipe, jps, jpe, kps, kpe,       &
                                grid%i_start(ij), grid%i_end(ij),   &
                                grid%j_start(ij), grid%j_end(ij),   &
                                k_start    , k_end                  )

    ENDDO !-------------------------------------------------------------
    !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate M11, M22 and M33
!_______________________________________________________________________

    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !----------------------------------------

      CALL calc_mii( nba_mij(ims,kms,jms,P_m11),       &
                     nba_mij(ims,kms,jms,P_m22),       &
                     nba_mij(ims,kms,jms,P_m33),       &
                     grid%defor11, grid%defor22,       &
                     grid%defor33, grid%defor12,       &
                     grid%defor13, grid%defor23,       &
                     nba_rij(ims,kms,jms,P_r12),       &
                     nba_rij(ims,kms,jms,P_r13),       &
                     nba_rij(ims,kms,jms,P_r23),       &
                     nba_rij(ims,kms,jms,P_smnsmn),    &
                     grid%tke_2,                       & 
                     grid%rdzw, grid%dx, grid%dy,      &
                     config_flags,                     &
                     ids, ide, jds, jde, kds, kde,     &
                     ims, ime, jms, jme, kms, kme,     &
                     ips, ipe, jps, jpe, kps, kpe,     &
                     grid%i_start(ij), grid%i_end(ij), &
                     grid%j_start(ij), grid%j_end(ij), &
                     k_start, k_end                    )

    ENDDO !-------------------------------------------------------------
    !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate M12
!_______________________________________________________________________

    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !----------------------------------------

      CALL calc_m12( nba_mij(ims,kms,jms,P_m12),       &
                     grid%defor11, grid%defor22,       &
                     grid%defor12, grid%defor13,       &
                     grid%defor23,                     &
                     nba_rij(ims,kms,jms,P_r12),       &
                     nba_rij(ims,kms,jms,P_r13),       &
                     nba_rij(ims,kms,jms,P_r23),       &
                     nba_rij(ims,kms,jms,P_smnsmn),    &
                     grid%tke_2,                       & 
                     grid%rdzw, grid%dx, grid%dy,      &
                     config_flags,                     &
                     ids, ide, jds, jde, kds, kde,     &
                     ims, ime, jms, jme, kms, kme,     &
                     ips, ipe, jps, jpe, kps, kpe,     &
                     grid%i_start(ij), grid%i_end(ij), &
                     grid%j_start(ij), grid%j_end(ij), &
                     k_start, k_end                    )

    ENDDO !-------------------------------------------------------------
    !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate M13
!_______________________________________________________________________

    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !----------------------------------------

      CALL calc_m13( nba_mij(ims,kms,jms,P_m13),       &
                     grid%defor11, grid%defor33,       &
                     grid%defor12, grid%defor13,       &
                     grid%defor23,                     &
                     nba_rij(ims,kms,jms,P_r12),       &
                     nba_rij(ims,kms,jms,P_r13),       &
                     nba_rij(ims,kms,jms,P_r23),       &
                     nba_rij(ims,kms,jms,P_smnsmn),    &
                     grid%tke_2,                       & 
                     grid%rdzw, grid%dx, grid%dy,      &
                     grid%fnm, grid%fnp,               &
                     config_flags,                     &
                     ids, ide, jds, jde, kds, kde,     &
                     ims, ime, jms, jme, kms, kme,     &
                     ips, ipe, jps, jpe, kps, kpe,     &
                     grid%i_start(ij), grid%i_end(ij), &
                     grid%j_start(ij), grid%j_end(ij), &
                     k_start, k_end                    )

    ENDDO !-------------------------------------------------------------
    !$OMP END PARALLEL DO
!_______________________________________________________________________
!
! Calculate M23
!_______________________________________________________________________

    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !----------------------------------------

      CALL calc_m23( nba_mij(ims,kms,jms,P_m23),       &
                     grid%defor22, grid%defor33,       &
                     grid%defor12, grid%defor13,       &
                     grid%defor23,                     &
                     nba_rij(ims,kms,jms,P_r12),       &
                     nba_rij(ims,kms,jms,P_r13),       &
                     nba_rij(ims,kms,jms,P_r23),       &
                     nba_rij(ims,kms,jms,P_smnsmn),    &
                     grid%tke_2,                       & 
                     grid%rdzw, grid%dx, grid%dy,      &
                     grid%fnm, grid%fnp,               &
                     config_flags,                     &
                     ids, ide, jds, jde, kds, kde,     &
                     ims, ime, jms, jme, kms, kme,     &
                     ips, ipe, jps, jpe, kps, kpe,     &
                     grid%i_start(ij), grid%i_end(ij), &
                     grid%j_start(ij), grid%j_end(ij), &
                     k_start, k_end                    )

    ENDDO !-------------------------------------------------------------
    !$OMP END PARALLEL DO
!_______________________________________________________________________
!
! Update boundary conditions and halos after calculating Mij
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_NBA_MIJ.inc"
#      include "PERIOD_EM_NBA_MIJ.inc"
#endif

    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !----------------------------------------

      CALL set_physical_bc3d( nba_mij(ims,kms,jms,P_m11), 'p',    &
                              config_flags,                     &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end                )
      
      CALL set_physical_bc3d( nba_mij(ims,kms,jms,P_m22), 'p',    &
                              config_flags,                     &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end                )
      
      CALL set_physical_bc3d( nba_mij(ims,kms,jms,P_m33), 'p',    &
                              config_flags,                     &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end                )

      CALL set_physical_bc3d( nba_mij(ims,kms,jms,P_m12), 'd',    &
                              config_flags,                     &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end                )

      CALL set_physical_bc3d( nba_mij(ims,kms,jms,P_m13), 'e',    &
                              config_flags,                     &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end                )

      CALL set_physical_bc3d( nba_mij(ims,kms,jms,P_m23), 'f',    &
                              config_flags,                     &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end                )

    ENDDO !-------------------------------------------------------------
    !$OMP END PARALLEL DO

!=======================================================================
!
!                                END NBA
!
!=======================================================================

  ENDIF !(config_flags%sfs_opt .EQ. 1) .OR. (config_flags%sfs_opt .EQ. 2)



IF ( config_flags%rsfs_opt.lt.6 .and. config_flags%rsfs_opt.ge.0 ) THEN

!=======================================================================
!
!                               BEGIN RSFS
!
!=======================================================================
!_______________________________________________________________________
!
! Computation of the RSFS stresses.
!
! RSFS stresses are computed as rtau_ij = <ur_i*ur_j> - <ur_i><ur_j>
! where <> denotes an explicit filter, and ur_i are the resolvable
! subfilter-scale (RSFS) velocities, computed up to level rsfs_opt. If
! rsfs_opt .EQ. 0 then ur are set equal to the velocities at the current
! time level contained in u_2, v_2 and w_2. This is the
! so-called "scale similarity model" of Bardina. If we are more ambitous
! we can reconstruct a higher order representation of the RSFS
! velocities instead. The higher-level RSFS velocities are obtained
! using a series expression in which each higher-order term in the
! expansion is provided by an additional application of explicit filters
! to the already-filtered velocities.
!
! In rec_0 the coefficients for the expansion are assigned to rcoef, and
! the first term in the series, which requires no filtering, is computed.
! If rsfs_opt .EQ. 0, then no reconstuction is performed and rcoef = 1.0,
! and the loop up to rsfs_opt is bypassed to **.
!
!
! Else, if rsfs_opt .EQ. 1, u_2, v_2, w_2 are filtered and
! returned as ufilt, vfilt, wfilt. The surface boundary conditions
! required for each variable prior to filtering are obtained from uvw_s
! which assigns 2D arrays with the appropriate surface values. These
! arrays are then passed along wioth the 3-D field to rsfs_filter where
! a tophat filter is applied. The upper boundary condition applied to
! each variable prior to filtering is is w = 0 and zero vertical
! gradient for u and v. After filtering, ufilt, vfilt, wfilt are sent to
! rec which computes the second term in the series for urec, vrec, wrec.
!
! If rsfs_opt .GT. 1, then for recn = 2:rsfs_opt the filtering is applied
! recursively, so rather than filtering u_2, v_2, w_2 we filter
! ufilt, vfilt, wfilt, which have already been filtered. Surface
! boundary conditions for the filtered variables are computed from the
! filtered fields in uvw_s.
!
! After each application of the filter the next term in the urec,
! vrec, wrec series is added in rec.
!
! **
! When urec, vrec, wrec have been computed up to reclevel their
! halos are updated and proj_c, etc. are used to project urec, etc. to
! the vorticity points where deformation and stress are computed:
!
! d_11, d_22, d_33 (du/dx, dv/dy, dw/dz) <-> 'c' cell center points
! d_12                 (du/dy, dv/dx)    <-> 'd' corner points
! d_13                 (du/dz, dw/dx)    <-> 'e' corner points
! d_23                 (dv/dz, dw/dy)    <-> 'f' corner points
!
! Next, urec, vrec, wrec are projected to c, d, e, f. Then the
! projected RSFS velocities are filtered. Finally the RSFS stresses
! are computed in rmij as:
!
! rm11 = <urecc*urecc>-<urecc>*<urecc>
! rm22 = <vrecc*vrecc>-<vrecc>*<vrecc>
! rm33 = <wrecc*wrecc>-<wrecc>*<wrecc>
! rm12 = <urecd*vrecd>-<urecd>*<vrecd>
! rm13 = <urece*wrece>-<urece>*<wrece>
! rm23 = <vrecf*wrecf>-<vrecf>*<wrecf>!
! rmij are then sent to the diffusion subroutines where they are
! added to the SGS stresses that have been computed elsewhere.
!_______________________________________________________________________

!_______________________________________________________________________
!
! Calculate zeroth term of reconstruction at u, v, w points.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rec_0( rsfs_rec(ims,kms,jms,P_urec),           &
                rsfs_rec(ims,kms,jms,P_vrec),           &
                rsfs_rec(ims,kms,jms,P_wrec),           &
                rcoef,                                  &
                grid%u_2, grid%v_2, grid%w_2,           &
                config_flags%rsfs_opt,                  &
                ids, ide, jds, jde, kds, kde,           &
                ims, ime, jms, jme, kms, kme,           &
                ips, ipe, jps, jpe, kps, kpe,           &
                grid%i_start(ij), grid%i_end(ij),       &
                grid%j_start(ij), grid%j_end(ij),       &
                k_start, k_end                          )

    CALL rec_0_pt( rsfs_rec(ims,kms,jms,P_trec),           &
                   rcoef,                                  &
                   grid%t_2,                        &
                   config_flags%rsfs_opt,                  &
                   ids, ide, jds, jde, kds, kde,           &
                   ims, ime, jms, jme, kms, kme,           &
                   ips, ipe, jps, jpe, kps, kpe,           &
                   grid%i_start(ij), grid%i_end(ij),       &
                   grid%j_start(ij), grid%j_end(ij),       &
                   k_start, k_end                          )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!    CALL writecheckm( grid%em_u_2, ims, ime, jms, jme, kms, kme, &
!                      config_flags%klevs, config_flags%kleve )
!    WRITE(*,*)'em_u_2'

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

  DO recn = 1, config_flags%rsfs_opt

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

    IF (recn .EQ. 1 ) THEN !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!_______________________________________________________________________
!
! First time through the reconstruction loop we filter u_2, v_2,
! w_2. These variables have full halos. We need to calculate surface
! values, us, vs, ws, and update their halos.
!_______________________________________________________________________

!_______________________________________________________________________
!
! Calculate u, v, w at the surface.
!_______________________________________________________________________

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL uvw_s( grid%us, grid%vs, grid%ws,                      &
                    grid%u_2, grid%v_2,                             &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    config_flags%rsfsvsbc_opt,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )

        CALL pt_s(  grid%ts,                                        &
                    grid%t_2,                                       &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    config_flags%rsfsvsbc_opt,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )

      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc2d_m( grid%us, 'u', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%vs, 'v', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%ws, 'w', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%ts, 't', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )


      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Filter u_2, v_2, w_2 using us, vs, ws to calculate ufilt,
! vfilt, wfilt.
!_______________________________________________________________________

 IF (config_flags%nested) THEN !use dummy variables to update b.c. for u_2, v_2 and w_2

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_rmij(i,k,j,P_rm11)=grid%u_2(i,k,j)
         rsfs_rmij(i,k,j,P_rm22)=grid%v_2(i,k,j)
         rsfs_rmij(i,k,j,P_rm33)=grid%w_2(i,k,j)
         rsfs_rmij(i,k,j,P_rm12)=grid%t_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_RMIJ.inc"
#      include "PERIOD_EM_RSFS_RMIJ.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm11), 'u', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm22), 'v', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )


        CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm33), 'w', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm12), 't', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_ufilt),   &
                          rsfs_rmij(ims,kms,jms,P_rm11),    &
                          grid%us, 'uu', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_vfilt),   &
                          rsfs_rmij(ims,kms,jms,P_rm22),    &
                          grid%vs, 'vv', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_wfilt),   &
                          rsfs_rmij(ims,kms,jms,P_rm33),    &
                          grid%ws, 'ww', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_tfilt),   &
                          rsfs_rmij(ims,kms,jms,P_rm12),     &
                          grid%ts, 'cc',  iexp,             &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO


 ELSE !not nested


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_ufilt),   &
                          grid%u_2,                         &
                          grid%us, 'uu', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_vfilt),   &
                          grid%v_2,                         &
                          grid%vs, 'vv', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_wfilt),   &
                          grid%w_2,                         &
                          grid%ws, 'ww', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_tfilt),   &
                          grid%t_2,                         &
                          grid%ts, 'cc',  iexp,             &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO

 ENDIF

  ELSE !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!_______________________________________________________________________
!
! If rsfs_opt > 1, then we filter the filtered velocities, rather than
! the velocities themselves. This requires:
! i)   update of halos for ufilt, vfilt, wfilt.
! ii)  calculation of consistent surface velocities, us, vs, ws using
!      ufilt, vfilt, wfilt.
! iii) update of halos for us, vs, ws.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_FILT.inc"
#      include "PERIOD_EM_RSFS_FILT.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_ufilt), 'u', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_vfilt), 'v', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_wfilt), 'w', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_tfilt), 't', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate u, v, w at the surface using ufilt, vfilt, wfilt.
!_______________________________________________________________________

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL uvw_s( grid%us, grid%vs, grid%ws,                      &
                    rsfs_filt(ims,kms,jms,P_ufilt),                 &
                    rsfs_filt(ims,kms,jms,P_vfilt),                 &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    config_flags%rsfsvsbc_opt,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )

        CALL pt_s(  grid%ts,                                        &
                    rsfs_filt(ims,kms,jms,P_tfilt),                 &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    config_flags%rsfsvsbc_opt,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )

      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halo for us, vs, ws.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc2d_m( grid%us, 'u', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%vs, 'v', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%ws, 'w', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%ts, 't', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

      ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate ufilt, vfilt, wfilt using ufilt, vfilt, wfilt.
!_______________________________________________________________________

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_ufilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_ufilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%us, 'uu', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

       ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_vfilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_vfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%vs, 'vv', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

      ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_wfilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_wfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%ws, 'ww', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

      ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_tfilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_tfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%ts, 'cc', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

      ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO


   ENDIF ! (recn .EQ. 1 ) ELSE !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!_______________________________________________________________________
!
! Add current reconstruction-level term to urec, vrec, wrec.
!_______________________________________________________________________

    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !------------------------------

      CALL rec( rsfs_rec(ims,kms,jms,P_urec),     &
                rsfs_rec(ims,kms,jms,P_vrec),     &
                rsfs_rec(ims,kms,jms,P_wrec),     &
                rsfs_filt(ims,kms,jms,P_ufilt),   &
                rsfs_filt(ims,kms,jms,P_vfilt),   &
                rsfs_filt(ims,kms,jms,P_wfilt),   &
                rcoef(recn),                      &
                ids, ide, jds, jde, kds, kde,     &
                ims, ime, jms, jme, kms, kme,     &
                ips, ipe, jps, jpe, kps, kpe,     &
                grid%i_start(ij), grid%i_end(ij), &
                grid%j_start(ij), grid%j_end(ij), &
                k_start, k_end                    )

      CALL rec_pt( rsfs_rec(ims,kms,jms,P_trec),     &
                   rsfs_filt(ims,kms,jms,P_tfilt),   &
                   rcoef(recn),                      &
                   ids, ide, jds, jde, kds, kde,     &
                   ims, ime, jms, jme, kms, kme,     &
                   ips, ipe, jps, jpe, kps, kpe,     &
                   grid%i_start(ij), grid%i_end(ij), &
                   grid%j_start(ij), grid%j_end(ij), &
                   k_start, k_end                    )

    ENDDO !---------------------------------------------------
    !$OMP END PARALLEL DO

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

  ENDDO !recn = 1, config_flags%rsfs_opt

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

!_______________________________________________________________________
!
! Now that reconstruction is over we need to project urec, vrec, wrec
! to vorticity points, then filter. This requires first updating halos
! for urec, vrec, wrec. wrec is returned from rec (or rec_0) starting at
! k=2, so first we need to assign wrec(*,1,*) its appropriate value,
! obtained from urec and vrec.
!_______________________________________________________________________

!_______________________________________________________________________
!
! Update halos for urec, vrec, wrec.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_REC.inc"
#      include "PERIOD_EM_RSFS_REC.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !----------------------------------------

    CALL set_physical_bc3d_m( rsfs_rec(ims,kms,jms,P_urec), 'u', &
                            config_flags,                      &
                            ids, ide, jds, jde, kds, kde,      &
                            ims, ime, jms, jme, kms, kme,      &
                            ips, ipe, jps, jpe, kps, kpe,      &
                            grid%i_start(ij), grid%i_end(ij),  &
                            grid%j_start(ij), grid%j_end(ij),  &
                            k_start, k_end                     )

    CALL set_physical_bc3d_m( rsfs_rec(ims,kms,jms,P_vrec), 'v', &
                            config_flags,                      &
                            ids, ide, jds, jde, kds, kde,      &
                            ims, ime, jms, jme, kms, kme,      &
                            ips, ipe, jps, jpe, kps, kpe,      &
                            grid%i_start(ij), grid%i_end(ij),  &
                            grid%j_start(ij), grid%j_end(ij),  &
                            k_start, k_end                     )

    CALL set_physical_bc3d_m( rsfs_rec(ims,kms,jms,P_wrec), 'w', &
                            config_flags,                      &
                            ids, ide, jds, jde, kds, kde,      &
                            ims, ime, jms, jme, kms, kme,      &
                            ips, ipe, jps, jpe, kps, kpe,      &
                            grid%i_start(ij), grid%i_end(ij),  &
                            grid%j_start(ij), grid%j_end(ij),  &
                            k_start, k_end                     )

    CALL set_physical_bc3d_m( rsfs_rec(ims,kms,jms,P_trec), 't', &
                            config_flags,                      &
                            ids, ide, jds, jde, kds, kde,      &
                            ims, ime, jms, jme, kms, kme,      &
                            ips, ipe, jps, jpe, kps, kpe,      &
                            grid%i_start(ij), grid%i_end(ij),  &
                            grid%j_start(ij), grid%j_end(ij),  &
                            k_start, k_end                     )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Compute us, vs, ws from urec, vrec.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

        CALL uvw_s( grid%us, grid%vs, grid%ws,                      &
                    grid%u_2, grid%v_2,                             &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    config_flags%rsfsvsbc_opt,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )

        CALL pt_s(  grid%ts,                                        &
                    grid%t_2,                                       &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    config_flags%rsfsvsbc_opt,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for us, vs, ws.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL set_physical_bc2d_m( grid%us, 'u', config_flags,       &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vs, 'v', config_flags,       &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%ws, 'w', config_flags,       &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%ts, 't', config_flags,       &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Assign ws to wrec(*,1,*) before projecting to c.
!_______________________________________________________________________


  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !----------------------------

    CALL update_w1( rsfs_rec(ims,kms,jms,P_wrec), grid%ws, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    ips, ipe, jps, jpe, kps, kpe,          &
                    grid%i_start(ij), grid%i_end(ij),      &
                    grid%j_start(ij), grid%j_end(ij),      &
                    k_start, k_end                         )

  ENDDO !-------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Now we have the RSFS velocities, urec, vrec, wrec at u, v, w points,
! as well as their halos and surface values.
! Next we project urec, vrec, wrec and their surface values.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL proj_c( rsfs_prec(ims,kms,jms,P_urecc),    &
                 rsfs_prec(ims,kms,jms,P_vrecc),    &
                 rsfs_prec(ims,kms,jms,P_wrecc),    &
                 rsfs_prec(ims,kms,jms,P_uurecc),   &
                 rsfs_prec(ims,kms,jms,P_vvrecc),   &
                 rsfs_prec(ims,kms,jms,P_wwrecc),   &
                 rsfs_rec(ims,kms,jms,P_urec),      &
                 rsfs_rec(ims,kms,jms,P_vrec),      &
                 rsfs_rec(ims,kms,jms,P_wrec), 'u', &
                 ids, ide, jds, jde, kds, kde,      &
                 ims, ime, jms, jme, kms, kme,      &
                 ips, ipe, jps, jpe, kps, kpe,      &
                 grid%i_start(ij), grid%i_end(ij),  &
                 grid%j_start(ij), grid%j_end(ij),  &
                 k_start, k_end                     )

    CALL proj_d( rsfs_prec(ims,kms,jms,P_urecd),   &
                 rsfs_prec(ims,kms,jms,P_vrecd),   &
                 rsfs_prec(ims,kms,jms,P_uvrecd),  &
                 rsfs_rec(ims,kms,jms,P_urec),     &
                 rsfs_rec(ims,kms,jms,P_vrec),     &
                 ids, ide, jds, jde, kds, kde,     &
                 ims, ime, jms, jme, kms, kme,     &
                 ips, ipe, jps, jpe, kps, kpe,     &
                 grid%i_start(ij), grid%i_end(ij), &
                 grid%j_start(ij), grid%j_end(ij), &
                 k_start, k_end                    )

    CALL proj_e( rsfs_prec(ims,kms,jms,P_urece),   &
                 rsfs_prec(ims,kms,jms,P_wrece),   &
                 rsfs_prec(ims,kms,jms,P_uwrece),  &
                 rsfs_rec(ims,kms,jms,P_urec),     &
                 rsfs_rec(ims,kms,jms,P_wrec),     &
                 grid%fnm, grid%fnp,               &
                 ids, ide, jds, jde, kds, kde,     &
                 ims, ime, jms, jme, kms, kme,     &
                 ips, ipe, jps, jpe, kps, kpe,     &
                 grid%i_start(ij), grid%i_end(ij), &
                 grid%j_start(ij), grid%j_end(ij), &
                 k_start, k_end                    )

    CALL proj_f( rsfs_prec(ims,kms,jms,P_vrecf),   &
                 rsfs_prec(ims,kms,jms,P_wrecf),   &
                 rsfs_prec(ims,kms,jms,P_vwrecf),  &
                 rsfs_rec(ims,kms,jms,P_vrec),     &
                 rsfs_rec(ims,kms,jms,P_wrec),     &
                 grid%fnm, grid%fnp,               &
                 ids, ide, jds, jde, kds, kde,     &
                 ims, ime, jms, jme, kms, kme,     &
                 ips, ipe, jps, jpe, kps, kpe,     &
                 grid%i_start(ij), grid%i_end(ij), &
                 grid%j_start(ij), grid%j_end(ij), &
                 k_start, k_end                    )

    CALL proj_tuvw( rsfs_prec(ims,kms,jms,P_trecu),   &
                   rsfs_prec(ims,kms,jms,P_trecv),   &
                   rsfs_prec(ims,kms,jms,P_trecw),   &
                   rsfs_prec(ims,kms,jms,P_utrecu),  &
                   rsfs_prec(ims,kms,jms,P_vtrecv),  &
                   rsfs_prec(ims,kms,jms,P_wtrecw),  &
                   rsfs_rec(ims,kms,jms,P_urec),     &
                   rsfs_rec(ims,kms,jms,P_vrec),     &
                   rsfs_rec(ims,kms,jms,P_wrec),     &
                   rsfs_rec(ims,kms,jms,P_trec),     &
                   grid%fnm, grid%fnp,               &
                   ids, ide, jds, jde, kds, kde,     &
                   ims, ime, jms, jme, kms, kme,     &
                   ips, ipe, jps, jpe, kps, kpe,     &
                   grid%i_start(ij), grid%i_end(ij), &
                   grid%j_start(ij), grid%j_end(ij), &
                   k_start, k_end                    )

    CALL proj_s( grid%usc, grid%vsc,                          &
                 grid%uusc, grid%vvsc, grid%wws,              &
                 grid%usd, grid%vsd, grid%uvsd,               &
                 grid%wse, grid%uwse,                         &
                 grid%wsf, grid%vwsf,                         &
                 grid%us, grid%vs, grid%ws,                   &
                 ids, ide, jds, jde,                          &
                 ims, ime, jms, jme,                          &
                 ips, ipe, jps, jpe,                          &
                 grid%i_start(ij), grid%i_end(ij),            &
                 grid%j_start(ij), grid%j_end(ij)             )

    CALL proj_s_pt( grid%tsu, grid%tsv,                          &
                    grid%utsu, grid%vtsv, grid%wtsw,             &
                    grid%us, grid%vs, grid%ws, grid%ts,          &
                    ids, ide, jds, jde,                          &
                    ims, ime, jms, jme,                          &
                    ips, ipe, jps, jpe,                          &
                    grid%i_start(ij), grid%i_end(ij),            &
                    grid%j_start(ij), grid%j_end(ij)             )

!if (config_flags%nested) then
!
!  print*,'vrec_f before bc update'
!  CALL writecheck_t2d( rsfs_prec(ims,kms,jms,P_vrecf),   &
!                       ims, ime, jms, jme, kms, kme,     &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )
!
!endif

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for variables before filtering.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_PREC.inc"
#      include "PERIOD_EM_RSFS_PREC.inc"
#      include "HALO_EM_RSFS_PUVWS.inc"
#      include "PERIOD_EM_RSFS_PUVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_urecc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uurecc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vrecc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vvrecc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wrecc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wwrecc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_urecd), 'd', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vrecd), 'd', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uvrecd), 'd', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_urece), 'e', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wrece), 'e', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uwrece), 'e', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vrecf), 'f', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wrecf), 'f', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vwrecf), 'f', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )


!VARIABLES FOR HEAT

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_trecu), 'u',  &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_trecv), 'v',  &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_trecw), 'w',  &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_utrecu), 'u', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vtrecv), 'v', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wtrecw), 'w', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

!2D VARIABLES

    CALL set_physical_bc2d_m( grid%usc, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uusc, 't', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vsc, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vvsc, 't', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wws, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%usd, 'd', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vsd, 'd', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uvsd, 'd', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wse, 'e', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uwse, 'e', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wsf, 'f', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vwsf, 'f', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

!2D VARIABLES FOR HEAT

    CALL set_physical_bc2d_m( grid%tsu, 'u', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%tsv, 'v', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%utsu, 'u', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vtsv, 'v', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wtsw, 'w', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Filter urecc, urecc*urecc, etc.
!_______________________________________________________________________

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_urecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_urecc),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%usc, 'cc', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uurecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uurecc),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uusc, 'cc', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vrecc),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vsc, 'cc', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vvrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vvrecc),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vvsc, 'cc', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wrecc),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%ws, 'cw', iexp,              &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wwrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wwrecc),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wws, 'cw', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_urecd)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_urecd),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%usd, 'dd', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vrecd)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

   CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vrecd),   &
                     rsfs_filt(ims,kms,jms,P_temp),    &
                     grid%vsd, 'dd', iexp,             &
                     ids, ide, jds, jde, kds, kde,     &
                     ims, ime, jms, jme, kms, kme,     &
                     ips, ipe, jps, jpe, kps, kpe,     &
                     grid%i_start(ij), grid%i_end(ij), &
                     grid%j_start(ij), grid%j_end(ij), &
                     k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uvrecd)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uvrecd),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uvsd, 'dd', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_urece)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_urece),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%us, 'ee', iexp,              &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wrece)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wrece),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wse, 'ew', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uwrece)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uwrece),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uwse, 'ew', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vrecf)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vrecf),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vs, 'ff', iexp,              &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wrecf)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wrecf),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wsf, 'fw', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vwrecf)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vwrecf),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vwsf, 'fw', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO


! VARIABLES FOR HEAT

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_rec(i,k,j,P_urec)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_rec(ims,kms,jms,P_urec),     &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%us, 'uu', iexp,                    &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_rec(i,k,j,P_vrec)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_rec(ims,kms,jms,P_vrec),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vs, 'vv', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_rec(i,k,j,P_wrec)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_rec(ims,kms,jms,P_wrec),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%ws, 'ww', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO


  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_trecu)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_trecu),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%tsu, 'uu', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_trecv)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_trecv),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%tsv, 'vv', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_trecw)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_trecw),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%ts, 'ww', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_utrecu)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_utrecu),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%utsu, 'uu', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vtrecv)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vtrecv),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vtsv, 'vv', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wtrecw)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wtrecw),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wtsw, 'ww', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate rtauij.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rmij( rsfs_rmij(ims,kms,jms,P_rm11),        &
               rsfs_prec(ims,kms,jms,P_urecc),       &
               rsfs_prec(ims,kms,jms,P_urecc),       &
               rsfs_prec(ims,kms,jms,P_uurecc), 'c', &
               ids, ide, jds, jde, kds, kde,         &
               ims, ime, jms, jme, kms, kme,         &
               ips, ipe, jps, jpe, kps, kpe,         &
               grid%i_start(ij), grid%i_end(ij),     &
               grid%j_start(ij), grid%j_end(ij),     &
               k_start, k_end                        )

     CALL rmij( rsfs_rmij(ims,kms,jms,P_rm22),        &
                rsfs_prec(ims,kms,jms,P_vrecc),       &
                rsfs_prec(ims,kms,jms,P_vrecc),       &
                rsfs_prec(ims,kms,jms,P_vvrecc), 'c', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( rsfs_rmij(ims,kms,jms,P_rm33),        &
                rsfs_prec(ims,kms,jms,P_wrecc),       &
                rsfs_prec(ims,kms,jms,P_wrecc),       &
                rsfs_prec(ims,kms,jms,P_wwrecc), 'c', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( rsfs_rmij(ims,kms,jms,P_rm12),        &
                rsfs_prec(ims,kms,jms,P_urecd),       &
                rsfs_prec(ims,kms,jms,P_vrecd),       &
                rsfs_prec(ims,kms,jms,P_uvrecd), 'd', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( rsfs_rmij(ims,kms,jms,P_rm13),        &
                rsfs_prec(ims,kms,jms,P_urece),       &
                rsfs_prec(ims,kms,jms,P_wrece),       &
                rsfs_prec(ims,kms,jms,P_uwrece), 'e', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( rsfs_rmij(ims,kms,jms,P_rm23),        &
                rsfs_prec(ims,kms,jms,P_vrecf),       &
                rsfs_prec(ims,kms,jms,P_wrecf),       &
                rsfs_prec(ims,kms,jms,P_vwrecf), 'f', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )


!VARIABLE FOR HEAT

     CALL rmij( rsfs_rmij(ims,kms,jms,P_rh1),        &
                rsfs_rec(ims,kms,jms,P_urec),       &
                rsfs_prec(ims,kms,jms,P_trecu),       &
                rsfs_prec(ims,kms,jms,P_utrecu), 'u', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( rsfs_rmij(ims,kms,jms,P_rh2),        &
                rsfs_rec(ims,kms,jms,P_vrec),       &
                rsfs_prec(ims,kms,jms,P_trecv),       &
                rsfs_prec(ims,kms,jms,P_vtrecv), 'v', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( rsfs_rmij(ims,kms,jms,P_rh3),        &
                rsfs_rec(ims,kms,jms,P_wrec),       &
                rsfs_prec(ims,kms,jms,P_trecw),       &
                rsfs_prec(ims,kms,jms,P_wtrecw), 'h', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for rtauij.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_RMIJ.inc"
#      include "PERIOD_EM_RSFS_RMIJ.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

     CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm11), 't', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm22), 't', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm33), 't', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm12), 'd', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm13), 'e', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rm23), 'f', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

!VARIABLE FOR HEAT

     CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rh1), 'u', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rh2), 'v', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     CALL set_physical_bc3d_m( rsfs_rmij(ims,kms,jms,P_rh3), 'w', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO



!---------------- RSFS for moist variables ---------------------
!+  XS added on 2019-9-5
IF ( n_moist .ge. PARAM_FIRST_SCALAR ) THEN

   do nq = PARAM_FIRST_SCALAR, n_moist
      ! use qvapor as a dummy array to store the moist variable
      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               qvapor(i,k,j) = moist(i,k,j,nq)
            ENDDO
         ENDDO
      ENDDO

      ! Calculate zeroth term of reconstruction at q points.

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rec_0_qmoist( rsfs_rec(ims,kms,jms,P_qrec),       &
              rcoef,                                  &
              qvapor,                                 &
              config_flags%rsfs_opt,                  &
              ids, ide, jds, jde, kds, kde,           &
              ims, ime, jms, jme, kms, kme,           &
              ips, ipe, jps, jpe, kps, kpe,           &
              grid%i_start(ij), grid%i_end(ij),       &
              grid%j_start(ij), grid%j_end(ij),       &
              k_start, k_end                          )

      ENDDO !-----------------------------------------------------
      !$OMP END PARALLEL DO

      !    CALL writecheckm( grid%em_u_2, ims, ime, jms, jme, kms, kme, &
      !                      config_flags%klevs, config_flags%kleve )
      !    WRITE(*,*)'em_u_2'

      !/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

      DO recn = 1, config_flags%rsfs_opt

         !/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

         IF (recn .EQ. 1 ) THEN !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------

              if (nq.le.2) then
               CALL qmoist_s(  grid%qss,                                        &
                    qvapor,                                       &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    1,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )
              else
               CALL qmoist_s(  grid%qss,                                        &
                    qvapor,                                       &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    2,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )
               end if
            ENDDO !---------------------------------------------------------
            !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------

               CALL set_physical_bc2d_m( grid%qss, 't', config_flags,       &
                    ids, ide, jds, jde,               &
                    ims, ime, jms, jme,               &
                    ips, ipe, jps, jpe,               &
                    grid%i_start(ij), grid%i_end(ij), &
                    grid%j_start(ij), grid%j_end(ij)  )


            ENDDO !---------------------------------------------------------
            !$OMP END PARALLEL DO


            IF (config_flags%nested) THEN !use dummy variables to update b.c.

               DO i = ims, ime
                  DO k = kms, kme
                     DO j = jms, jme
                        rsfs_filt(i,k,j,P_temp)=qvapor(i,k,j)
                     ENDDO
                  ENDDO
               ENDDO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_FILT.inc"
#      include "PERIOD_EM_RSFS_FILT.inc"
#endif

               !$OMP PARALLEL DO   &
               !$OMP PRIVATE ( ij )
               DO ij = 1 , grid%num_tiles !------------------------------------

                  CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_temp), 't', &
                       config_flags,                        &
                       ids, ide, jds, jde, kds, kde,        &
                       ims, ime, jms, jme, kms, kme,        &
                       ips, ipe, jps, jpe, kps, kpe,        &
                       grid%i_start(ij), grid%i_end(ij),    &
                       grid%j_start(ij), grid%j_end(ij),    &
                       k_start, k_end                       )

               ENDDO !---------------------------------------------------------
               !$OMP END PARALLEL DO


               !$OMP PARALLEL DO   &
               !$OMP PRIVATE ( ij )
               DO ij = 1 , grid%num_tiles !------------------------------------

                  CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_qfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),     &
                          grid%qss, 'cc',  iexp,             &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

               ENDDO !---------------------------------------------------------
               !$OMP END PARALLEL DO


            ELSE !not nested


               !$OMP PARALLEL DO   &
               !$OMP PRIVATE ( ij )
               DO ij = 1 , grid%num_tiles !------------------------------------

                  CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_qfilt),   &
                          qvapor,                         &
                          grid%qss, 'cc',  iexp,             &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

               ENDDO !---------------------------------------------------------
               !$OMP END PARALLEL DO

            ENDIF

         ELSE !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_FILT.inc"
#      include "PERIOD_EM_RSFS_FILT.inc"
#endif

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------

               CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_qfilt), 't', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

            ENDDO !---------------------------------------------------------
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------
              if (nq.eq.1) then
               CALL qmoist_s(  grid%qss,                                        &
                    rsfs_filt(ims,kms,jms,P_qfilt),                 &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    1,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )
             else
              CALL qmoist_s(  grid%qss,                                        &
                   rsfs_filt(ims,kms,jms,P_qfilt),                 &
                   grid%ht, grid%rdx, grid%rdy,                    &
                   grid%cf1, grid%cf2, grid%cf3,                   &
                   2,                      &
                   ids, ide, jds, jde, kds, kde,                   &
                   ims, ime, jms, jme, kms, kme,                   &
                   ips, ipe, jps, jpe, kps, kpe,                   &
                   grid%i_start(ij), grid%i_end(ij),               &
                   grid%j_start(ij), grid%j_end(ij),               &
                   k_start, k_end                                  )
              end if
            ENDDO !---------------------------------------------------------
            !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------

               CALL set_physical_bc2d_m( grid%qss, 't', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

            ENDDO !-------------------------------------------------
            !$OMP END PARALLEL DO

            DO i = ims, ime
               DO k = kms, kme
                  DO j = jms, jme
                     rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_qfilt)
                  ENDDO
               ENDDO
            ENDDO

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !----------------------------

               CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_qfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%qss, 'cc', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

            ENDDO !-------------------------------------------------
            !$OMP END PARALLEL DO


         ENDIF ! (recn .EQ. 1 ) ELSE !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles !------------------------------

            CALL rec_qmoist( rsfs_rec(ims,kms,jms,P_qrec),     &
                   rsfs_filt(ims,kms,jms,P_qfilt),   &
                   rcoef(recn),                      &
                   ids, ide, jds, jde, kds, kde,     &
                   ims, ime, jms, jme, kms, kme,     &
                   ips, ipe, jps, jpe, kps, kpe,     &
                   grid%i_start(ij), grid%i_end(ij), &
                   grid%j_start(ij), grid%j_end(ij), &
                   k_start, k_end                    )

         ENDDO !---------------------------------------------------
         !$OMP END PARALLEL DO

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

      ENDDO !recn = 1, config_flags%rsfs_opt

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/


#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_REC.inc"
#      include "PERIOD_EM_RSFS_REC.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------------------

         CALL set_physical_bc3d_m( rsfs_rec(ims,kms,jms,P_qrec), 't', &
                            config_flags,                      &
                            ids, ide, jds, jde, kds, kde,      &
                            ims, ime, jms, jme, kms, kme,      &
                            ips, ipe, jps, jpe, kps, kpe,      &
                            grid%i_start(ij), grid%i_end(ij),  &
                            grid%j_start(ij), grid%j_end(ij),  &
                            k_start, k_end                     )

      ENDDO !-----------------------------------------------------
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------
        if (nq.le.2) then
        CALL qmoist_s(  grid%qss,                                        &
                   qvapor,                                       &
                   grid%ht, grid%rdx, grid%rdy,                    &
                   grid%cf1, grid%cf2, grid%cf3,                   &
                   1,                      &
                   ids, ide, jds, jde, kds, kde,                   &
                   ims, ime, jms, jme, kms, kme,                   &
                   ips, ipe, jps, jpe, kps, kpe,                   &
                   grid%i_start(ij), grid%i_end(ij),               &
                   grid%j_start(ij), grid%j_end(ij),               &
                   k_start, k_end                                  )

        else
         CALL qmoist_s(  grid%qss,                                        &
                    qvapor,                                       &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    2,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )
        end if
      ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL set_physical_bc2d_m( grid%qss, 't', config_flags,       &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

      ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL proj_quvw( rsfs_prec(ims,kms,jms,P_qrecu),   &
                   rsfs_prec(ims,kms,jms,P_qrecv),   &
                   rsfs_prec(ims,kms,jms,P_qrecw),   &
                   rsfs_prec(ims,kms,jms,P_uqrecu),  &
                   rsfs_prec(ims,kms,jms,P_vqrecv),  &
                   rsfs_prec(ims,kms,jms,P_wqrecw),  &
                   rsfs_rec(ims,kms,jms,P_urec),     &
                   rsfs_rec(ims,kms,jms,P_vrec),     &
                   rsfs_rec(ims,kms,jms,P_wrec),     &
                   rsfs_rec(ims,kms,jms,P_qrec),     &
                   grid%fnm, grid%fnp,               &
                   ids, ide, jds, jde, kds, kde,     &
                   ims, ime, jms, jme, kms, kme,     &
                   ips, ipe, jps, jpe, kps, kpe,     &
                   grid%i_start(ij), grid%i_end(ij), &
                   grid%j_start(ij), grid%j_end(ij), &
                   k_start, k_end                    )


         CALL proj_s_qmoist( grid%qsu, grid%qsv,                      &
                    grid%uqsu, grid%vqsv, grid%wqsw,             &
                    grid%us, grid%vs, grid%ws, grid%qss,          &
                    ids, ide, jds, jde,                          &
                    ims, ime, jms, jme,                          &
                    ips, ipe, jps, jpe,                          &
                    grid%i_start(ij), grid%i_end(ij),            &
                    grid%j_start(ij), grid%j_end(ij)             )

!if (config_flags%nested) then
!
!  print*,'vrec_f before bc update'
!  CALL writecheck_t2d( rsfs_prec(ims,kms,jms,P_vrecf),   &
!                       ims, ime, jms, jme, kms, kme,     &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )
!
!endif

      ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for variables before filtering.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_PREC.inc"
#      include "PERIOD_EM_RSFS_PREC.inc"
#      include "HALO_EM_RSFS_PUVWS.inc"
#      include "PERIOD_EM_RSFS_PUVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------


!VARIABLES FOR MOISTURE

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_qrecu), 'u',  &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_qrecv), 'v',  &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_qrecw), 'w',  &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uqrecu), 'u', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vqrecv), 'v', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wqrecw), 'w', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

!2D VARIABLES FOR MOISTURE

         CALL set_physical_bc2d_m( grid%qsu, 'u', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

         CALL set_physical_bc2d_m( grid%qsv, 'v', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

         CALL set_physical_bc2d_m( grid%uqsu, 'u', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

         CALL set_physical_bc2d_m( grid%vqsv, 'v', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

         CALL set_physical_bc2d_m( grid%wqsw, 'w', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

      ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO



      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_qrecu)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_qrecu),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%qsu, 'uu', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_qrecv)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_qrecv),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%qsv, 'vv', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_qrecw)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_qrecw),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%qss, 'ww', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uqrecu)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uqrecu),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uqsu, 'uu', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vqrecv)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vqrecv),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vqsv, 'vv', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wqrecw)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wqrecw),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wqsw, 'ww', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate rtauij.
!_______________________________________________________________________

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         !VARIABLE FOR MOISTURE


         CALL rmij( rsfs_qmij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+1),        &
                rsfs_rec(ims,kms,jms,P_urec),       &
                rsfs_prec(ims,kms,jms,P_qrecu),       &
                rsfs_prec(ims,kms,jms,P_uqrecu), 'u', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

         CALL rmij( rsfs_qmij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+2),        &
                rsfs_rec(ims,kms,jms,P_vrec),       &
                rsfs_prec(ims,kms,jms,P_qrecv),       &
                rsfs_prec(ims,kms,jms,P_vqrecv), 'v', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

         CALL rmij( rsfs_qmij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+3),        &
                rsfs_rec(ims,kms,jms,P_wrec),       &
                rsfs_prec(ims,kms,jms,P_qrecw),       &
                rsfs_prec(ims,kms,jms,P_wqrecw), 'h', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for rtauij.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_QMIJ.inc"
#      include "PERIOD_EM_RSFS_QMIJ.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------



         CALL set_physical_bc3d_m( rsfs_qmij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+1), 'u', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

         CALL set_physical_bc3d_m( rsfs_qmij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+2), 'v', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

         CALL set_physical_bc3d_m( rsfs_qmij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+3), 'w', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

      ENDDO !-----------------------------------------------------
      !$OMP END PARALLEL DO
   enddo  ! do nq = 1, n_moist
ENDIF ! if (n_moist >= 1)





IF ( n_scalar .ge. PARAM_FIRST_SCALAR ) THEN

   do nq = PARAM_FIRST_SCALAR, n_scalar
      ! use qvapor as a dummy array to store the moist variable
      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               qvapor(i,k,j) = scalar(i,k,j,nq)
            ENDDO
         ENDDO
      ENDDO

      ! Calculate zeroth term of reconstruction at q points.

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rec_0_qmoist( rsfs_rec(ims,kms,jms,P_qrec),       &
              rcoef,                                  &
              qvapor,                                 &
              config_flags%rsfs_opt,                  &
              ids, ide, jds, jde, kds, kde,           &
              ims, ime, jms, jme, kms, kme,           &
              ips, ipe, jps, jpe, kps, kpe,           &
              grid%i_start(ij), grid%i_end(ij),       &
              grid%j_start(ij), grid%j_end(ij),       &
              k_start, k_end                          )

      ENDDO !-----------------------------------------------------
      !$OMP END PARALLEL DO

      !    CALL writecheckm( grid%em_u_2, ims, ime, jms, jme, kms, kme, &
      !                      config_flags%klevs, config_flags%kleve )
      !    WRITE(*,*)'em_u_2'

      !/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

      DO recn = 1, config_flags%rsfs_opt

         !/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

         IF (recn .EQ. 1 ) THEN !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------

               CALL qmoist_s(  grid%qss,                                        &
                    qvapor,                                       &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    2,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )

            ENDDO !---------------------------------------------------------
            !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------

               CALL set_physical_bc2d_m( grid%qss, 't', config_flags,       &
                    ids, ide, jds, jde,               &
                    ims, ime, jms, jme,               &
                    ips, ipe, jps, jpe,               &
                    grid%i_start(ij), grid%i_end(ij), &
                    grid%j_start(ij), grid%j_end(ij)  )


            ENDDO !---------------------------------------------------------
            !$OMP END PARALLEL DO


            IF (config_flags%nested) THEN !use dummy variables to update b.c.

               DO i = ims, ime
                  DO k = kms, kme
                     DO j = jms, jme
                        rsfs_filt(i,k,j,P_temp)=qvapor(i,k,j)
                     ENDDO
                  ENDDO
               ENDDO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_FILT.inc"
#      include "PERIOD_EM_RSFS_FILT.inc"
#endif

               !$OMP PARALLEL DO   &
               !$OMP PRIVATE ( ij )
               DO ij = 1 , grid%num_tiles !------------------------------------

                  CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_temp), 't', &
                       config_flags,                        &
                       ids, ide, jds, jde, kds, kde,        &
                       ims, ime, jms, jme, kms, kme,        &
                       ips, ipe, jps, jpe, kps, kpe,        &
                       grid%i_start(ij), grid%i_end(ij),    &
                       grid%j_start(ij), grid%j_end(ij),    &
                       k_start, k_end                       )

               ENDDO !---------------------------------------------------------
               !$OMP END PARALLEL DO


               !$OMP PARALLEL DO   &
               !$OMP PRIVATE ( ij )
               DO ij = 1 , grid%num_tiles !------------------------------------

                  CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_qfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),     &
                          grid%qss, 'cc',  iexp,             &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

               ENDDO !---------------------------------------------------------
               !$OMP END PARALLEL DO


            ELSE !not nested


               !$OMP PARALLEL DO   &
               !$OMP PRIVATE ( ij )
               DO ij = 1 , grid%num_tiles !------------------------------------

                  CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_qfilt),   &
                          qvapor,                         &
                          grid%qss, 'cc',  iexp,             &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

               ENDDO !---------------------------------------------------------
               !$OMP END PARALLEL DO

            ENDIF

         ELSE !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_FILT.inc"
#      include "PERIOD_EM_RSFS_FILT.inc"
#endif

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------

               CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_qfilt), 't', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

            ENDDO !---------------------------------------------------------
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------

               CALL qmoist_s(  grid%qss,                                        &
                    rsfs_filt(ims,kms,jms,P_qfilt),                 &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    2,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )

            ENDDO !---------------------------------------------------------
            !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------------

               CALL set_physical_bc2d_m( grid%qss, 't', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

            ENDDO !-------------------------------------------------
            !$OMP END PARALLEL DO

            DO i = ims, ime
               DO k = kms, kme
                  DO j = jms, jme
                     rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_qfilt)
                  ENDDO
               ENDDO
            ENDDO

            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !----------------------------

               CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_qfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%qss, 'cc', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

            ENDDO !-------------------------------------------------
            !$OMP END PARALLEL DO


         ENDIF ! (recn .EQ. 1 ) ELSE !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles !------------------------------

            CALL rec_qmoist( rsfs_rec(ims,kms,jms,P_qrec),     &
                   rsfs_filt(ims,kms,jms,P_qfilt),   &
                   rcoef(recn),                      &
                   ids, ide, jds, jde, kds, kde,     &
                   ims, ime, jms, jme, kms, kme,     &
                   ips, ipe, jps, jpe, kps, kpe,     &
                   grid%i_start(ij), grid%i_end(ij), &
                   grid%j_start(ij), grid%j_end(ij), &
                   k_start, k_end                    )

         ENDDO !---------------------------------------------------
         !$OMP END PARALLEL DO

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

      ENDDO !recn = 1, config_flags%rsfs_opt

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/


#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_REC.inc"
#      include "PERIOD_EM_RSFS_REC.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------------------

         CALL set_physical_bc3d_m( rsfs_rec(ims,kms,jms,P_qrec), 't', &
                            config_flags,                      &
                            ids, ide, jds, jde, kds, kde,      &
                            ims, ime, jms, jme, kms, kme,      &
                            ips, ipe, jps, jpe, kps, kpe,      &
                            grid%i_start(ij), grid%i_end(ij),  &
                            grid%j_start(ij), grid%j_end(ij),  &
                            k_start, k_end                     )

      ENDDO !-----------------------------------------------------
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL qmoist_s(  grid%qss,                                        &
                    qvapor,                                       &
                    grid%ht, grid%rdx, grid%rdy,                    &
                    grid%cf1, grid%cf2, grid%cf3,                   &
                    2,                      &
                    ids, ide, jds, jde, kds, kde,                   &
                    ims, ime, jms, jme, kms, kme,                   &
                    ips, ipe, jps, jpe, kps, kpe,                   &
                    grid%i_start(ij), grid%i_end(ij),               &
                    grid%j_start(ij), grid%j_end(ij),               &
                    k_start, k_end                                  )

      ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL set_physical_bc2d_m( grid%qss, 't', config_flags,       &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

      ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL proj_quvw( rsfs_prec(ims,kms,jms,P_qrecu),   &
                   rsfs_prec(ims,kms,jms,P_qrecv),   &
                   rsfs_prec(ims,kms,jms,P_qrecw),   &
                   rsfs_prec(ims,kms,jms,P_uqrecu),  &
                   rsfs_prec(ims,kms,jms,P_vqrecv),  &
                   rsfs_prec(ims,kms,jms,P_wqrecw),  &
                   rsfs_rec(ims,kms,jms,P_urec),     &
                   rsfs_rec(ims,kms,jms,P_vrec),     &
                   rsfs_rec(ims,kms,jms,P_wrec),     &
                   rsfs_rec(ims,kms,jms,P_qrec),     &
                   grid%fnm, grid%fnp,               &
                   ids, ide, jds, jde, kds, kde,     &
                   ims, ime, jms, jme, kms, kme,     &
                   ips, ipe, jps, jpe, kps, kpe,     &
                   grid%i_start(ij), grid%i_end(ij), &
                   grid%j_start(ij), grid%j_end(ij), &
                   k_start, k_end                    )


         CALL proj_s_qmoist( grid%qsu, grid%qsv,                      &
                    grid%uqsu, grid%vqsv, grid%wqsw,             &
                    grid%us, grid%vs, grid%ws, grid%qss,          &
                    ids, ide, jds, jde,                          &
                    ims, ime, jms, jme,                          &
                    ips, ipe, jps, jpe,                          &
                    grid%i_start(ij), grid%i_end(ij),            &
                    grid%j_start(ij), grid%j_end(ij)             )

!if (config_flags%nested) then
!
!  print*,'vrec_f before bc update'
!  CALL writecheck_t2d( rsfs_prec(ims,kms,jms,P_vrecf),   &
!                       ims, ime, jms, jme, kms, kme,     &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )
!
!endif

      ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for variables before filtering.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_PREC.inc"
#      include "PERIOD_EM_RSFS_PREC.inc"
#      include "HALO_EM_RSFS_PUVWS.inc"
#      include "PERIOD_EM_RSFS_PUVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------


!VARIABLES FOR MOISTURE

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_qrecu), 'u',  &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_qrecv), 'v',  &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_qrecw), 'w',  &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uqrecu), 'u', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vqrecv), 'v', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

         CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wqrecw), 'w', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

!2D VARIABLES FOR MOISTURE

         CALL set_physical_bc2d_m( grid%qsu, 'u', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

         CALL set_physical_bc2d_m( grid%qsv, 'v', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

         CALL set_physical_bc2d_m( grid%uqsu, 'u', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

         CALL set_physical_bc2d_m( grid%vqsv, 'v', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

         CALL set_physical_bc2d_m( grid%wqsw, 'w', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

      ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO



      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_qrecu)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_qrecu),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%qsu, 'uu', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_qrecv)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_qrecv),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%qsv, 'vv', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_qrecw)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_qrecw),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%qss, 'ww', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uqrecu)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uqrecu),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uqsu, 'uu', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vqrecv)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vqrecv),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vqsv, 'vv', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
         DO k = kms, kme
            DO j = jms, jme
               rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wqrecw)
            ENDDO
         ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wqrecw),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wqsw, 'ww', iexp,                  &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate rtauij.
!_______________________________________________________________________

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------

         !VARIABLE FOR MOISTURE


         CALL rmij( rsfs_qsij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+1),        &
                rsfs_rec(ims,kms,jms,P_urec),       &
                rsfs_prec(ims,kms,jms,P_qrecu),       &
                rsfs_prec(ims,kms,jms,P_uqrecu), 'u', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

         CALL rmij( rsfs_qsij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+2),        &
                rsfs_rec(ims,kms,jms,P_vrec),       &
                rsfs_prec(ims,kms,jms,P_qrecv),       &
                rsfs_prec(ims,kms,jms,P_vqrecv), 'v', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

         CALL rmij( rsfs_qsij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+3),        &
                rsfs_rec(ims,kms,jms,P_wrec),       &
                rsfs_prec(ims,kms,jms,P_qrecw),       &
                rsfs_prec(ims,kms,jms,P_wqrecw), 'h', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for rtauij.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_QSIJ.inc"
#      include "PERIOD_EM_RSFS_QSIJ.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !--------------------------------



         CALL set_physical_bc3d_m( rsfs_qsij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+1), 'u', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

         CALL set_physical_bc3d_m( rsfs_qsij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+2), 'v', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

         CALL set_physical_bc3d_m( rsfs_qsij(ims,kms,jms,(nq-PARAM_FIRST_SCALAR)*3+3), 'w', &
                             config_flags,                       &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

      ENDDO !-----------------------------------------------------
      !$OMP END PARALLEL DO
   enddo  ! do nq = 1, n_scalar
ENDIF ! if (n_scalar >= 1)



ENDIF !rsfs_opt .LE. 6

!+
!---------------------------------------------------------------




IF ( config_flags%dwl_opt .GE. 1 ) THEN

!=======================================================================
!
! START DWL
!
!=======================================================================
!
! Compute the Dynamic Wong-Lilly eddy viscosity coefficient (DWL_K):
!
! DWL_K = {(L_ij-H_ij)<S_ij>}/(b*{<S_ij><S_ij>})
!
! {} denotes local averaging, <> denotes a test filter, b is a constant
! that depends on the ratio of the explicit grid to the test-grid, S_ij
! is the deformation tensor and L_ij and H_ij will be discussed below.
! All terms on the RHS are projected to the cell center before
! computation of DWL_K.
!
! The dynamic SGS model is based on the use of two grids, one being
! the "explicit" grid, given by dx, dy and dz, and another, the
! so-called "test grid" that is coarser.  Here the test grid has a
! grid spacing twice that of the explicit grid.
!
! L_ij is straightforwardly computed using the velocities on the
! explicit grid at the current time level (em_u_2, em_v_2, em_w_2)
! and the test filter. The velocities are first projected to the
! cell center (uc, vc, wc), various products of the velocities are
! computed (uuc, uvc, etc.) and then they are passed through the test
! filter (tfuc, tfuuc, tfvc, etc.) This gives
! L_ij=<uc_i*uc_j>-<uc_i>*<uc_j>. E.g., L_11= tfuuc - tfu*tfuc
! where uc_i is the velocity at the
!
! H_ij is the difference between the RSFS stresses computed on the
! test grid (trtau_ij), and the RSFS stresses computed on the
! explicit grid (rtau_ij) then passed through test filter (tfrtau_ij):
! H_ij = tfrtau_ij - trtau_ij
! If there is no reconstruction (taurec .EQ. 0) then H_ij = 0.
!
! Otherwise the RSFS stresses must be computed on the coarse grid as
! well. This involves first projecting the velocities on the explicit
! grid at the current time-level onto the test grid, which is achieved
! by passing em_u_2, em_v_2, em_w_2 through the test filter. Next we
! follow the same steps as detailed in the text describing the
! reconstruction loop above, however here every application of the
! explicit filter is followed by application of the test filter. This
! means that halos must be updated between the two filtering steps.
! Once the reconstruction is complete we again have urec, vrec and
! wrec. These are again projected and filtered, and again the test
! filter is applied after the explicit filter. Finally the RSFS stresses
! on the test grid, trtau_ij, are computed. The first term in H_ij,
! tfrtau_ij, are obtained by passing rtauij through the test filter.
!
! We now have all of the terms needed to compute DWL_K. We compute the
! numerator and denominator separately and filter each once with the
! explicit filter and twice with the test filter. We then combine the
! numerator and denominator and we are done.
!
!=======================================================================

!_______________________________________________________________________

IF ( config_flags%dwl_opt .EQ. 1 ) THEN

!_______________________________________________________________________
!
! dwl_opt 1 computes the RSFS stresses on the test grid (below),
! and also filtersvthe RSFS stresses on the explicit grid
! (computed already if rsfs_opt .lt. 6).
!
! Warning! if dwl_opt = 1, rsfs_opt cannot be 6
!_______________________________________________________________________

!_______________________________________________________________________
!
! Compute the RSFS stresses on the test grid.
!_______________________________________________________________________

!_______________________________________________________________________
!
! First compute kinematic surface values of w_2 using u_2 and v_2.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------

                CALL uvw_s( grid%us, grid%vs, grid%ws,        &
                            grid%u_2, grid%v_2,               &
                            grid%ht, grid%rdx, grid%rdy,      &
                            grid%cf1, grid%cf2, grid%cf3,     &
                            config_flags%rsfsvsbc_opt,        &
                            ids, ide, jds, jde, kds, kde,     &
                            ims, ime, jms, jme, kms, kme,     &
                            ips, ipe, jps, jpe, kps, kpe,     &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij), &
                            k_start, k_end                    )

            ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for us, vs, ws
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc2d_m( grid%us, 'u', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%vs, 'v', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%ws, 'w', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

      ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Compute tfu, tfv, tfw by test-filtering u_2, v_2, w_2.
!_______________________________________________________________________

! ASW debug
!  print*,'compute tfu'
!  CALL writecheckm_drm( tfuc_drm(ims,kms,jms,P_tfuc),       &
!                        ims, ime, jms, jme, kms, kme,     &
!                        grid%i_start(ij), grid%i_end(ij), &
!                        grid%j_start(ij), grid%j_end(ij), &
!                        k_start, k_end                    )


 IF (config_flags%nested) THEN !use dummy variables to update b.c. for u_2, v_2 and w_2

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         sijc_drm(ims,kms,jms,P_s11c)=grid%u_2(i,k,j)
         sijc_drm(ims,kms,jms,P_s22c)=grid%v_2(i,k,j)
         sijc_drm(ims,kms,jms,P_s33c)=grid%w_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO

#ifdef DM_PARALLEL
#      include "HALO_EM_SIJC_DRM.inc"
#      include "PERIOD_EM_SIJC_DRM.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc3d_m( sijc_drm(ims,kms,jms,P_s11c), 'u', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( sijc_drm(ims,kms,jms,P_s22c), 'v', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )


        CALL set_physical_bc3d_m( sijc_drm(ims,kms,jms,P_s33c), 'w', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )


      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL rsfs_filter( tfu_drm(ims,kms,jms,P_tfu),       &
                          sijc_drm(ims,kms,jms,P_s11c),     &
                        !  grid%us, 'uu', iexp,              &     seems to be a bug. XS
                          grid%us, 'uu', itst,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( tfu_drm(ims,kms,jms,P_tfv),       &
                          sijc_drm(ims,kms,jms,P_s22c),     &
                        !  grid%vs, 'vv', iexp,              &  seems to be a bug. XS
                          grid%vs, 'vv', itst,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( tfu_drm(ims,kms,jms,P_tfw),       &
                          sijc_drm(ims,kms,jms,P_s33c),     &
                        !  grid%ws, 'ww', iexp,              &  seems to be a bug. XS
                          grid%ws, 'ww', itst,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )


      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO


 ELSE !not nested

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !------------------------------

              CALL rsfs_filter( tfu_drm(ims,kms,jms,P_tfu),       &
                                grid%u_2, grid%us, 'uu', itst,    &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start, k_end                    )

              CALL rsfs_filter( tfu_drm(ims,kms,jms,P_tfv),       &
                                grid%v_2, grid%vs, 'vv', itst,    &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start, k_end                    )

              CALL rsfs_filter( tfu_drm(ims,kms,jms,P_tfw),       &
                                grid%w_2, grid%ws, 'ww', itst,    &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

 ENDIF

! ASW debug
!   print*,'end of compute tfu'
!   CALL writecheckm_drm( tfuc_drm(ims,kms,jms,P_tfuc),       &
!                         ims, ime, jms, jme, kms, kme,     &
!                         grid%i_start(ij), grid%i_end(ij), &
!                         grid%j_start(ij), grid%j_end(ij), &
!                         k_start, k_end                    )


!_______________________________________________________________________
!
! Now we have the test filtered u,v, and w; tfu. tfv. tfw. tfw(*,kts,*)
! is not correct, but it is not needed.
!
! Now compute T_RSFS, the RSFS stresses computed from the test filtered variables.
! Computation should be identical to that for t_ij, with different arguments.
!_______________________________________________________________________

IF ( config_flags%rsfs_opt .GE. 0 ) THEN ! ASW only enter this for DRM, not for DWL
!print*,'should not be here'


IF ( config_flags%rsfs_opt .EQ. 6 ) THEN
print*,' must have rsfs_opt < 6 to use dwl_opt 1'
STOP
ENDIF

!_______________________________________________________________________
!
! Calculate zeroth term of reconstruction at u, v, w points.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rec_0( rsfs_rec(ims,kms,jms,P_urec),     &
                rsfs_rec(ims,kms,jms,P_vrec),     &
                rsfs_rec(ims,kms,jms,P_wrec),     &
                rcoef,                            &
                tfu_drm(ims,kms,jms,P_tfu),       &
                tfu_drm(ims,kms,jms,P_tfv),       &
                tfu_drm(ims,kms,jms,P_tfw),       &
                config_flags%rsfs_opt,            &
                ids, ide, jds, jde, kds, kde,     &
                ims, ime, jms, jme, kms, kme,     &
                ips, ipe, jps, jpe, kps, kpe,     &
                grid%i_start(ij), grid%i_end(ij), &
                grid%j_start(ij), grid%j_end(ij), &
                k_start, k_end                    )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! wrec(*,kts,*) is not correct, but it is not needed.
!_______________________________________________________________________

!_______________________________________________________________________
!
! Loop up to level of reconstruction if reclevel > 0
!_______________________________________________________________________

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

  DO recn = 1, config_flags%rsfs_opt

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

    IF (recn .EQ. 1 ) THEN !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!_______________________________________________________________________
!
! First time through the reconstruction loop we filter tfu, tfv,
! tfw. Update halos amd calculate corresponding surface
! values, us, vs, ws.
!_______________________________________________________________________


#ifdef DM_PARALLEL
#      include "HALO_EM_TFU_DRM.inc"
#      include "PERIOD_EM_TFU_DRM.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc3d_m( tfu_drm(ims,kms,jms,P_tfu),         &
                                'u', config_flags,                  &
                                ids, ide, jds, jde, kds, kde,       &
                                ims, ime, jms, jme, kms, kme,       &
                                ips, ipe, jps, jpe, kps, kpe,       &
                                grid%i_start(ij), grid%i_end(ij),   &
                                grid%j_start(ij), grid%j_end(ij),   &
                                k_start, k_end                      )

        CALL set_physical_bc3d_m( tfu_drm(ims,kms,jms,P_tfv),         &
                                'v', config_flags,                  &
                                ids, ide, jds, jde, kds, kde,       &
                                ims, ime, jms, jme, kms, kme,       &
                                ips, ipe, jps, jpe, kps, kpe,       &
                                grid%i_start(ij), grid%i_end(ij),   &
                                grid%j_start(ij), grid%j_end(ij),   &
                                k_start, k_end                      )

        CALL set_physical_bc3d_m( tfu_drm(ims,kms,jms,P_tfw),         &
                                'w', config_flags,                  &
                                ids, ide, jds, jde, kds, kde,       &
                                ims, ime, jms, jme, kms, kme,       &
                                ips, ipe, jps, jpe, kps, kpe,       &
                                grid%i_start(ij), grid%i_end(ij),   &
                                grid%j_start(ij), grid%j_end(ij),   &
                                k_start, k_end                      )

      ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO
!_______________________________________________________________________
!
! Compute kinematic surface values of w_2 using tfu and tfv.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------

                CALL uvw_s( grid%us, grid%vs, grid%ws,        &
                            tfu_drm(ims,kms,jms,P_tfu),       &
                            tfu_drm(ims,kms,jms,P_tfv),       &
                            grid%ht, grid%rdx, grid%rdy,      &
                            grid%cf1, grid%cf2, grid%cf3,     &
                            config_flags%rsfsvsbc_opt,        &
                            ids, ide, jds, jde, kds, kde,     &
                            ims, ime, jms, jme, kms, kme,     &
                            ips, ipe, jps, jpe, kps, kpe,     &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij), &
                            k_start, k_end                    )

            ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for us, vs, ws before filtering
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc2d_m( grid%us, 'u', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%vs, 'v', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%ws, 'w', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

      ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! (explicit) Filter tfu, tfv, tfw to calculate ufilt, vfilt, wfilt.
!_______________________________________________________________________

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_ufilt),   &
                          tfu_drm(ims,kms,jms,P_tfu),       &
                          grid%us, 'uu', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_vfilt),   &
                          tfu_drm(ims,kms,jms,P_tfv),       &
                          grid%vs, 'vv', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_wfilt),   &
                          tfu_drm(ims,kms,jms,P_tfw),       &
                          grid%ws, 'ww', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

       ENDDO !---------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
!  wfilt(*,kts,*) is not correct, but it is not needed.
!_______________________________________________________________________


                ELSE ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!_______________________________________________________________________
!
! Update halos for ufilt, vfilt, wfilt before (explicit) filtering
! ufilt, vfilt, wfilt again
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_FILT.inc"
#      include "PERIOD_EM_RSFS_FILT.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_ufilt), 'u', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_vfilt), 'v', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_wfilt), 'w', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Compute us, vs, ws using ufilt and vfilt.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------

                CALL uvw_s( grid%us, grid%vs, grid%ws,        &
                            rsfs_filt(ims,kms,jms,P_ufilt),   &
                            rsfs_filt(ims,kms,jms,P_vfilt),   &
                            grid%ht, grid%rdx, grid%rdy,      &
                            grid%cf1, grid%cf2, grid%cf3,     &
                            config_flags%rsfsvsbc_opt,        &
                            ids, ide, jds, jde, kds, kde,     &
                            ims, ime, jms, jme, kms, kme,     &
                            ips, ipe, jps, jpe, kps, kpe,     &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij), &
                            k_start, k_end                    )

            ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for us, vs, ws  before (explicit) filtering ufilt, vfilt,
! wfilt again
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc2d_m( grid%us, 'u', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%vs, 'v', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%ws, 'w', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

      ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! (explicit) Filter ufilt, vfilt, wfilt to calculate ufilt, vfilt, wfilt
!_______________________________________________________________________


      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_ufilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_ufilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%us, 'uu', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

       ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_vfilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_vfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%vs, 'vv', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

       ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_wfilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------


        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_wfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%ws, 'ww', iexp,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

! print*,'ufilt'
! CALL writecheckm_drm( rsfs_filt(ims,kms,jms,P_ufilt),     &
!                       ims, ime, jms, jme, kms, kme,     &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )
! print*,'vfilt'
! CALL writecheckm_drm( rsfs_filt(ims,kms,jms,P_vfilt),   &
!                       ims, ime, jms, jme, kms, kme ,    &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )

       ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

  ENDIF ! (recn .eq. 1 ) ELSE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!_______________________________________________________________________
!
! Now we apply the test filter on top of the explicit filter.
!_______________________________________________________________________
!_______________________________________________________________________
!
! Update halos for ufilt, vfilt, wfilt before computing us, vs, ws
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_FILT.inc"
#      include "PERIOD_EM_RSFS_FILT.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_ufilt), 'u', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_vfilt), 'v', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

        CALL set_physical_bc3d_m( rsfs_filt(ims,kms,jms,P_wfilt), 'w', &
                                config_flags,                        &
                                ids, ide, jds, jde, kds, kde,        &
                                ims, ime, jms, jme, kms, kme,        &
                                ips, ipe, jps, jpe, kps, kpe,        &
                                grid%i_start(ij), grid%i_end(ij),    &
                                grid%j_start(ij), grid%j_end(ij),    &
                                k_start, k_end                       )

      ENDDO !---------------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Compute us, vs, ws using ufilt and vfilt.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------

                CALL uvw_s( grid%us, grid%vs, grid%ws,        &
                            rsfs_filt(ims,kms,jms,P_ufilt),   &
                            rsfs_filt(ims,kms,jms,P_vfilt),   &
                            grid%ht, grid%rdx, grid%rdy,      &
                            grid%cf1, grid%cf2, grid%cf3,     &
                            config_flags%rsfsvsbc_opt,        &
                            ids, ide, jds, jde, kds, kde,     &
                            ims, ime, jms, jme, kms, kme,     &
                            ips, ipe, jps, jpe, kps, kpe,     &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij), &
                            k_start, k_end                    )

            ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for us, vs, ws before filtering
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !------------------------------------

        CALL set_physical_bc2d_m( grid%us, 'u', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%vs, 'v', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

        CALL set_physical_bc2d_m( grid%ws, 'w', config_flags,       &
                                ids, ide, jds, jde,               &
                                ims, ime, jms, jme,               &
                                ips, ipe, jps, jpe,               &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij)  )

      ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Test-filter tufilt, tvfilt, twfilt
!_______________________________________________________________________


      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_ufilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_ufilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%us, 'uu', itst,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

       ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_vfilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_vfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%vs, 'vv', itst,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

       ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

      DO i = ims, ime
      DO k = kms, kme
      DO j = jms, jme
         rsfs_filt(i,k,j,P_temp)=rsfs_filt(i,k,j,P_wfilt)
      ENDDO
      ENDDO
      ENDDO

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles !----------------------------

        CALL rsfs_filter( rsfs_filt(ims,kms,jms,P_wfilt),   &
                          rsfs_filt(ims,kms,jms,P_temp),    &
                          grid%ws, 'ww', itst,              &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

       ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate urec, vrec, wrec.
!_______________________________________________________________________


    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij )
    DO ij = 1 , grid%num_tiles !------------------------------

      CALL rec( rsfs_rec(ims,kms,jms,P_urec), &
                rsfs_rec(ims,kms,jms,P_vrec),     &
                rsfs_rec(ims,kms,jms,P_wrec),     &
                rsfs_filt(ims,kms,jms,P_ufilt),   &
                rsfs_filt(ims,kms,jms,P_vfilt),   &
                rsfs_filt(ims,kms,jms,P_wfilt),   &
                rcoef(recn),                      &
                ids, ide, jds, jde, kds, kde,     &
                ims, ime, jms, jme, kms, kme,     &
                ips, ipe, jps, jpe, kps, kpe,     &
                grid%i_start(ij), grid%i_end(ij), &
                grid%j_start(ij), grid%j_end(ij), &
                k_start, k_end                    )

    ENDDO !---------------------------------------------------
    !$OMP END PARALLEL DO

!_______________________________________________________________________
!
!  wrec(*,kts,*) is not correct, but it is not needed.
!_______________________________________________________________________


!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

              ENDDO !recn = 1, config_flags%reclevel

!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/


!_______________________________________________________________________
!
! Now that reconstruction is over we will need to project urec, vrec,
! wrec to vorticity points.
!_______________________________________________________________________

!_______________________________________________________________________
!
! First update halos for urec, vrec, wrec, to compute us, vs, ws
!_______________________________________________________________________


#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_REC.inc"
#      include "PERIOD_EM_RSFS_REC.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !------------------------------

     CALL set_physical_bc3d_m( rsfs_rec(ims,kms,jms,P_urec),       &
                             'u', config_flags,                  &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     CALL set_physical_bc3d_m( rsfs_rec(ims,kms,jms,P_vrec),       &
                             'v', config_flags,                  &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     ENDDO !---------------------------------------------
     !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Compute us, vs, ws consistent with urec and vrec.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------

                CALL uvw_s( grid%us, grid%vs, grid%ws,        &
                            rsfs_rec(ims,kms,jms,P_urec),     &
                            rsfs_rec(ims,kms,jms,P_vrec),     &
                            grid%ht, grid%rdx, grid%rdy,      &
                            grid%cf1, grid%cf2, grid%cf3,     &
                            config_flags%rsfsvsbc_opt,        &
                            ids, ide, jds, jde, kds, kde,     &
                            ims, ime, jms, jme, kms, kme,     &
                            ips, ipe, jps, jpe, kps, kpe,     &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij), &
                            k_start, k_end                    )

            ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! replace w1 in wrec with ws before projecting to cell center
!_______________________________________________________________________


  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------

                CALL update_w1( rsfs_rec(ims,kms,jms,P_wrec),     &
                                grid%ws,                          &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start, k_end                    )

   ENDDO !-----------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
!  Update b.c. on surface variables and wrec. wrec(*,kts,*) is now correct.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_REC.inc"
#      include "PERIOD_EM_RSFS_REC.inc"
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !------------------------------

     CALL set_physical_bc3d_m( rsfs_rec(ims,kms,jms,P_wrec),       &
                             'w', config_flags,                  &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             grid%i_start(ij), grid%i_end(ij),   &
                             grid%j_start(ij), grid%j_end(ij),   &
                             k_start, k_end                      )

     CALL set_physical_bc2d_m( grid%us, 'u', config_flags,       &
                             ids, ide, jds, jde,               &
                             ims, ime, jms, jme,               &
                             ips, ipe, jps, jpe,               &
                             grid%i_start(ij), grid%i_end(ij), &
                             grid%j_start(ij), grid%j_end(ij)  )

     CALL set_physical_bc2d_m( grid%vs, 'v', config_flags,       &
                             ids, ide, jds, jde,               &
                             ims, ime, jms, jme,               &
                             ips, ipe, jps, jpe,               &
                             grid%i_start(ij), grid%i_end(ij), &
                             grid%j_start(ij), grid%j_end(ij)  )

     CALL set_physical_bc2d_m( grid%ws, 'w', config_flags,       &
                             ids, ide, jds, jde,               &
                             ims, ime, jms, jme,               &
                             ips, ipe, jps, jpe,               &
                             grid%i_start(ij), grid%i_end(ij), &
                             grid%j_start(ij), grid%j_end(ij)  )

     ENDDO !---------------------------------------------
     !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Now we have the RSFS velocities, urec, vrec, wrec at u, v, w points,
! as well as their halos and surface values.
! Next we project urec, vrec, wrec and their surface values.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL proj_c( rsfs_prec(ims,kms,jms,P_urecc),    &
                 rsfs_prec(ims,kms,jms,P_vrecc),    &
                 rsfs_prec(ims,kms,jms,P_wrecc),    &
                 rsfs_prec(ims,kms,jms,P_uurecc),   &
                 rsfs_prec(ims,kms,jms,P_vvrecc),   &
                 rsfs_prec(ims,kms,jms,P_wwrecc),   &
                 rsfs_rec(ims,kms,jms,P_urec),      &
                 rsfs_rec(ims,kms,jms,P_vrec),      &
                 rsfs_rec(ims,kms,jms,P_wrec), 'u', &
                 ids, ide, jds, jde, kds, kde,      &
                 ims, ime, jms, jme, kms, kme,      &
                 ips, ipe, jps, jpe, kps, kpe,      &
                 grid%i_start(ij), grid%i_end(ij),  &
                 grid%j_start(ij), grid%j_end(ij),  &
                 k_start, k_end                     )

    CALL proj_d( rsfs_prec(ims,kms,jms,P_urecd),   &
                 rsfs_prec(ims,kms,jms,P_vrecd),   &
                 rsfs_prec(ims,kms,jms,P_uvrecd),  &
                 rsfs_rec(ims,kms,jms,P_urec),     &
                 rsfs_rec(ims,kms,jms,P_vrec),     &
                 ids, ide, jds, jde, kds, kde,     &
                 ims, ime, jms, jme, kms, kme,     &
                 ips, ipe, jps, jpe, kps, kpe,     &
                 grid%i_start(ij), grid%i_end(ij), &
                 grid%j_start(ij), grid%j_end(ij), &
                 k_start, k_end                    )

    CALL proj_e( rsfs_prec(ims,kms,jms,P_urece),   &
                 rsfs_prec(ims,kms,jms,P_wrece),   &
                 rsfs_prec(ims,kms,jms,P_uwrece),  &
                 rsfs_rec(ims,kms,jms,P_urec),     &
                 rsfs_rec(ims,kms,jms,P_wrec),     &
                 grid%fnm, grid%fnp,               &
                 ids, ide, jds, jde, kds, kde,     &
                 ims, ime, jms, jme, kms, kme,     &
                 ips, ipe, jps, jpe, kps, kpe,     &
                 grid%i_start(ij), grid%i_end(ij), &
                 grid%j_start(ij), grid%j_end(ij), &
                 k_start, k_end                    )

    CALL proj_f( rsfs_prec(ims,kms,jms,P_vrecf),   &
                 rsfs_prec(ims,kms,jms,P_wrecf),   &
                 rsfs_prec(ims,kms,jms,P_vwrecf),  &
                 rsfs_rec(ims,kms,jms,P_vrec),     &
                 rsfs_rec(ims,kms,jms,P_wrec),     &
                 grid%fnm, grid%fnp,               &
                 ids, ide, jds, jde, kds, kde,     &
                 ims, ime, jms, jme, kms, kme,     &
                 ips, ipe, jps, jpe, kps, kpe,     &
                 grid%i_start(ij), grid%i_end(ij), &
                 grid%j_start(ij), grid%j_end(ij), &
                 k_start, k_end                    )

    CALL proj_s( grid%usc, grid%vsc,               &
                 grid%uusc, grid%vvsc, grid%wws,   &
                 grid%usd, grid%vsd, grid%uvsd,    &
                 grid%wse, grid%uwse,              &
                 grid%wsf, grid%vwsf,              &
                 grid%us, grid%vs, grid%ws,        &
                 ids, ide, jds, jde,               &
                 ims, ime, jms, jme,               &
                 ips, ipe, jps, jpe,               &
                 grid%i_start(ij), grid%i_end(ij), &
                 grid%j_start(ij), grid%j_end(ij)  )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for variables before filtering.
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_PREC.inc"
#      include "PERIOD_EM_RSFS_PREC.inc"
#      include "HALO_EM_RSFS_PUVWS.inc"
#      include "PERIOD_EM_RSFS_PUVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_urecc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uurecc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vrecc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vvrecc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wrecc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wwrecc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_urecd), 'd', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vrecd), 'd', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uvrecd), 'd', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_urece), 'e', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wrece), 'e', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uwrece), 'e', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vrecf), 'f', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wrecf), 'f', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vwrecf), 'f', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc2d_m( grid%usc, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uusc, 't', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vsc, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vvsc, 't', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wws, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%usd, 'd', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vsd, 'd', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uvsd, 'd', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wse, 'e', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uwse, 'e', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wsf, 'f', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vwsf, 'f', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Explicit filter urecc, urecc*urecc, etc.
!_______________________________________________________________________

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_urecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_urecc),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%usc, 'cc', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uurecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uurecc),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uusc, 'cc', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vrecc),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vsc, 'cc', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vvrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vvrecc),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vvsc, 'cc', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wrecc),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%ws, 'cw', iexp,              &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wwrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wwrecc),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wws, 'cw', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_urecd)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_urecd),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%usd, 'dd', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vrecd)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

   CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vrecd),   &
                     rsfs_filt(ims,kms,jms,P_temp),    &
                     grid%vsd, 'dd', iexp,             &
                     ids, ide, jds, jde, kds, kde,     &
                     ims, ime, jms, jme, kms, kme,     &
                     ips, ipe, jps, jpe, kps, kpe,     &
                     grid%i_start(ij), grid%i_end(ij), &
                     grid%j_start(ij), grid%j_end(ij), &
                     k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uvrecd)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uvrecd),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uvsd, 'dd', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_urece)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_urece),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%us, 'ee', iexp,              &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wrece)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wrece),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wse, 'ew', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uwrece)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uwrece),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uwse, 'ew', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vrecf)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vrecf),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vs, 'ff', iexp,              &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wrecf)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wrecf),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wsf, 'fw', iexp,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vwrecf)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vwrecf),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vwsf, 'fw', iexp,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Done explicit filtering. Now test filter. First find surface us, vs
! relative to filtered urecc, vrecc; update halos, then project those
! surface values to c, d, e, f points. Then assign wsc to wrecc. Then
! update halos for urecc, uurecc, etc. Then test filter.

!_______________________________________________________________________
!
! Compute surface values relative to urecc, vrecc before filtering again
! urecc, etc. are already projected to stress nodes. Therefore update
! surface values at stress nodes as well. No need to project surface
! variables. But do need to update halos.
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL uvw_s_c( grid%usc, grid%vsc, grid%ws,      &
                  grid%uusc, grid%vvsc, grid%wws,   &
                  grid%usd, grid%vsd, grid%uvsd,    &
                  grid%wse, grid%uwse,              &
                  grid%wsf, grid%vwsf,              &
                  rsfs_prec(ims,kms,jms,P_urecc),   &
                  rsfs_prec(ims,kms,jms,P_vrecc),   &
                  rsfs_prec(ims,kms,jms,P_wrecc),   &
                  rsfs_prec(ims,kms,jms,P_uurecc),  &
                  rsfs_prec(ims,kms,jms,P_vvrecc),  &
                  rsfs_prec(ims,kms,jms,P_wwrecc),  &
                  rsfs_prec(ims,kms,jms,P_urecd),   &
                  rsfs_prec(ims,kms,jms,P_vrecd),   &
                  rsfs_prec(ims,kms,jms,P_uvrecd),  &
                  rsfs_prec(ims,kms,jms,P_wrece),   &
                  rsfs_prec(ims,kms,jms,P_uwrece),  &
                  rsfs_prec(ims,kms,jms,P_wrecf),   &
                  rsfs_prec(ims,kms,jms,P_vwrecf),  &
                  grid%ht, grid%rdx, grid%rdy,      &
                  grid%cf1, grid%cf2, grid%cf3,     &
                  config_flags%rsfsvsbc_opt,        &
                  ids, ide, jds, jde, kds, kde,     &
                  ims, ime, jms, jme, kms, kme,     &
                  ips, ipe, jps, jpe, kps, kpe,     &
                  grid%i_start(ij), grid%i_end(ij), &
                  grid%j_start(ij), grid%j_end(ij), &
                  k_start, k_end                    )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos before filtering again urecc, vrecc, usc, etc.
!_______________________________________________________________________


#ifdef DM_PARALLEL
#      include "HALO_EM_RSFS_PREC.inc"
#      include "PERIOD_EM_RSFS_PREC.inc"
#      include "HALO_EM_RSFS_UVWS.inc"
#      include "PERIOD_EM_RSFS_UVWS.inc"
#      include "HALO_EM_RSFS_PUVWS.inc"
#      include "PERIOD_EM_RSFS_PUVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !----------------------------

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_urecc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uurecc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vrecc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vvrecc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wrecc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wwrecc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_urecd), 'd', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vrecd), 'd', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uvrecd), 'd', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_urece), 'e', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wrece), 'e', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_uwrece), 'e', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vrecf), 'f', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_wrecf), 'f', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( rsfs_prec(ims,kms,jms,P_vwrecf), 'f', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc2d_m( grid%usc, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uusc, 't', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vsc, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vvsc, 't', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%ws, 't', config_flags,       &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wws, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%usd, 'd', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vsd, 'd', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uvsd, 'd', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wse, 'e', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uwse, 'e', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wsf, 'f', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vwsf, 'f', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

  ENDDO !-------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! TEST filter urecc, urecc*urecc, etc.
!_______________________________________________________________________

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_urecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_urecc),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%usc, 'cc', itst,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uurecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uurecc),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uusc, 'cc', itst,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vrecc),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vsc, 'cc', itst,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vvrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vvrecc),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vvsc, 'cc', itst,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wrecc),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%ws, 'cw', itst,              &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wwrecc)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wwrecc),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wws, 'cw', itst,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_urecd)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_urecd),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%usd, 'dd', itst,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vrecd)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

   CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vrecd),   &
                     rsfs_filt(ims,kms,jms,P_temp),    &
                     grid%vsd, 'dd', itst,             &
                     ids, ide, jds, jde, kds, kde,     &
                     ims, ime, jms, jme, kms, kme,     &
                     ips, ipe, jps, jpe, kps, kpe,     &
                     grid%i_start(ij), grid%i_end(ij), &
                     grid%j_start(ij), grid%j_end(ij), &
                     k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uvrecd)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uvrecd),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uvsd, 'dd', itst,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_urece)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_urece),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%us, 'ee', itst,              &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wrece)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wrece),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wse, 'ew', itst,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_uwrece)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_uwrece),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%uwse, 'ew', itst,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vrecf)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vrecf),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vs, 'ff', itst,              &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_wrecf)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_wrecf),   &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%wsf, 'fw', itst,             &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

  DO i = ims, ime
  DO k = kms, kme
  DO j = jms, jme
     rsfs_filt(i,k,j,P_temp)=rsfs_prec(i,k,j,P_vwrecf)
  ENDDO
  ENDDO
  ENDDO

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rsfs_filter( rsfs_prec(ims,kms,jms,P_vwrecf),  &
                      rsfs_filt(ims,kms,jms,P_temp),    &
                      grid%vwsf, 'fw', itst,            &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      ips, ipe, jps, jpe, kps, kpe,     &
                      grid%i_start(ij), grid%i_end(ij), &
                      grid%j_start(ij), grid%j_end(ij), &
                      k_start, k_end                    )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO
!_______________________________________________________________________
!
! Done test filtering.
!_______________________________________________________________________
!

!_______________________________________________________________________
!
! Calculate tgrmij = RSFS stresses on the test grid
!_______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL rmij( tgrmij_drm(ims,kms,jms,P_tgrm11),     &
               rsfs_prec(ims,kms,jms,P_urecc),       &
               rsfs_prec(ims,kms,jms,P_urecc),       &
               rsfs_prec(ims,kms,jms,P_uurecc), 'c', &
               ids, ide, jds, jde, kds, kde,         &
               ims, ime, jms, jme, kms, kme,         &
               ips, ipe, jps, jpe, kps, kpe,         &
               grid%i_start(ij), grid%i_end(ij),     &
               grid%j_start(ij), grid%j_end(ij),     &
               k_start, k_end                        )

     CALL rmij( tgrmij_drm(ims,kms,jms,P_tgrm22),     &
                rsfs_prec(ims,kms,jms,P_vrecc),       &
                rsfs_prec(ims,kms,jms,P_vrecc),       &
                rsfs_prec(ims,kms,jms,P_vvrecc), 'c', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( tgrmij_drm(ims,kms,jms,P_tgrm33),     &
                rsfs_prec(ims,kms,jms,P_wrecc),       &
                rsfs_prec(ims,kms,jms,P_wrecc),       &
                rsfs_prec(ims,kms,jms,P_wwrecc), 'c', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( tgrmij_drm(ims,kms,jms,P_tgrm12),     &
                rsfs_prec(ims,kms,jms,P_urecd),       &
                rsfs_prec(ims,kms,jms,P_vrecd),       &
                rsfs_prec(ims,kms,jms,P_uvrecd), 'd', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( tgrmij_drm(ims,kms,jms,P_tgrm13),     &
                rsfs_prec(ims,kms,jms,P_urece),       &
                rsfs_prec(ims,kms,jms,P_wrece),       &
                rsfs_prec(ims,kms,jms,P_uwrece), 'e', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

     CALL rmij( tgrmij_drm(ims,kms,jms,P_tgrm23),     &
                rsfs_prec(ims,kms,jms,P_vrecf),       &
                rsfs_prec(ims,kms,jms,P_wrecf),       &
                rsfs_prec(ims,kms,jms,P_vwrecf), 'f', &
                ids, ide, jds, jde, kds, kde,         &
                ims, ime, jms, jme, kms, kme,         &
                ips, ipe, jps, jpe, kps, kpe,         &
                grid%i_start(ij), grid%i_end(ij),     &
                grid%j_start(ij), grid%j_end(ij),     &
                k_start, k_end                        )

  ENDDO !---------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Test-filter the explicit RSFS stresses (rtau) -> tfrtau
!_______________________________________________________________________

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( tfrmij_drm(ims,kms,jms,P_tfrm11), &
                          rsfs_rmij(ims,kms,jms,P_rm11),    &
                          grid%us, 'tc', itst,              & !us is placeholder
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( tfrmij_drm(ims,kms,jms,P_tfrm22), &
                          rsfs_rmij(ims,kms,jms,P_rm22),    &
                          grid%us, 'tc', itst,              & !us is placeholder
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( tfrmij_drm(ims,kms,jms,P_tfrm33), &
                          rsfs_rmij(ims,kms,jms,P_rm33),    &
                          grid%us, 'tc', itst,              & !us is placeholder
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( tfrmij_drm(ims,kms,jms,P_tfrm12), &
                          rsfs_rmij(ims,kms,jms,P_rm12),    &
                          grid%us, 'td', itst,              & !us is placeholder
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( tfrmij_drm(ims,kms,jms,P_tfrm13), &
                          rsfs_rmij(ims,kms,jms,P_rm13),    &
                          grid%us, 'te', itst,              & !us is placeholder
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

        CALL rsfs_filter( tfrmij_drm(ims,kms,jms,P_tfrm23), &
                          rsfs_rmij(ims,kms,jms,P_rm23),    &
                          grid%us, 'tf', itst,              & !us is placeholder
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

       ENDDO !---------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate Hij = tfrtauij - trtauij
!_______________________________________________________________________

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !----------------------------

                CALL hij( hij_drm(ims,kms,jms,P_h11),       &
                          hij_drm(ims,kms,jms,P_h22),       &
                          hij_drm(ims,kms,jms,P_h33),       &
                          hij_drm(ims,kms,jms,P_h12),       &
                          hij_drm(ims,kms,jms,P_h13),       &
                          hij_drm(ims,kms,jms,P_h23),       &
                          tgrmij_drm(ims,kms,jms,P_tgrm11), &
                          tgrmij_drm(ims,kms,jms,P_tgrm22), &
                          tgrmij_drm(ims,kms,jms,P_tgrm33), &
                          tgrmij_drm(ims,kms,jms,P_tgrm12), &
                          tgrmij_drm(ims,kms,jms,P_tgrm13), &
                          tgrmij_drm(ims,kms,jms,P_tgrm23), &
                          tfrmij_drm(ims,kms,jms,P_tfrm11), &
                          tfrmij_drm(ims,kms,jms,P_tfrm22), &
                          tfrmij_drm(ims,kms,jms,P_tfrm33), &
                          tfrmij_drm(ims,kms,jms,P_tfrm12), &
                          tfrmij_drm(ims,kms,jms,P_tfrm13), &
                          tfrmij_drm(ims,kms,jms,P_tfrm23), &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

       ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO
!_______________________________________________________________________
!
! Update halos for Hij
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_HIJ_DRM.inc"
#      include "PERIOD_EM_HIJ_DRM.inc"
#endif

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !----------------------------

    CALL set_physical_bc3d_m( hij_drm(ims,kms,jms,P_h11), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( hij_drm(ims,kms,jms,P_h22), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( hij_drm(ims,kms,jms,P_h33), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( hij_drm(ims,kms,jms,P_h12), 'd', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( hij_drm(ims,kms,jms,P_h13), 'e', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( hij_drm(ims,kms,jms,P_h23), 'f', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

! print*,'h13'
! CALL writecheckm_drm( hij_drm(ims,kms,jms,P_h13),       &
!                       ims, ime, jms, jme, kms, kme,     &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )
! print*,'h23'
! CALL writecheckm_drm( hij_drm(ims,kms,jms,P_h23),       &
!                       ims, ime, jms, jme, kms, kme ,    &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )

       ENDDO !-------------------------------------------------
      !$OMP END PARALLEL DO

ENDIF ! ASW make sure we don't enter RSFS for DWL

ENDIF ! config_flags%dwl_opt .EQ. 1

!_______________________________________________________________________
!
! If dwl_opt = 1, then we are done with the RSFS stresses and Hij.
! If dwl_opt = 2, this is where DRM model computation begins.
!
! Now compute Lij.
!
! Project variables to cell center
!_______________________________________________________________________


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------------

              CALL proj_c( uc_drm(ims,kms,jms,P_uc),          &
                           uc_drm(ims,kms,jms,P_vc),          &
                           uc_drm(ims,kms,jms,P_wc),          &
                           uc_drm(ims,kms,jms,P_uuc),         &
                           uc_drm(ims,kms,jms,P_vvc),         &
                           uc_drm(ims,kms,jms,P_wwc),         &
                           grid%u_2, grid%v_2, grid%w_2, 'u', &
                           ids, ide, jds, jde, kds, kde,      &
                           ims, ime, jms, jme, kms, kme,      &
                           ips, ipe, jps, jpe, kps, kpe,      &
                           grid%i_start(ij), grid%i_end(ij),  &
                           grid%j_start(ij), grid%j_end(ij),  &
                           k_start, k_end                     )


              CALL proj_c( sijc_drm(ims,kms,jms,P_s12c),      &
                           sijc_drm(ims,kms,jms,P_s13c),      &
                           sijc_drm(ims,kms,jms,P_s23c),      &
                           tfu_drm(ims,kms,jms,P_tfu),        & !placeholder
                           tfu_drm(ims,kms,jms,P_tfv),        & !placeholder
                           tfu_drm(ims,kms,jms,P_tfw),        & !placeholder
                           grid%defor12, grid%defor13,        &
                           grid%defor23, 'd',                 &
                           ids, ide, jds, jde, kds, kde,      &
                           ims, ime, jms, jme, kms, kme,      &
                           ips, ipe, jps, jpe, kps, kpe,      &
                           grid%i_start(ij), grid%i_end(ij),  &
                           grid%j_start(ij), grid%j_end(ij),  &
                           k_start, k_end                     )

! Do not project these deformations (already at cell center)
! but divide by two

              CALL proj_c( sijc_drm(ims,kms,jms,P_s11c),      &
                           sijc_drm(ims,kms,jms,P_s22c),      &
                           sijc_drm(ims,kms,jms,P_s33c),      &
                           tfu_drm(ims,kms,jms,P_tfu),        & !placeholder
                           tfu_drm(ims,kms,jms,P_tfv),        & !placeholder
                           tfu_drm(ims,kms,jms,P_tfw),        & !placeholder
                           grid%defor11, grid%defor22,        &
                           grid%defor33, 'c',                 &
                           ids, ide, jds, jde, kds, kde,      &
                           ims, ime, jms, jme, kms, kme,      &
                           ips, ipe, jps, jpe, kps, kpe,      &
                           grid%i_start(ij), grid%i_end(ij),  &
                           grid%j_start(ij), grid%j_end(ij),  &
                           k_start, k_end                     )

              IF ( config_flags%rsfs_opt .GE. 0 ) THEN  ! USE HIJ !Julia

              CALL proj_c( hijc_drm(ims,kms,jms,P_h12c),      &
                           hijc_drm(ims,kms,jms,P_h13c),      &
                           hijc_drm(ims,kms,jms,P_h23c),      &
                           tfu_drm(ims,kms,jms,P_tfu),        & !placeholder
                           tfu_drm(ims,kms,jms,P_tfv),        & !placeholder
                           tfu_drm(ims,kms,jms,P_tfw),        & !placeholder
                           hij_drm(ims,kms,jms,P_h12),        &
                           hij_drm(ims,kms,jms,P_h13),        &
                           hij_drm(ims,kms,jms,P_h23), 'h',   &
                           ids, ide, jds, jde, kds, kde,      &
                           ims, ime, jms, jme, kms, kme,      &
                           ips, ipe, jps, jpe, kps, kpe,      &
                           grid%i_start(ij), grid%i_end(ij),  &
                           grid%j_start(ij), grid%j_end(ij),  &
                           k_start, k_end                     )

              ENDIF


      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate u*v, u*w and v*w at cell center nodes above the surface:
! uvc = uc*vc, uwc = uc*wc, vwc = vc*wc
!_______________________________________________________________________

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles !------------------------------

              CALL product( uc_drm(ims,kms,jms,P_uvc),        &
                            uc_drm(ims,kms,jms,P_uc),         &
                            uc_drm(ims,kms,jms,P_vc),         &
                            ids, ide, jds, jde, kds, kde,     &
                            ims, ime, jms, jme, kms, kme,     &
                            ips, ipe, jps, jpe, kps, kpe,     &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij), &
                            k_start, k_end                    )

              CALL product( uc_drm(ims,kms,jms,P_uwc),        &
                            uc_drm(ims,kms,jms,P_uc),         &
                            uc_drm(ims,kms,jms,P_wc),         &
                            ids, ide, jds, jde, kds, kde,     &
                            ims, ime, jms, jme, kms, kme,     &
                            ips, ipe, jps, jpe, kps, kpe,     &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij), &
                            k_start, k_end                    )

              CALL product( uc_drm(ims,kms,jms,P_vwc),        &
                            uc_drm(ims,kms,jms,P_vc),         &
                            uc_drm(ims,kms,jms,P_wc),         &
                            ids, ide, jds, jde, kds, kde,     &
                            ims, ime, jms, jme, kms, kme,     &
                            ips, ipe, jps, jpe, kps, kpe,     &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij), &
                            k_start, k_end                    )


            ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!______________________________________________________________________
!
! Here project variables straight down to the surface--no horizontal
! projection. Since all input variables are at cell center, all
! surface variables will be as well.
!______________________________________________________________________

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

    CALL uvw_s_c( grid%usc, grid%vsc, grid%ws,      &
                  grid%uusc, grid%vvsc, grid%wws,   &
                  grid%usd, grid%vsd, grid%uvsd,    &
                  grid%wse, grid%uwse,              &
                  grid%wsf, grid%vwsf,              &
                  uc_drm(ims,kms,jms,P_uc),         &
                  uc_drm(ims,kms,jms,P_vc),         &
                  uc_drm(ims,kms,jms,P_wc),         &
                  uc_drm(ims,kms,jms,P_uuc),        &
                  uc_drm(ims,kms,jms,P_vvc),        &
                  uc_drm(ims,kms,jms,P_wwc),        &
                  uc_drm(ims,kms,jms,P_uc),         & !uc is placeholder
                  uc_drm(ims,kms,jms,P_uc),         & !uc is placeholder
                  uc_drm(ims,kms,jms,P_uvc),        &
                  uc_drm(ims,kms,jms,P_uc),         & !uc is placeholder
                  uc_drm(ims,kms,jms,P_uwc),        &
                  uc_drm(ims,kms,jms,P_uc),         & !uc is placeholder
                  uc_drm(ims,kms,jms,P_vwc),        &
                  grid%ht, grid%rdx, grid%rdy,      &
                  grid%cf1, grid%cf2, grid%cf3,     &
                  config_flags%rsfsvsbc_opt,        &
                  ids, ide, jds, jde, kds, kde,     &
                  ims, ime, jms, jme, kms, kme,     &
                  ips, ipe, jps, jpe, kps, kpe,     &
                  grid%i_start(ij), grid%i_end(ij), &
                  grid%j_start(ij), grid%j_end(ij), &
                  k_start, k_end                    )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos for all variables that will be filtered before filtering.
!_______________________________________________________________________


#ifdef DM_PARALLEL
#      include "HALO_EM_UC_DRM.inc"
#      include "PERIOD_EM_UC_DRM.inc"
#      include "HALO_EM_SIJC_DRM.inc"
#      include "PERIOD_EM_SIJC_DRM.inc"
#      include "HALO_EM_RSFS_PUVWS.inc"
#      include "PERIOD_EM_RSFS_PUVWS.inc"
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !----------------------------

    CALL set_physical_bc3d_m( uc_drm(ims,kms,jms,P_uc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( uc_drm(ims,kms,jms,P_vc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( uc_drm(ims,kms,jms,P_wc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( uc_drm(ims,kms,jms,P_uuc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( uc_drm(ims,kms,jms,P_vvc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( uc_drm(ims,kms,jms,P_wwc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( uc_drm(ims,kms,jms,P_uvc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( uc_drm(ims,kms,jms,P_uwc), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( uc_drm(ims,kms,jms,P_vwc), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( sijc_drm(ims,kms,jms,P_s11c), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( sijc_drm(ims,kms,jms,P_s22c), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( sijc_drm(ims,kms,jms,P_s33c), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc3d_m( sijc_drm(ims,kms,jms,P_s12c), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( sijc_drm(ims,kms,jms,P_s13c), 't', &
                            config_flags,                        &
                            ids, ide, jds, jde, kds, kde,        &
                            ims, ime, jms, jme, kms, kme,        &
                            ips, ipe, jps, jpe, kps, kpe,        &
                            grid%i_start(ij), grid%i_end(ij),    &
                            grid%j_start(ij), grid%j_end(ij),    &
                            k_start, k_end                       )

    CALL set_physical_bc3d_m( sijc_drm(ims,kms,jms,P_s23c), 't', &
                            config_flags,                         &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                            grid%i_start(ij), grid%i_end(ij),     &
                            grid%j_start(ij), grid%j_end(ij),     &
                            k_start, k_end                        )

    CALL set_physical_bc2d_m( grid%usc, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uusc, 't', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vsc, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vvsc, 't', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%ws, 't', config_flags,       &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%wws, 't', config_flags,      &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uvsd, 'd', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%uwse, 'e', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

    CALL set_physical_bc2d_m( grid%vwsf, 'f', config_flags,     &
                            ids, ide, jds, jde,               &
                            ims, ime, jms, jme,               &
                            ips, ipe, jps, jpe,               &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij)  )

  ENDDO !-------------------------------------------------
  !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Test-filter uc -> tfuc, etc. There is no need to explicitly filter
! these variables first as the prognostic variables are already
! considered to have been explicitly filtered.
!_______________________________________________________________________

! ASW debug
! print*,'Test-filter uc -> tfuc'
! CALL writecheckm_drm( tfuc_drm(ims,kms,jms,P_tfuc),       &
!                       ims, ime, jms, jme, kms, kme,     &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL drm_filter( tfuc_drm(ims,kms,jms,P_tfuc),     &
                         uc_drm(ims,kms,jms,P_uc),         &
                         grid%usc, 'cc', itst,             &
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfuc_drm(ims,kms,jms,P_tfvc),     &
                         uc_drm(ims,kms,jms,P_vc),         &
                         grid%vsc, 'cc', itst,             &
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfuc_drm(ims,kms,jms,P_tfwc),     &
                         uc_drm(ims,kms,jms,P_wc),         &
                         grid%ws, 'cw', itst,              &
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfuc_drm(ims,kms,jms,P_tfuuc),    &
                         uc_drm(ims,kms,jms,P_uuc),        &
                         grid%uusc, 'cc', itst,            &
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfuc_drm(ims,kms,jms,P_tfvvc),    &
                         uc_drm(ims,kms,jms,P_vvc),        &
                         grid%vvsc, 'cc', itst,            &
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfuc_drm(ims,kms,jms,P_tfwwc),    &
                         uc_drm(ims,kms,jms,P_wwc),        &
                         grid%wws, 'cw', itst,             &
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfuc_drm(ims,kms,jms,P_tfuvc),    &
                         uc_drm(ims,kms,jms,P_uvc),        &
                         grid%uvsd, 'cc', itst,            &
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfuc_drm(ims,kms,jms,P_tfuwc),    &
                         uc_drm(ims,kms,jms,P_uwc),        &
                         grid%uwse, 'cw', itst,            &
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfuc_drm(ims,kms,jms,P_tfvwc),    &
                         uc_drm(ims,kms,jms,P_vwc),        &
                         grid%vwsf, 'cw', itst,            &
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfsijc_drm(ims,kms,jms,P_tfs11c), &
                         sijc_drm(ims,kms,jms,P_s11c),     &
                         grid%us, 'sh', itst,              & !us is placeholder.
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfsijc_drm(ims,kms,jms,P_tfs22c), &
                         sijc_drm(ims,kms,jms,P_s22c),     &
                         grid%us, 'sh', itst,              & !us is placeholder.
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfsijc_drm(ims,kms,jms,P_tfs33c), &
                         sijc_drm(ims,kms,jms,P_s33c),     &
                         grid%us, 'sh', itst,              & !us is placeholder.
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfsijc_drm(ims,kms,jms,P_tfs12c), &
                         sijc_drm(ims,kms,jms,P_s12c),     &
                         grid%us, 'sh', itst,              & !us is placeholder.
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfsijc_drm(ims,kms,jms,P_tfs13c), &
                         sijc_drm(ims,kms,jms,P_s13c),     &
                         grid%us, 'sh', itst,              & !us is placeholder.
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )

        CALL drm_filter( tfsijc_drm(ims,kms,jms,P_tfs23c), &
                         sijc_drm(ims,kms,jms,P_s23c),     &
                         grid%us, 'sh', itst,              & !us is placeholder.
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         grid%i_start(ij), grid%i_end(ij), &
                         grid%j_start(ij), grid%j_end(ij), &
                         k_start, k_end                    )


       ENDDO !---------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Calculate numerator and denominator of the DWL_K
!_______________________________________________________________________

! ASW debug
! print*,'before drmnd tfu'
! CALL writecheckm_drm( tfu_drm(ims,kms,jms,P_tfu),       &
!                       ims, ime, jms, jme, kms, kme,     &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )

! ASW debug
! print*,'before drmnd tfuc'
! CALL writecheckm_drm( tfuc_drm(ims,kms,jms,P_tfuc),       &
!                       ims, ime, jms, jme, kms, kme,     &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------------

        tmp = 1.0 / 3.0
        
              CALL drmnd( evc_drm(ims,kms,jms,P_drmn),      &
                          evc_drm(ims,kms,jms,P_drmd),      &
                          hij_drm(ims,kms,jms,P_h11),       &
                          hij_drm(ims,kms,jms,P_h22),       &
                          hij_drm(ims,kms,jms,P_h33),       &
                          hijc_drm(ims,kms,jms,P_h12c),     &
                          hijc_drm(ims,kms,jms,P_h13c),     &
                          hijc_drm(ims,kms,jms,P_h23c),     &
                          tfsijc_drm(ims,kms,jms,P_tfs11c), &
                          tfsijc_drm(ims,kms,jms,P_tfs22c), &
                          tfsijc_drm(ims,kms,jms,P_tfs33c), &
                          tfsijc_drm(ims,kms,jms,P_tfs12c), &
                          tfsijc_drm(ims,kms,jms,P_tfs13c), &
                          tfsijc_drm(ims,kms,jms,P_tfs23c), &
                          tfuc_drm(ims,kms,jms,P_tfuc),     &
                          tfuc_drm(ims,kms,jms,P_tfvc),     &
                          tfuc_drm(ims,kms,jms,P_tfwc),     &
                          tfuc_drm(ims,kms,jms,P_tfuuc),    &
                          tfuc_drm(ims,kms,jms,P_tfvvc),    &
                          tfuc_drm(ims,kms,jms,P_tfwwc),    &
                          tfuc_drm(ims,kms,jms,P_tfuvc),    &
                          tfuc_drm(ims,kms,jms,P_tfuwc),    &
                          tfuc_drm(ims,kms,jms,P_tfvwc),    &
                          config_flags%dwl_opt,             &
						  grid%BN2,             & !jss buoy vars
				          tmp,                & !jss buoy vars
						  grid%defor11,grid%defor22,grid%defor33,grid%defor12,                       &
						  grid%defor13,grid%defor23,                                           &              
						  ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos and filter
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_EVC_DRM.inc"
#      include "PERIOD_EM_EVC_DRM.inc"
#endif

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------------

         CALL set_physical_bc3d_m( evc_drm(ims,kms,jms,P_drmn),      &
                                 't', config_flags,                &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start, k_end                    )

         CALL set_physical_bc3d_m( evc_drm(ims,kms,jms,P_drmd),      &
                                 't', config_flags,                &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! First apply the explicit filter
!_______________________________________________________________________


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( evc_drm(ims,kms,jms,P_fdrmn),     &
                          evc_drm(ims,kms,jms,P_drmn),      &
                          grid%us, 'kc', iexp,              & !us is placeholder.
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )
     ENDDO !---------------------------------------------------
     !$OMP END PARALLEL DO


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( evc_drm(ims,kms,jms,P_fdrmd),     &
                          evc_drm(ims,kms,jms,P_drmd),      &
                          grid%us, 'kc', iexp,              & !us is placeholder.
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )
     ENDDO !---------------------------------------------------
     !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos and test filter
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_EVC_DRM.inc"
#      include "PERIOD_EM_EVC_DRM.inc"
#endif

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------------

         CALL set_physical_bc3d_m( evc_drm(ims,kms,jms,P_fdrmn),     &
                                 't', config_flags,                &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start, k_end                    )

         CALL set_physical_bc3d_m( evc_drm(ims,kms,jms,P_fdrmd),     &
                                 't',  config_flags,               &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! test filter
!_______________________________________________________________________


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( evc_drm(ims,kms,jms,P_drmn),      &
                          evc_drm(ims,kms,jms,P_fdrmn),     &
                          grid%us, 'kc', itst,              & !us is placeholder.
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )
     ENDDO !---------------------------------------------------
     !$OMP END PARALLEL DO

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( evc_drm(ims,kms,jms,P_drmd),      &
                          evc_drm(ims,kms,jms,P_fdrmd),     &
                          grid%us, 'kc', itst,              & !us is placeholder.
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )
     ENDDO !---------------------------------------------------
     !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! REPEAT filtering of ev coefficients
!_______________________________________________________________________


!_______________________________________________________________________
!
! Update halos and filter
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_EVC_DRM.inc"
#      include "PERIOD_EM_EVC_DRM.inc"
#endif

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------------

         CALL set_physical_bc3d_m( evc_drm(ims,kms,jms,P_drmn),      &
                                 't', config_flags,                &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start, k_end                    )

         CALL set_physical_bc3d_m( evc_drm(ims,kms,jms,P_drmd),      &
                                 't', config_flags,                &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! First apply the explicit filter
!_______________________________________________________________________


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( evc_drm(ims,kms,jms,P_fdrmn),     &
                          evc_drm(ims,kms,jms,P_drmn),      &
                          grid%us, 'kc', iexp,              & !us is placeholder.
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )
     ENDDO !---------------------------------------------------
     !$OMP END PARALLEL DO


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( evc_drm(ims,kms,jms,P_fdrmd),     &
                          evc_drm(ims,kms,jms,P_drmd),      &
                          grid%us, 'kc', iexp,              &  !us is placeholder.
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )
     ENDDO !---------------------------------------------------
     !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Update halos and test filter
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_EVC_DRM.inc"
#      include "PERIOD_EM_EVC_DRM.inc"
#endif

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------------

         CALL set_physical_bc3d_m( evc_drm(ims,kms,jms,P_fdrmn),     &
                                 't', config_flags,                &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start, k_end                    )

         CALL set_physical_bc3d_m( evc_drm(ims,kms,jms,P_fdrmd),     &
                                 't',  config_flags,               &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start, k_end                    )

      ENDDO !---------------------------------------------------
      !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! test filter
!_______________________________________________________________________


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( evc_drm(ims,kms,jms,P_drmn),      &
                          evc_drm(ims,kms,jms,P_fdrmn),     &
                          grid%us, 'kc', itst,              & !us is placeholder.
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )
     ENDDO !---------------------------------------------------
     !$OMP END PARALLEL DO

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------

        CALL rsfs_filter( evc_drm(ims,kms,jms,P_drmd),      &
                          evc_drm(ims,kms,jms,P_fdrmd),     &
                          grid%us, 'kc', itst,              & !us is placeholder.
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          ips, ipe, jps, jpe, kps, kpe,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                    )
     ENDDO !---------------------------------------------------
     !$OMP END PARALLEL DO

!_______________________________________________________________________
!
! Compute DWL eddy viscosity coefficient
!_______________________________________________________________________

!!! ASW: potential halos bug?
!_______________________________________________________________________
!
! Update halos and Compute DWL eddy viscosity coefficient
!_______________________________________________________________________

#ifdef DM_PARALLEL
#      include "HALO_EM_EVC_DRM.inc"
#      include "PERIOD_EM_EVC_DRM.inc"
#endif
!!! ASW


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !------------------------------

             CALL drmevc( grid%xkmv, grid%xkhv, grid%xkmh, grid%xkhh, &
                          evc_drm(ims,kms,jms,P_drmn),                &
                          evc_drm(ims,kms,jms,P_drmd),                &
                          grid%bn2,                                   &
						  config_flags%dwl_pr_opt,                    & ! ASW
					      config_flags%dwl_stab_opt,                  & ! ASW
                          ids, ide, jds, jde, kds, kde,               &
                          ims, ime, jms, jme, kms, kme,               &
                          ips, ipe, jps, jpe, kps, kpe,               &
                          grid%i_start(ij), grid%i_end(ij),           &
                          grid%j_start(ij), grid%j_end(ij),           &
                          k_start, k_end                              )

      ENDDO !---------------------------------------------------
     !$OMP END PARALLEL DO


! print*,'h13'
! CALL writecheckm_drm( hij_drm(ims,kms,jms,P_h13),       &
!                       ims, ime, jms, jme, kms, kme,     &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )
! print*,'h23'
! CALL writecheckm_drm( hij_drm(ims,kms,jms,P_h23),       &
!                       ims, ime, jms, jme, kms, kme ,    &
!                       grid%i_start(ij), grid%i_end(ij), &
!                       grid%j_start(ij), grid%j_end(ij), &
!                       k_start, k_end                    )

!_______________________________________________________________________
!
! BCs and Halos for XKMH, XKMV, XKHH, XKHV
!_______________________________________________________________________

!_______________________________________________________________________
!
! Update halos for XKMH, XKMV, XKHH, XKHV, NW_TAU13 and NW_TAU23
!_______________________________________________________________________

! Julia
#ifdef DM_PARALLEL
#      include "HALO_EM_EVC_DRM.inc"
#      include "PERIOD_EM_EVC_DRM.inc"
#      include "HALO_EM_EV.inc" 
#      include "PERIOD_BDY_EM_PHY_BC.inc" 
#endif

  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij )
  DO ij = 1 , grid%num_tiles !--------------------------------

        CALL set_physical_bc3d_m( grid%xkmh, 't', config_flags,     &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start, k_end                    )

        CALL set_physical_bc3d_m( grid%xkhh, 't', config_flags,     &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start, k_end                    )

        CALL set_physical_bc3d_m( grid%xkmv, 't', config_flags,     &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start, k_end                    )

        CALL set_physical_bc3d_m( grid%xkhv, 't', config_flags,     &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start, k_end                    )

  ENDDO !-----------------------------------------------------
  !$OMP END PARALLEL DO

  IF ( config_flags%can_opt .EQ. 1 ) THEN

 !_______________________________________________________________________
 !
 ! Calculate canopy stresses near the surface ...
 !_______________________________________________________________________

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !--------------------------------

           CALL drm_nw_canopy( grid%nwtau13, grid%nwtau23,       &
                               grid%u_2, grid%v_2,               &
                               grid%w_2, grid%rdz,               &
                               grid%rdzw, grid%dx,               &
                               grid%fnm, grid%fnp,               &
                               !config_flags%z_rough,       & !ASW
                               config_flags%canfact,             &
                               grid%ustm,                        &
                               ids, ide, jds, jde, kds, kde,     &
                               ims, ime, jms, jme, kms, kme,     &
                               ips, ipe, jps, jpe, kps, kpe,     &
                               grid%i_start(ij), grid%i_end(ij), &
                               grid%j_start(ij), grid%j_end(ij), &
                               k_start, k_end                    )

     ENDDO !-----------------------------------------------------
     !$OMP END PARALLEL DO

 !_______________________________________________________________________
 !
 ! BCs and halos for NWTAU13 and NWTAU23
 !_______________________________________________________________________

 #ifdef DM_PARALLEL
 #      include "HALO_EM_NWTAU.inc"
 #      include "PERIOD_EM_NWTAU.inc"
 #endif

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles !--------------------------------

           CALL set_physical_bc3d_m( grid%nwtau13, 'e', config_flags,  &
                                   ids, ide, jds, jde, kds, kde,     &
                                   ims, ime, jms, jme, kms, kme,     &
                                   ips, ipe, jps, jpe, kps, kpe,     &
                                   grid%i_start(ij), grid%i_end(ij), &
                                   grid%j_start(ij), grid%j_end(ij), &
                                   k_start, k_end                    )

           CALL set_physical_bc3d_m( grid%nwtau23, 'f', config_flags,  &
                                   ids, ide, jds, jde, kds, kde,     &
                                   ims, ime, jms, jme, kms, kme,     &
                                   ips, ipe, jps, jpe, kps, kpe,     &
                                   grid%i_start(ij), grid%i_end(ij), &
                                   grid%j_start(ij), grid%j_end(ij), &
                                   k_start, k_end                    )

     ENDDO !-----------------------------------------------------
     !$OMP END PARALLEL DO

  ENDIF

!=======================================================================
!
!                     END DRM IMPLEMENTATION
!
!=======================================================================

ENDIF !dwl_opt .GT. 0

END SUBROUTINE sfs_driver

END MODULE module_sfs_driver
