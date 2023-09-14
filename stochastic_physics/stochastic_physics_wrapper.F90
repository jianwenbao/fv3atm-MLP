module stochastic_physics_wrapper_mod

  use machine, only: kind_phys

  implicit none

  ! For stochastic physics pattern generation
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlat
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlon
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sppt_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: shum_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebu_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebv_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sfc_wts
  real(kind=kind_phys), dimension(:,:,:,:), allocatable, save :: spp_wts

!mlp
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_tcnv_cpl,tcnvtend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_qcnv_cpl,qcnvtend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_ucnv_cpl,ucnvtend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_vcnv_cpl,vcnvtend_cpl
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_tcnv_diag,tcnvtend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_qcnv_diag,qcnvtend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_ucnv_diag,ucnvtend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_vcnv_diag,vcnvtend_diag


  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_tmp_cpl,tmptend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_qmp_cpl,qmptend_cpl
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_tmp_diag,tmptend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_qmp_diag,qmptend_diag


  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_tpbl_cpl,tpbltend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_qpbl_cpl,qpbltend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_upbl_cpl,upbltend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_vpbl_cpl,vpbltend_cpl
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_tpbl_diag,tpbltend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_qpbl_diag,qpbltend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_upbl_diag,upbltend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_vpbl_diag,vpbltend_diag

  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_tshalcnv_cpl,tshalcnvtend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_qshalcnv_cpl,qshalcnvtend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_ushalcnv_cpl,ushalcnvtend_cpl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_vshalcnv_cpl,vshalcnvtend_cpl
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_tshalcnv_diag,tshalcnvtend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_qshalcnv_diag,qshalcnvtend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_ushalcnv_diag,ushalcnvtend_diag
! real(kind=kind_phys), dimension(:,:,:), allocatable, save :: mlp_pert_vshalcnv_diag,vshalcnvtend_diag

!mlp

  logical, save :: is_initialized = .false.
  integer, save :: lsoil = -999
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: smc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: stc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: slc
  !
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: vfrac
  !albedo
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: snoalb
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnsf
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnwf
  !emissivity
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: semis
  !roughness length for land
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: zorll

  !real(kind=kind_phys), dimension(:,:),   allocatable, save :: stype
  integer, dimension(:,:),   allocatable, save :: stype

  ! For cellular automata
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: sst
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lmsk
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lake
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: condition
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca_deep_cpl, ca_turb_cpl, ca_shal_cpl
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca1_cpl, ca2_cpl, ca3_cpl


!----------------
! Public Entities
!----------------
! functions
  public stochastic_physics_wrapper

  contains

!-------------------------------
!  CCPP step
!-------------------------------
  subroutine stochastic_physics_wrapper (GFS_Control, GFS_Data, Atm_block, ierr)

#ifdef _OPENMP
    use omp_lib
#endif

    use GFS_typedefs,       only: GFS_control_type, GFS_data_type
    use mpp_mod,            only: FATAL, mpp_error
    use block_control_mod,  only: block_control_type
    use atmosphere_mod,     only: Atm, mygrid

    use stochastic_physics,           only: init_stochastic_physics, run_stochastic_physics
    use cellular_automata_global_mod, only: cellular_automata_global
    use cellular_automata_sgs_mod,    only: cellular_automata_sgs
    use lndp_apply_perts_mod,         only: lndp_apply_perts

    implicit none

    type(GFS_control_type),   intent(inout) :: GFS_Control
    type(GFS_data_type),      intent(inout) :: GFS_Data(:)
    type(block_control_type), intent(inout) :: Atm_block
    integer,                  intent(out)   :: ierr

    integer :: nthreads, nb, levs, maxblk, nblks, n, v
    logical :: param_update_flag

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif
    ierr = 0

    levs   = GFS_Control%levs
    maxblk = maxval(GFS_Control%blksz)
    nblks  = Atm_block%nblks

    ! Initialize

    initalize_stochastic_physics: if (.not. is_initialized) then

      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb  .OR. GFS_Control%do_mlp .OR. (GFS_Control%lndp_type > 0) .OR. GFS_Control%do_spp) then
         allocate(xlat(1:nblks,maxblk))
         allocate(xlon(1:nblks,maxblk))
         do nb=1,nblks
            xlat(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlat(:)
            xlon(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlon(:)
         end do
!mlp  
      !print*,'allocate do_mlp not initialized',GFS_Control%do_mlp
      if(GFS_Control%do_mlp)then
         allocate(mlp_pert_tcnv_cpl (1:nblks,maxblk,1:levs))
         allocate(tcnvtend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_qcnv_cpl (1:nblks,maxblk,1:levs))
         allocate(qcnvtend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_ucnv_cpl (1:nblks,maxblk,1:levs))
         allocate(ucnvtend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_vcnv_cpl (1:nblks,maxblk,1:levs))
         allocate(vcnvtend_cpl (1:nblks,maxblk,1:levs))

         allocate(mlp_pert_tmp_cpl (1:nblks,maxblk,1:levs))
         allocate(tmptend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_qmp_cpl (1:nblks,maxblk,1:levs))
         allocate(qmptend_cpl (1:nblks,maxblk,1:levs))

         allocate(mlp_pert_tpbl_cpl (1:nblks,maxblk,1:levs))
         allocate(tpbltend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_qpbl_cpl (1:nblks,maxblk,1:levs))
         allocate(qpbltend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_upbl_cpl (1:nblks,maxblk,1:levs))
         allocate(upbltend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_vpbl_cpl (1:nblks,maxblk,1:levs))
         allocate(vpbltend_cpl (1:nblks,maxblk,1:levs))

         allocate(mlp_pert_tshalcnv_cpl (1:nblks,maxblk,1:levs))
         allocate(tshalcnvtend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_qshalcnv_cpl (1:nblks,maxblk,1:levs))
         allocate(qshalcnvtend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_ushalcnv_cpl (1:nblks,maxblk,1:levs))
         allocate(ushalcnvtend_cpl (1:nblks,maxblk,1:levs))
         allocate(mlp_pert_vshalcnv_cpl (1:nblks,maxblk,1:levs))
         allocate(vshalcnvtend_cpl (1:nblks,maxblk,1:levs))


      endif 

! end mlp
        ! Initialize stochastic physics
        call init_stochastic_physics(levs, GFS_Control%blksz, GFS_Control%dtp, GFS_Control%sppt_amp,                           &
            GFS_Control%input_nml_file, GFS_Control%fn_nml, GFS_Control%nlunit, xlon, xlat, GFS_Control%do_sppt, GFS_Control%do_shum, &
!mlp
            GFS_Control%do_skeb, GFS_Control%lndp_type,GFS_Control%do_mlp, GFS_Control%do_mlp_cnv,GFS_Control%do_mlp_shalcnv,   &
            GFS_Control%do_mlp_mp,GFS_Control%do_mlp_pbl, mlp_pert_tcnv_cpl,   &
                      mlp_pert_qcnv_cpl, mlp_pert_ucnv_cpl, mlp_pert_vcnv_cpl,  mlp_pert_tmp_cpl, mlp_pert_qmp_cpl, &
                      mlp_pert_tpbl_cpl, mlp_pert_qpbl_cpl, mlp_pert_upbl_cpl,  mlp_pert_vpbl_cpl, &
                      mlp_pert_tshalcnv_cpl,mlp_pert_qshalcnv_cpl, mlp_pert_ushalcnv_cpl,  mlp_pert_vshalcnv_cpl, &
!mlp
             GFS_Control%n_var_lndp, GFS_Control%use_zmtnblck, GFS_Control%skeb_npass,     &
            GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,    &
            GFS_Control%n_var_spp, GFS_Control%spp_var_list, GFS_Control%spp_prt_list, GFS_Control%spp_stddev_cutoff, GFS_Control%do_spp,                            &
            GFS_Control%ak, GFS_Control%bk, nthreads, GFS_Control%master, GFS_Control%communicator, ierr)
            if (ierr/=0)  then
                    write(6,*) 'call to init_stochastic_physics failed'
                    write(6,*) 'call to init_stochastic_physics domlp',GFS_Control%do_mlp
                    write(6,*) 'call to init_stochastic_physics vsha',mlp_pert_vshalcnv_cpl
                    write(6,*) 'call to init_stochastic_physics domlpcnv',GFS_Control%do_mlp_cnv
                    write(6,*) 'call to init_stochastic_physics domlpshal',GFS_Control%do_mlp_shalcnv
                    write(6,*) 'call to init_stochastic_physics domlpmp',GFS_Control%do_mlp_mp
                    write(6,*) 'call to init_stochastic_physics domlppbl',GFS_Control%do_mlp_pbl
                    return
            endif
      end if
      if (GFS_Control%do_sppt) then
         allocate(sppt_wts(1:nblks,maxblk,1:levs))
      end if
      if (GFS_Control%do_shum) then
         allocate(shum_wts(1:nblks,maxblk,1:levs))
      end if
      if (GFS_Control%do_skeb) then
         allocate(skebu_wts(1:nblks,maxblk,1:levs))
         allocate(skebv_wts(1:nblks,maxblk,1:levs))
      end if
      if ( GFS_Control%do_spp ) then
         allocate(spp_wts(1:nblks,maxblk,1:levs,1:GFS_Control%n_var_spp))
         do n=1,GFS_Control%n_var_spp
           select case (trim(GFS_Control%spp_var_list(n)))
           case('pbl')
             GFS_Control%spp_pbl = 1
           case('sfc')
             GFS_Control%spp_sfc = 1
           case('mp')
             GFS_Control%spp_mp = 7
           case('rad')
             GFS_Control%spp_rad = 1
           case('gwd')
             GFS_Control%spp_gwd = 1
           end select
         end do
      end if
      if ( GFS_Control%lndp_type == 2 ) then
          allocate(sfc_wts(1:nblks,maxblk,1:GFS_Control%n_var_lndp))
          if ( (GFS_Control%lsm == GFS_Control%lsm_noah) .or. (GFS_Control%lsm == GFS_Control%lsm_noahmp)) then
            lsoil = GFS_Control%lsoil
          elseif (GFS_Control%lsm == GFS_Control%lsm_ruc) then
            lsoil = GFS_Control%lsoil_lsm
          endif
          allocate(smc   (1:nblks, maxblk, lsoil))
          do v = 1,GFS_Control%n_var_lndp
            select case (trim(GFS_Control%lndp_var_list(v)))
            case('smc')
              allocate(slc   (1:nblks, maxblk, lsoil))
              allocate(stype (1:nblks, maxblk))
            case('stc')
              allocate(stc   (1:nblks, maxblk, lsoil))
            case('vgf') 
              allocate(vfrac (1:nblks, maxblk))
            case('alb') 
              allocate(alnsf (1:nblks, maxblk))
              allocate(alnwf (1:nblks, maxblk))
            case('sal') 
              allocate(snoalb(1:nblks, maxblk))
            case('emi') 
              allocate(semis (1:nblks, maxblk))
            case('zol') 
              allocate(zorll (1:nblks, maxblk))
            endselect 
          enddo 
      endif


      if ( GFS_Control%lndp_type == 1 ) then ! this scheme sets perts once
         allocate(sfc_wts(1:nblks, maxblk, GFS_Control%n_var_lndp))
         call run_stochastic_physics(levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz,       &
                                     sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts,         &
                                     skebv_wts=skebv_wts, sfc_wts=sfc_wts,                              &
!mlp
                                     spp_wts=spp_wts,&
                      mlp_pert_tcnv_cpl=mlp_pert_tcnv_cpl,tcnvtend_cpl=tcnvtend_cpl,         &
                      mlp_pert_qcnv_cpl=mlp_pert_qcnv_cpl,qcnvtend_cpl=qcnvtend_cpl,         &
                      mlp_pert_ucnv_cpl=mlp_pert_ucnv_cpl,ucnvtend_cpl=ucnvtend_cpl,             &
                      mlp_pert_vcnv_cpl=mlp_pert_vcnv_cpl,vcnvtend_cpl=vcnvtend_cpl,             &
                      mlp_pert_tmp_cpl=mlp_pert_tmp_cpl,tmptend_cpl=tmptend_cpl,               &
                      mlp_pert_qmp_cpl=mlp_pert_qmp_cpl,qmptend_cpl=qmptend_cpl,               &
                      mlp_pert_tpbl_cpl=mlp_pert_tpbl_cpl,tpbltend_cpl=tpbltend_cpl,         &
                      mlp_pert_qpbl_cpl=mlp_pert_qpbl_cpl,qpbltend_cpl=qpbltend_cpl,         &
                      mlp_pert_upbl_cpl=mlp_pert_upbl_cpl,upbltend_cpl=upbltend_cpl,             &
                      mlp_pert_vpbl_cpl=mlp_pert_vpbl_cpl,vpbltend_cpl=vpbltend_cpl,             &
                      mlp_pert_tshalcnv_cpl=mlp_pert_tshalcnv_cpl,tshalcnvtend_cpl=tshalcnvtend_cpl,    &
                      mlp_pert_qshalcnv_cpl=mlp_pert_qshalcnv_cpl,qshalcnvtend_cpl=qshalcnvtend_cpl,     &
                      mlp_pert_ushalcnv_cpl=mlp_pert_ushalcnv_cpl,ushalcnvtend_cpl=ushalcnvtend_cpl,      &
                      mlp_pert_vshalcnv_cpl=mlp_pert_vshalcnv_cpl,vshalcnvtend_cpl=vshalcnvtend_cpl,      &
                      nthreads=nthreads)

         ! Copy contiguous data back
         do nb=1,nblks
            GFS_Data(nb)%Coupling%sfc_wts(:,:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
         end do
         deallocate(sfc_wts)
      end if
      ! Consistency check for cellular automata
      if(GFS_Control%do_ca)then
        ! DH* The current implementation of cellular_automata assumes that all blocksizes are the
        ! same - abort if this is not the case, otherwise proceed with Atm_block%blksz(1) below
        if (.not. minval(Atm_block%blksz) == maxblk) then
           call mpp_error(FATAL, 'Logic errror: cellular_automata not compatible with non-uniform blocksizes')
        end if
        if(GFS_Control%ca_sgs)then
           allocate(sst         (1:nblks, maxblk))
           allocate(lmsk        (1:nblks, maxblk))
           allocate(lake        (1:nblks, maxblk))
           allocate(condition   (1:nblks, maxblk))
           allocate(ca_deep_cpl (1:nblks, maxblk))
           allocate(ca_turb_cpl (1:nblks, maxblk))
           allocate(ca_shal_cpl (1:nblks, maxblk))
        endif
        if(GFS_Control%ca_global)then
          ! Allocate contiguous arrays; no need to copy in (intent out)
          allocate(ca1_cpl (1:nblks, maxblk))
          allocate(ca2_cpl (1:nblks, maxblk))
          allocate(ca3_cpl (1:nblks, maxblk))
        endif
      endif

      is_initialized = .true.

    else initalize_stochastic_physics
!mlp  

! previous pert ?
      if (GFS_Control%do_mlp) then
         do nb=1,nblks
              if(nb.eq.1)then
                !print*,"in wrapper tpbl tend is ", GFS_Data(nb)%Coupling%tpbltend(:,:)
              !  print*,"in wrapper tpbl tend is ", GFS_Data(nb)%Coupling%tpbltend(:,:)
              !  print*,"in wrapper tmp tend is ", GFS_Data(nb)%Coupling%tmptend(:,:)
              endif
             mlp_pert_tcnv_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_tcnv(:,:)
             tcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%tcnvtend(:,:)
             mlp_pert_qcnv_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_qcnv(:,:)
             qcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%qcnvtend(:,:)
             mlp_pert_ucnv_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_ucnv(:,:)
             ucnvtend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%ucnvtend(:,:)
             mlp_pert_vcnv_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_vcnv(:,:)
             vcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%vcnvtend(:,:)

             mlp_pert_tmp_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_tmp(:,:)
             tmptend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%tmptend(:,:)
             mlp_pert_qmp_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_qmp(:,:)
             qmptend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%qmptend(:,:)

             mlp_pert_tpbl_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_tpbl(:,:)
             tpbltend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%tpbltend(:,:)
             mlp_pert_qpbl_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_qpbl(:,:)
             qpbltend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%qpbltend(:,:)
             mlp_pert_upbl_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_upbl(:,:)
             upbltend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%upbltend(:,:)
             mlp_pert_vpbl_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_vpbl(:,:)
             vpbltend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%vpbltend(:,:)

             mlp_pert_tshalcnv_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_tshalcnv(:,:)
             tshalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%tshalcnvtend(:,:)
             mlp_pert_qshalcnv_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_qshalcnv(:,:)
             qshalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%qshalcnvtend(:,:)
             mlp_pert_ushalcnv_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_ushalcnv(:,:)
             ushalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%ushalcnvtend(:,:)
             mlp_pert_vshalcnv_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%mlp_pert_vshalcnv(:,:)
             vshalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Coupling%vshalcnvtend(:,:)
         enddo
      endif 
! end mlp
      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. GFS_Control%do_mlp .OR. (GFS_Control%lndp_type == 2) .OR. GFS_Control%do_spp) then
         call run_stochastic_physics(levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz, &
                                 sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
! mlp
                                 spp_wts=spp_wts,      &   
                 mlp_pert_tcnv_cpl=mlp_pert_tcnv_cpl,tcnvtend_cpl=tcnvtend_cpl, &   
                 mlp_pert_qcnv_cpl=mlp_pert_qcnv_cpl,qcnvtend_cpl=qcnvtend_cpl, &   
                 mlp_pert_ucnv_cpl=mlp_pert_ucnv_cpl,ucnvtend_cpl=ucnvtend_cpl, &   
                 mlp_pert_vcnv_cpl=mlp_pert_vcnv_cpl,vcnvtend_cpl=vcnvtend_cpl, &   
                 mlp_pert_tmp_cpl=mlp_pert_tmp_cpl,tmptend_cpl=tmptend_cpl, &   
                 mlp_pert_qmp_cpl=mlp_pert_qmp_cpl,qmptend_cpl=qmptend_cpl, &   
                 mlp_pert_tpbl_cpl=mlp_pert_tpbl_cpl,tpbltend_cpl=tpbltend_cpl, &   
                 mlp_pert_qpbl_cpl=mlp_pert_qpbl_cpl,qpbltend_cpl=qpbltend_cpl, &   
                 mlp_pert_upbl_cpl=mlp_pert_upbl_cpl,upbltend_cpl=upbltend_cpl, &   
                 mlp_pert_vpbl_cpl=mlp_pert_vpbl_cpl,vpbltend_cpl=vpbltend_cpl, &   
                 mlp_pert_tshalcnv_cpl=mlp_pert_tshalcnv_cpl,tshalcnvtend_cpl=tshalcnvtend_cpl, &
                 mlp_pert_qshalcnv_cpl=mlp_pert_qshalcnv_cpl,qshalcnvtend_cpl=qshalcnvtend_cpl, &
                 mlp_pert_ushalcnv_cpl=mlp_pert_ushalcnv_cpl,ushalcnvtend_cpl=ushalcnvtend_cpl, &
                 mlp_pert_vshalcnv_cpl=mlp_pert_vshalcnv_cpl,vshalcnvtend_cpl=vshalcnvtend_cpl, &
                  nthreads=nthreads)
         ! Copy contiguous data back
         if (GFS_Control%do_sppt) then
            do nb=1,nblks
                GFS_Data(nb)%Coupling%sppt_wts(:,:) = sppt_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
! mlp
! need to call new stochastic for mlp?
          if (GFS_Control%do_mlp) then
            do nb=1,nblks
                GFS_Data(nb)%Coupling%mlp_pert_tcnv(:,:) = mlp_pert_tcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%tcnvtend(:,:) = tcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_qcnv(:,:) = mlp_pert_qcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%qcnvtend(:,:) = tcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_ucnv(:,:) = mlp_pert_ucnv_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%ucnvtend(:,:) = ucnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_vcnv(:,:) = mlp_pert_vcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%vcnvtend(:,:) = vcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
!diag
                GFS_Data(nb)%Intdiag%mlp_pert_tcnv(:,:) = mlp_pert_tcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%tcnvtend(:,:) = tcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_qcnv(:,:) = mlp_pert_qcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%qcnvtend(:,:) = tcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_ucnv(:,:) = mlp_pert_ucnv_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%ucnvtend(:,:) = ucnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_vcnv(:,:) =mlp_pert_vcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%vcnvtend(:,:) = vcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)

!mp
                GFS_Data(nb)%Coupling%mlp_pert_tmp(:,:) = mlp_pert_tmp_cpl(nb,1:GFS_Control%blksz(nb),:)
           !    GFS_Data(nb)%Coupling%tmptend(:,:) = tmptend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_qmp(:,:) = mlp_pert_qmp_cpl(nb,1:GFS_Control%blksz(nb),:)
           !    GFS_Data(nb)%Coupling%qmptend(:,:) = qmptend_cpl(nb,1:GFS_Control%blksz(nb),:)
! diag
                GFS_Data(nb)%Intdiag%mlp_pert_tmp(:,:) = mlp_pert_tmp_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%tmptend(:,:) = tmptend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_qmp(:,:) = mlp_pert_qmp_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%qmptend(:,:) = qmptend_cpl(nb,1:GFS_Control%blksz(nb),:)
!pbl
                GFS_Data(nb)%Coupling%mlp_pert_tpbl(:,:) = mlp_pert_tpbl_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%tpbltend(:,:) = tpbltend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_qpbl(:,:) = mlp_pert_qpbl_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%qpbltend(:,:) = tpbltend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_upbl(:,:) = mlp_pert_upbl_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%upbltend(:,:) = upbltend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_vpbl(:,:) = mlp_pert_vpbl_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%vpbltend(:,:) = vpbltend_cpl(nb,1:GFS_Control%blksz(nb),:)
    
               !print*,"max min qpbltend_cpl in wrapper are " ,minval(qpbltend_cpl),maxval(qpbltend_cpl)
              ! print*,"max min qpbltend in wrapper are " ,minval(GFS_Data(nb)%Coupling%qpbltend(:,:)),maxval(GFS_Data(nb)%Coupling%qpbltend(:,:))
!diag
                GFS_Data(nb)%Intdiag%mlp_pert_tpbl(:,:) = mlp_pert_tpbl_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%tpbltend(:,:) = tpbltend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_qpbl(:,:) = mlp_pert_qpbl_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%qpbltend(:,:) = tpbltend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_upbl(:,:) = mlp_pert_upbl_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%upbltend(:,:) = upbltend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_vpbl(:,:) = mlp_pert_vpbl_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%vpbltend(:,:) = vpbltend_cpl(nb,1:GFS_Control%blksz(nb),:)
!shalcnv
                GFS_Data(nb)%Coupling%mlp_pert_tshalcnv(:,:) = mlp_pert_tshalcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%tshalcnvtend(:,:) = tshalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_qshalcnv(:,:) = mlp_pert_qshalcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%qshalcnvtend(:,:) = tshalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_ushalcnv(:,:) = mlp_pert_ushalcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%ushalcnvtend(:,:) = ushalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%mlp_pert_vshalcnv(:,:) = mlp_pert_vshalcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
            !   GFS_Data(nb)%Coupling%vshalcnvtend(:,:) = vshalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
!diag
                GFS_Data(nb)%Intdiag%mlp_pert_tshalcnv(:,:) = mlp_pert_tshalcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%tshalcnvtend(:,:) = tshalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_qshalcnv(:,:) = mlp_pert_qshalcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%qshalcnvtend(:,:) = tshalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_ushalcnv(:,:) = mlp_pert_ushalcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%ushalcnvtend(:,:) = ushalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%mlp_pert_vshalcnv(:,:) = mlp_pert_vshalcnv_cpl(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Intdiag%vshalcnvtend(:,:) = vshalcnvtend_cpl(nb,1:GFS_Control%blksz(nb),:)



            end do
          end if ! mlp
         end if
             print*,"wrapper max min cplqshalcnvtend  " ,minval(qshalcnvtend_cpl),maxval(qshalcnvtend_cpl)
         if (GFS_Control%do_shum) then
            do nb=1,nblks
                GFS_Data(nb)%Coupling%shum_wts(:,:) = shum_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_skeb) then
            do nb=1,nblks
                GFS_Data(nb)%Coupling%skebu_wts(:,:) = skebu_wts(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%skebv_wts(:,:) = skebv_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_spp) then
            do n=1,GFS_Control%n_var_spp
               select case (trim(GFS_Control%spp_var_list(n)))
               case('pbl')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_pbl(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('sfc')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_sfc(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('mp')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_mp(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('gwd')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_gwd(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('rad')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_rad(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               end select
            end do
         end if

         if (GFS_Control%lndp_type == 2) then ! save wts, and apply lndp scheme
             do nb=1,nblks
                GFS_Data(nb)%Coupling%sfc_wts(:,:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
             end do
 
             do nb=1,nblks
                do v = 1,GFS_Control%n_var_lndp
                  ! used to identify locations with land model (=soil) 
                  if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                     smc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%smois(1:GFS_Control%blksz(nb),1:lsoil)
                  else  ! noah or noah-MP
                     smc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%smc(1:GFS_Control%blksz(nb),1:lsoil)
                  endif

                  select case (trim(GFS_Control%lndp_var_list(v)))
                  case('smc')
                      ! stype used to fetch soil params
                      stype(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%stype(1:GFS_Control%blksz(nb))
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                         slc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%sh2o(1:GFS_Control%blksz(nb),1:lsoil)
                      else  ! noah or noah-MP
                         slc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%slc(1:GFS_Control%blksz(nb),1:lsoil)
                      endif
                  case('stc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                         stc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%tslb(1:GFS_Control%blksz(nb),1:lsoil)
                      else ! noah or noah-MP 
                         stc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%stc(1:GFS_Control%blksz(nb),1:lsoil)
                      endif 
                  case('vgf')
                      if ( (GFS_Control%lsm == GFS_Control%lsm_noahmp) ) then 
                         ! assumes iopt_dveg = 4 (will be checked later)
                         vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%shdmax(1:GFS_Control%blksz(nb))
                      else ! ruc or noah-MP
                         vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%vfrac(1:GFS_Control%blksz(nb))
                      endif
                  case('alb')
                      alnsf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%alnsf(1:GFS_Control%blksz(nb))
                      alnwf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%alnwf(1:GFS_Control%blksz(nb))
                  case('sal')
                      snoalb(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Sfcprop%snoalb(1:GFS_Control%blksz(nb))
                  case('emi')
                      semis(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Radtend%semis(1:GFS_Control%blksz(nb))
                  case('zol')
                      zorll(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%zorll(1:GFS_Control%blksz(nb))
                  endselect
              enddo
             enddo


             param_update_flag = .false.
             ! noah and noah-MP treated differently, as global cycle doesn't overwrite shdmax for Noah-MP 
             ! determine whether land paramaters have been over-written to
             ! trigger applying perturbations (logic copied from GFS_driver),
             if ( (GFS_Control%lsm == GFS_Control%lsm_noah) .and. GFS_Control%nscyc >  0) then
                 if (mod(GFS_Control%kdt,GFS_Control%nscyc) == 1 ) then
                   param_update_flag = .true.
                 endif
             endif
             if ( ( GFS_Control%nscyc ==  0 .or. GFS_Control%lsm == GFS_Control%lsm_noahmp) .and. GFS_Control%first_time_step ) then
             ! call once at start of the forecast.
                    param_update_flag = .true.
             endif
              
             call lndp_apply_perts(GFS_Control%blksz, GFS_Control%lsm, GFS_Control%lsm_noah, GFS_Control%lsm_ruc,             &
                               GFS_Control%lsm_noahmp, GFS_Control%iopt_dveg, lsoil, GFS_Control%dtp, GFS_Control%kdt,        &
                               GFS_Control%n_var_lndp, GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,                  &
                               sfc_wts, xlon, xlat, stype, GFS_Control%pores, GFS_Control%resid,param_update_flag,            &
                               smc, slc, stc, vfrac, alnsf, alnwf, snoalb, semis, zorll, ierr)

             if (ierr/=0)  then
                    write(6,*) 'call to GFS_apply_lndp failed'
                    return
             endif

             do nb=1,nblks
                do v = 1,GFS_Control%n_var_lndp

                  select case (trim(GFS_Control%lndp_var_list(v)))
                  case('smc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                           GFS_Data(nb)%Sfcprop%smois(1:GFS_Control%blksz(nb),1:lsoil) = smc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                           GFS_Data(nb)%Sfcprop%sh2o(1:GFS_Control%blksz(nb),1:lsoil)  = slc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      else  ! noah or noah-MP
                           GFS_Data(nb)%Sfcprop%smc(1:GFS_Control%blksz(nb),1:lsoil) = smc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                           GFS_Data(nb)%Sfcprop%slc(1:GFS_Control%blksz(nb),1:lsoil) = slc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      endif
                  case('stc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                           GFS_Data(nb)%Sfcprop%tslb(1:GFS_Control%blksz(nb),1:lsoil)  = stc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      else ! noah or noah-MP 
                           GFS_Data(nb)%Sfcprop%stc(1:GFS_Control%blksz(nb),1:lsoil) = stc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      endif 
                  case('vgf')
                      if ( (GFS_Control%lsm == GFS_Control%lsm_noahmp) ) then 
                        GFS_Data(nb)%Sfcprop%shdmax(1:GFS_Control%blksz(nb))  = vfrac(nb,1:GFS_Control%blksz(nb))
                      else 
                        GFS_Data(nb)%Sfcprop%vfrac(1:GFS_Control%blksz(nb))  = vfrac(nb,1:GFS_Control%blksz(nb))
                      endif
                  case('alb')
                       GFS_Data(nb)%Sfcprop%alnsf(1:GFS_Control%blksz(nb))  = alnsf(nb,1:GFS_Control%blksz(nb))
                       GFS_Data(nb)%Sfcprop%alnwf(1:GFS_Control%blksz(nb))  = alnwf(nb,1:GFS_Control%blksz(nb))
                  case('sal')
                        GFS_Data(nb)%Sfcprop%snoalb(1:GFS_Control%blksz(nb)) = snoalb(nb,1:GFS_Control%blksz(nb))
                  case('emi')
                        GFS_Data(nb)%Radtend%semis(1:GFS_Control%blksz(nb))  = semis(nb,1:GFS_Control%blksz(nb))
                  case('zol')
                        GFS_Data(nb)%Sfcprop%zorll(1:GFS_Control%blksz(nb))  = zorll(nb,1:GFS_Control%blksz(nb))
                  end select   
                enddo 
            enddo
         endif ! lndp block
      endif ! if do* block

      if (GFS_Control%do_ca) then

       if(GFS_Control%ca_sgs)then
         do nb=1,nblks
             sst        (nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Sfcprop%tsfco(:)
             lmsk       (nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Sfcprop%slmsk(:)
             lake       (nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Sfcprop%lakefrac(:)
             condition  (nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Coupling%condition(:)
             ca_deep_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Coupling%ca_deep(:)
             ca_turb_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Coupling%ca_turb(:)
             ca_shal_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Coupling%ca_shal(:)
         enddo
         call cellular_automata_sgs(GFS_Control%kdt,GFS_control%dtp,GFS_control%restart,GFS_Control%first_time_step,              &
            sst,lmsk,lake,condition,ca_deep_cpl,ca_turb_cpl,ca_shal_cpl, Atm(mygrid)%domain_for_coupler,nblks,                    &
            Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy, levs,                        &
            GFS_Control%nthresh,GFS_Control%tile_num,GFS_Control%nca,GFS_Control%ncells,GFS_Control%nlives,                       &
            GFS_Control%nfracseed, GFS_Control%nseed,GFS_Control%iseed_ca,                                                        &
            GFS_Control%nspinup,GFS_Control%ca_trigger,Atm_block%blksz(1),GFS_Control%master,GFS_Control%communicator)
         ! Copy contiguous data back as needed
         do nb=1,nblks
             GFS_Data(nb)%Coupling%ca_deep(:) = ca_deep_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca_turb(:) = ca_turb_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca_shal(:) = ca_shal_cpl (nb,1:GFS_Control%blksz(nb))
         enddo
       endif
       if(GFS_Control%ca_global)then
          call cellular_automata_global(GFS_Control%kdt,GFS_control%restart,GFS_Control%first_time_step,ca1_cpl,ca2_cpl,ca3_cpl,                &
            Atm(mygrid)%domain_for_coupler, nblks,Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy,levs, &
            GFS_Control%nca_g,GFS_Control%ncells_g,GFS_Control%nlives_g,GFS_Control%nfracseed,GFS_Control%nseed_g,                              &
            GFS_Control%iseed_ca,GFS_control%tile_num,GFS_Control%ca_smooth,GFS_Control%nspinup,Atm_block%blksz(1),                             &
            GFS_Control%nsmooth,GFS_Control%ca_amplitude,GFS_Control%master,GFS_Control%communicator)
          ! Copy contiguous data back
          do nb=1,nblks
             GFS_Data(nb)%Coupling%ca1(:) = ca1_cpl(nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca2(:) = ca2_cpl(nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca3(:) = ca3_cpl(nb,1:GFS_Control%blksz(nb))
          enddo
       endif

      endif !do_ca

    endif initalize_stochastic_physics

  end subroutine stochastic_physics_wrapper


  subroutine stochastic_physics_wrapper_end (GFS_Control)

  use GFS_typedefs,       only: GFS_control_type, GFS_data_type
  use stochastic_physics, only: finalize_stochastic_physics

  implicit none

  type(GFS_control_type),   intent(inout) :: GFS_Control

  if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. GFS_Control%do_mlp .OR. (GFS_Control%lndp_type > 0) .OR. GFS_Control%do_spp) then
      if (allocated(xlat)) deallocate(xlat)
      if (allocated(xlon)) deallocate(xlon)
      if (GFS_Control%do_sppt) then
         if (allocated(sppt_wts)) deallocate(sppt_wts)
!mlp
      if (GFS_Control%do_mlp) then
          if (allocated ( mlp_pert_tcnv_cpl )) deallocate(mlp_pert_tcnv_cpl)
          if (allocated ( tcnvtend_cpl )) deallocate(tcnvtend_cpl)
          if (allocated ( mlp_pert_qcnv_cpl )) deallocate(mlp_pert_qcnv_cpl)
          if (allocated ( qcnvtend_cpl )) deallocate(qcnvtend_cpl)
          if (allocated ( mlp_pert_ucnv_cpl )) deallocate(mlp_pert_ucnv_cpl)
          if (allocated ( ucnvtend_cpl )) deallocate(ucnvtend_cpl)
          if (allocated ( mlp_pert_vcnv_cpl )) deallocate(mlp_pert_vcnv_cpl)
          if (allocated ( vcnvtend_cpl )) deallocate(vcnvtend_cpl)
!       ! if (allocated ( mlp_pert_tcnv_diag )) deallocate(mlp_pert_tcnv_diag)
!       ! if (allocated ( tcnvtend_diag )) deallocate(tcnvtend_diag)
!       ! if (allocated ( mlp_pert_qcnv_diag )) deallocate(mlp_pert_qcnv_diag)
!       ! if (allocated ( qcnvtend_diag )) deallocate(qcnvtend_diag)
!       ! if (allocated ( mlp_pert_ucnv_diag )) deallocate(mlp_pert_ucnv_diag)
!       ! if (allocated ( ucnvtend_diag )) deallocate(ucnvtend_diag)
!       ! if (allocated ( mlp_pert_vcnv_diag )) deallocate(mlp_pert_vcnv_diag)
!       ! if (allocated ( vcnvtend_diag )) deallocate(vcnvtend_diag)
!
!
          if (allocated ( mlp_pert_tmp_cpl )) deallocate(mlp_pert_tmp_cpl)
          if (allocated ( tmptend_cpl )) deallocate(tmptend_cpl)
          if (allocated ( mlp_pert_qmp_cpl )) deallocate(mlp_pert_qmp_cpl)
          if (allocated ( qmptend_cpl )) deallocate(qmptend_cpl)
!!        if (allocated ( mlp_pert_tmp_diag )) deallocate(mlp_pert_tmp_diag)
!!        if (allocated ( tmptend_diag )) deallocate(tmptend_diag)
!!        if (allocated ( mlp_pert_qmp_diag )) deallocate(mlp_pert_qmp_diag)
!!        if (allocated ( qmptend_diag )) deallocate(qmptend_diag)
!
          if (allocated ( mlp_pert_tpbl_cpl )) deallocate(mlp_pert_tpbl_cpl)
          if (allocated ( tpbltend_cpl )) deallocate(tpbltend_cpl)
          if (allocated ( mlp_pert_qpbl_cpl )) deallocate(mlp_pert_qpbl_cpl)
          if (allocated ( qpbltend_cpl )) deallocate(qpbltend_cpl)
          if (allocated ( mlp_pert_upbl_cpl )) deallocate(mlp_pert_upbl_cpl)
          if (allocated ( upbltend_cpl )) deallocate(upbltend_cpl)
          if (allocated ( mlp_pert_vpbl_cpl )) deallocate(mlp_pert_vpbl_cpl)
          if (allocated ( vpbltend_cpl )) deallocate(vpbltend_cpl)
!!        if (allocated ( mlp_pert_tpbl_diag )) deallocate(mlp_pert_tpbl_diag)
!!        if (allocated ( tpbltend_diag )) deallocate(tpbltend_diag)
!!        if (allocated ( mlp_pert_qpbl_diag )) deallocate(mlp_pert_qpbl_diag)
!!        if (allocated ( qpbltend_diag )) deallocate(qpbltend_diag)
!!        if (allocated ( mlp_pert_upbl_diag )) deallocate(mlp_pert_upbl_diag)
!!        if (allocated ( upbltend_diag )) deallocate(upbltend_diag)
!!        if (allocated ( mlp_pert_vpbl_diag )) deallocate(mlp_pert_vpbl_diag)
!!        if (allocated ( vpbltend_diag )) deallocate(vpbltend_diag)
!
          if (allocated ( mlp_pert_tshalcnv_cpl )) deallocate(mlp_pert_tshalcnv_cpl)
          if (allocated ( tshalcnvtend_cpl )) deallocate(tshalcnvtend_cpl)
          if (allocated ( mlp_pert_qshalcnv_cpl )) deallocate(mlp_pert_qshalcnv_cpl)
          if (allocated ( qshalcnvtend_cpl )) deallocate(qshalcnvtend_cpl)
          if (allocated ( mlp_pert_ushalcnv_cpl )) deallocate(mlp_pert_ushalcnv_cpl)
          if (allocated ( ushalcnvtend_cpl )) deallocate(ushalcnvtend_cpl)
          if (allocated ( mlp_pert_vshalcnv_cpl )) deallocate(mlp_pert_vshalcnv_cpl)
          if (allocated ( vshalcnvtend_cpl )) deallocate(vshalcnvtend_cpl)
!!        if (allocated ( mlp_pert_tshalcnv_diag )) deallocate(mlp_pert_tshalcnv_diag)
!!        if (allocated ( tshalcnvtend_diag )) deallocate(tshalcnvtend_diag)
!!        if (allocated ( mlp_pert_qshalcnv_diag )) deallocate(mlp_pert_qshalcnv_diag)
!!        if (allocated ( qshalcnvtend_diag )) deallocate(qshalcnvtend_diag)
!!        if (allocated ( mlp_pert_ushalcnv_diag )) deallocate(mlp_pert_ushalcnv_diag)
!!        if (allocated ( ushalcnvtend_diag )) deallocate(ushalcnvtend_diag)
!!        if (allocated ( mlp_pert_vshalcnv_diag )) deallocate(mlp_pert_vshalcnv_diag)
!!        if (allocated ( vshalcnvtend_diag )) deallocate(vshalcnvtend_diag)
!mlp pert?
      end if ! mlp
      end if ! sppt
      if (GFS_Control%do_shum) then
         if (allocated(shum_wts)) deallocate(shum_wts)
      end if
      if (GFS_Control%do_skeb) then
         if (allocated(skebu_wts)) deallocate(skebu_wts)
         if (allocated(skebv_wts)) deallocate(skebv_wts)
      end if
      if (GFS_Control%do_spp) then
         if (allocated(spp_wts)) deallocate(spp_wts)
      end if
      if ( GFS_Control%lndp_type == 2 ) then
         lsoil = -999
         if (allocated(sfc_wts)) deallocate(sfc_wts)
      end if
      if (GFS_Control%lndp_type == 2) then
          if (allocated(smc))    deallocate(smc)
          if (allocated(slc))    deallocate(slc)
          if (allocated(stc))    deallocate(stc)
          if (allocated(stype))  deallocate(stype)
          if (allocated(vfrac))  deallocate(vfrac)
          if (allocated(snoalb)) deallocate(snoalb)
          if (allocated(alnsf))  deallocate(alnsf)
          if (allocated(alnwf))  deallocate(alnwf)
          if (allocated(semis))  deallocate(semis)
          if (allocated(zorll))  deallocate(zorll)
      endif
      call finalize_stochastic_physics()
   endif
   if(GFS_Control%do_ca)then
        if(GFS_Control%ca_sgs)then
           deallocate(sst         )
           deallocate(lmsk        )
           deallocate(lake        )
           deallocate(condition   )
           deallocate(ca_deep_cpl )
           deallocate(ca_turb_cpl )
           deallocate(ca_shal_cpl )
        endif
        if(GFS_Control%ca_global)then
            deallocate(ca1_cpl )
            deallocate(ca2_cpl )
            deallocate(ca3_cpl )
        endif
   endif
  end subroutine stochastic_physics_wrapper_end

end module stochastic_physics_wrapper_mod
