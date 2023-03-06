subroutine read_mtg_irs(mype,val_irs,ithin,isfcalc,rmesh,jsatid,gstime,&
     infile,lunout,obstype,nread,ndata,nodata,twind,sis,&
     mype_root,mype_sub,npe_sub,mpi_comm_sub,nobs, &
     nrec_start,nrec_start_ears,nrec_start_db,dval_use)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_mtg_irs                read bufr format mtg_irs data
! prgmmr :   derber          org: np20                date: 2022-02-33
!
! abstract:  This routine reads BUFR format radiance 
!            files.  Based on read_iasi Optionally, the data are thinned to 
!            a specified resolution using simple quality control checks.
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!   2022-02-33  derber  - read mtg_irs data in bufr format
!
!   input argument list:
!     mype     - mpi task id
!     val_irs  - weighting factor applied to super obs
!     ithin    - flag to thin data
!     isfcalc  - when set to one, calculate surface characteristics using
!                method that accounts for the size/shape of the fov. 
!                when not one, calculate surface characteristics using
!                bilinear interpolation.
!     rmesh    - thinning mesh size (km)
!     jsatid   - satellite id
!     gstime   - analysis time in minutes from reference date
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window (hours)
!     sis      - sensor/instrument/satellite indicator
!     mype_root - "root" task for sub-communicator
!     mype_sub - mpi task id within sub-communicator
!     npe_sub  - number of data read tasks
!     mpi_comm_sub - sub-communicator for data read
!     nrec_start - first subset with useful information
!     nrec_start_ears - first ears subset with useful information
!     nrec_start_db - first db subset with useful information
!
!   output argument list:
!     nread    - number of BUFR MTG_IRS observations read
!     ndata    - number of BUFR MTG_IRS profiles retained for further processing
!     nodata   - number of BUFR MTG_IRS observations retained for further processing
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
! Use modules
  use kinds, only: r_kind,r_double,i_kind
  use satthin, only: super_val,itxmax,makegrids,map2tgrid,destroygrids, &
      finalcheck,checkob,score_crit
  use satthin, only: radthin_time_info,tdiff2crit
  use obsmod,  only: time_window_max
  use radinfo, only:iuse_rad,nuchan,nusis,jpch_rad,crtm_coeffs_path,use_edges, &
      radedge1,radedge2,radstart,radstep
  use crtm_module, only: success, &
      crtm_kind => fp
  use crtm_planck_functions, only: crtm_planck_temperature
  use crtm_spccoeff, only: sc,crtm_spccoeff_load,crtm_spccoeff_destroy
  use gridmod, only: diagnostic_reg,regional,nlat,nlon,&
      tll2xy,txy2ll,rlats,rlons
  use constants, only: zero,deg2rad,rad2deg,r60inv,one,ten,r100
  use gsi_4dvar, only: l4dvar,l4densvar,iwinbgn,winlen
  use calc_fov_crosstrk, only: instrument_init, fov_check, fov_cleanup
  use deter_sfc_mod, only: deter_sfc,deter_sfc_fov
  use obsmod, only: bmiss
  use gsi_nstcouplermod, only:nst_gsi,nstinfo
  use gsi_nstcouplermod, only: gsi_nstcoupler_skindepth, gsi_nstcoupler_deter
  use mpimod, only: npe
  use gsi_io, only: verbose
! use radiance_mod, only: rad_obs_type

  implicit none

! BUFR format for MTG_IRS SPOT 
! Input variables
  integer(i_kind)  ,intent(in   ) :: mype,nrec_start,nrec_start_ears,nrec_start_db
  integer(i_kind)  ,intent(in   ) :: ithin
  integer(i_kind)  ,intent(inout) :: isfcalc
  integer(i_kind)  ,intent(in   ) :: lunout
  integer(i_kind)  ,intent(in   ) :: mype_root
  integer(i_kind)  ,intent(in   ) :: mype_sub
  integer(i_kind)  ,intent(in   ) :: npe_sub
  integer(i_kind)  ,intent(in   ) :: mpi_comm_sub  
  character(len=*), intent(in   ) :: infile
  character(len=10),intent(in   ) :: jsatid
  character(len=*), intent(in   ) :: obstype
  character(len=20),intent(in   ) :: sis
  real(r_kind)     ,intent(in   ) :: twind
  real(r_kind)     ,intent(inout) :: val_irs
  real(r_kind)     ,intent(in   ) :: gstime
  real(r_kind)     ,intent(in   ) :: rmesh
  logical          ,intent(in   ) :: dval_use

! Output variables
  integer(i_kind)  ,intent(inout) :: nread
  integer(i_kind),dimension(npe)  ,intent(inout) :: nobs
  integer(i_kind)  ,intent(  out) :: ndata,nodata
  

! BUFR file sequencial number
!  character(len=512)  :: table_file
  integer(i_kind)     :: lnbufr = 10

! Variables for BUFR IO    
  real(r_double) :: crchn_reps
  real(r_double),dimension(13) :: allspot
  real(r_double),allocatable,dimension(:,:) :: allchan
  real(r_double),dimension(3,10):: cscale 
  real(r_double),dimension(7):: cloud_frac
  integer(i_kind) :: bufr_size
  
  real(r_kind)      :: step, start,step_adjust
  character(len=8)  :: subset
  character(len=4)  :: senname
  character(len=80) :: allspotlist
  character(len=40) :: infile2
  integer(i_kind)   :: iret,ireadsb,ireadmg,irec,next, nrec_startx
  integer(i_kind),allocatable,dimension(:) :: nrec


! Work variables for time
  integer(i_kind)   :: idate
  integer(i_kind)   :: idate5(5)
  real(r_kind)      :: sstime, tdiff, t4dv
  integer(i_kind)   :: nmind


! Other work variables
  real(r_kind)     :: piece
  real(r_kind)     :: rsat, dlon, dlat
  real(r_kind)     :: dlon_earth,dlat_earth,dlon_earth_deg,dlat_earth_deg
  real(r_kind)     :: lza, lzaest,sat_height_ratio
  real(r_kind)     :: pred, crit1, dist1
  real(r_kind)     :: sat_zenang
  real(crtm_kind)  :: radiance
  real(r_kind)     :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr
  real(r_kind)     :: zob,tref,dtw,dtc,tz_tr
  real(r_kind),dimension(0:4) :: rlndsea
  real(r_kind),dimension(0:3) :: sfcpct
  real(r_kind),dimension(0:3) :: ts
  real(r_kind),dimension(10) :: sscale
  real(crtm_kind),allocatable,dimension(:) :: temperature
  real(r_kind),allocatable,dimension(:) :: scalef
  real(r_kind),allocatable,dimension(:,:):: data_all
  real(r_kind) cdist,disterr,disterrmax,dlon00,dlat00

  logical          :: outside,iuse,assim,valid
  logical          :: mtg_irs,quiet,cloud_info

  integer(i_kind)  :: ifov, instr, iscn, ioff, sensorindex
  integer(i_kind)  :: i, j, l, iskip, ifovn, bad_line, ksatid, kidsat
  integer(i_kind)  :: nreal, isflg
  integer(i_kind)  :: itx, k, nele, itt, n
  integer(i_kind):: iexponent,maxinfo, bufr_nchan
  integer(i_kind):: idomsfc(1)
  integer(i_kind):: ntest
  integer(i_kind):: error_status, irecx,ierr
  integer(i_kind):: radedge_min, radedge_max
  integer(i_kind)   :: subset_start, subset_end, satinfo_nchan, sc_chan, bufr_chan
  integer(i_kind)   :: sfc_channel_index
  integer(i_kind),allocatable, dimension(:) :: channel_number, sc_index, bufr_index
  integer(i_kind),allocatable, dimension(:) :: bufr_chan_test
  character(len=20),dimension(1):: sensorlist


! Set standard parameters
  character(8),parameter:: fov_flag="crosstrk"
  integer(i_kind),parameter:: sfc_channel=1271
  integer(i_kind),parameter:: ichan=-999  ! fov-based surface code is not channel specific for mtg_irs
  real(r_kind),parameter:: expansion=one         ! exansion factor for fov-based surface code.
                                                 ! use one for ir sensors.
  real(r_kind),parameter:: R90    =  90._r_kind
  real(r_kind),parameter:: R360   = 360._r_kind
  real(r_kind),parameter:: tbmin  = 50._r_kind
  real(r_kind),parameter:: tbmax  = 550._r_kind
  real(r_kind),parameter:: earth_radius = 6371000._r_kind
  integer(i_kind),parameter :: ilon = 3
  integer(i_kind),parameter :: ilat = 4
  real(r_kind)    :: ptime,timeinflat,crit0
  real(r_kind),dimension(8,2) :: data2
  integer(i_kind) :: ithin_time,n_tbin,it_mesh,jstart
  character(80):: datastr2='CHSF SMRA SSDR SQFA RRIB CONFLG CSSQ CDSQ'
  logical print_verbose

  print_verbose=.false.
  if(verbose)print_verbose=.true.

! Initialize variables
  maxinfo    =  47
  disterrmax=zero
  ntest=0
  if(dval_use) maxinfo=maxinfo+2
  nreal  = maxinfo + nstinfo

  ndata = 0
  nodata = 0
  mtg_irs = obstype == 'mtg_irs'

  bad_line=-1

  if (nst_gsi > 0 ) then
    call gsi_nstcoupler_skindepth(obstype, zob)         ! get penetration depth (zob) for the obstype
  endif

  if(jsatid == 'mtg1')kidsat=72
 
!  write(6,*)'READ_MTG_IRS: mype, mype_root,mype_sub, npe_sub,mpi_comm_sub', &
!          mype, mype_root,mype_sub,mpi_comm_sub

  radedge_min = 0
  radedge_max = 100000

! Find the mtg_irs offset in the jpch_rad list.  This is for the iuse flag
! and count the number of cahnnels in the satinfo file 
  ioff=jpch_rad
  subset_start = 0
  subset_end = 0
  assim = .false.
  do i=1,jpch_rad
     if (trim(nusis(i))==trim(sis)) then
        ioff = min(ioff,i)   ! mtg_irs offset
        if (subset_start == 0) then
          step  = radstep(i)
          start = radstart(i)
          if (radedge1(i)/=-1 .and. radedge2(i)/=-1) then
             radedge_min=min(radedge1(i),radedge_min)
             radedge_max=max(radedge2(i),radedge_max)
          end if
          subset_start = i
        endif
        if (iuse_rad(i) > 0) assim = .true.  ! Are any of the MTG_IRS channels being used?
        subset_end = i 
     endif
  end do 
  satinfo_nchan = subset_end - subset_start + 1
  allocate(channel_number(satinfo_nchan))
  allocate(sc_index(satinfo_nchan))
  allocate(bufr_index(satinfo_nchan)) 
  ioff = ioff - 1 

  step_adjust = 0.625_r_kind
! If all channels of a given sensor are set to monitor or not
! assimilate mode (iuse_rad<1), reset relative weight to zero.
! We do not want such observations affecting the relative
! weighting between observations within a given thinning group.
  if (.not.assim) val_irs=zero

  if (mype_sub==mype_root)write(6,*)'READ_MTG_IRS:  mtg_irs offset ',ioff

  senname = 'MTG_IRS'
  
  allspotlist= &
   'SAID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH SAZA BEARAZ SOZA SOLAZI'

! load spectral coefficient structure  
  quiet=.not. verbose
  sensorlist(1)=sis
  if( crtm_coeffs_path /= "" ) then
     if(mype_sub==mype_root .and. print_verbose) write(6,*)'READ_MTG_IRS: crtm_spccoeff_load() on path "'//trim(crtm_coeffs_path)//'"'
     error_status = crtm_spccoeff_load(sensorlist,&
        File_Path = crtm_coeffs_path,quiet=quiet )
  else
     error_status = crtm_spccoeff_load(sensorlist,quiet=quiet)
  endif

  if (error_status /= success) then
     write(6,*)'READ_MTG_IRS:  ***ERROR*** crtm_spccoeff_load error_status=',error_status,&
        '   TERMINATE PROGRAM EXECUTION'
     call stop2(71)
  endif

! Find the channels being used (from satinfo file) in the spectral coef. structure.
  do i=subset_start,subset_end
     channel_number(i -subset_start +1) = nuchan(i)
  end do
  sc_index(:) = 0
  satinfo_chan: do i=1,satinfo_nchan
     spec_coef: do l=1,sc(1)%n_channels
        if ( channel_number(i) == sc(1)%sensor_channel(l) ) then
           sc_index(i) = l
           exit spec_coef
        endif
     end do spec_coef
  end do  satinfo_chan

!  find MTG_IRS sensorindex
  sensorindex = 0
  if ( sc(1)%sensor_id(1:4) == 'mtg_irs' .or. sc(1)%sensor_id(1:4) == 'MTG_IRS') then
     sensorindex = 1
  else
     write(6,*)'READ_MTG_IRS: sensorindex not set  NO MTG_IRS DATA USED'
     write(6,*)'READ_MTG_IRS: We are looking for ', sc(1)%sensor_id, '   TERMINATE PROGRAM EXECUTION'
     call stop2(71)
  end if

! Calculate parameters needed for FOV-based surface calculation.
  if (isfcalc==1)then
     instr=18  !*****NEED TO FIND CORRECT VALUE FOR THIS******
     call instrument_init(instr, jsatid, expansion, valid)
     if (.not. valid) then
        if (assim) then 
           write(6,*)'READ_MTG_IRS:  ***ERROR*** IN SETUP OF FOV-SFC CODE. STOP'
           call stop2(71)
        else
           call fov_cleanup
           isfcalc = 0
           write(6,*)'READ_MTG_IRS:  ***ERROR*** IN SETUP OF FOV-SFC CODE'
        endif
     endif
  endif

  if (isfcalc==1)then
     rlndsea = zero
  else
     rlndsea(0) = zero                       
     rlndsea(1) = 10._r_kind
     rlndsea(2) = 15._r_kind
     rlndsea(3) = 10._r_kind
     rlndsea(4) = 30._r_kind
  endif

  call radthin_time_info(obstype, jsatid, sis, ptime, ithin_time)
  if( ptime > 0.0_r_kind) then
     n_tbin=nint(2*time_window_max/ptime)
  else
     n_tbin=1
  endif
! Make thinning grids
  call makegrids(rmesh,ithin,n_tbin=n_tbin)

! Allocate arrays to hold data
! The number of channels in obtained from the satinfo file being used.
  nele=nreal+satinfo_nchan
  allocate(data_all(nele,itxmax),nrec(itxmax))
  allocate(temperature(1))   ! dependent on # of channels in the bufr file
  allocate(allchan(3,1))     ! actual values set after ireadsb
  allocate(bufr_chan_test(1))! actual values set after ireadsb
  allocate(scalef(1))

! Big loop to read data file
  next=0
  irec=0
  nrec=999999

  nrec_startx=nrec_start
  infile2=trim(infile)         ! Set bufr subset names based on type of data to read

!    Open BUFR file
  open(lnbufr,file=trim(infile2),form='unformatted',status='old',iostat=ierr)

  if(ierr /= 0) then 
    if(mype == 0) then
       write(6,*) ' MTG IRS file ',infile2,' not available '
    end if
    return
  end if
! Open BUFR table
  call openbf(lnbufr,'IN',lnbufr)
  call datelen(10)

  irecx = 0
  read_subset: do while(ireadmg(lnbufr,subset,idate)>=0)
     irecx = irecx + 1
     if(irecx < nrec_startx) cycle read_subset
     irec = irec + 1
     next=next+1
     if(next == npe_sub)next=0
     if(next /= mype_sub)cycle read_subset

     read_loop: do while (ireadsb(lnbufr)==0)

!       Get the size of the channels and radiance (allchan) array
        call ufbint(lnbufr,crchn_reps,1,1,iret,'(CHNM)')
        bufr_nchan = int(crchn_reps)

        bufr_size = size(temperature,1)
        if ( bufr_size /= bufr_nchan ) then ! Re-allocation if number of channels has changed
!          Allocate the arrays needed for the channel and radiance array
           deallocate(temperature,allchan,bufr_chan_test,scalef)
           allocate(temperature(bufr_nchan))   ! dependent on # of channels in the bufr file
           allocate(allchan(3,bufr_nchan))
           allocate(bufr_chan_test(bufr_nchan))
           allocate(scalef(bufr_nchan))
           bufr_chan_test(:)=0
        endif       !  allocation if

        call ufbint(lnbufr,allspot,13,1,iret,allspotlist)
        if(iret /= 1) cycle read_loop

!       Extract satellite id.  If not the one we want, read next subset
        ksatid=nint(allspot(1))
        if(ksatid /= kidsat) cycle read_loop

!       Check observing position
        dlat_earth = allspot(8)   ! latitude
        dlon_earth = allspot(9)   ! longitude
        if( abs(dlat_earth) > R90  .or. abs(dlon_earth) > R360 .or. &
           (abs(dlat_earth) == R90 .and. dlon_earth /= ZERO) )then
           write(6,*)'READ_MTG_IRS:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
              ' STRANGE OBS POINT (LAT,LON):', dlat_earth, dlon_earth
           cycle read_loop
        endif

!       Retrieve observing position
        if(dlon_earth >= R360)then
           dlon_earth = dlon_earth - R360
        else if(dlon_earth < ZERO)then
           dlon_earth = dlon_earth + R360
        endif

        dlat_earth_deg = dlat_earth
        dlon_earth_deg = dlon_earth
        dlat_earth = dlat_earth * deg2rad
        dlon_earth = dlon_earth * deg2rad

!       If regional, map obs lat,lon to rotated grid.
        if(regional)then

!          Convert to rotated coordinate.  dlon centered on 180 (pi),
!          so always positive for limited area
           call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
           if(diagnostic_reg) then
              call txy2ll(dlon,dlat,dlon00,dlat00)
              ntest=ntest+1
              cdist=sin(dlat_earth)*sin(dlat00)+cos(dlat_earth)*cos(dlat00)* &
                   (sin(dlon_earth)*sin(dlon00)+cos(dlon_earth)*cos(dlon00))
              cdist=max(-one,min(cdist,one))
              disterr=acos(cdist)*rad2deg
              disterrmax=max(disterrmax,disterr)
           end if

!          Check to see if in domain.  outside=.true. if dlon_earth,
!          dlat_earth outside domain, =.false. if inside
           if(outside) cycle read_loop

!       Global case 
        else
           dlat = dlat_earth
           dlon = dlon_earth
           call grdcrd1(dlat,rlats,nlat,1)
           call grdcrd1(dlon,rlons,nlon,1)
        endif

!       Check obs time
        idate5(1) = nint(allspot(2)) ! year
        idate5(2) = nint(allspot(3)) ! month
        idate5(3) = nint(allspot(4)) ! day
        idate5(4) = nint(allspot(5)) ! hour
        idate5(5) = nint(allspot(6)) ! minute

        if( idate5(1) < 1900 .or. idate5(1) > 3000 .or. &
            idate5(2) < 1    .or. idate5(2) >   12 .or. &
            idate5(3) < 1    .or. idate5(3) >   31 .or. &
            idate5(4) <0     .or. idate5(4) >   24 .or. &
            idate5(5) <0     .or. idate5(5) >   60 )then

            write(6,*)'READ_MTG_IRS:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
                 ' STRANGE OBS TIME (YMDHM):', idate5(1:5)
            cycle read_loop

        endif

!       Retrieve obs time
        call w3fs21(idate5,nmind)
        t4dv = (real(nmind-iwinbgn,r_kind) + real(allspot(7),r_kind)*r60inv)*r60inv ! add in seconds
        sstime = real(nmind,r_kind) + real(allspot(7),r_kind)*r60inv ! add in seconds
        tdiff = (sstime - gstime)*r60inv

        if (l4dvar.or.l4densvar) then
           if (t4dv<zero .OR. t4dv>winlen) cycle read_loop
        else
           if (abs(tdiff)>twind) cycle read_loop
        endif

!       Increment nread counter by satinfo_nchan
        nread = nread + satinfo_nchan

        crit0 = 0.01_r_kind
        timeinflat=6.0_r_kind
        call tdiff2crit(tdiff,ptime,ithin_time,timeinflat,crit0,crit1,it_mesh)
        call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis,it_mesh=it_mesh)

        if(.not. iuse)cycle read_loop

!       Observational info
        sat_zenang  = allspot(10)-90._r_kind            ! satellite zenith angle

!       Check  satellite zenith angle (SAZA)
        if(abs(sat_zenang) > 90._r_kind) then
           write(6,*)'READ_MTG_IRS:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
              ' STRANGE OBS INFO(FOVN,SLNM,SAZA,BEARAZ):', ifov, iscn, allspot(10),allspot(11)
           cycle read_loop
        endif

!       "Score" observation.  We use this information to identify "best" obs
!       Locate the observation on the analysis grid.  Get sst and land/sea/ice
!       mask.  
!       isflg    - surface flag
!             0 sea
!             1 land
!             2 sea ice
!             3 snow
!             4 mixed 

!       When using FOV-based surface code, must screen out obs with bad fov numbers.
        if (isfcalc == 1) then
           call fov_check(ifov,instr,ichan,valid)
           if (.not. valid) cycle read_loop

!       When isfcalc is set to one, calculate surface fields using size/shape of fov.
!       Otherwise, use bilinear interpolation.

           call deter_sfc_fov(fov_flag,ifov,instr,ichan,real(allspot(11),r_kind),dlat_earth_deg, &
                           dlon_earth_deg,expansion,t4dv,isflg,idomsfc(1), &
                           sfcpct,vfr,sty,vty,stp,sm,ff10,sfcr,zz,sn,ts,tsavg)
        else
           call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg,idomsfc(1),sfcpct, &
              ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)
        endif

!       Set common predictor parameters
        crit1 = crit1 + rlndsea(isflg)
 
        call checkob(one,crit1,itx,iuse)
        if(.not. iuse)cycle read_loop

!  Read diagnostic information
!        (*,1) and (*,2) represent band 1 (longwave) and band 2(midwave)
!        respectively
!        (1,*) scale factor for imagery ::  rvalues = rvalues/10.0**iscale
!        (2,*) mean value for imagery
!        (2,*) standard deviation for imagery
!        (4,*) scoreQuantization Factor (same for bands 1 and 2)
!        (5,*) global_pcr_scores
!        (6,*) global_pcrs_quality
!              IAND(global_pcrs_quality_lw(columns_array(i),rows_array(i)),1)
!              (0=valid, 1=invalid, 15=missing)
!        (7,*) spatial_sample_quality
!              spatial_sample_quality, 13-bit (bit 1 defined here as most significant)
!              Bit  Meaning
!              1-4 reserved
!              5   solar_straylight_correction_warning
!              6   solar_straylight_warning
!              7   noisy_detector_sample_warning
!              8   undersaturated_detector_sample_warning
!              9   saturated_detector_sample_warning
!              10   dust
!              11   cloudy
!              12   limb_view
!              13   space_view

!        (8,*) detector_sample quality
!              detector_sample_quality, 4-bit (bit 1 defined here as most significant)
!              Bit  Meaning
!              1   excluded_detector_sample
!              2   noisy_detector_sample
!              3   undersaturated_detector_sample
!              4   saturated_detector_sample
 
        call ufbseq(lnbufr,data2,8,2,iret,datastr2)

!  Add qc based on data2 here!!!

!       Read MTG_IRS channel number(CHNM) and radiance (SRAD) and band (1
!       mid-wave, 2 long wave)
        call ufbseq(lnbufr,allchan,3,bufr_nchan,iret,'CHNM SRAD TOBD')

!       Coordinate bufr channels with satinfo file channels
!       If this is the first time or a change in the bufr channels is detected, sync with satinfo file
        if (ANY(int(allchan(1,:)) /= bufr_chan_test(:))) then
           sfc_channel_index = 0
           bufr_index(:) = 0
           bufr_chans: do l=1,bufr_nchan
              bufr_chan_test(l) = int(allchan(1,l))                      ! Copy this bufr channel selection into array for comparison to next profile
              satinfo_chans: do i=1,satinfo_nchan                        ! Loop through sensor (mtg_irs) channels in the satinfo file
                 if ( channel_number(i) == int(allchan(1,l)) ) then      ! Channel found in both bufr and satinfo file
                    bufr_index(i) = l
                    if ( channel_number(i) == sfc_channel) sfc_channel_index = l
                    exit satinfo_chans                                   ! go to next bufr channel
                 endif
              end do  satinfo_chans
           end do bufr_chans
        endif

!       if (sfc_channel_index == 0) then
!         write(6,*)'READ_MTG_IRS: ***ERROR*** SURFACE CHANNEL USED FOR QC WAS NOT FOUND'
!         cycle read_loop
!       endif

!$omp parallel do schedule(dynamic,1) private(i,sc_chan,bufr_chan,radiance)
        channel_loop: do i=1,satinfo_nchan
           bufr_chan = bufr_index(i)
           if (bufr_chan /= 0 ) then
!          check that channel number is within reason
             if (( allchan(2,bufr_chan) > zero .and. allchan(2,bufr_chan) < 99999._r_kind)) then  ! radiance bounds
               radiance = allchan(2,bufr_chan)*scalef(bufr_chan)
               sc_chan = sc_index(i)
               call crtm_planck_temperature(sensorindex,sc_chan,radiance,temperature(bufr_chan))
             else
                temperature(bufr_chan) = tbmin
             endif
           else
              temperature(bufr_chan) = tbmin
           end if
        end do channel_loop

!       Check for reasonable temperature values
        iskip = 0
        skip_loop: do i=1,satinfo_nchan
           if ( bufr_index(i) == 0 ) cycle skip_loop
           bufr_chan = bufr_index(i)
           if(temperature(bufr_chan) <= tbmin .or. temperature(bufr_chan) > tbmax ) then
              temperature(bufr_chan) = min(tbmax,max(tbmin,temperature(bufr_chan)))
              if(iuse_rad(ioff+i) >= 0)iskip = iskip + 1
           endif
        end do skip_loop

        if(iskip > 0)then
           if(print_verbose)write(6,*) ' READ_MTG_IRS : iskip > 0 ',iskip
           cycle read_loop 
        end if

!       crit1=crit1 + ten*float(iskip)

!       If the surface channel exists (~960.0 cm-1) and the AVHRR cloud information is missing, use an
!       estimate of the surface temperature to determine if the profile may be clear.
!       if (.not. cloud_info) then
!          pred = tsavg*0.98_r_kind - temperature(sfc_channel_index)
!          pred = max(pred,zero)
!          crit1=crit1 + pred
!       endif

!       Map obs to grids
        if (pred == zero) then
           call finalcheck(dist1,crit1,itx,iuse)
        else
           call finalcheck(one,crit1,itx,iuse)
        endif
        if(.not. iuse)cycle read_loop

!
!       interpolate NSST variables to Obs. location and get dtw, dtc, tz_tr
!
        if ( nst_gsi > 0 ) then
           tref  = ts(0)
           dtw   = zero
           dtc   = zero
           tz_tr = one
           if ( sfcpct(0) > zero ) then
              call gsi_nstcoupler_deter(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
           endif
        endif

        rsat=allspot(1) 
        data_all(1,itx) = rsat                      ! satellite ID 
        data_all(2,itx) = t4dv                      ! time diff (obs-anal) (hrs)
        data_all(3,itx) = dlon                      ! grid relative longitude
        data_all(4,itx) = dlat                      ! grid relative latitude
        data_all(5,itx) = sat_zenang*deg2rad        ! satellite zenith angle (rad)
        data_all(6,itx) = allspot(11)               ! satellite azimuth angle (deg)
        data_all(7,itx) = lza                       ! look angle (rad)
        data_all(8,itx) = ifovn                     ! fov number
        data_all(9,itx) = allspot(12)               ! solar zenith angle (deg)
        data_all(10,itx)= allspot(13)               ! solar azimuth angle (deg)
        data_all(11,itx) = sfcpct(0)                ! sea percentage of
        data_all(12,itx) = sfcpct(1)                ! land percentage
        data_all(13,itx) = sfcpct(2)                ! sea ice percentage
        data_all(14,itx) = sfcpct(3)                ! snow percentage
        data_all(15,itx)= ts(0)                     ! ocean skin temperature
        data_all(16,itx)= ts(1)                     ! land skin temperature
        data_all(17,itx)= ts(2)                     ! ice skin temperature
        data_all(18,itx)= ts(3)                     ! snow skin temperature
        data_all(19,itx)= tsavg                     ! average skin temperature
        data_all(20,itx)= vty                       ! vegetation type
        data_all(21,itx)= vfr                       ! vegetation fraction
        data_all(22,itx)= sty                       ! soil type
        data_all(23,itx)= stp                       ! soil temperature
        data_all(24,itx)= sm                        ! soil moisture
        data_all(25,itx)= sn                        ! snow depth
        data_all(26,itx)= zz                        ! surface height
        data_all(27,itx)= idomsfc(1) + 0.001_r_kind ! dominate surface type
        data_all(28,itx)= sfcr                      ! surface roughness
        data_all(29,itx)= ff10                      ! ten meter wind factor
        data_all(30,itx)= dlon_earth_deg            ! earth relative longitude (degrees)
        data_all(31,itx)= dlat_earth_deg            ! earth relative latitude (degrees)
        data_all(32,itx)= data2(1,1)                  ! qc variable CHSF channnel group 1
        data_all(33,itx)= data2(1,2)                  ! qc variable CHSF channnel group 2
        data_all(34,itx)= data2(2,1)                  ! qc variable SMRA channnel group 1
        data_all(35,itx)= data2(2,2)                  ! qc variable SMRA channnel group 2
        data_all(36,itx)= data2(3,1)                  ! qc variable SSDR channnel group 1
        data_all(37,itx)= data2(3,2)                  ! qc variable SSDR channnel group 2
        data_all(38,itx)= data2(4,1)                  ! qc variable SQFA channnel group 1
        data_all(39,itx)= data2(4,2)                  ! qc variable SQFA channnel group 2
        data_all(40,itx)= data2(5,1)                  ! qc variable RRIB channnel group 1
        data_all(41,itx)= data2(5,2)                  ! qc variable RRIB channnel group 2
        data_all(42,itx)= data2(6,1)                  ! qc variable CONFLG channnel group 1
        data_all(43,itx)= data2(6,2)                  ! qc variable CONFLG channnel group 2
        data_all(44,itx)= data2(7,1)                  ! qc variable CSSQ channnel group 1
        data_all(45,itx)= data2(7,2)                  ! qc variable CSSQ channnel group 2
        data_all(46,itx)= data2(8,1)                  ! qc variable CDSQ channnel group 1
        data_all(47,itx)= data2(8,2)                  ! qc variable CDSQ channnel group 2

        if(dval_use) then
           data_all(maxinfo-1,itx)= val_irs
           data_all(maxinfo,itx)= itt
        end if


        if ( nst_gsi > 0 ) then
           data_all(maxinfo+1,itx) = tref         ! foundation temperature
           data_all(maxinfo+2,itx) = dtw          ! dt_warm at zob
           data_all(maxinfo+3,itx) = dtc          ! dt_cool at zob
           data_all(maxinfo+4,itx) = tz_tr        ! d(Tz)/d(Tr)
        endif

!       Put satinfo defined channel temperatures into data array
        do l=1,satinfo_nchan
           i = bufr_index(l)
           data_all(l+nreal,itx) = temperature(i)   ! brightness temerature
        end do
        nrec(itx)=irec

     enddo read_loop

  enddo read_subset

  call closbf(lnbufr)
  close(lnbufr)

  deallocate(temperature, allchan, bufr_chan_test,scalef)
  deallocate(channel_number,sc_index)
  deallocate(bufr_index)
! deallocate crtm info
  error_status = crtm_spccoeff_destroy()
  if (error_status /= success) &
    write(6,*)'OBSERVER:  ***ERROR*** crtm_destroy error_status=',error_status

! If multiple tasks read input bufr file, allow each tasks to write out
! information it retained and then let single task merge files together

  call combine_radobs(mype_sub,mype_root,npe_sub,mpi_comm_sub,&
     nele,itxmax,nread,ndata,data_all,score_crit,nrec)

! Allow single task to check for bad obs, update superobs sum,
! and write out data to scratch file for further processing.
  if (mype_sub==mype_root.and.ndata>0) then

!    Identify "bad" observation (unreasonable brightness temperatures).
!    Update superobs sum according to observation location

     do n=1,ndata
        do i=1,satinfo_nchan
           if(data_all(i+nreal,n) > tbmin .and. &
              data_all(i+nreal,n) < tbmax)nodata=nodata+1
        end do
     end do

     if(dval_use .and. assim)then
        do n=1,ndata
          itt=nint(data_all(33,n))
          super_val(itt)=super_val(itt)+val_irs
        end do
     end if

!    Write final set of "best" observations to output file
     call count_obs(ndata,nele,ilat,ilon,data_all,nobs)
     write(lunout) obstype,sis,nreal,satinfo_nchan,ilat,ilon
     write(lunout) ((data_all(k,n),k=1,nele),n=1,ndata)
  
  endif


  deallocate(data_all,nrec) ! Deallocate data arrays
  call destroygrids    ! Deallocate satthin arrays

! Deallocate arrays and nullify pointers.
  if(isfcalc == 1) call fov_cleanup

  if(diagnostic_reg .and. ntest > 0 .and. mype_sub==mype_root) &
     write(6,*)'READ_MTG_IRS:  mype,ntest,disterrmax=',&
        mype,ntest,disterrmax
  
  return
end subroutine read_mtg_irs
