module m_namelists

  implicit none

  ! Namelist limits
  integer, parameter :: rowmax = 200, colmax = 3, lenmax = 200
  integer, parameter :: slenmax = 1024, smax = 10

  ! System namelist
  character(len=slenmax), save :: ibasedir, obasedir, tabledir, griddata
  logical, save :: createsubdirs, forcefilescan, verbose
  namelist /sys/ ibasedir, obasedir, tabledir, griddata, createsubdirs, &
    forcefilescan, verbose

  ! Model namelist
  character(len=slenmax), save :: model_id, institute_id
  character(len=slenmax), dimension(smax), save :: institution, source, &
    references, contact
  character(len=(slenmax+1)*smax), save :: institution1, source1, references1, contact1
  character(len=slenmax), save :: tagoyr, tagoyrbgc, tagomon, tagomonbgc, tagoday, &
    tagodaybgc, tagimon, tagiday, tagamon, tagaday, taga6hr, taga6hri, taga3hr, &
    taga3hri, taglmon, taglday, tagl3hr, tagl3hri, taglyr
  character(len=slenmax), save :: secindexfile, ocngridfile, ocninitfile, ocnmertfile, &
    rhotablesuff, atmgridfile, ocnregnfile
  logical, save :: linebreaks
  character(len=slenmax), save :: parent_source_id, coordtable, namelist_file_json, &
    atmgrid, atmgrid_label, atmgrid_resolution, ocngrid, ocngrid_label, ocngrid_resolution, &
    icegrid, icegrid_label, icegrid_resolution, lndgrid, lndgrid_label, lndgrid_resolution
  namelist /model/ parent_source_id, atmgrid, atmgrid_label, atmgrid_resolution, &
    ocngrid, ocngrid_label, ocngrid_resolution, icegrid, icegrid_label, icegrid_resolution, &
    lndgrid, lndgrid_label, lndgrid_resolution, model_id, source, institution, institute_id, &
    references, contact, tagoyr, tagoyrbgc, tagomon, tagomonbgc, tagoday, tagodaybgc, tagimon, &
    tagiday, tagamon, tagaday, taga6hr, taga6hri, taga3hr, taga3hri, taglyr, taglmon, taglday, &
    tagl3hr, tagl3hri, rhotablesuff, coordtable, secindexfile, atmgridfile, ocngridfile, &
    ocninitfile, ocnmertfile, ocnregnfile, linebreaks

  ! Experiment namelist
  character(len=slenmax), save :: casename, experiment_id, parent_experiment_id, &
    parent_experiment_rip, isubdir, osubdir, membertag
  character(len=slenmax), dimension(smax), save :: history, comment, forcing
  character(len=(slenmax+1)*smax), save :: history1, comment1, forcing1
  integer, save :: realization, exprefyear, year1, yearn, month1, monthn
  real(kind=8), save :: branch_time
  logical, save :: dry_run, plevdummy, readdummy, add_fill_day, do_fx, do_oyr, do_oyrbgc, &
    do_amon, do_omon, do_omonbgc, do_oimon, do_lmon, do_limon, do_aero, do_day, do_6hrlev, &
    do_6hrlevi, do_6hrplev, do_3hr, do_3hri, do_2d, do_3d, do_xd, do_bgc, newcolumnorder, &
    scanallfiles
  integer, save :: physics_version = 1, initialization_method = 1, ivnmpos = 1, ovnmpos = 2
  logical, save :: do_ofx, do_6hrPlevPt, do_AERday, do_AERhr, do_AERmonZ, do_CF3hr, do_CFday, &
    do_CFmon, do_CFsubhr, do_E1hrClimMon, do_E1hr, do_E3hr, do_E3hrPt, do_E6hrZ, do_Eday, &
    do_EdayZ, do_Efx, do_Emon, do_EmonZ, do_Esubhr, do_Eyr, do_Oclim, do_Oday, do_Odaybgc, &
    do_Odec, do_COfx, do_SIday
  character(len=slenmax), save :: activity_id, parent_variant_label, parent_mip_era, mip_era, &
    sub_experiment_id, parent_sub_experiment, parent_activity_id, branch_method, &
    parent_time_units, tracking_prefix, variant_label, source_type
  real(kind=8), save :: branch_time_in_child, branch_time_in_parent
  character(len=slenmax), save :: forcing_index, physics_index, realization_index, &
    initialization_index
  namelist /experiment/ activity_id, variant_label, parent_variant_label, parent_mip_era, &
    mip_era, sub_experiment_id, parent_sub_experiment, parent_activity_id, &
    branch_time_in_child, branch_time_in_parent, branch_method, parent_time_units, &
    tracking_prefix, physics_index, initialization_index, realization_index, forcing_index, &
    source_type, casename, experiment_id, history, comment, forcing, realization, branch_time, &
    exprefyear, parent_experiment_id, parent_experiment_rip, year1, month1, yearn, monthn, &
    do_fx, do_amon, do_oyr, do_oyrbgc, do_omon, do_omonbgc, do_oimon, do_lmon, do_limon, &
    do_aero, do_day, do_6hrlev, do_6hrlevi, do_6hrplev, do_3hr, do_3hri, do_2d, do_3d, do_xd, &
    do_bgc, do_ofx, do_6hrPlevPt, do_AERday, do_AERhr, do_AERmonZ, do_CF3hr, do_CFday, do_CFmon, &
    do_CFsubhr, do_E1hrClimMon, do_E1hr, do_E3hr, do_E3hrPt, do_E6hrZ, do_Eday, do_EdayZ, do_Efx, &
    do_Eyr, do_Emon, do_EmonZ, do_Esubhr, do_Oclim, do_Oday, do_Odaybgc, do_Odec, do_COfx, &
    do_SIday, dry_run, plevdummy, readdummy, add_fill_day, newcolumnorder, isubdir, osubdir, &
    scanallfiles, membertag

  ! Tables
  logical, save :: dfx, damon, domon, doimon, daero, dday, d6hrlev, d6hrlevi, d6hrplev, &
    d3hr, d3hri, dlmon, dlimon, doyr, doyrbgc, domonbgc
  character(len=lenmax), dimension(colmax, rowmax), save :: vfx, vamon, vomon, voimon, vaero, &
    vday, v6hrlev, v6hrlevi, v6hrplev, v3hr, v3hri, vlmon, vlimon, voyr, voyrbgc, vomonbgc
  character(len=slenmax), save :: pfx, pamon, pomon, poimon, paero, pday, p6hrlev, p6hrlevi, &
    p6hrplev, p3hr, p3hri, plmon, plimon, poyr, poyrbgc, pomonbgc
  character(len=slenmax), save :: tfx, tamon, tomon, toimon, taero, tday, t6hrlev, t6hrlevi, &
    t6hrplev, t3hr, t3hri, tlmon, tlimon, toyr, toyrbgc, tomonbgc, tgrids
  integer, save :: nfx, namon, nomon, noimon, naero, nday, n6hrlev, n6hrlevi, n6hrplev, &
    n3hr, n3hri, nlmon, nlimon, noyr, noyrbgc, nomonbgc
  integer, save :: rfx, ramon, romon, roimon, raero, rday, r6hrlev, r6hrlevi, r6hrplev, &
    r3hr, r3hri, rlmon, rlimon, royr, royrbgc, romonbgc
  namelist /table_grids/ tgrids
  namelist /table_fx/ dfx, pfx, tfx, vfx
  namelist /table_amon/ damon, pamon, tamon, ramon, vamon
  namelist /table_aero/ daero, paero, taero, raero, vaero
  namelist /table_oyr/ doyr, poyr, toyr, royr, voyr
  namelist /table_oyrbgc/ doyrbgc, poyrbgc, toyrbgc, royrbgc, voyrbgc
  namelist /table_omon/ domon, pomon, tomon, romon, vomon
  namelist /table_omonbgc/ domonbgc, pomonbgc, tomonbgc, romonbgc, vomonbgc
  namelist /table_oimon/ doimon, poimon, toimon, roimon, voimon
  namelist /table_lmon/ dlmon, plmon, tlmon, rlmon, vlmon
  namelist /table_limon/ dlimon, plimon, tlimon, rlimon, vlimon
  namelist /table_day/ dday, pday, tday, rday, vday
  namelist /table_6hrlev/ d6hrlev, p6hrlev, t6hrlev, r6hrlev, v6hrlev
  namelist /table_6hrlevi/ d6hrlevi, p6hrlevi, t6hrlevi, r6hrlevi, v6hrlevi
  namelist /table_6hrplev/ d6hrplev, p6hrplev, t6hrplev, r6hrplev, v6hrplev
  namelist /table_3hr/ d3hr, p3hr, t3hr, r3hr, v3hr
  namelist /table_3hri/ d3hri, p3hri, t3hri, r3hri, v3hri
  logical, save :: dofx, d6hrPlevPt, dAERday, dAERhr, dAERmonZ, dCF3hr, dCFday, dCFmon, &
    dCFsubhr, dE1hrClimMon, dE1hr, dE3hr, dE3hrPt, dE6hrZ, dEday, dEdayZ, dEfx, dEyr, dEmon, &
    dEmonZ, dEsubhr, dOclim, dOday, dOdaybgc, dOdec, dCOfx, dSIday
  character(len=lenmax), dimension(colmax, rowmax), save :: vofx, v6hrPlevPt, vAERday, vAERhr, &
    vAERmonZ, vCF3hr, vCFday, vCFmon, vCFsubhr, vE1hrClimMon, vE1hr, vE3hr, vE3hrPt, vE6hrZ, &
    vEday, vEdayZ, vEfx, vEyr, vEmon, vEmonZ, vEsubhr, vOclim, vOday, vOdaybgc, vOdec, vCOfx, &
    vSIday
  character(len=slenmax), save :: pofx, p6hrPlevPt, pAERday, pAERhr, pAERmonZ, pCF3hr, pCFday, &
    pCFmon, pCFsubhr, pE1hrClimMon, pE1hr, pE3hr, pE3hrPt, pE6hrZ, pEday, pEdayZ, pEfx, pEyr, &
    pEmon, pEmonZ, pEsubhr, pOclim, pOday, pOdaybgc, pOdec, pCOfx, pSIday
  character(len=slenmax), save :: tofx, t6hrPlevPt, tAERday, tAERhr, tAERmonZ, tCF3hr, tCFday, &
    tCFmon, tCFsubhr, tE1hrClimMon, tE1hr, tE3hr, tE3hrPt, tE6hrZ, tEday, tEdayZ, tEfx, tEyr, &
    tEmon, tEmonZ, tEsubhr, tOclim, tOday, tOdaybgc, tOdec, tCOfx, tSIday
  integer, save :: nofx, n6hrPlevPt, nAERday, nAERhr, nAERmonZ, nCF3hr, nCFday, nCFmon, &
    nCFsubhr, nE1hrClimMon, nE1hr, nE3hr, nE3hrPt, nE6hrZ, nEday, nEdayZ, nEfx, nEyr, nEmon, &
    nEmonZ, nEsubhr, nOclim, nOday, nOdaybgc, nOdec, nCOfx, nSIday
  integer, save :: rofx, r6hrPlevPt, rAERday, rAERhr, rAERmonZ, rCF3hr, rCFday, rCFmon, &
    rCFsubhr, rE1hrClimMon, rE1hr, rE3hr, rE3hrPt, rE6hrZ, rEday, rEdayZ, rEfx, rEyr, rEmon, &
    rEmonZ, rEsubhr, rOclim, rOday, rOdaybgc, rOdec, rCOfx, rSIday
  namelist /table_ofx/ dofx, pofx, tofx, vofx, rofx
  namelist /table_6hrPlevPt/ d6hrPlevPt, p6hrPlevPt, t6hrPlevPt, v6hrPlevPt, r6hrPlevPt
  namelist /table_AERday/ dAERday, pAERday, tAERday, vAERday, rAERday
  namelist /table_AERhr/ dAERhr, pAERhr, tAERhr, vAERhr, rAERhr
  namelist /table_AERmonZ/ dAERmonZ, pAERmonZ, tAERmonZ, vAERmonZ, rAERmonZ
  namelist /table_CF3hr/ dCF3hr, pCF3hr, tCF3hr, vCF3hr, rCF3hr
  namelist /table_CFday/ dCFday, pCFday, tCFday, vCFday, rCFday
  namelist /table_CFmon/ dCFmon, pCFmon, tCFmon, vCFmon, rCFmon
  namelist /table_CFsubhr/ dCFsubhr, pCFsubhr, tCFsubhr, vCFsubhr, rCFsubhr
  namelist /table_E1hrClimMon/ dE1hrClimMon, pE1hrClimMon, tE1hrClimMon, vE1hrClimMon, &
    rE1hrClimMon
  namelist /table_E1hr/ dE1hr, pE1hr, tE1hr, vE1hr, rE1hr
  namelist /table_E3hr/ dE3hr, pE3hr, tE3hr, vE3hr, rE3hr
  namelist /table_E3hrPt/ dE3hrPt, pE3hrPt, tE3hrPt, vE3hrPt, rE3hrPt
  namelist /table_E6hrZ/ dE6hrZ, pE6hrZ, tE6hrZ, vE6hrZ, rE6hrZ
  namelist /table_Eday/ dEday, pEday, tEday, vEday, rEday
  namelist /table_EdayZ/ dEdayZ, pEdayZ, tEdayZ, vEdayZ, rEdayZ
  namelist /table_Efx/ dEfx, pEfx, tEfx, vEfx, rEfx
  namelist /table_Emon/ dEmon, pEmon, tEmon, vEmon, rEmon
  namelist /table_EmonZ/ dEmonZ, pEmonZ, tEmonZ, vEmonZ, rEmonZ
  namelist /table_Esubhr/ dEsubhr, pEsubhr, tEsubhr, vEsubhr, rEsubhr
  namelist /table_Eyr/ dEyr, pEyr, tEyr, vEyr, rEyr
  namelist /table_Oclim/ dOclim, pOclim, tOclim, vOclim, rOclim
  namelist /table_Oday/ dOday, pOday, tOday, vOday, rOday
  namelist /table_Odaybgc/ dOdaybgc, pOdaybgc, tOdaybgc, vOdaybgc, rOdaybgc
  namelist /table_Odec/ dOdec, pOdec, tOdec, vOdec, rOdec
  namelist /table_COfx/ dCOfx, pCOfx, tCOfx, vCOfx, rCOfx
  namelist /table_SIday/ dSIday, pSIday, tSIday, vSIday, rSIday

  ! Misc
  integer :: istatus, funit
  character(len=slenmax), save :: fnm, fnm2, itag, nmlfp, vsingle = ' '
  character(len=slenmax), save :: nmlfpsys, nmlfpmod, nmlfpexp, nmlfpvar

  ! Time related variables
  logical, save :: linstant
  integer, save :: year, month, rec
  real(kind=8), save :: tval(1), tval2(2), tbnd(2), mbnd(2), tbnds(2,1), mbnds(2,1)
  character(len=slenmax), save :: calendar = 'noleap', &
    calunits = 'days since 1850-01-01 00:00:00'

contains

  ! -----------------------------------------------------------------

  subroutine read_namelists

    implicit none

    integer :: n, pi, po
    logical :: fexist

    ! Initialise namelist variables
    ibasedir      = ' '
    obasedir      = ' '
    tabledir      = ' '
    griddata      = ' '
    tagoyr        = ' '
    tagoyrbgc     = ' '
    tagomon       = ' '
    tagomonbgc    = ' '
    tagoday       = ' '
    tagodaybgc    = ' '
    tagimon       = ' '
    tagiday       = ' '
    tagamon       = ' '
    tagaday       = ' '
    taga6hr       = ' '
    taga6hri      = ' '
    taga3hr       = ' '
    taga3hri      = ' '
    taglyr        = ' '
    taglmon       = ' '
    taglday       = ' '
    tagl3hr       = ' '
    tagl3hri      = ' '
    atmgridfile   = ' '
    ocngridfile   = ' '
    ocninitfile   = ' '
    ocnmertfile   = ' '
    ocnregnfile   = ' '
    secindexfile  = ' '
    rhotablesuff  = 'OnRho'
    coordtable    = 'CMIP7_coordinate.json'
    year1         = 0
    month1        = 1
    yearn         = 0
    monthn        = 12
    createsubdirs = .true.
    forcefilescan = .true.
    verbose       = .true.
    do_fx         = .true.
    do_oyr        = .true.
    do_oyrbgc     = .true.
    do_omon       = .true.
    do_omonbgc    = .true.
    do_oimon      = .true.
    do_amon       = .true.
    do_aero       = .true.
    do_lmon       = .true.
    do_limon      = .true.
    do_day        = .true.
    do_6hrlev     = .true.
    do_6hrlevi    = .true.
    do_6hrplev    = .true.
    do_3hr        = .true.
    do_3hri       = .true.
    do_2d         = .true.
    do_3d         = .true.
    do_xd         = .true.
    do_bgc        = .true.
    do_ofx        = .true.
    do_6hrPlevPt  = .true.
    do_AERday     = .true.
    do_AERhr      = .true.
    do_AERmonZ    = .true.
    do_CF3hr      = .true.
    do_CFday      = .true.
    do_CFmon      = .true.
    do_CFsubhr    = .true.
    do_E1hrClimMon= .true.
    do_E1hr       = .true.
    do_E3hr       = .true.
    do_E3hrPt     = .true.
    do_E6hrZ      = .true.
    do_Eday       = .true.
    do_EdayZ      = .true.
    do_Efx        = .true.
    do_Emon       = .true.
    do_EmonZ      = .true.
    do_Esubhr     = .true.
    do_Eyr        = .true.
    do_Oclim      = .true.
    do_Oday       = .true.
    do_Odaybgc    = .true.
    do_Odec       = .true.
    do_COfx       = .true.
    do_SIday      = .true.
    dry_run       = .false.
    plevdummy     = .false.
    readdummy     = .false.
    add_fill_day  = .false.
    newcolumnorder= .true.
    scanallfiles  = .true.

    casename      = ' '
    experiment_id = ' '
    institute_id  = ' '
    institution   = ' '
    source        = ' '
    contact       = ' '
    history       = ' '
    comment       = ' '
    references    = ' '
    model_id      = ' '
    forcing       = ' '
    realization   = 1
    branch_time   = 0.0
    parent_experiment_id = ' '
    parent_experiment_rip = ' '
    isubdir       = ' '
    osubdir       = ' '
    membertag     = ' '

    dfx           = .false.
    doyr          = .false.
    doyrbgc       = .false.
    domon         = .false.
    domonbgc      = .false.
    doimon        = .false.
    damon         = .false.
    daero         = .false.
    dlmon         = .false.
    dlimon        = .false.
    dday          = .false.
    d6hrlev       = .false.
    d6hrlevi      = .false.
    d6hrplev      = .false.
    d3hr          = .false.
    d3hri         = .false.
    dofx          = .false.
    d6hrPlevPt    = .false.
    dAERday       = .false.
    dAERhr        = .false.
    dAERmonZ      = .false.
    dCF3hr        = .false.
    dCFday        = .false.
    dCFmon        = .false.
    dCFsubhr      = .false.
    dE1hrClimMon  = .false.
    dE1hr         = .false.
    dE3hr         = .false.
    dE3hrPt       = .false.
    dE6hrZ        = .false.
    dEday         = .false.
    dEdayZ        = .false.
    dEfx          = .false.
    dEmon         = .false.
    dEmonZ        = .false.
    dEsubhr       = .false.
    dEyr          = .false.
    dOclim        = .false.
    dOday         = .false.
    dOdaybgc      = .false.
    dOdec         = .false.
    dCOfx         = .false.
    dSIday        = .false.

    vfx           = ' '
    voyr          = ' '
    voyrbgc       = ' '
    vomon         = ' '
    vomonbgc      = ' '
    voimon        = ' '
    vamon         = ' '
    vaero         = ' '
    vlmon         = ' '
    vlimon        = ' '
    vday          = ' '
    v6hrlev       = ' '
    v6hrlevi      = ' '
    v6hrplev      = ' '
    v3hr          = ' '
    v3hri         = ' '
    vofx          = ' '
    v6hrPlevPt    = ' '
    vAERday       = ' '
    vAERhr        = ' '
    vAERmonZ      = ' '
    vCF3hr        = ' '
    vCFday        = ' '
    vCFmon        = ' '
    vCFsubhr      = ' '
    vE1hrClimMon  = ' '
    vE1hr         = ' '
    vE3hr         = ' '
    vE3hrPt       = ' '
    vE6hrZ        = ' '
    vEday         = ' '
    vEdayZ        = ' '
    vEfx          = ' '
    vEmon         = ' '
    vEmonZ        = ' '
    vEsubhr       = ' '
    vEyr          = ' '
    vOclim        = ' '
    vOday         = ' '
    vOdaybgc      = ' '
    vOdec         = ' '
    vCOfx         = ' '
    vSIday        = ' '

    pfx           = ' '
    poyr          = ' '
    poyrbgc       = ' '
    pomon         = ' '
    pomonbgc      = ' '
    poimon        = ' '
    pamon         = ' '
    paero         = ' '
    plmon         = ' '
    plimon        = ' '
    pday          = ' '
    p6hrlev       = ' '
    p6hrlevi      = ' '
    p6hrplev      = ' '
    p3hr          = ' '
    p3hri         = ' '
    pofx          = ' '
    p6hrPlevPt    = ' '
    pAERday       = ' '
    pAERhr        = ' '
    pAERmonZ      = ' '
    pCF3hr        = ' '
    pCFday        = ' '
    pCFmon        = ' '
    pCFsubhr      = ' '
    pE1hrClimMon  = ' '
    pE1hr         = ' '
    pE3hr         = ' '
    pE3hrPt       = ' '
    pE6hrZ        = ' '
    pEday         = ' '
    pEdayZ        = ' '
    pEfx          = ' '
    pEmon         = ' '
    pEmonZ        = ' '
    pEsubhr       = ' '
    pEyr          = ' '
    pOclim        = ' '
    pOday         = ' '
    pOdaybgc      = ' '
    pOdec         = ' '
    pCOfx         = ' '
    pSIday        = ' '

    royr          = 1000000
    royrbgc       = 1000000
    romon         = 1000000
    romonbgc      = 1000000
    roimon        = 1000000
    ramon         = 1000000
    raero         = 1000000
    rlmon         = 1000000
    rlimon        = 1000000
    rday          = 1000000
    r6hrlev       = 1000000
    r6hrlevi      = 1000000
    r6hrplev      = 1000000
    r3hr          = 1000000
    r3hri         = 1000000
    rofx          = 1000000
    r6hrPlevPt    = 1000000
    rAERday       = 1000000
    rAERhr        = 1000000
    rAERmonZ      = 1000000
    rCF3hr        = 1000000
    rCFday        = 1000000
    rCFmon        = 1000000
    rCFsubhr      = 1000000
    rE1hrClimMon  = 1000000
    rE1hr         = 1000000
    rE3hr         = 1000000
    rE3hrPt       = 1000000
    rE6hrZ        = 1000000
    rEday         = 1000000
    rEdayZ        = 1000000
    rEfx          = 1000000
    rEmon         = 1000000
    rEmonZ        = 1000000
    rEsubhr       = 1000000
    rEyr          = 1000000
    rOclim        = 1000000
    rOday         = 1000000
    rOdaybgc      = 1000000
    rOdec         = 1000000
    rCOfx         = 1000000
    rSIday        = 1000000

    tgrids        = 'CMIP7_grids.json'
    tfx           = 'CMIP6_fx.json'
    toyr          = 'CMIP6_Oyr.json'
    toyrbgc       = 'CMIP6_Oyr.json'
    tomon         = 'CMIP6_Omon.json'
    tomonbgc      = 'CMIP6_Omon.json'
    toimon        = 'CMIP6_OImon.json'
    tamon         = 'CMIP6_Amon.json'
    taero         = 'CMIP6_aero.json'
    tlmon         = 'CMIP6_Lmon.json'
    tlimon        = 'CMIP6_LImon.json'
    tday          = 'CMIP6_da.json'
    t6hrlev       = 'CMIP6_6hrLev.json'
    t6hrlevi      = 'CMIP6_6hrLev.json'
    t6hrplev      = 'CMIP6_6hrPlev.json'
    t3hr          = 'CMIP6_3hr.json'
    t3hri         = 'CMIP6_3hr.json'

    tofx          = 'CMIP6_ofx.json'
    t6hrPlevPt    = 'CMIP6_6hrPlevPt.json'
    tAERday       = 'CMIP6_AERday.json'
    tAERhr        = 'CMIP6_AERhr.json'
    tAERmonZ      = 'CMIP6_AERmonZ.json'
    tCF3hr        = 'CMIP6_CF3hr.json'
    tCFday        = 'CMIP6_CFday.json'
    tCFmon        = 'CMIP6_CFmon.json'
    tCFsubhr      = 'CMIP6_CFsubhr.json'
    tE1hrClimMon  = 'CMIP6_E1hrClimMon.json'
    tE1hr         = 'CMIP6_E1hr.json'
    tE3hr         = 'CMIP6_E3hr.json'
    tE3hrPt       = 'CMIP6_E3hrPt.json'
    tE6hrZ        = 'CMIP6_E6hrZ.json'
    tEday         = 'CMIP6_Eday.json'
    tEdayZ        = 'CMIP6_EdayZ.json'
    tEfx          = 'CMIP6_Efx.json'
    tEmon         = 'CMIP6_Emon.json'
    tEmonZ        = 'CMIP6_EmonZ.json'
    tEsubhr       = 'CMIP6_Esubhr.json'
    tEyr          = 'CMIP6_Eyr.json'
    tOclim        = 'CMIP6_Oclim.json'
    tOday         = 'CMIP6_Oday.json'
    tOdaybgc      = 'CMIP6_Oday.json'
    tOdec         = 'CMIP6_Odec.json'
    tCOfx         = 'CMIP6_COfx.json'
    tSIday        = 'CMIP6_SIday.json'

    ! Read namelists
    if (iargc() /= 1 .and. iargc() /= 2 .and. iargc() /= 4 .and. iargc() /= 5) then
      write(*, *) 'Usage: noresm2cmor <path to master namelist file> [<cmor variable to process>]'
      write(*, *) '       or'
      write(*, *) '       noresm2cmor <system nml-file> <model nml-file> <exp nml-file> ' // &
        '<variable nml-file> [<cmor variable to process>]'
      stop
    else if (iargc() == 1 .or. iargc() == 2) then
      nmlfp = ' '
      call getarg(1, nmlfp)
      inquire(file=trim(nmlfp), exist=fexist)
      if (.not. fexist) then
        write(*, *) 'cannot find namelist file'//trim(nmlfp)
        stop
      end if
      nmlfpsys = nmlfp
      nmlfpmod = nmlfp
      nmlfpexp = nmlfp
      nmlfpvar = nmlfp
      if (iargc() == 2) then
        call getarg(2, vsingle)
      end if
    else if (iargc() == 4 .or. iargc() == 5) then
      nmlfpsys = ' '
      call getarg(1, nmlfpsys)
      inquire(file=trim(nmlfpsys), exist=fexist)
      if (.not. fexist) then
        write(*, *) 'cannot find namelist file'//trim(nmlfpsys)
        stop
      end if
      nmlfpmod = ' '
      call getarg(2, nmlfpmod)
      inquire(file=trim(nmlfpmod), exist=fexist)
      if (.not. fexist) then
        write(*, *) 'cannot find namelist file'//trim(nmlfpmod)
        stop
      end if
      nmlfpexp = ' '
      call getarg(3, nmlfpexp)
      inquire(file=trim(nmlfpexp), exist=fexist)
      if (.not. fexist) then
        write(*, *) 'cannot find namelist file'//trim(nmlfpexp)
        stop
      end if
      nmlfpvar = ' '
      call getarg(4, nmlfpvar)
      inquire(file=trim(nmlfpvar), exist=fexist)
      if (.not. fexist) then
        write(*, *) 'cannot find namelist file'//trim(nmlfpvar)
        stop
      end if
      if (iargc() == 5) then
        call getarg(5, vsingle)
      end if
    end if

    funit = get_free_unit()
    open(funit, file=trim(nmlfpsys), status='old', action='read', recl=200)
    read(funit, nml=sys, iostat=istatus)
    close(funit)
    if (istatus /= 0) then
      write(*, *) 'Problem reading namelist system in file '//trim(nmlfpsys)
      stop
    end if
    griddata = trim(griddata)//'/'
    tabledir = trim(tabledir)//'/'

    open(funit, file=trim(nmlfpmod), status='old', action='read', recl=200)
    read(funit, nml=model, iostat=istatus)
    close(funit)
    if (istatus /= 0) then
      write(*, *) 'Problem reading namelist model in file '//trim(nmlfpmod)
      stop
    end if

    open(funit, file=trim(nmlfpexp), status='old', action='read', recl=200)
    read(funit, nml=experiment, iostat=istatus)
    close(funit)
    if (istatus /= 0) then
      write(*, *) 'Problem reading namelist experiment in file '//trim(nmlfpexp)
      stop
    end if

    open(funit, file=trim(nmlfpvar), status='old', action='read', recl=200)
    read(funit, nml=table_grids, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table grids not in namelist file. Using CMIP5 default'

    read(funit, nml=table_fx, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table FX not in namelist file. Skipping table...'

    read(funit, nml=table_oyr, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Oyr not in namelist file. Skipping table...'

    read(funit, nml=table_oyrbgc, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Oyrbgc not in namelist file. Skipping table...'

    read(funit, nml=table_omon, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Omon not in namelist file. Skipping table...'

    read(funit, nml=table_omonbgc, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Omonbgc not in namelist file. Skipping table...'

    read(funit, nml=table_oimon, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table OImon not in namelist file. Skipping table...'

    read(funit, nml=table_amon, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Amon not in namelist file. Skipping table...'

    read(funit, nml=table_aero, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Aero not in namelist file. Skipping table...'

    read(funit, nml=table_lmon, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Lmon not in namelist file. Skipping table...'

    read(funit, nml=table_limon, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table LImon not in namelist file. Skipping table...'

    read(funit, nml=table_day, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Day not in namelist file. Skipping table...'

    read(funit, nml=table_6hrlev, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table 6hrLev not in namelist file. Skipping table...'

    read(funit, nml=table_6hrlevi, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table 6hrLevi not in namelist file. Skipping table...'

    read(funit, nml=table_6hrplev, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table 6hrPLev not in namelist file. Skipping table...'

    read(funit, nml=table_3hr, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table 3hr not in namelist file. Skipping table...'

    read(funit, nml=table_3hri, iostat=istatus)
    rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table 3hri not in namelist file. Skipping table...'

    read(funit, nml=table_ofx, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Ofx not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_6hrPlevPt, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table 6hrPlevPt not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_AERday, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table AERday not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_AERhr, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table AERhr not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_AERmonZ, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table AERmonZ not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_CF3hr, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table CF3hr not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_CFday, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table CFday not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_CFmon, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table CFmon not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_CFsubhr, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table CFsubhr not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_E1hrClimMon, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table E1hrClimMon not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_E1hr, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table E1hr not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_E3hr, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table E3hr not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_E3hrPt, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table E3hrPt not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_E6hrZ, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table E6hrZ not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_Eday, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Eday not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_EdayZ, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table EdayZ not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_Efx, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Efx not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_Emon, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Emon not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_EmonZ, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table EmonZ not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_Esubhr, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Esubhr not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_Eyr, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Eyr not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_Oclim, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Oclim not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_Oday, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Oday not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_Odaybgc, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Odaybgc not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_Odec, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table Odec not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_COfx, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table COfx not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    read(funit, nml=table_SIday, iostat=istatus)
    if (istatus /= 0) write(*, *) &
      'WARNING: Table SIday not in namelist file. Skipping table...'
    rewind(funit, iostat=istatus)
    if (istatus /= 0) write(*, *) 'cannot rewind namelist file'

    close(funit)

    ! Map column order of variable definitions
    if (newcolumnorder) then
      ivnmpos = 2
      ovnmpos = 1
    end if

    ! Merge global string arrays
    call merge_strarr(slenmax, smax, source, source1, linebreaks)
    call merge_strarr(slenmax, smax, history, history1, linebreaks)
    call merge_strarr(slenmax, smax, comment, comment1, linebreaks)
    call merge_strarr(slenmax, smax, references, references1, linebreaks)
    call merge_strarr(slenmax, smax, forcing, forcing1, .false.)
    call merge_strarr(slenmax, smax, contact, contact1, .false.)
    call merge_strarr(slenmax, smax, institution, institution1, .false.)

    ! Count number of table entries
    nfx = 0
    noyr = 0
    noyrbgc = 0
    nomon = 0
    nomonbgc = 0
    noimon = 0
    namon = 0
    naero = 0
    nlmon = 0
    nlimon = 0
    nday = 0
    n6hrlev = 0
    n6hrlevi = 0
    n6hrplev = 0
    n3hr = 0
    n3hri = 0

    do n = 1, rowmax
      if (len_trim(vfx(1, n)) /= 0) nfx = nfx + 1
      if (len_trim(voyr(1, n)) /= 0) noyr = noyr + 1
      if (len_trim(voyrbgc(1, n)) /= 0) noyrbgc = noyrbgc + 1
      if (len_trim(vomon(1, n)) /= 0) nomon = nomon + 1
      if (len_trim(vomonbgc(1, n)) /= 0) nomonbgc = nomonbgc + 1
      if (len_trim(voimon(1, n)) /= 0) noimon = noimon + 1
      if (len_trim(vamon(1, n)) /= 0) namon = namon + 1
      if (len_trim(vaero(1, n)) /= 0) naero = naero + 1
      if (len_trim(vlmon(1, n)) /= 0) nlmon = nlmon + 1
      if (len_trim(vlimon(1, n)) /= 0) nlimon = nlimon + 1
      if (len_trim(vday(1, n)) /= 0) nday = nday + 1
      if (len_trim(v6hrlev(1, n)) /= 0) n6hrlev = n6hrlev + 1
      if (len_trim(v6hrlevi(1, n)) /= 0) n6hrlevi = n6hrlevi + 1
      if (len_trim(v6hrplev(1, n)) /= 0) n6hrplev = n6hrplev + 1
      if (len_trim(v3hr(1, n)) /= 0) n3hr = n3hr + 1
      if (len_trim(v3hri(1, n)) /= 0) n3hri = n3hri + 1
    end do

    if (len_trim(vsingle) /= 0) then
      pi = ivnmpos
      po = ovnmpos
      do n = 1, rowmax
        if (trim(vsingle) /= trim(vfx(po, n))) vfx(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(voyr(po, n))) voyr(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(voyrbgc(po, n))) voyrbgc(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vomon(po, n))) vomon(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vomonbgc(po, n))) vomonbgc(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(voimon(po, n))) voimon(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vamon(po, n))) vamon(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vaero(po, n))) vaero(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vlmon(po, n))) vlmon(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vlimon(po, n))) vlimon(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vday(po, n))) vday(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(v6hrlev(po, n))) v6hrlev(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(v6hrlevi(po, n))) v6hrlevi(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(v6hrplev(po, n))) v6hrplev(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(v3hr(po, n))) v3hr(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(v3hri(po, n))) v3hri(pi, n) = 'SKIP'
      end do
    end if

    nofx = 0
    n6hrPlevPt = 0
    nAERday = 0
    nAERhr = 0
    nAERmonZ = 0
    nCF3hr = 0
    nCFday = 0
    nCFmon = 0
    nCFsubhr = 0
    nE1hrClimMon = 0
    nE1hr = 0
    nE3hr = 0
    nE3hrPt = 0
    nE6hrZ = 0
    nEday = 0
    nEdayZ = 0
    nEfx = 0
    nEmon = 0
    nEmonZ = 0
    nEsubhr = 0
    nEyr = 0
    nOclim = 0
    nOday = 0
    nOdaybgc = 0
    nOdec = 0
    nCOfx = 0
    nSIday = 0

    do n = 1, rowmax
      if (len_trim(vofx(1, n)) /= 0) nofx = nofx + 1
      if (len_trim(v6hrPlevPt(1, n)) /= 0) n6hrPlevPt = n6hrPlevPt + 1
      if (len_trim(vAERday(1, n)) /= 0) nAERday = nAERday + 1
      if (len_trim(vAERhr(1, n)) /= 0) nAERhr = nAERhr + 1
      if (len_trim(vAERmonZ(1, n)) /= 0) nAERmonZ = nAERmonZ + 1
      if (len_trim(vCF3hr(1, n)) /= 0) nCF3hr = nCF3hr + 1
      if (len_trim(vCFday(1, n)) /= 0) nCFday = nCFday + 1
      if (len_trim(vCFmon(1, n)) /= 0) nCFmon = nCFmon + 1
      if (len_trim(vCFsubhr(1, n)) /= 0) nCFsubhr = nCFsubhr + 1
      if (len_trim(vE1hrClimMon(1, n)) /= 0) nE1hrClimMon = nE1hrClimMon + 1
      if (len_trim(vE1hr(1, n)) /= 0) nE1hr = nE1hr + 1
      if (len_trim(vE3hr(1, n)) /= 0) nE3hr = nE3hr + 1
      if (len_trim(vE3hrPt(1, n)) /= 0) nE3hrPt = nE3hrPt + 1
      if (len_trim(vE6hrZ(1, n)) /= 0) nE6hrZ = nE6hrZ + 1
      if (len_trim(vEday(1, n)) /= 0) nEday = nEday + 1
      if (len_trim(vEdayZ(1, n)) /= 0) nEdayZ = nEdayZ + 1
      if (len_trim(vEfx(1, n)) /= 0) nEfx = nEfx + 1
      if (len_trim(vEmon(1, n)) /= 0) nEmon = nEmon + 1
      if (len_trim(vEmonZ(1, n)) /= 0) nEmonZ = nEmonZ + 1
      if (len_trim(vEsubhr(1, n)) /= 0) nEsubhr = nEsubhr + 1
      if (len_trim(vEyr(1, n)) /= 0) nEyr = nEyr + 1
      if (len_trim(vOclim(1, n)) /= 0) nOclim = nOclim + 1
      if (len_trim(vOday(1, n)) /= 0) nOday = nOday + 1
      if (len_trim(vOdaybgc(1, n)) /= 0) nOdaybgc = nOdaybgc + 1
      if (len_trim(vOdec(1, n)) /= 0) nOdec = nOdec + 1
      if (len_trim(vCOfx(1, n)) /= 0) nCOfx = nCOfx + 1
      if (len_trim(vSIday(1, n)) /= 0) nSIday = nSIday + 1
    end do

    if (len_trim(vsingle) /= 0) then
      pi = ivnmpos
      po = ovnmpos
      do n = 1, rowmax
        if (trim(vsingle) /= trim(vofx(po, n))) vofx(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(v6hrPlevPt(po, n))) v6hrPlevPt(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vAERday(po, n))) vAERday(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vAERhr(po, n))) vAERhr(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vAERmonZ(po, n))) vAERmonZ(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vCF3hr(po, n))) vCF3hr(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vCFday(po, n))) vCFday(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vCFmon(po, n))) vCFmon(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vCFsubhr(po, n))) vCFsubhr(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vE1hrClimMon(po, n))) vE1hrClimMon(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vE1hr(po, n))) vE1hr(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vE3hr(po, n))) vE3hr(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vE3hrPt(po, n))) vE3hrPt(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vE6hrZ(po, n))) vE6hrZ(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vEday(po, n))) vEday(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vEdayZ(po, n))) vEdayZ(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vEfx(po, n))) vEfx(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vEmon(po, n))) vEmon(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vEmonZ(po, n))) vEmonZ(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vEsubhr(po, n))) vEsubhr(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vEyr(po, n))) vEyr(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vOclim(po, n))) vOclim(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vOday(po, n))) vOday(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vOdaybgc(po, n))) vOdaybgc(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vOdec(po, n))) vOdec(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vCOfx(po, n))) vCOfx(pi, n) = 'SKIP'
        if (trim(vsingle) /= trim(vSIday(po, n))) vSIday(pi, n) = 'SKIP'
      end do
    end if

    ! Skip deselected namelists
    if (.not. do_fx)      nfx = 0
    if (.not. do_oyr)     noyr = 0
    if (.not. do_oyrbgc .or. .not. do_bgc)  noyrbgc = 0
    if (.not. do_omon)    nomon = 0
    if (.not. do_omonbgc .or. .not. do_bgc) nomonbgc = 0
    if (.not. do_oimon)   noimon = 0
    if (.not. do_amon)    namon = 0
    if (.not. do_aero)    naero = 0
    if (.not. do_lmon)    nlmon = 0
    if (.not. do_limon)   nlimon = 0
    if (.not. do_day)     nday = 0
    if (.not. do_6hrlev)  n6hrlev = 0
    if (.not. do_6hrlevi) n6hrlevi = 0
    if (.not. do_6hrplev) n6hrplev = 0
    if (.not. do_3hr)     n3hr = 0
    if (.not. do_3hri)    n3hri = 0
    if (.not. do_ofx)         nofx = 0
    if (.not. do_6hrPlevPt)   n6hrPlevPt = 0
    if (.not. do_AERday)      nAERday = 0
    if (.not. do_AERhr)       nAERhr = 0
    if (.not. do_AERmonZ)     nAERmonZ = 0
    if (.not. do_CF3hr)       nCF3hr = 0
    if (.not. do_CFday)       nCFday = 0
    if (.not. do_CFmon)       nCFmon = 0
    if (.not. do_CFsubhr)     nCFsubhr = 0
    if (.not. do_E1hrClimMon) nE1hrClimMon = 0
    if (.not. do_E1hr)        nE1hr = 0
    if (.not. do_E3hr)        nE3hr = 0
    if (.not. do_E3hrPt)      nE3hrPt = 0
    if (.not. do_E6hrZ)       nE6hrZ = 0
    if (.not. do_Eday)        nEday = 0
    if (.not. do_EdayZ)       nEdayZ = 0
    if (.not. do_Efx)         nEfx = 0
    if (.not. do_Emon)        nEmon = 0
    if (.not. do_EmonZ)       nEmonZ = 0
    if (.not. do_Esubhr)      nEsubhr = 0
    if (.not. do_Eyr)         nEyr = 0
    if (.not. do_Oclim)       nOclim = 0
    if (.not. do_Oday)        nOday = 0
    if (.not. do_Odaybgc .or. .not. do_bgc) nOdaybgc = 0
    if (.not. do_Odec)        nOdec = 0
    if (.not. do_COfx)        nCOfx = 0
    if (.not. do_SIday)       nSIday = 0

    ! Extend input path if necessary
    if (len_trim(isubdir) > 0) &
      ibasedir = trim(ibasedir)//'/'//trim(isubdir)

    ! Modify output path and create output folder
    obasedir = trim(obasedir)//'/'//trim(osubdir)
    call system('mkdir -p '//trim(obasedir))

  end subroutine read_namelists

  ! -----------------------------------------------------------------

  subroutine print_namelists

    implicit none

    integer :: n

    write(*, *)
    write(*, *) 'System namelist:'
    write(*, *) ' input directory  = ', trim(ibasedir)
    write(*, *) ' output directory = ', trim(obasedir)
    write(*, *) ' table directory  = ', trim(tabledir)
    write(*, *) ' grid data dir.   = ', trim(griddata)
    write(*, *) ' create sub-dirs  = ', createsubdirs
    write(*, *) ' verbose          = ', verbose
    write(*, *)
    write(*, *) 'Model namelist:'
    write(*, *) ' institution      = ', trim(institution1)
    write(*, *) ' model id         = ', trim(model_id)
    write(*, *) ' source           = ', trim(source1)
    write(*, *) ' references       = ', trim(references1)
    write(*, *) ' contact          = ', trim(contact1)
    write(*, *) ' tag annual ocn   = ', trim(tagoyr)
    write(*, *) ' tag annual bgc   = ', trim(tagoyrbgc)
    write(*, *) ' tag monthly ocn  = ', trim(tagomon)
    write(*, *) ' tag monthly bgc  = ', trim(tagomonbgc)
    write(*, *) ' tag daily ocn    = ', trim(tagoday)
    write(*, *) ' tag daily bgc    = ', trim(tagodaybgc)
    write(*, *) ' tag monthly ice  = ', trim(tagimon)
    write(*, *) ' tag daily ice    = ', trim(tagiday)
    write(*, *) ' tag monthly atm  = ', trim(tagamon)
    write(*, *) ' tag daily atm    = ', trim(tagaday)
    write(*, *) ' tag 6hourly atm  = ', trim(taga6hr)
    write(*, *) ' tag 6hourly insa = ', trim(taga6hri)
    write(*, *) ' tag 3hourly atm  = ', trim(taga3hr)
    write(*, *) ' tag 3hourly insa = ', trim(taga3hri)
    write(*, *) ' tag yearly lnd   = ', trim(taglyr)
    write(*, *) ' tag monthly lnd  = ', trim(taglmon)
    write(*, *) ' tag daily lnd    = ', trim(taglday)
    write(*, *) ' tag 3hourly lnd  = ', trim(tagl3hr)
    write(*, *) ' tag 3hourly insl = ', trim(tagl3hri)
    write(*, *) ' atmos grid file  = ', trim(atmgridfile)
    write(*, *) ' ocean grid file  = ', trim(ocngridfile)
    write(*, *) ' ocean ini file   = ', trim(ocninitfile)
    write(*, *) ' ocean sec file   = ', trim(secindexfile)
    write(*, *) ' ocean moc file   = ', trim(ocnmertfile)
    write(*, *) ' ocean reg file   = ', trim(ocnregnfile)
    write(*, *) ' allow line break = ', linebreaks

    write(*, *)
    write(*, *) 'Experiment namelist:'
    write(*, *) ' case name        = ', trim(casename)
    write(*, *) ' experiment id    = ', trim(experiment_id)
    write(*, *) ' history          = ', trim(history1)
    write(*, *) ' comment          = ', trim(comment1)
    write(*, *) ' forcing          = ', trim(forcing1)
    write(*, *) ' realization      = ', realization
    write(*, *) ' start year       = ', year1
    write(*, *) ' end year         = ', yearn
    write(*, *) ' start month      = ', month1
    write(*, *) ' end month        = ', monthn
    write(*, *) ' add dummy day    = ', add_fill_day
    write(*, *) ' dry run          = ', dry_run

    write(*, *)
    write(*, *) 'Table grids:'
    write(*, *) ' grid table file  = ', trim(tgrids)

    write(*, *)
    write(*, *) 'Table fx:'
    do n = 1, nfx
      write(*, '(1X,A20,A20,A20)') vfx(:, n)
    end do

    write(*, *)
    write(*, *) 'Table yr:'
    do n = 1, noyr
      write(*, '(1X,A20,A20,A20)') voyr(:, n)
    end do

    write(*, *)
    write(*, *) 'Table yr (bgc):'
    do n = 1, noyrbgc
      write(*, '(1X,A20,A20,A20)') voyrbgc(:, n)
    end do

    write(*, *)
    write(*, *) 'Table omon:'
    do n = 1, nomon
      write(*, '(1X,A20,A20,A20)') vomon(:, n)
    end do

    write(*, *)
    write(*, *) 'Table omon (bgc):'
    do n = 1, nomonbgc
      write(*, '(1X,A20,A20,A20)') vomonbgc(:, n)
    end do

    write(*, *)
    write(*, *) 'Table oimon:'
    do n = 1, noimon
      write(*, '(1X,A20,A20,A20)') voimon(:, n)
    end do

    write(*, *)
    write(*, *) 'Table amon:'
    do n = 1, namon
      write(*, '(1X,A20,A20,A20)') vamon(:, n)
    end do

    write(*, *)
    write(*, *) 'Table aero:'
    do n = 1, naero
      write(*, '(1X,A20,A20,A20)') vaero(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Lmon:'
    do n = 1, nlmon
      write(*, '(1X,A20,A20,A20)') vlmon(:, n)
    end do

    write(*, *)
    write(*, *) 'Table LImon:'
    do n = 1, nlimon
      write(*, '(1X,A20,A20,A20)') vlimon(:, n)
    end do

    write(*, *)
    write(*, *) 'Table day:'
    do n = 1, nday
      write(*, '(1X,A20,A20,A20)') vday(:, n)
    end do

    write(*, *)
    write(*, *) 'Table 6hrlev:'
    do n = 1, n6hrlev
      write(*, '(1X,A20,A20,A20)') v6hrlev(:, n)
    end do

    write(*, *)
    write(*, *) 'Table 6hrlevi:'
    do n = 1, n6hrlevi
      write(*, '(1X,A20,A20,A20)') v6hrlevi(:, n)
    end do

    write(*, *)
    write(*, *) 'Table 6hrplev:'
    do n = 1, n6hrplev
      write(*, '(1X,A20,A20,A20)') v6hrplev(:, n)
    end do

    write(*, *)
    write(*, *) 'Table 3hr:'
    do n = 1, n3hr
      write(*, '(1X,A20,A20,A20)') v3hr(:, n)
    end do

    write(*, *)
    write(*, *) 'Table 3hri:'
    do n = 1, n3hri
      write(*, '(1X,A20,A20,A20)') v3hri(:, n)
    end do

    write(*, *)
    write(*, *) 'Table ofx:'
    do n = 1, nofx
      write(*, '(1X,A20,A20,A20)') vofx(:, n)
    end do

    write(*, *)
    write(*, *) 'Table 6hrPlevPt:'
    do n = 1, n6hrPlevPt
      write(*, '(1X,A20,A20,A20)') v6hrPlevPt(:, n)
    end do

    write(*, *)
    write(*, *) 'Table AERday:'
    do n = 1, nAERday
      write(*, '(1X,A20,A20,A20)') vAERday(:, n)
    end do

    write(*, *)
    write(*, *) 'Table AERhr:'
    do n = 1, nAERhr
      write(*, '(1X,A20,A20,A20)') vAERhr(:, n)
    end do

    write(*, *)
    write(*, *) 'Table AERmonZ:'
    do n = 1, nAERmonZ
      write(*, '(1X,A20,A20,A20)') vAERmonZ(:, n)
    end do

    write(*, *)
    write(*, *) 'Table CF3hr:'
    do n = 1, nCF3hr
      write(*, '(1X,A20,A20,A20)') vCF3hr(:, n)
    end do

    write(*, *)
    write(*, *) 'Table CFday:'
    do n = 1, nCFday
      write(*, '(1X,A20,A20,A20)') vCFday(:, n)
    end do

    write(*, *)
    write(*, *) 'Table CFmon:'
    do n = 1, nCFmon
      write(*, '(1X,A20,A20,A20)') vCFmon(:, n)
    end do

    write(*, *)
    write(*, *) 'Table CFsubhr:'
    do n = 1, nCFsubhr
      write(*, '(1X,A20,A20,A20)') vCFsubhr(:, n)
    end do

    write(*, *)
    write(*, *) 'Table E1hrClimMon:'
    do n = 1, nE1hrClimMon
      write(*, '(1X,A20,A20,A20)') vE1hrClimMon(:, n)
    end do

    write(*, *)
    write(*, *) 'Table E1hr:'
    do n = 1, nE1hr
      write(*, '(1X,A20,A20,A20)') vE1hr(:, n)
    end do

    write(*, *)
    write(*, *) 'Table E3hr:'
    do n = 1, nE3hr
      write(*, '(1X,A20,A20,A20)') vE3hr(:, n)
    end do

    write(*, *)
    write(*, *) 'Table E3hrPt:'
    do n = 1, nE3hrPt
      write(*, '(1X,A20,A20,A20)') vE3hrPt(:, n)
    end do

    write(*, *)
    write(*, *) 'Table E6hrZ:'
    do n = 1, nE6hrZ
      write(*, '(1X,A20,A20,A20)') vE6hrZ(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Eday:'
    do n = 1, nEday
      write(*, '(1X,A20,A20,A20)') vEday(:, n)
    end do

    write(*, *)
    write(*, *) 'Table EdayZ:'
    do n = 1, nEdayZ
      write(*, '(1X,A20,A20,A20)') vEdayZ(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Efx:'
    do n = 1, nEfx
      write(*, '(1X,A20,A20,A20)') vEfx(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Emon:'
    do n = 1, nEmon
      write(*, '(1X,A20,A20,A20)') vEmon(:, n)
    end do

    write(*, *)
    write(*, *) 'Table EmonZ:'
    do n = 1, nEmonZ
      write(*, '(1X,A20,A20,A20)') vEmonZ(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Esubhr:'
    do n = 1, nEsubhr
      write(*, '(1X,A20,A20,A20)') vEsubhr(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Eyr:'
    do n = 1, nEyr
      write(*, '(1X,A20,A20,A20)') vEyr(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Oclim:'
    do n = 1, nOclim
      write(*, '(1X,A20,A20,A20)') vOclim(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Oday:'
    do n = 1, nOday
      write(*, '(1X,A20,A20,A20)') vOday(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Odaybgc:'
    do n = 1, nOdaybgc
      write(*, '(1X,A20,A20,A20)') vOdaybgc(:, n)
    end do

    write(*, *)
    write(*, *) 'Table Odec:'
    do n = 1, nOdec
      write(*, '(1X,A20,A20,A20)') vOdec(:, n)
    end do

    write(*, *)
    write(*, *) 'Table COfx:'
    do n = 1, nCOfx
      write(*, '(1X,A20,A20,A20)') vCOfx(:, n)
    end do

    write(*, *)
    write(*, *) 'Table SIday:'
    do n = 1, nSIday
      write(*, '(1X,A20,A20,A20)') vSIday(:, n)
    end do

  end subroutine print_namelists

  ! -----------------------------------------------------------------

  subroutine merge_strarr(slen, sdm, strin, strout, lb)

    implicit none

    integer, intent(in) :: sdm, slen
    character(len=slen), dimension(sdm), intent(in) :: strin
    character(len=(slen+1)*sdm), intent(out) :: strout
    logical, intent(in) :: lb

    integer :: n, pos

    strout = ' '
    pos = 0
    do n = 1, sdm
      if (len_trim(strin(n)) > 0) then
        if (pos /= 0) then
          pos = pos + 1
          if (lb) then
            strout(pos:pos) = achar(10)
          else
            strout(pos:pos) = ' '
          end if
        end if
        strout(pos+1:pos+len_trim(strin(n))) = trim(strin(n))
        pos = pos + len_trim(strin(n))
      end if
    end do
    pos = 1
    do
      if (index(strout(pos:), '\n') > 0) then
        pos = pos + index(strout(pos:), '\n')
        strout(pos-1:pos) = ' '//achar(10)
        pos = pos + 2
        if (pos >= len(strout)) exit
      else
        exit
      end if
    end do

  end subroutine merge_strarr

  ! -----------------------------------------------------------------

  function get_free_unit() result(free_unit)

    implicit none

    integer :: free_unit
    logical :: in_use

    do free_unit = 10, 99
      inquire(free_unit, opened=in_use)
      if (.not. in_use) exit
    end do

  end function get_free_unit

  ! -----------------------------------------------------------------

  subroutine write_namelist_json(grid, grid_label, grid_resolution, varname)

    use json_module

    implicit none

#ifdef MPI
    include 'mpif.h'
    integer :: mpirank, mpisize, mpierror
#endif

    character(len=*), intent(in) :: grid, grid_label, grid_resolution, varname
    character :: yyyymm1*6, yyyymm2*6, c2*2, r3*3
    type(json_core) :: json
    type(json_value), pointer :: p

    call json%initialize()
    call json%create_object(p, '')

    ! mapped from cmor2
    call json%add(p, 'outpath', trim(obasedir))
    call json%add(p, 'experiment_id', trim(experiment_id))
    call json%add(p, 'institution_id', trim(institute_id))
    call json%add(p, 'institution', trim(institution1))
    call json%add(p, 'source_id', trim(model_id))
    call json%add(p, 'source', trim(source1))
    call json%add(p, 'calendar', trim(calendar))
    call json%add(p, 'realization_index', trim(realization_index))
    call json%add(p, 'physics_index', trim(physics_index))
    call json%add(p, 'initialization_index', trim(initialization_index))
    call json%add(p, 'contact', trim(contact1))
    call json%add(p, 'history', trim(history1))
    call json%add(p, 'comment', trim(comment1))
    call json%add(p, 'references', trim(references1))
    call json%add(p, 'model_id', trim(model_id))
    call json%add(p, 'run_variant', trim(forcing1))
    call json%add(p, 'branch_time', branch_time)
    call json%add(p, 'parent_experiment_id', trim(parent_experiment_id))
    ! new for cmor3
    call json%add(p, 'forcing_index', trim(forcing_index))
    call json%add(p, 'parent_variant_label', trim(parent_variant_label))
    call json%add(p, '_controlled_vocabulary_file', 'cmor-cvs.json')
    call json%add(p, '_AXIS_ENTRY_FILE', 'CMIP7_coordinate.json')
    call json%add(p, '_FORMULA_VAR_FILE', 'CMIP7_formula_terms.json')
    call json%add(p, '_cmip7_option', 1)
    call json%add(p, 'activity_id', trim(activity_id))
    call json%add(p, 'source_type', trim(source_type))
    call json%add(p, 'sub_experiment_id', trim(sub_experiment_id))
    call json%add(p, 'parent_sub_experiment_id', trim(parent_sub_experiment))
    call json%add(p, 'parent_mip_era', trim(parent_mip_era))
    call json%add(p, 'mip_era', trim(mip_era))
    call json%add(p, 'parent_activity_id', trim(parent_activity_id))
    call json%add(p, 'parent_source_id', trim(parent_source_id))
    call json%add(p, 'grid', trim(grid))
    call json%add(p, 'grid_label', trim(grid_label))
    call json%add(p, 'nominal_resolution', trim(grid_resolution))
    call json%add(p, 'branch_method', trim(branch_method))
    call json%add(p, 'branch_time_in_child', branch_time_in_child)
    call json%add(p, 'branch_time_in_parent', branch_time_in_parent)
    call json%add(p, 'parent_time_units', trim(parent_time_units))
    call json%add(p, 'tracking_prefix', trim(tracking_prefix))
    call json%add(p, 'output_path_template', &
      '<activity_id><institution_id><source_id><experiment_id><_member_id><table><variable_id><grid_label><version>')
    call json%add(p, 'output_file_template', &
      '<variable_id><table><source_id><experiment_id><_member_id><grid_label>')
    call json%add(p, 'license_id', 'CC-BY-4-0')
    call json%add(p, 'archive_id', 'WCRP')
    call json%add(p, 'frequency', 'mon')
    call json%add(p, 'region', 'glb')
    call json%add(p, 'branded_variable', trim(varname))
    call json%add(p, 'drs_specs', 'MIP-DRS7')

    write(yyyymm1, '(I4.4,I2.2)') year1, month1
    write(yyyymm2, '(I4.4,I2.2)') yearn, monthn
    write(r3, '(I3.3)') realization
#ifdef MPI
    call mpi_comm_rank(mpi_comm_world, mpirank, mpierror)
    write(c2, '(I2.2)') mpirank
#else
    c2 = '00'
#endif

    namelist_file_json = 'namelist_'//trim(casename)//'_'// &
      trim(varname)//'_'//yyyymm1//'-'//yyyymm2//'_'//r3//'_'//c2//'.json'
    call json%print(p, namelist_file_json)
    call json%destroy(p)            ! cleanup
    write(*, *) 'json-namelist file: ', trim(namelist_file_json)

  end subroutine write_namelist_json

end module m_namelists