!-------------------------------------------------------------------------------
! NNLOJET
!-------------------------------------------------------------------------------
!
! MODULE: Observables [external members: `_obs` and `_jet`]
!
! DESCRIPTION:
!> Module to facilitate the evaluation of observables and their management.
!>
!> EvalObs_mod is a collection of all evaluation functions.
!>
!> Observables_mod registers the observables for different process types
!> and supplies a function to compute the whole set and store the result in a
!> cache. Each Observable_t has a name and an identifier for it's
!> evaluation. Initialization allows to pick the relevant observables and
!> associate each one with a unique name. Here, we can also pick more
!> suitable names for the same evaluation function such as:
!> -----
!>   eval_minv34 : invariant of particles 3 & 4
!>     > for process "ZJ" map it to the name `mll` [lepton-pair]
!>     > for process "WJ" map it to the name `mln` [lepton & neutrino]
!>   eval_pT4 : transverse momentum of 4
!>     > for process "ZJ" map it to the name `pTl+` [lepton^+]
!>     > for process "WJ" map it to the name `ETmiss` [neutrino: missing ET]
!> -----
!> Any attempt to evaluate an unregistered observable will error out, which
!> is a desirable feature.
!>
!> Module also performs the successive jet recombination algorithm.
!> The interface is intentionally chosen very similar to FastJet to allow for
!> an easy transition to the external library if desired.
!> In order to use FastJet compile with `make jet=fastjet`
!> Native implementation was validated against FastJet and performs better
!
!> We also accommodate certain jet observables that cannot be used in the
!> histogramming but can be used in selectors and are only active during
!> the jet recombination. These observables receive negative identifiers
!> and are treated separately in Selectors_mod and directly passed to
!> the Jet algorithm
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
module EvalObs_mod
!-------------------------------------------------------------------------------
  use KinData_mod
  implicit none
  public

  private :: list_cache
  
  !--------------------------!
  !  Observable identifiers  !
  !--------------------------!

  ! number of reconstructed jets
  integer, parameter :: id_njets                 = 1
  ! number of FS partons = `nproto`
  integer, parameter :: id_npartons              = 2
  ! 1st jet (leading)
  integer, parameter :: id_pt_j1                 = 3
  integer, parameter :: id_y_j1                  = 4
  integer, parameter :: id_abs_y_j1              = 5
  integer, parameter :: id_eta_j1                = 6
  integer, parameter :: id_abs_eta_j1            = 7
  ! 2nd jet (sub-leading)
  integer, parameter :: id_pt_j2                 = 8
  integer, parameter :: id_y_j2                  = 9
  integer, parameter :: id_abs_y_j2              = 10
  integer, parameter :: id_eta_j2                = 11
  integer, parameter :: id_abs_eta_j2            = 12
  ! 3rd jet
  integer, parameter :: id_pt_j3                 = 13
  integer, parameter :: id_y_j3                  = 14
  integer, parameter :: id_abs_y_j3              = 15
  integer, parameter :: id_eta_j3                = 16
  integer, parameter :: id_abs_eta_j3            = 17
  ! 4th jet
  integer, parameter :: id_pt_j4                 = 18
  integer, parameter :: id_y_j4                  = 19
  integer, parameter :: id_abs_y_j4              = 20
  integer, parameter :: id_eta_j4                = 21
  integer, parameter :: id_abs_eta_j4            = 22
  ! sum pt
  integer, parameter :: id_sum_pt_jets           = 23
  ! di-jet invariant mass (j1+j2)
  integer, parameter :: id_minv_j12              = 24
  ! tri-jet invariant mass (j1+j2+j3)
  integer, parameter :: id_minv_j123             = 25
  ! more di-jet observables
  integer, parameter :: id_deltapt_j12           = 26
  integer, parameter :: id_deltaeta_j12          = 27
  integer, parameter :: id_deltay_j12            = 28
  integer, parameter :: id_deltaphi_j12          = 29
  integer, parameter :: id_chi_j12               = 30
  integer, parameter :: id_ystar_j12             = 31
  integer, parameter :: id_yboost_j12            = 32
  integer, parameter :: id_max_y_j12             = 33
  ! avergae of pt's
  integer, parameter :: id_ptavg_j12             = 34
  integer, parameter :: id_ptavg_j123            = 35
  integer, parameter :: id_ptavg_jall            = 36
  integer, parameter :: id_ptavg_geom_jall       = 37
  integer, parameter :: id_etavg_j12             = 38
  integer, parameter :: id_etavg_j123            = 39
  integer, parameter :: id_etavg_jall            = 40
  integer, parameter :: id_etavg_geom_jall       = 41
  integer, parameter :: id_et2_j12               = 42
  ! H_T (only jet momenta)
  integer, parameter :: id_ht_jets               = 43
  integer, parameter :: id_ht_part               = 44
  ! tau Jet observables (used in H+J)
  integer, parameter :: id_tau_j1                = 45
  integer, parameter :: id_sum_tau_jets          = 46
  integer, parameter :: id_max_tau_jet           = 47
  !>----- "perfect" jet observables (no cuts applied)
  integer, parameter :: id_pt_j1_nocut           = 48
  integer, parameter :: id_pt_j2_nocut           = 49
  integer, parameter :: id_ht_jall_nocut         = 50
  integer, parameter :: id_ptavg_j12_nocut       = 51
  integer, parameter :: id_ptavg_jall_nocut      = 52
  integer, parameter :: id_ptavg_geom_j12_nocut  = 53
  integer, parameter :: id_ptavg_geom_jall_nocut = 54
  !>----- jet observables using R=1 jets
  integer, parameter :: id_pt_j1_R1              = 55
  integer, parameter :: id_pt_j2_R1              = 56
  integer, parameter :: id_ht_jall_R1            = 57
  integer, parameter :: id_ptavg_j12_R1          = 58
  integer, parameter :: id_ptavg_jall_R1         = 59
  integer, parameter :: id_ptavg_geom_j12_R1     = 60
  integer, parameter :: id_ptavg_geom_jall_R1    = 61
  !>----- gauge-boson production
  ! 1st lepton [npar-1]
  integer, parameter :: id_pt_l1                 = 62
  integer, parameter :: id_y_l1                  = 63
  integer, parameter :: id_abs_y_l1              = 64
  ! 2nd lepton [npar]
  integer, parameter :: id_pt_l2                 = 65
  integer, parameter :: id_y_l2                  = 66
  integer, parameter :: id_abs_y_l2              = 67
  ! 3rd lepton [npar-3]
  integer, parameter :: id_pt_l3                 = 68
  integer, parameter :: id_y_l3                  = 69
  integer, parameter :: id_abs_y_l3              = 70
  ! 4th lepton [npar-2]
  integer, parameter :: id_pt_l4                 = 71
  integer, parameter :: id_y_l4                  = 72
  integer, parameter :: id_abs_y_l4              = 73
  ! gauge boson: l1+l2
  integer, parameter :: id_pt_V                  = 74
  integer, parameter :: id_minv_V                = 75
  integer, parameter :: id_y_V                   = 76
  integer, parameter :: id_abs_y_V               = 77
  integer, parameter :: id_eta_V                 = 78
  integer, parameter :: id_abs_eta_V             = 79
  integer, parameter :: id_mt_V                  = 80
  integer, parameter :: id_et_V                  = 81
  ! gauge boson: l3+l4
  integer, parameter :: id_pt_V2                 = 82
  integer, parameter :: id_minv_V2               = 83
  integer, parameter :: id_y_V2                  = 84
  integer, parameter :: id_abs_y_V2              = 85
  integer, parameter :: id_eta_V2                = 86
  integer, parameter :: id_abs_eta_V2            = 87
  integer, parameter :: id_mt_V2                 = 88
  integer, parameter :: id_et_V2                 = 89
  ! gauge boson: l1+l2+l3+l4
  integer, parameter :: id_pt_V4l                = 90
  integer, parameter :: id_minv_V4l              = 91
  integer, parameter :: id_y_V4l                 = 92
  integer, parameter :: id_minv_V23              = 93
  integer, parameter :: id_minv_V14              = 94
  integer, parameter :: id_abs_y_V4l             = 95
  integer, parameter :: id_eta_V4l               = 96
  integer, parameter :: id_abs_eta_V4l           = 97
  integer, parameter :: id_minv_mZ1              = 98
  integer, parameter :: id_minv_mZ2              = 99
  integer, parameter :: id_minv_mZSFOS1st        = 100
  integer, parameter :: id_minv_mZSFOS2nd        = 101
  integer, parameter :: id_minv_mZSFOS3rd        = 102
  integer, parameter :: id_minv_mZSFOS4th        = 103
  ! ordering switch of the SFOS lepton pairs close to mZ
  integer, parameter :: id_minvSFOS_vs_mZ        = 104
  ! Collins-Soper frame (angles wrt l1)
  integer, parameter :: id_cos_theta_CS          = 105
  integer, parameter :: id_cos_theta_CS_0PT      = 106
  integer, parameter :: id_theta_CS              = 107
  integer, parameter :: id_phi_CS                = 108
  ! projectors for the angular coefficients A_i [1606.00689]
  integer, parameter :: id_proj_A0               = 109
  integer, parameter :: id_proj_A1               = 110
  integer, parameter :: id_proj_A2               = 111
  integer, parameter :: id_proj_A3               = 112
  integer, parameter :: id_proj_A4               = 113
  integer, parameter :: id_proj_A5               = 114
  integer, parameter :: id_proj_A6               = 115
  integer, parameter :: id_proj_A7               = 116
  integer, parameter :: id_proj_A0mA2            = 117
  integer, parameter :: id_proj_A0mA2alt         = 118
  ! Vjj (leading + sub-leading jet)
  integer, parameter :: id_pt_Vjj                = 119
  ! Vj ( V + leading jet)
  integer, parameter :: id_pt_Vj                 = 120
  integer, parameter :: id_ystar_Vj              = 121
  integer, parameter :: id_yboost_Vj             = 122
  !>----- H -> gamma gamma decay
  ! leading gamma [npar, npar-1]
  integer, parameter :: id_idx_g1                = 123
  integer, parameter :: id_pt_g1                 = 124
  integer, parameter :: id_y_g1                  = 125
  integer, parameter :: id_abs_y_g1              = 126
  ! sub-leading gamma [npar, npar-1]
  integer, parameter :: id_idx_g2                = 127
  integer, parameter :: id_pt_g2                 = 128
  integer, parameter :: id_y_g2                  = 129
  integer, parameter :: id_abs_y_g2              = 130
  ! gamma-gamma separation [npar, npar-1]
  integer, parameter :: id_deltaphi_g1g2         = 131
  integer, parameter :: id_deltay_g1g2           = 132
  integer, parameter :: id_abs_deltay_g1g2       = 133
  integer, parameter :: id_ptt_g1g2              = 134
  !>----- H -> ZZ --> 4 l decay
  ! leading lepton [npar, npar-1, npar-2, npar-3]
  integer, parameter :: id_pt_l1st               = 135
  integer, parameter :: id_y_l1st                = 136
  integer, parameter :: id_abs_y_l1st            = 137
  integer, parameter :: id_R_l1st_j1             = 138
  ! sub-leading lepton [npar, npar-1, npar-2, npar-3]
  integer, parameter :: id_pt_l2nd               = 139
  integer, parameter :: id_y_l2nd                = 140
  integer, parameter :: id_abs_y_l2nd            = 141
  integer, parameter :: id_R_l2nd_j1             = 142
  ! third-leading lepton [npar, npar-1, npar-2, npar-3]
  integer, parameter :: id_pt_l3rd               = 143
  integer, parameter :: id_y_l3rd                = 144
  integer, parameter :: id_abs_y_l3rd            = 145
  integer, parameter :: id_R_l3rd_j1             = 146
  ! least-leading lepton [npar, npar-1, npar-2, npar-3]
  integer, parameter :: id_pt_l4th               = 147
  integer, parameter :: id_y_l4th                = 148
  integer, parameter :: id_abs_y_l4th            = 149
  integer, parameter :: id_R_l4th_j1             = 150
  !>----- more
  integer, parameter :: id_ht_full               = 151
  integer, parameter :: id_min_dR_l1j            = 152
  integer, parameter :: id_min_dR_l2j            = 153
  integer, parameter :: id_min_dR_l3j            = 154
  integer, parameter :: id_min_dR_l4j            = 155
  integer, parameter :: id_min_dR_l12j           = 156
  integer, parameter :: id_min_dR_l34j           = 157
  integer, parameter :: id_min_dR_l1234j         = 158
  integer, parameter :: id_dR_l12                = 159
  integer, parameter :: id_dR_l34                = 160
  integer, parameter :: id_min_dR_l_sf           = 161
  integer, parameter :: id_min_dR_l_df           = 162
  integer, parameter :: id_phi_star_eta          = 163
  integer, parameter :: id_costhetastar          = 164
  integer, parameter :: id_deltaphi_Vj1          = 165
  integer, parameter :: id_deltay_Vj1            = 166
  integer, parameter :: id_ydiff_Vj1             = 167
  integer, parameter :: id_ysum_Vj1              = 168
  !>----- dynamical scales for Z+jets
  integer, parameter :: id_ZJ_dscale1            = 169
  integer, parameter :: id_ZJ_dscale3            = 170
  integer, parameter :: id_ZJ_HT                 = 171
  integer, parameter :: id_ZJ_HTshape05          = 172
  integer, parameter :: id_ZJ_HTshape20          = 173
  !>----- dynamical scales for Hto4l+jets
  integer, parameter :: id_H4lJ_dscale1          = 174
  integer, parameter :: id_H4lJ_dscale2          = 175
  !>----- VFH_cuts observables
  integer, parameter :: id_VFH_deltay            = 176
  integer, parameter :: id_VFH_y1xy2             = 177
  !>----- dynamical scale  for VFH
  integer, parameter :: id_VFH_dscale            = 178
  !>----- observables for VFH
  integer, parameter :: id_VFH_z3                = 179
  !>----- observables for DIS
  !> global observables
  integer, parameter :: id_DIS_Q2                = 180
  integer, parameter :: id_DIS_x                 = 181
  integer, parameter :: id_DIS_y                 = 182
  integer, parameter :: id_DIS_cosgammah         = 183
  integer, parameter :: id_DIS_sqrtQ2            = 184
  integer, parameter :: id_DIS_W2                = 185
  integer, parameter :: id_DIS_W                 = 186
  integer, parameter :: id_DIS_nu                = 187
  
  !> DIS scale choices
  integer, parameter :: id_DIS_dscl              = 188
  integer, parameter :: id_DIS_dsclj1            = 189
  integer, parameter :: id_DIS_dsclj2            = 190
  integer, parameter :: id_DIS_dsclj3            = 191
  integer, parameter :: id_DIS_dsclj4            = 192
  integer, parameter :: id_DIS_dsclj1l           = 193
  integer, parameter :: id_DIS_dsclj2l           = 194
  integer, parameter :: id_DIS_dsclj3l           = 195
  integer, parameter :: id_DIS_dsclj4l           = 196
  integer, parameter :: id_DIS_dsclZEUS          = 197
  integer, parameter :: id_DIS_dsclZEUS2         = 198
  integer, parameter :: id_DIS_dscl3j            = 199

  !> observables in the Breit frame
  !> assume above observable ID's for Breit frame
  integer, parameter :: id_DIS_visbyW            = 200
  integer, parameter :: id_DIS_etas              = 201
  integer, parameter :: id_DIS_Lgxi2             = 202
  integer, parameter :: id_DIS_xi2               = 203
  integer, parameter :: id_DIS_xi3               = 204
  integer, parameter :: id_DIS_Thrust            = 205
  integer, parameter :: id_DIS_Thrust_c          = 206
  integer, parameter :: id_DIS_JB                = 207
  integer, parameter :: id_DIS_JM2               = 208
  integer, parameter :: id_DIS_C                 = 209
  integer, parameter :: id_DIS_resrho            = 210
  !> Observables in proton gamma frame
  integer, parameter :: id_gammap_eta_j1         = 211
  integer, parameter :: id_gammap_eta_j2         = 212
  integer, parameter :: id_DIS_delta_phis        = 213
  !> Observables in the HERA frame 
  integer, parameter :: id_DIS_hera_deltaeta_j12 = 214
  integer, parameter :: id_DIS_hera_etaavg_j12   = 215
  integer, parameter :: id_DIS_hera_etj1         = 216
  integer, parameter :: id_DIS_hera_etj2         = 217
  integer, parameter :: id_DIS_hera_etj3         = 218
  integer, parameter :: id_DIS_hera_etj4         = 219
  integer, parameter :: id_DIS_hera_etaj1        = 220
  integer, parameter :: id_DIS_hera_etaj2        = 221
  integer, parameter :: id_DIS_hera_etaj3        = 222
  integer, parameter :: id_DIS_hera_etaj4        = 223
  integer, parameter :: id_DIS_hera_et2j1        = 224
  integer, parameter :: id_DIS_hera_et2j2        = 225
  integer, parameter :: id_DIS_hera_et2j3        = 226
  integer, parameter :: id_DIS_hera_et2j4        = 227
  integer, parameter :: id_DIS_hera_et2byQ2j1    = 228
  integer, parameter :: id_DIS_hera_et2byQ2j2    = 229
  integer, parameter :: id_DIS_hera_et2byQ2j3    = 230
  integer, parameter :: id_DIS_hera_et2byQ2j4    = 231
  integer, parameter :: id_DIS_hera_deta1        = 232
  integer, parameter :: id_DIS_hera_deta2        = 233

  !> some special observables
  integer, parameter :: id_DIS_hera_eta_j1_0     = 234
  integer, parameter :: id_DIS_hera_eta_j2_0     = 235
  integer, parameter :: id_DIS_hera_xgamma       = 236

  !> DIS reweighting
  integer, parameter :: id_DIS_hera_FJ_reweight  = 237

  !> epem Observables
  integer, parameter :: id_epem_C                = 238
  integer, parameter :: id_epem_D                = 239
  integer, parameter :: id_epem_WJB              = 240
  integer, parameter :: id_epem_TJB              = 241
  integer, parameter :: id_epem_T                = 242
  integer, parameter :: id_epem_1mT              = 243
  integer, parameter :: id_epem_cosT             = 244
  integer, parameter :: id_epem_HJM              = 245
  integer, parameter :: id_epem_SHM              = 246
  integer, parameter :: id_epem_costh_ej1        = 247
  integer, parameter :: id_epem_costh_en3j       = 248
  integer, parameter :: id_epem_chi              = 249
  !> ycut transitions in the JADE algorithm (for DIS & epem)
  integer, parameter :: id_y45                   = 250
  integer, parameter :: id_y34                   = 251
  integer, parameter :: id_y23                   = 252
  integer, parameter :: id_y12                   = 253
  integer, parameter :: id_y01                   = 254

  !>----- EXTRA observables
  !> One cannot fill histograms with extra obs

  !> DIS Jet-Observable identifiers in HERA frame
  integer, parameter :: id_hera_ptj              = 255
  integer, parameter :: id_hera_yj               = 256
  integer, parameter :: id_hera_abs_yj           = 257
  integer, parameter :: id_hera_etaj             = 258
  integer, parameter :: id_hera_abs_etaj         = 259
  integer, parameter :: id_hera_etj              = 260
  integer, parameter :: id_hera_xjets            = 261
  integer, parameter :: id_hera_H1xjets          = 262
  integer, parameter :: id_DIS_ptl               = 263

  !> photon cuts
  integer, parameter :: id_photon_pt             = 264
  integer, parameter :: id_photon_y              = 265
 !> MiNLO 
  integer, parameter :: id_sudakov               = 266
  integer, parameter :: id_minlo                 = 267

  !> Spike plots
  integer, parameter :: id_s13                   = 268
  integer, parameter :: id_s14                   = 269
  integer, parameter :: id_s15                   = 270
  integer, parameter :: id_s16                   = 271
  integer, parameter :: id_s23                   = 272
  integer, parameter :: id_s24                   = 273
  integer, parameter :: id_s25                   = 274
  integer, parameter :: id_s26                   = 275
  integer, parameter :: id_s34                   = 276
  integer, parameter :: id_s35                   = 277
  integer, parameter :: id_s36                   = 278
  integer, parameter :: id_s45                   = 279
  integer, parameter :: id_s46                   = 280
  integer, parameter :: id_s56                   = 281

  !> Dijet scale choices
  integer, parameter :: id_pt_j1expYstar         = 282

  ! total number of observables
  integer, parameter :: n_obs = 282
  
  integer, parameter :: max_obs_name_length = 32
  integer, parameter :: max_obs_desc_length = 256

  !------------------------------!
  !  Jet-Observable identifiers  !
  !------------------------------!

  integer, parameter :: idj_pt           = 1
  integer, parameter :: idj_y            = 2
  integer, parameter :: idj_abs_y        = 3
  integer, parameter :: idj_eta          = 4
  integer, parameter :: idj_abs_eta      = 5
  integer, parameter :: idj_et           = 6
  !> process-specific jet selectors:
  integer, parameter :: idj_dR_l1        = 7
  integer, parameter :: idj_dR_l2        = 8
  integer, parameter :: idj_dR_l3        = 9
  integer, parameter :: idj_dR_l4        = 10
  integer, parameter :: idj_min_dR_l12   = 11
  integer, parameter :: idj_min_dR_l34   = 12
  integer, parameter :: idj_min_dR_l1234 = 13
  !> number of jet-observables
  integer, parameter :: n_obsj = 13


  !-------------------!
  !  Jet Observables  !
  !-------------------!

  public :: setEtrecom
  public :: init_jet, algo_jet, setRcone_jet
  public :: setSelector_jet
  public :: digest_jet
  public :: min_shat_jet
  public :: initBuffer_jet, destroyBuffer_jet
#ifdef USEFASTJET
  public :: cluster_fastjet
#endif

  logical :: Et_recom =.false.

  !> flag choosing the algorithm
  !>   0: none
  !>   1: anti-kt
  !>   2: C/A
  !>   3: kt
  !>   4: jade
  integer :: jetalg = -1
  integer :: ip_jet ! power of the generalized kt algorithm
  !> R parameter of the algorithm set in the runcard
  double precision :: inputRcone = -1d0
  !> ycut in the jade algorithm
  double precision :: jade_ycut = -1d0
  !> threadlocal copies for R [& squared]
  !> typically will be set to inputRcone but multi-runs can 
  !> have their own values
  double precision :: Rcone   = -1d0
  double precision :: RconeSq = -1d0
!$OMP THREADPRIVATE(Rcone, RconeSq)

  !> type to collect information needed in the recombination
  type Jet_t
    !> flag to identify status
    !>   < 0: ProtoJet (subject to further recombinations)
    !>   = 0: Nothing (merged into another Protojet)
    !>   > 0: Jet (identified as a hard jet)
    integer :: isJet
    double precision, dimension(4) :: p
    double precision :: pt2, y, phi
  end type Jet_t

  !> threadprivate buffer to minimize re-allocation.
  type(Jet_t), allocatable, dimension(:) :: jets
  double precision, allocatable, dimension(:) :: diB  ! "B" = beams
  double precision, allocatable, dimension(:,:) :: dij
  !> information on currently stored data
  integer :: npar_jet, nproto_jet, njets_jet
  !> store the clustering history
  !> if parton `i` is merged into `j`:
  !>  clusterHist_jet(i) = j
  !> final jet will have
  !>  clusterHist_jet(i) = i
  integer, allocatable, dimension(:) :: clusterHist_jet
  !> ordering according to pt (descending)
  !>  leading jet: pt_sort_jet(1)
  !>  => ptj1 = jets(pt_sort_jet(1))%pt
  !>  etc...
  integer, allocatable, dimension(:) :: pt_sort_jet
!$OMP THREADPRIVATE(jets, diB, dij)
!$OMP THREADPRIVATE(npar_jet, nproto_jet, njets_jet)
!$OMP THREADPRIVATE(clusterHist_jet, pt_sort_jet)

  type JetSelector_t
    logical :: qactive = .false.
    character(len=max_obs_name_length) :: name
    double precision :: lower = 0d0
    double precision :: upper = 0d0
    logical :: qlower = .false.
    logical :: qupper = .false.
    logical :: qinvert = .false.
  end type JetSelector_t

  type(JetSelector_t), dimension(n_obsj) :: list_obsj


  !----------------------------------!
  !  Cache for computed observables  !
  !----------------------------------!

  type EvalCache_t
    logical :: qcached = .false.
    double precision :: value
  end type EvalCache_t

  type(EvalCache_t), dimension(n_obs) :: list_cache
!$OMP THREADPRIVATE(list_cache)


contains
    

  !-----------------------------------------------------------------------------
  !> @brief
  !> determine if the associated `eval_obs` call would be valid
  !> manual (+global) observables always return true
  !> the check if manual observables have been set is done in getManual_obs
  !
  !> @param[in] obsId unique identifier for the observabes
  !> @param[in] npar kinematics specified by the number of particles
  !-----------------------------------------------------------------------------
  pure logical function isValid_obs(obsId, npar)
    use DIS_mod
    integer, intent(in) :: obsId, npar
    isValid_obs = .true.
    
    !> Dont cut on extra jet cuts
    if(isExtraj_obs(obsId)) isValid_obs = .false.
    
    select case (obsId)
      ! 1st jet (leading)
      case ( id_pt_j1, id_y_j1, id_abs_y_j1, id_eta_j1, id_abs_eta_j1,  &
           & id_dis_hera_etaj1,id_dis_hera_etj1,id_dis_hera_et2j1,  &
           & id_dis_hera_et2byQ2j1)
        if ( kin(npar)%njets < 1 ) isValid_obs = .false.
      case (id_pt_Vj)
        if ( kin(npar)%njets < 1 ) isValid_obs = .false.
      ! Added by KR
      case (id_ystar_Vj)
        if ( kin(npar)%njets < 1 ) isValid_obs = .false.
      case (id_yboost_Vj)
        if ( kin(npar)%njets < 1 ) isValid_obs = .false.
      ! End: Added by KR
      case (id_deltaphi_Vj1, id_deltay_Vj1, id_ydiff_Vj1, id_ysum_Vj1)
        if ( kin(npar)%njets < 1 ) isValid_obs = .false.
      ! 2nd jet (sub-leading)
      case ( id_pt_j2, id_y_j2, id_abs_y_j2, id_eta_j2, id_abs_eta_j2,  &
           & id_dis_hera_etaj2,id_dis_hera_etj2,id_dis_hera_et2j2,  &
           & id_dis_hera_et2byQ2j2)
        if ( kin(npar)%njets < 2 ) isValid_obs = .false.
      ! 3rd jet
      case ( id_pt_j3, id_y_j3, id_abs_y_j3, id_eta_j3, id_abs_eta_j3,  &
           & id_dis_hera_etaj3,id_dis_hera_etj3,id_dis_hera_et2j3,  &
           & id_dis_hera_et2byQ2j3,id_dis_hera_deta1,id_dis_hera_deta2)
        if ( kin(npar)%njets < 3 ) isValid_obs = .false.
      ! 4th jet
      case ( id_pt_j4, id_y_j4, id_abs_y_j4, id_eta_j4, id_abs_eta_j4,  &
           & id_dis_hera_etaj4,id_dis_hera_etj4,id_dis_hera_et2j4,  &
           & id_dis_hera_et2byQ2j4)
        if ( kin(npar)%njets < 4 ) isValid_obs = .false.
        !> exclusive 3 jet oriented epem event shapes
     case(id_epem_costh_en3j,id_epem_chi)
        if (kin(npar)%njets .ne. 3) isValid_obs = .false.
      ! di-jet observables
      case ( id_minv_j12, id_deltapt_j12, id_deltaeta_j12, id_deltay_j12,  &
           & id_deltaphi_j12, id_chi_j12, id_ystar_j12, id_yboost_j12, id_max_y_j12 )
        if ( kin(npar)%njets < 2 ) isValid_obs = .false.
      case (id_pt_Vjj)
        if ( kin(npar)%njets < 2 ) isValid_obs = .false.
      ! VFH 
      case (id_VFH_deltay, id_VFH_y1xy2)
        if ( kin(npar)%njets < 2 ) isValid_obs = .false.
      case (id_VFH_z3)
        if ( kin(npar)%njets < 3) isValid_obs = .false.
      ! DIS
      case ( id_ptavg_jall, id_ptavg_geom_jall, id_etavg_jall,  &
           & id_etavg_geom_jall,id_gammap_eta_j1)
        if ( kin(npar)%njets < 1 ) isValid_obs = .false.
      case ( id_DIS_etas,id_DIS_xi2, id_ptavg_j12, id_DIS_hera_deltaeta_j12,  &
           & id_DIS_Lgxi2,id_DIS_hera_etaavg_j12, id_etavg_j12, id_et2_j12,  &
           & id_gammap_eta_j2, id_DIS_delta_phis, id_DIS_hera_xgamma)
        if ( kin(npar)%njets < 2 ) isValid_obs = .false.
      case (id_DIS_xi3, id_ptavg_j123, id_etavg_j123)
        if ( kin(npar)%njets < 3 ) isValid_obs = .false.
      case (id_DIS_hera_eta_j1_0)
        if ( dis%rawjets<1) isValid_obs = .false.
      case (id_DIS_hera_eta_j2_0)
        if ( dis%rawjets<2) isValid_obs = .false.       
    end select
  end function isValid_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> determine if the associated observable is a global one
  !> global variables are the same for the whole event, i.e. also for the
  !> counter-events and should be set manually by the user
  !> and we therefore use the invalid flag to check if the varibale has been set
  !> already.
  !
  !> @param[in] obsId unique identifier for the observabes
  !-----------------------------------------------------------------------------
  pure logical function isGlobal_obs(obsId)
    integer, intent(in) :: obsId
    isGlobal_obs = .false.
    select case (obsId)
      ! DIS
      case ( id_DIS_sqrtQ2, id_DIS_Q2, id_DIS_x, id_DIS_y,  &
      &      id_DIS_W2, id_DIS_W, id_DIS_cosgammah, id_DIS_nu, id_DIS_resrho )
        isGlobal_obs = .true.
    end select
  end function isGlobal_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> determine if the associated observable is a manual one
  !> manual variables have to be set by the user but, in contrast to a 
  !> global observable, will change for each event. 
  !> isCached_obs is used to test for validity
  !
  !> @param[in] obsId unique identifier for the observabes
  !-----------------------------------------------------------------------------
  pure logical function isManual_obs(obsId)
    integer, intent(in) :: obsId
    isManual_obs = .false.
    !> all global observables are "manual": have to be set manually (at least once)
    if ( isGlobal_obs(obsId) ) then
      isManual_obs = .true.
      return
    end if
    !> transitions are no longer manual observables
    ! select case (obsId)
    !   !> JADE observables
    !   !> manual because set in the clustering algorithm
    !   case (id_y45, id_y34, id_y23, id_y12, id_y01)
    !     isManual_obs = .true.
    ! end select
  end function isManual_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> determine if the associated observable is an extra one
  !> extra variables are only used to facilitate the infrastructure of
  !> registering new observables / selectors... evaluation of these observables
  !> is left to the user
  !
  !> @param[in] obsId unique identifier for the observabes
  !-----------------------------------------------------------------------------
  pure logical function isExtraj_obs(obsId)
    integer, intent(in) :: obsId
    isExtraj_obs = .false.
    select case (obsId)
      ! DIS
      case(id_hera_abs_etaj, id_hera_etj, id_hera_ptj, id_hera_yj,id_hera_abs_yj, &
           & id_hera_etaj,id_hera_xjets,id_hera_H1xjets)
        isExtraj_obs = .true.
    end select
  end function isExtraj_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> set the value for a manual observable
  !
  !> @param[in] obsId unique identifier for the observabe
  !> @param[in] value value of the global observable
  !-----------------------------------------------------------------------------
  subroutine setManual_obs(obsId, value)
    integer, intent(in) :: obsId
    double precision, intent(in) :: value
    if ( .not. isManual_obs(obsId) ) then
      print*, "setManual_obs: tried to set a non-manual observable: ", obsId
      stop
    end if
    list_cache(obsId)%value   = value
    list_cache(obsId)%qcached = .true.
  end subroutine setManual_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> get the value for a manual observable
  !
  !> @param[in] obsId unique identifier for the observabe
  !> @return value of the observable
  !-----------------------------------------------------------------------------
  double precision function getManual_obs(obsId)
    integer, intent(in) :: obsId
    if ( .not. isManual_obs(obsId) ) then
      print*, "getManual_obs: tried to access a non-manual observable!"
      stop
    end if
    if (.not.list_cache(obsId)%qcached) then
      print*, "getManual_obs: manual observale not set:", obsId
      stop
    end if
    getManual_obs = list_cache(obsId)%value
  end function getManual_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> see whether an observables was cached already from outside
  !>
  !> @param[in] obsId unique identifier for the observable
  !-----------------------------------------------------------------------------
  pure logical function isCached_obs(obsID)
    integer, intent(in) :: obsId
    isCached_obs = list_cache(obsId)%qcached
    return
  end function isCached_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> determine if the observable is an integer type
  !
  !> @param[in] obsId unique identifier for the observabes
  !-----------------------------------------------------------------------------
  pure logical function isInteger_obs(obsId)
    integer, intent(in) :: obsId
    isInteger_obs = .false.
    select case (obsId)
      ! number of reconstructed jets
      case (id_njets, id_npartons)
        isInteger_obs = .true.
      ! index of leading / sub-leading gamma(leptons)
      case (id_idx_g1, id_idx_g2)
        isInteger_obs = .true.
    end select
  end function isInteger_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> reset the cache (do it before calling the eval_obs routines)
  !-----------------------------------------------------------------------------
  subroutine clearCache_obs()
    integer :: i
    do i=1,n_obs
      !> global observables should persist beyond different phase-space points
      !> or between real ME & subtraction terms
      !> manual observables are re-set for each kinematics
      if (.not.isGlobal_obs(i)) list_cache(i)%qcached = .false.
    end do
  end subroutine clearCache_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> Wrapper function for observables
  !
  !> @param[in] obsId unique identifier for the observabes
  !> @param[in] npar kinematics specified by the number of particles
  !-----------------------------------------------------------------------------
  double precision function eval_obs(obsId, npar)
    integer, intent(in) :: obsId, npar
    double precision :: eval_sudakov_minlo, eval_rewgt_minlo
    !> already computed this observable: pull from the cache
    if ( list_cache(obsId)%qcached ) then
      eval_obs = list_cache(obsId)%value
      return
    end if
    !> actually compute the observable
    select case (obsId)
      ! number of reconstructed jets (integer)
      case (id_njets)
        eval_obs = DBLE( kin(npar)%njets ) ! back-conversion to INT is exact
      case (id_npartons)
        eval_obs = DBLE( kin(npar)%nproto )  ! # of FS partons
      ! 1st lepton
      case (id_pt_l1)
        eval_obs = eval_pt_l1(npar)
      case (id_y_l1)
        eval_obs = eval_y_l1(npar)
      case (id_abs_y_l1)
        eval_obs = abs(eval_y_l1(npar))
      ! 2nd lepton
      case (id_pt_l2)
        eval_obs = eval_pt_l2(npar)
      case (id_y_l2)
        eval_obs = eval_y_l2(npar)
      case (id_abs_y_l2)
        eval_obs = abs(eval_y_l2(npar))
      ! 3st lepton
      case (id_pt_l3)
        eval_obs = eval_pt_l3(npar)
      case (id_y_l3)
        eval_obs = eval_y_l3(npar)
      case (id_abs_y_l3)
        eval_obs = abs(eval_y_l3(npar))
      ! 4nd lepton
      case (id_pt_l4)
        eval_obs = eval_pt_l4(npar)
      case (id_y_l4)
        eval_obs = eval_y_l4(npar)
      case (id_abs_y_l4)
        eval_obs = abs(eval_y_l4(npar))
      ! leading gamma / lepton in Z
      case (id_idx_g1)
        call set_idx_g12(npar)
        eval_obs = list_cache(id_idx_g1)%value
        return
      case (id_pt_g1)
        eval_obs = eval_pt_g1(npar)
      case (id_y_g1)
        eval_obs = eval_y_g1(npar)
      case (id_abs_y_g1)
        eval_obs = abs(eval_y_g1(npar))
      ! sub-leading gamma / lepton in Z
      case (id_idx_g2)
        call set_idx_g12(npar)
        eval_obs = list_cache(id_idx_g2)%value
        return
      case (id_pt_g2)
        eval_obs = eval_pt_g2(npar)
      case (id_y_g2)
        eval_obs = eval_y_g2(npar)
      case (id_abs_y_g2)
        eval_obs = abs(eval_y_g2(npar))
      ! gamma-gamma/lepton-lepton separation
      case (id_deltaphi_g1g2)
        eval_obs = eval_deltaphi_g1g2(npar)
      case (id_deltay_g1g2)
        eval_obs = eval_deltay_g1g2(npar)
      case (id_abs_deltay_g1g2)
        eval_obs = abs(eval_deltay_g1g2(npar))
      case (id_ptt_g1g2)
        eval_obs = eval_ptt_g1g2(npar)
      ! boson - leading jet
      case (id_deltaphi_Vj1)
        eval_obs = eval_deltaphi_Vj1(npar)
      case (id_deltay_Vj1)
        eval_obs = eval_deltay_Vj1(npar)
      case (id_ydiff_Vj1)
        eval_obs = eval_ydiff_Vj1(npar)
      case (id_ysum_Vj1)
        eval_obs = eval_ysum_Vj1(npar)
      !>----- H -> ZZ --> 4 l decay
      ! leading lepton
      case (id_pt_l1st)
        eval_obs = eval_pt_4l(npar,1)
      case (id_y_l1st)
        eval_obs = eval_y_4l(npar,1)
      case (id_abs_y_l1st)
        eval_obs = abs(eval_y_4l(npar,1))
      case (id_R_l1st_j1)
        eval_obs = eval_dR_lj1(npar,1)
      ! sub-leading lepton
      case (id_pt_l2nd)
        eval_obs = eval_pt_4l(npar,2)
      case (id_y_l2nd)
        eval_obs = eval_y_4l(npar,2)
      case (id_abs_y_l2nd)
        eval_obs = abs(eval_y_4l(npar,2))
      case (id_R_l2nd_j1)
        eval_obs = eval_dR_lj1(npar,2)
      ! third-leading lepton
      case (id_pt_l3rd)
        eval_obs = eval_pt_4l(npar,3)
      case (id_y_l3rd)
        eval_obs = eval_y_4l(npar,3)
      case (id_abs_y_l3rd)
        eval_obs = abs(eval_y_4l(npar,3))
      case (id_R_l3rd_j1)
        eval_obs = eval_dR_lj1(npar,3)
      ! least-leading lepton
      case (id_pt_l4th)
        eval_obs = eval_pt_4l(npar,4)
      case (id_y_l4th)
        eval_obs = eval_y_4l(npar,4)
      case (id_abs_y_l4th)
        eval_obs = abs(eval_y_4l(npar,4))
      case (id_R_l4th_j1)
        eval_obs = eval_dR_lj1(npar,4)
      ! gauge boson = l1+l2
      case (id_pt_V)
        eval_obs = eval_pt_V(npar)
      case (id_minv_V)
        eval_obs = eval_minv_V(npar)
      case (id_y_V)
        eval_obs = eval_y_V(npar)
      case (id_abs_y_V)
        eval_obs = abs(eval_y_V(npar))
      case (id_eta_V)
        eval_obs = eval_eta_V(npar)
      case (id_abs_eta_V)
        eval_obs = abs(eval_eta_V(npar))
      case (id_mt_V)
        eval_obs = eval_mt_V(npar)
      case (id_et_V)
        eval_obs = eval_et_V(npar)
      ! gauge boson = l3+l4
      case (id_pt_V2)
        eval_obs = eval_pt_V2(npar)
      case (id_minv_V2)
        eval_obs = eval_minv_V2(npar)
      case (id_y_V2)
        eval_obs = eval_y_V2(npar)
      case (id_abs_y_V2)
        eval_obs = abs(eval_y_V2(npar))
      case (id_eta_V2)
        eval_obs = eval_eta_V2(npar)
      case (id_abs_eta_V2)
        eval_obs = abs(eval_eta_V2(npar))
      case (id_mt_V2)
        eval_obs = eval_mt_V2(npar)
      case (id_et_V2)
        eval_obs = eval_et_V2(npar)
      ! gauge boson = l1+l2+l3+l4
      case (id_pt_V4l)
        eval_obs = eval_pt_V4l(npar)
      case (id_minv_V4l)
        eval_obs = eval_minv_V4l(npar)
      case (id_minv_V23)
        eval_obs = eval_minv_V23(npar)
      case (id_minv_V14)
        eval_obs = eval_minv_V14(npar)
      case (id_y_V4l)
        eval_obs = eval_y_V4l(npar)
      case (id_abs_y_V4l)
        eval_obs = abs(eval_y_V4l(npar))
      case (id_eta_V4l)
        eval_obs = eval_eta_V4l(npar)
      case (id_abs_eta_V4l)
        eval_obs = abs(eval_eta_V4l(npar))
      case (id_minv_mZ1)
        eval_obs = eval_minv_mZ1(npar)
      case (id_minvSFOS_vs_mZ)
        eval_obs = eval_minvSFOS_vs_mZ(npar)
      case (id_minv_mZ2)
        eval_obs = eval_minv_mZ2(npar)
      case (id_minv_mZSFOS1st)
        eval_obs = eval_minv_mZSFOS(npar,1)
      case (id_minv_mZSFOS2nd)
        eval_obs = eval_minv_mZSFOS(npar,2)
      case (id_minv_mZSFOS3rd)
        eval_obs = eval_minv_mZSFOS(npar,3)
      case (id_minv_mZSFOS4th)
        eval_obs = eval_minv_mZSFOS(npar,4)
      ! thrust dependend variables
      case ( id_epem_T, id_epem_1mT, id_epem_cosT,  &
           & id_epem_WJB, id_epem_TJB, id_epem_HJM, id_epem_SHM)
         call set_Taxis(npar)
         eval_obs = list_cache(obsId)%value
         return 
      !> Theta depende
      case (id_epem_costh_ej1)
         eval_obs = eval_epem_costh_ej1(npar)
      case (id_epem_costh_en3j)
         eval_obs = eval_epem_costh_en3j(npar)
      case (id_epem_chi)
         eval_obs = eval_epem_chi(npar)
      case (id_epem_C)
         eval_obs = eval_epem_C(npar) 
      case (id_epem_D)
         eval_obs = eval_epem_D(npar)  
             
      ! JADE ycut transitions
      case (id_y45, id_y34, id_y23, id_y12, id_y01)
        call set_jade_ytrans(npar)
        eval_obs = list_cache(obsId)%value
        return

      ! Collins-Soper frame (angles wrt l1 = lm)
      case(id_cos_theta_CS, id_theta_CS, id_phi_CS)
        call set_Collins_Soper(npar)
        eval_obs = list_cache(obsId)%value
        return

     case(id_cos_theta_CS_0PT)
        eval_obs = eval_cos_theta_CS_0PT(npar)
        return

      ! projectors for the angular coefficients A_i [1606.00689]
      case ( id_proj_A0, id_proj_A1, id_proj_A2,  &
           & id_proj_A3, id_proj_A4, id_proj_A5,  &
           & id_proj_A6, id_proj_A7, id_proj_A0mA2, id_proj_A0mA2alt)
        call set_proj_Ai(npar)
        eval_obs = list_cache(obsId)%value
        return

      ! VJ... observables
      case (id_pt_Vjj)
        eval_obs = eval_pt_Vjj(npar)
      case (id_pt_Vj)
        eval_obs = eval_pt_Vj(npar)
      ! Added by KR
      case (id_ystar_Vj)
        eval_obs = eval_ystar_Vj(npar)
      case (id_yboost_Vj)
        eval_obs = eval_yboost_Vj(npar)
      ! End: Added by KR
      ! 1st jet (leading)
      case (id_pt_j1)
        eval_obs = eval_pt_j1(npar)
      case (id_y_j1)
        eval_obs = eval_y_j1(npar)
      case (id_abs_y_j1)
        eval_obs = abs(eval_y_j1(npar))
      case (id_eta_j1)
        eval_obs = eval_eta_j1(npar)
      case (id_abs_eta_j1)
        eval_obs = abs(eval_eta_j1(npar))
      ! 2nd jet (sub-leading)
      case (id_pt_j2)
        eval_obs = eval_pt_j2(npar)
      case (id_y_j2)
        eval_obs = eval_y_j2(npar)
      case (id_abs_y_j2)
        eval_obs = abs(eval_y_j2(npar))
      case (id_eta_j2)
        eval_obs = eval_eta_j2(npar)
      case (id_abs_eta_j2)
        eval_obs = abs(eval_eta_j2(npar))
      ! 3rd jet
      case (id_pt_j3)
        eval_obs = eval_pt_j3(npar)
      case (id_y_j3)
        eval_obs = eval_y_j3(npar)
      case (id_abs_y_j3)
        eval_obs = abs(eval_y_j3(npar))
      case (id_eta_j3)
        eval_obs = eval_eta_j3(npar)
      case (id_abs_eta_j3)
        eval_obs = abs(eval_eta_j3(npar))
      ! 4th jet
      case (id_pt_j4)
        eval_obs = eval_pt_j4(npar)
      case (id_y_j4)
        eval_obs = eval_y_j4(npar)
      case (id_abs_y_j4)
        eval_obs = abs(eval_y_j4(npar))
      case (id_eta_j4)
        eval_obs = eval_eta_j4(npar)
      case (id_abs_eta_j4)
        eval_obs = abs(eval_eta_j4(npar))
      case (id_sum_pt_jets)
        eval_obs = eval_sum_pt_jets(npar)
      ! di-jet invariant mass (j1+j2)
      case (id_minv_j12)
        eval_obs = eval_minv_j12(npar)
      ! more di-jet observables
      case(id_deltapt_j12)
        eval_obs = eval_deltapt_j12(npar)
      case(id_deltaeta_j12)
        eval_obs = eval_deltaeta_j12(npar)
      case(id_deltay_j12)
        eval_obs = eval_deltay_j12(npar)
      case(id_deltaphi_j12)
        eval_obs = eval_deltaphi_j12(npar)
      case(id_chi_j12)
        eval_obs = eval_chi_j12(npar)
      case(id_ystar_j12)
        eval_obs = eval_ystar_j12(npar)
      case(id_yboost_j12)
        eval_obs = eval_yboost_j12(npar)
      case(id_max_y_j12)
        eval_obs = eval_max_y_j12(npar)
      case (id_pt_j1expYstar)
        eval_obs = eval_pt_j1expYstar(npar)

      ! tau-jet observables (used in H+J)
      case (id_tau_j1)
        eval_obs = eval_tau_j1(npar)
      case (id_sum_tau_jets)
        eval_obs = eval_sum_tau_jets(npar)
      case (id_max_tau_jet)
        eval_obs = eval_max_tau_jet(npar)

      ! tri-jet invariant mass (j1+j2+j3)
      case (id_minv_j123)
        eval_obs = eval_minv_j123(npar)
      ! avergae of pt's
      case (id_ptavg_j12)
        eval_obs = eval_ptavg_j12(npar)
      case (id_ptavg_j123)
        eval_obs = eval_ptavg_j123(npar)
      case (id_ptavg_jall)
        eval_obs = eval_ptavg_jall(npar)
      case (id_ptavg_geom_jall)
        eval_obs = eval_ptavg_geom_jall(npar)
      ! avergae of et's
      case (id_etavg_j12)
        eval_obs = eval_etavg_j12(npar)
      case (id_etavg_j123)
        eval_obs = eval_etavg_j123(npar)
      case (id_etavg_jall)
        eval_obs = eval_etavg_jall(npar)
      case (id_etavg_geom_jall)
        eval_obs = eval_etavg_geom_jall(npar)
      case (id_et2_j12)
        eval_obs = eval_et2_j12(npar)

      !> "perfect" jet observables
      case ( id_pt_j1_nocut, id_pt_j2_nocut, id_ht_jall_nocut,  &
           & id_ptavg_j12_nocut, id_ptavg_jall_nocut,  &
           & id_ptavg_geom_j12_nocut, id_ptavg_geom_jall_nocut)
        call set_jets_nocut(npar)
        eval_obs = list_cache(obsId)%value
        return

      !> jet observables using R=1
      case ( id_pt_j1_R1, id_pt_j2_R1, id_ht_jall_R1,  &
           & id_ptavg_j12_R1, id_ptavg_jall_R1,  &
           & id_ptavg_geom_j12_R1, id_ptavg_geom_jall_R1)
        call set_jets_R1(npar)
        eval_obs = list_cache(obsId)%value
        return

      ! photon observables
      case (id_photon_pt)
        eval_obs = eval_photon_pt(npar)
      case (id_photon_y)
        eval_obs = eval_photon_y(npar)
      ! DIS observables
      case (id_DIS_Lgxi2)
         eval_obs = eval_dis_Lgxi2(npar)
      case (id_DIS_visbyW)
         eval_obs = eval_dis_visbyW(npar)
      case (id_DIS_etas)
         eval_obs = eval_dis_etas(npar)
      case (id_DIS_xi2)
         eval_obs = eval_dis_xi2(npar)
      case (id_DIS_Thrust)
         eval_obs = eval_dis_thrust(npar)
      case (id_DIS_Thrust_c)
         eval_obs = eval_dis_thrust_c(npar)
      case (id_DIS_JB)
         eval_obs = eval_dis_JB(npar)
      case (id_DIS_JM2)
         eval_obs = eval_dis_JM2(npar)
      case (id_DIS_C)
         eval_obs = eval_dis_C(npar)
      case (id_dis_hera_xgamma)
         eval_obs = eval_dis_hera_xgamma(npar)
      case (id_gammap_eta_j1)
         eval_obs = eval_gammap_eta_j1()
      case (id_gammap_eta_j2)
         eval_obs = eval_gammap_eta_j2()
      case (id_DIS_delta_phis)
         eval_obs = eval_dis_delta_phis()
      case (id_DIS_hera_etaj1)
         eval_obs=eval_dis_hera_etaj1(npar)
      case (id_DIS_hera_etaj2)
         eval_obs=eval_dis_hera_etaj2(npar)
      case (id_DIS_hera_etaj3)
         eval_obs=eval_dis_hera_etaj3(npar)
      case (id_DIS_hera_etaj4)
         eval_obs=eval_dis_hera_etaj4(npar)
      case (id_DIS_hera_etj1)
         eval_obs=eval_dis_hera_etj1(npar)
      case (id_DIS_hera_etj2)
         eval_obs=eval_dis_hera_etj2(npar)
      case (id_DIS_hera_etj3)
         eval_obs=eval_dis_hera_etj3(npar)
      case (id_DIS_hera_etj4)
         eval_obs=eval_dis_hera_etj4(npar)
      case (id_DIS_hera_et2j1)
         eval_obs=eval_dis_hera_et2j1(npar)
      case (id_DIS_hera_et2j2)
         eval_obs=eval_dis_hera_et2j2(npar)
      case (id_DIS_hera_et2j3)
         eval_obs=eval_dis_hera_et2j3(npar)
      case (id_DIS_hera_et2j4)
         eval_obs=eval_dis_hera_et2j4(npar)
      case (id_DIS_hera_et2byq2j1)
         eval_obs=eval_dis_hera_et2byq2j1(npar)
      case (id_DIS_hera_et2byq2j2)
         eval_obs=eval_dis_hera_et2byq2j2(npar)
      case (id_DIS_hera_et2byq2j3)
         eval_obs=eval_dis_hera_et2byq2j3(npar)
      case (id_DIS_hera_et2byq2j4)
         eval_obs=eval_dis_hera_et2byq2j4(npar)
      case (id_DIS_hera_eta_j1_0)
         eval_obs = eval_dis_hera_eta_j1_0()
      case (id_DIS_hera_eta_j2_0)
         eval_obs = eval_dis_hera_eta_j2_0()
      case (id_DIS_hera_deta1)
         eval_obs = eval_dis_hera_deta1(npar)
      case (id_DIS_hera_deta2)
         eval_obs = eval_dis_hera_deta2(npar)
      case(id_DIS_hera_FJ_reweight)
         eval_obs = eval_DIS_hera_FJ_reweight(npar)
      case (id_DIS_hera_deltaeta_j12)
         eval_obs = eval_dis_hera_deltaeta_j12()
      case (id_DIS_hera_etaavg_j12)
         eval_obs = eval_dis_hera_etaavg_j12()
      case (id_DIS_xi3)
         eval_obs = eval_dis_xi3(npar)
      case (id_DIS_dscl)
         eval_obs = eval_dis_dscl(npar)
      case (id_DIS_dscl3j)
         eval_obs = eval_dis_dscl3j(npar)
      case (id_DIS_dsclj1)
         eval_obs = eval_dis_dsclj1(npar)
      case (id_DIS_dsclj2)
         eval_obs = eval_dis_dsclj2(npar)
      case (id_DIS_dsclj3)
         eval_obs = eval_dis_dsclj3(npar)
      case (id_DIS_dsclj4)
         eval_obs = eval_dis_dsclj4(npar)
      case (id_DIS_dsclZEUS)
         eval_obs = eval_dis_dsclZEUS(npar)
      case (id_DIS_dsclj1l)
         eval_obs = eval_dis_dsclj1l(npar)
      case (id_DIS_dsclj2l)
         eval_obs = eval_dis_dsclj2l(npar)
      case (id_DIS_dsclj3l)
         eval_obs = eval_dis_dsclj3l(npar)
      case (id_DIS_dsclj4l)
         eval_obs = eval_dis_dsclj4l(npar)
      case (id_DIS_dsclZEUS2)
         eval_obs = eval_dis_dsclZEUS2(npar)
      case (id_DIS_ptl)
         eval_obs = eval_dis_ptl(npar)
      ! more
      case (id_ht_jets)
        eval_obs = eval_ht_jets(npar)
      case (id_ht_part)
        eval_obs = eval_ht_part(npar)
      case (id_ht_full)
        eval_obs = eval_ht_full(npar)
      case (id_min_dR_l1j)
        eval_obs = eval_min_dR_l1j(npar)
      case (id_min_dR_l2j)
        eval_obs = eval_min_dR_l2j(npar)
      case (id_min_dR_l3j)
        eval_obs = eval_min_dR_l3j(npar)
      case (id_min_dR_l4j)
        eval_obs = eval_min_dR_l4j(npar)
      case (id_min_dR_l1234j)
        eval_obs = eval_min_dR_l1234j(npar)
      case (id_min_dR_l12j)
        eval_obs = eval_min_dR_l12j(npar)
      case (id_min_dR_l34j)
        eval_obs = eval_min_dR_l34j(npar)
      case (id_dR_l12)
        eval_obs = eval_dR_l12(npar)
      case (id_dR_l34)
        eval_obs = eval_dR_l34(npar)
      case (id_min_dR_l_sf)
        eval_obs = eval_min_dR_l_sf(npar)
      case (id_min_dR_l_df)
        eval_obs = eval_min_dR_l_df(npar)
      case (id_phi_star_eta)
        eval_obs = eval_phi_star_eta(npar)
      case (id_ZJ_dscale1)
        eval_obs = eval_ZJ_dscale1(npar)
      case (id_ZJ_dscale3)
        eval_obs = eval_ZJ_dscale3(npar)
      case (id_ZJ_HT)
        eval_obs = eval_ZJ_HT(npar)
      case (id_ZJ_HTshape05)
        eval_obs = eval_ZJ_HT_shape(npar, 0.5d0)
      case (id_ZJ_HTshape20)
        eval_obs = eval_ZJ_HT_shape(npar, 2d0)
      case (id_H4lJ_dscale1)
        eval_obs = eval_H4lJ_dscale1(npar)
      case (id_H4lJ_dscale2)
        eval_obs = eval_H4lJ_dscale2(npar)
      ! photon decay angles in the Collins-Soper reference frame
      case (id_costhetastar)
        eval_obs = eval_costhetastar(npar)
      ! VFH
      case (id_VFH_deltay)
        eval_obs = eval_VFH_deltay(npar)
      case (id_VFH_y1xy2)
        eval_obs = eval_VFH_y1xy2(npar)
      case (id_VFH_dscale)
        eval_obs = eval_VFH_dscale(npar)
      case (id_VFH_z3)
        eval_obs = eval_VFH_z3(npar)
      ! MiNLO observables
      case (id_sudakov)
         eval_obs = eval_sudakov_minlo(npar)
      case (id_minlo)
         eval_obs = eval_rewgt_minlo(npar)
      ! Spikes
      case (id_s13)
        eval_obs = eval_sij(npar, 1, 3)
      case (id_s14)
        eval_obs = eval_sij(npar, 1, 4)
      case (id_s15)
        eval_obs = eval_sij(npar, 1, 5)
      case (id_s16)
        eval_obs = eval_sij(npar, 1, 6)
      case (id_s23)
        eval_obs = eval_sij(npar, 2, 3)
      case (id_s24)
        eval_obs = eval_sij(npar, 2, 4)
      case (id_s25)
        eval_obs = eval_sij(npar, 2, 5)
      case (id_s26)
        eval_obs = eval_sij(npar, 2, 6)
      case (id_s34)
        eval_obs = eval_sij(npar, 3, 4)
      case (id_s35)
        eval_obs = eval_sij(npar, 3, 5)
      case (id_s36)
        eval_obs = eval_sij(npar, 3, 6)
      case (id_s45)
        eval_obs = eval_sij(npar, 4, 5)
      case (id_s46)
        eval_obs = eval_sij(npar, 4, 6)
      case (id_s56)
        eval_obs = eval_sij(npar, 5, 6)
      ! escape
      case default
        eval_obs = 0d0
        print*, "eval_obs: unknown obsId or unset manual obs: ", obsId
        stop
    end select
    !> save into cache
    list_cache(obsId)%value   = eval_obs
    list_cache(obsId)%qcached = .true.
  end function eval_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> Wrapper function for jet observables.
  !> These are evaluated during the recombination algorithm in order to
  !> evaluate jet-veto conditions and cannot be used for histograms!
  !
  !> @param[in] obsIdj unique identifier for the jet observabes
  !> @param[in] npar kinematics specified by the number of particles
  !> @param[in] v four-vector of the jet
  !-----------------------------------------------------------------------------
  pure double precision function eval_obsj(obsIdj, npar, v)
    integer, intent(in) :: obsIdj, npar
    double precision, intent(in), dimension(4) :: v
    select case (obsIdj)
      !> jet properties
      case (idj_pt)
        eval_obsj = v1_pt(v)
      case (idj_y)
        eval_obsj = v1_y(v)
      case (idj_abs_y)
        eval_obsj = abs(v1_y(v))
      case (idj_eta)
        eval_obsj = v1_eta(v)
      case (idj_abs_eta)
        eval_obsj = abs(v1_eta(v))
      case (idj_et)
        eval_obsj = v1_et(v)
      !> properties dependent on npar kinematics
      case(idj_dR_l1)
        eval_obsj = v2_delta_R(kin(npar)%p(:,npar-1), v)
      case(idj_dR_l2)
        eval_obsj = v2_delta_R(kin(npar)%p(:,npar),   v)
      case(idj_dR_l3)
        eval_obsj = v2_delta_R(kin(npar)%p(:,npar-3), v)
      case(idj_dR_l4)
        eval_obsj = v2_delta_R(kin(npar)%p(:,npar-2), v)
      case(idj_min_dR_l12)
        eval_obsj = min( v2_delta_R(kin(npar)%p(:,npar-1), v),  &
                       & v2_delta_R(kin(npar)%p(:,npar),   v) )
      case(idj_min_dR_l34)
        eval_obsj = min( v2_delta_R(kin(npar)%p(:,npar-3), v),  &
                       & v2_delta_R(kin(npar)%p(:,npar-2), v) )
      case(idj_min_dR_l1234)
        eval_obsj = min( v2_delta_R(kin(npar)%p(:,npar-1), v),  &
                       & v2_delta_R(kin(npar)%p(:,npar),   v),  &
                       & v2_delta_R(kin(npar)%p(:,npar-3), v),  &
                       & v2_delta_R(kin(npar)%p(:,npar-2), v) )
      !> escape
      case default
        eval_obsj = 0d0
        ! stop "eval_obsj: unknown obsIdj"
      end select
  end function eval_obsj


  !------------------------------------------------------------------------------
  !> @brief
  !> Setter routine to determine whether Et weighted recombination scheme should be
  !> used or not
  !> @param[in] true or false
  !-------------------------------------------------------------------------------
  subroutine setEtrecom(val)
    logical, intent(in) :: val
    Et_recom=val
  end subroutine setEtrecom


  !-----------------------------------------------------------------------------
  !> @brief
  !> set the algorithm and cone size
  !
  !> @param[in] algorithm string identifying the algorithm
  !> @param[in] param parameter of the jet algorithm
  !>              * ycut for JADE
  !>              * R else
  !-----------------------------------------------------------------------------
  subroutine init_jet(algorithm, param)
    use StrHelper_mod
    use Process_mod
    character(len=*), intent(in) :: algorithm
    double precision, intent(in) :: param
    !> global setting for the algorithm
    select case (to_lower_str(algorithm))
      case ("none")
        jetalg = 0
        ip_jet = 0
      case ("antikt")
        jetalg = 1
        ip_jet = -1
      case ("cam")
        jetalg = 2
        ip_jet = 0
      case ("kt")
        jetalg = 3
        ip_jet = +1
      case ("jade")
        jetalg = 4
        ip_jet = 0
        jade_ycut = param
      case default
        print*, "unknown jet algorithm: ", algorithm
        stop
    end select
    inputRcone = param
    call setRcone_jet()
    !> register observables
    list_obsj(idj_pt)      = JetSelector_t(.true., "jets_pt")
    list_obsj(idj_y)       = JetSelector_t(.true., "jets_y")
    list_obsj(idj_abs_y)   = JetSelector_t(.true., "jets_abs_y")
    list_obsj(idj_eta)     = JetSelector_t(.true., "jets_eta")
    list_obsj(idj_abs_eta) = JetSelector_t(.true., "jets_abs_eta")
    list_obsj(idj_et)      = JetSelector_t(.true., "jets_et")

    !> process-specific jet-observables (default init is qactive=.false.)
    select case (name_proc)

      case ("Z","ZJ","ZJJ","ZJJJ")
        list_obsj(idj_min_dR_l12)   = JetSelector_t(.true., "jets_min_dr_lj")

      case ("WM","WMJ","WMJJ","WMJJJ")
        list_obsj(idj_dR_l1)        = JetSelector_t(.true., "jets_min_dr_lj")

      case ("WP","WPJ","WPJJ","WPJJJ")
        list_obsj(idj_dR_l2)        = JetSelector_t(.true., "jets_min_dr_lj")

      case ("HTO2PM","HTO2PJM","HTO2P","HTO2PJ","HTO2PJJ","HTO2PJJJ")
        list_obsj(idj_min_dR_l12)   = JetSelector_t(.true., "jets_min_dr_gj")

      case ("HTO2TAUM","HTO2TAUJM","HTO2TAU","HTO2TAUJ","HTO2TAUJJ","HTO2TAUJJJ")
        list_obsj(idj_min_dR_l12)   = JetSelector_t(.true., "jets_min_dr_tauj")

      case ("HTO4EM","HTO4EJM","HTO4E","HTO4EJ","HTO4EJJ","HTO4EJJJ")
        list_obsj(idj_min_dR_l1234) = JetSelector_t(.true., "jets_min_dr_lj")

      case ("HTO2E2MUM","HTO2E2MUJM","HTO2E2MU","HTO2E2MUJ","HTO2E2MUJJ","HTO2E2MUJJJ")
        list_obsj(idj_min_dR_l12)   = JetSelector_t(.true., "jets_min_dr_muj")
        list_obsj(idj_min_dR_l34)   = JetSelector_t(.true., "jets_min_dr_ej")
        list_obsj(idj_min_dR_l1234) = JetSelector_t(.true., "jets_min_dr_lj")

        !> @todo For H --> Wp Wm: only two charged leptons
        !> ATTENTION: depending on ordering (Wp,Wm) vs. (Wm,Wp),
        !> the charged leptons sit in (l2,l3) or (l1,l4)

    end select

  end subroutine init_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> get the algorithm id
  !-----------------------------------------------------------------------------
  integer function algo_jet()
    algo_jet = jetalg
  end function algo_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> set the cone size
  !
  !> @param[in] R optional argument for the R-parameter in the jet algorithm
  !-----------------------------------------------------------------------------
  subroutine setRcone_jet(R)
    double precision, optional, intent(in) :: R
    double precision :: Rcone_val
    if (present(R)) then
      !> used in multi-runs to set different thread-local values
      Rcone_val = R
    else
      !> default to using the input Rcone value
      Rcone_val = inputRcone
    end if
    Rcone   = Rcone_val
    RconeSq = Rcone_val**2
  end subroutine setRcone_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> set the selector for obsId: lower < obs < upper
  !> @note omit both `lower` and `upper` to disable the selector
  !
  !> @param[in] obsId unique identifyer of the jet observable
  !> @param[in] lower optional lower cut
  !> @param[in] upper optional upper cut
  !> @param[in] invert optional flag controlling whether we "accept" or "reject"
  !-----------------------------------------------------------------------------
  subroutine setSelector_jet(obsId, lower, upper, invert)
    integer, intent(in) :: obsId
    double precision, optional, intent(in) :: lower, upper
    logical, optional, intent(in) :: invert

    if (.not.list_obsj(-obsId)%qactive) then
      print*, "setSelector_jet: inactive jet observable:", obsId
      stop
    end if

    if (present(lower)) list_obsj(-obsId)%lower = lower
    if (present(upper)) list_obsj(-obsId)%upper = upper
    list_obsj(-obsId)%qlower = present(lower)
    list_obsj(-obsId)%qupper = present(upper)
    if (present(invert)) then
      list_obsj(-obsId)%qinvert = invert
    end if ! default init is .false. (see above)

  end subroutine setSelector_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> digest information on the registered selectors
  !
  !> @param[in] unit optional logical unit connected for formatted output (file, screen)
  !-----------------------------------------------------------------------------
  subroutine digest_jet(unit)
    integer, optional, intent(in) :: unit
    integer :: i, iout
    character(len=*), parameter :: fmt_algo    = '(1x,"*",3x,"jet algorithm:",A10,4x,"R = ",1f4.2,61x,"*")'
    character(len=*), parameter :: fmt_jade    = '(1x,"*",3x,"jet algorithm:",A10,4x,"ycut = ",1f4.2,61x,"*")'
    character(len=*), parameter :: fmt_min_max = '(1x,"*",3x,A6,":",1f12.2," <= ",A15," <= ",1f12.2,43x,"*")'
    character(len=*), parameter :: fmt_min     = '(1x,"*",3x,A6,":",1f12.2," <= ",A15,59x,"*")'
    character(len=*), parameter :: fmt_max     = '(1x,"*",3x,A6,":",16x,A15," <= ",1f12.2,43x,"*")'
    character(len=6) :: seltype

    if (present(unit)) then
      iout = unit
    else
      iout = 6 ! output to screen
    endif

    select case (jetalg)
    case (0)
      write (iout,fmt_algo) "none", Rcone
    case (1)
      write (iout,fmt_algo) "anti-kt", Rcone
    case (2)
      write(iout,fmt_algo) "C/A", Rcone
    case (3)
      write(iout,fmt_algo) "kT", Rcone
    case (4)
      write(iout,fmt_jade) "JADE", jade_ycut
    end select

    do i = 1,n_obsj
        if (.not.list_obsj(i)%qactive) cycle
        seltype = ''
        if ( list_obsj(i)%qinvert ) then
          seltype = 'reject'
        else
          seltype = 'accept'
        endif
        if (list_obsj(i)%qlower .and. list_obsj(i)%qupper) then
          write(iout,fmt_min_max) seltype, list_obsj(i)%lower, list_obsj(i)%name, list_obsj(i)%upper
        else if (list_obsj(i)%qlower .and. .not.list_obsj(i)%qupper) then
          write(iout,fmt_min) seltype, list_obsj(i)%lower, list_obsj(i)%name
        else if (.not.list_obsj(i)%qlower .and. list_obsj(i)%qupper) then
          write(iout,fmt_max) seltype, list_obsj(i)%name, list_obsj(i)%upper
        endif
    end do

  end subroutine digest_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> determine from the jet selectors the minimal partonic 
  !> centre-of-mass energy (squared)
  !
  !> @return minimal value for shat [GeV^2]
  !-----------------------------------------------------------------------------
  double precision function min_shat_jet()
    double precision :: min_valpt, min_valet
    min_shat_jet = 0d0
    !> only pt and et of the jets are dimensionful variables I have
    !> note: if minnjets = 0 this does not impose a shat_min cut!!!
    min_valpt = 0d0
    if ( .not.list_obsj(idj_pt)%qinvert .and. list_obsj(idj_pt)%qlower ) then
      min_valpt = list_obsj(idj_pt)%lower
    end if
    min_valet = 0d0
    if ( .not.list_obsj(idj_et)%qinvert .and. list_obsj(idj_et)%qlower ) then
      min_valet = list_obsj(idj_et)%lower
    end if
    !>---------------------------------------------------
    !> !!! this is WRONG !!!
    !> never uncomment! Left in so noone get's the idea of implementing such a line
    !> see comments in min_shat_sel
    ! if (      list_obsj(idj_pt)%qinvert .and. list_obsj(idj_pt)%qupper ) then
    !   min_valpt = list_obsj(idj_pt)%upper
    ! end if 
    !>---------------------------------------------------
    min_shat_jet = (2d0*max(min_valpt,min_valet))**2
    return
  end function min_shat_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> routine to initialize the (threadprivate) buffer for the jet-algorithm
  !
  !> @param[in] maxnproto largest array size needed for the jet-algorithm
  !-----------------------------------------------------------------------------
  subroutine initBuffer_jet(maxnproto)
    integer, intent(in) :: maxnproto
    !> DIS using JADE needs another pseudo-particle
    integer :: nonpis
    common /nonpisflag/ nonpis
    integer :: size

    !> no jet algo
    if (jetalg == 0) return

    size = maxnproto
    if (nonpis == 1 .and. jetalg == 4) then
      size = maxnproto + 1
    end if

    allocate( jets(size) )
    allocate( diB(size) )
    allocate( dij(size,size) )
    allocate( clusterHist_jet(size) )
    allocate( pt_sort_jet(size) )
  end subroutine initBuffer_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> routine to deallocate the (threadprivate) buffer for the jet-algorithm
  !-----------------------------------------------------------------------------
  subroutine destroyBuffer_jet()

    !> no jet algo
    if (jetalg == 0) return

    deallocate(jets)
    deallocate(diB)
    deallocate(dij)
    deallocate(clusterHist_jet)
    deallocate(pt_sort_jet)
  end subroutine destroyBuffer_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> Initialise & prepare the local arrays for the jet recombination algorithm
  !
  !> @param[in] npar kinematics specified by the number of particles
  !-----------------------------------------------------------------------------
  subroutine initRecomb_jet(npar)
    integer, intent(in) :: npar
    type(Kin_t), pointer :: kd
    integer :: i

    !> prettify code
    kd => kin(npar)

    !> copy over the momenta and initialize
    npar_jet   = npar
    njets_jet  = -1
    nproto_jet = 0
    do i = kd%ijets_lower, kd%ijets_upper
      nproto_jet = nproto_jet + 1
      jets(nproto_jet)%isJet = -1 ! all proto-jets
      jets(nproto_jet)%p(:)  =         kd%p(:,i)
      jets(nproto_jet)%pt2   = v1_pt2( kd%p(:,i) )
      jets(nproto_jet)%y     = v1_y  ( kd%p(:,i) )
      jets(nproto_jet)%phi   = v1_phi( kd%p(:,i) )
    end do

    if (nproto_jet .ne. kd%nproto) then
      print*, "initRecomb_jet: nproto don't match"
      stop
    end if

  end subroutine initRecomb_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> Apply the jet recombination
  !
  !> @param[in] RSq optional argument to override the Rcone value (squared!)
  !-----------------------------------------------------------------------------
  subroutine cluster_jet(RSq)
    double precision, optional, intent(in) :: RSq
    double precision :: RSq_val
    !>-----
    double precision :: diB_min, dij_min
    integer :: iB_min, ij_min(2)
    integer :: i,j,k

    if (present(RSq)) then
      RSq_val = RSq
    else
      RSq_val = RconeSq
    end if

    !> init distance measures & clustering history
    do i=1,nproto_jet
      clusterHist_jet(i) = i
      diB(i) = del_iB(jets(i))
      do j=i+1,nproto_jet
        dij(i,j) = del_ij(jets(i),jets(j))
      end do
    end do

    !> perform the recombination algorithm
    njets_jet = 0
    do k=1,nproto_jet
      ! find the minimum
      diB_min = -1d0
      iB_min  = -1
      dij_min = -1d0
      ij_min  = -1
      do i=1,nproto_jet
        if ( .not. jets(i)%isJet < 0 ) cycle ! not a proto-jet
        if ( (diB(i) < diB_min) .or. (diB_min < 0d0) ) then
          diB_min = diB(i)
          iB_min  = i
        endif
        do j=i+1,nproto_jet
          if ( .not. jets(j)%isJet < 0 ) cycle ! not a proto-jet
          if ( (dij(i,j) < dij_min) .or. (dij_min < 0d0) ) then
            dij_min = dij(i,j)
            ij_min  = (/i,j/)
          endif
        end do
      end do
      ! merge with beam or two proto-jets
      if ( (dij_min > 0d0) .and. (dij_min < diB_min) ) then
        ! dij smallest -> merge the two jets (`j` into `i`)
        i = ij_min(1)
        j = ij_min(2)
        call merge(jets(i),jets(j))
        ! store clustering
        clusterHist_jet(j) = i  ! 2nd argument merged into 1st
        ! update distance measures
        diB(i) = del_iB(jets(i))
        do j=1,nproto_jet
          if ( .not. jets(j)%isJet < 0 ) cycle ! not a proto-jet
          if (j==i) cycle
          if (i<j) dij(i,j) = del_ij(jets(i),jets(j))
          if (i>j) dij(j,i) = del_ij(jets(i),jets(j))
        end do
      else
        ! diB smallest -> make it a jet
        jets(iB_min)%isJet = 1
        njets_jet = njets_jet + 1
      end if
    end do

    return

  contains ! functions local to `cluster_jet`

    pure double precision function del_iB(j)
      type(Jet_t), intent(in) :: j
      del_iB = j%pt2**(ip_jet)
    end function del_iB

    pure double precision function del_ij(j1,j2)
      use Constants_mod
      type(Jet_t), intent(in) :: j1,j2
      double precision :: del_y, del_phi, deltaRSq
      del_y = j1%y - j2%y
      del_phi = abs(j1%phi - j2%phi)
      if(del_phi.gt.pi) del_phi = 2d0*pi-del_phi
      deltaRSq = del_y**2 + del_phi**2
      del_ij = min(j1%pt2**(ip_jet),j2%pt2**(ip_jet))*deltaRSq/RSq_val
    end function del_ij

    pure subroutine merge(j1,j2)
      type(Jet_t), intent(inout) :: j1,j2

      if(.not.Et_recom) then
        j1%p   = j1%p + j2%p ! add four-momenta
      else
        j1%p   = Etrecom(j1,j2) ! Use Et-weighted recombination
      endif
      j1%pt2 = v1_pt2(j1%p)
      j1%y   = v1_y  (j1%p)
      j1%phi = v1_phi(j1%p)
      ! second jet eliminated
      j2%isJet = 0
    end subroutine merge

    pure function Etrecom(j1,j2)
      use Constants_mod
      type(Jet_t), intent(in) :: j1,j2
      double precision, dimension (4) :: Etrecom
      double precision :: Et1,Et2,Etsum,eta,phi,paux1,paux2

      Et1=v1_et(j1%p)
      Et2=v1_et(j2%p)
      Etsum= Et1+Et2
      paux1=v1_phi(j1%p)
      paux2=v1_phi(j2%p)
      eta=(Et1*v1_eta(j1%p)+Et2*v1_eta(j2%p))/Etsum
      if(abs(paux1-paux2).gt.pi) then
        if(paux1.lt.paux2) then
          paux1=2d0*pi+paux1
        else
          paux2=2d0*pi+paux2
        endif
      endif

      phi=(Et1*paux1+Et2*paux2)/Etsum
      if(phi.gt.2d0*pi) phi=phi-2d0*pi
      
      !> combination such that final state is massless
      Etrecom(1)=Etsum*dcos(phi)
      Etrecom(2)=Etsum*dsin(phi)
      Etrecom(3)=Etsum*dsinh(eta)
      Etrecom(4)=Etsum*dcosh(eta)
    end function Etrecom

  end subroutine cluster_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> Apply the JADE recombination
  !
  !> @param[in] ycut optional argument to override the ycut value
  !> @param[in] compute_ytrans optional argument to switch to the computation
  !>            of the ycut transition values y_{n n+1}
  !>            the reconstructed jet momenta correspond to ycut --> infty
  !-----------------------------------------------------------------------------
  subroutine cluster_jade_jet(ycut, compute_ytrans)
    use DIS_mod
    double precision, optional, intent(in) :: ycut
    logical, optional, intent(in) :: compute_ytrans
    double precision :: ycut_val
    logical :: qytrans, qdone
    !>-----
    double precision :: dij_min, dij_last
    integer :: ij_min(2)
    integer :: i,j,k
    !> sum up all 4-momenta of particles involved in the clustering to
    !> determine the normalisation in the distance measure ( = minv)
    double precision, dimension(4) :: psum 
    double precision :: M2_norm 
    !> DIS using JADE needs another pseudo-particle
    integer :: nonpis
    common /nonpisflag/ nonpis
    type(Kin_t), pointer :: kd
    double precision :: pL_miss
    double precision, dimension(4) :: prem
    !> in DIS JADE we want to ignore "+1" (remnant) in the counting
    integer :: count_offset

    if (present(ycut)) then
      ycut_val = ycut
    else
      ycut_val = jade_ycut
    end if

    if (present(compute_ytrans)) then
      qytrans = compute_ytrans
    else
      qytrans = .false.  ! default: normal clustering till ycut_val
    end if

    count_offset = 0
    !> DIS with the JADE algorithm registers a pseudo-particle
    !> which goes along the beam direction (proton remnant)
    if (nonpis == 1) then
      !> don't count the remnant in the observables (y45,...)
      count_offset = 1
      !> prettify code
      kd => kin(npar_jet)
      !> determine missing long. momentum of the event
      pL_miss = getEp() - getEe()  ! incoming proton & electrom
      do i = 3,npar_jet
        pL_miss = pL_miss - kd%p(3,i)  ! all outgoing particles
      end do
      !> set up the momentum
      prem(:) = 0d0
      prem(3) = pL_miss       ! z direction
      prem(4) = abs(pL_miss)  ! energy
      !> register as a proto-jet
      nproto_jet = nproto_jet + 1
      jets(nproto_jet)%isJet = -1 ! init as proto-jet
      jets(nproto_jet)%p(:)  =         prem(:)
      jets(nproto_jet)%pt2   = v1_pt2( prem(:) )
      jets(nproto_jet)%y     = v1_y  ( prem(:) )
      jets(nproto_jet)%phi   = v1_phi( prem(:) )

      ! print*, "check momentum"
      ! call print_momenta(npar_jet)
      ! print*, "x1_kin   ", x1_kin
      ! print*, "parton   ", kd%p(:,1)
      ! print*, "proton   ", 0d0, 0d0,  getEp(), getEp()
      ! print*, "         ", kd%p(:,1) / x1_kin
      ! print*, "electron ", 0d0, 0d0, -getEe(), getEe()
      ! print*, "         ", kd%p(:,2)
      ! print*, "remnant  ", prem
      ! print*, '# partons', nproto_jet
      ! do k=1,nproto_jet
      !   print*, ">", k, jets(k)%isJet, jets(k)%p(:)
      ! end do
    end if

    !> get the invariant mass of the particles involved in the clustering
    psum = 0d0
    do i=1,nproto_jet
      psum(:) = psum(:) + jets(i)%p(:)
    end do
    M2_norm = v1_minv2(psum)
    ! !> for DIS, this is the same as W^2
    ! !> for epem, this is the same as roots^2
    ! print*, "M2_norm", M2_norm
    ! print*, "DIS_W2 ", getManual_obs(id_DIS_W2)
    ! read*,

    !> init the ycut transitions 
    !> everything must be valid all the time... default transition = 0
    !> better solution requires some major rewrite in how Observables are treated
    if (qytrans) then
      ! print*,
      ! print*, "------------------------ cluster_jade_jet", M2_norm
      list_cache(id_y45)%value   = 0d0
      list_cache(id_y45)%qcached = .true.
      list_cache(id_y34)%value   = 0d0
      list_cache(id_y34)%qcached = .true.
      list_cache(id_y23)%value   = 0d0
      list_cache(id_y23)%qcached = .true.
      list_cache(id_y12)%value   = 0d0
      list_cache(id_y12)%qcached = .true.
      list_cache(id_y01)%value   = 0d0
      list_cache(id_y01)%qcached = .true.
    end if

    !> init distance measures & clustering history
    do i=1,nproto_jet
      clusterHist_jet(i) = i
      ! print*, ">", i, jets(i)%isJet, jets(i)%p(:)
      do j=i+1,nproto_jet
        dij(i,j) = del_ij(jets(i),jets(j))
        ! print*, "# dij", i, j, dij(i,j)
      end do
    end do

    !> perform the recombination algorithm
    dij_last = -1d0
    qdone = .false.
    do k=1,nproto_jet

      !> find the minimum distance
      dij_min = -1d0
      ij_min  = -1
      do i=1,nproto_jet
        if ( .not. jets(i)%isJet < 0 ) cycle ! not a proto-jet
        do j=i+1,nproto_jet
          if ( .not. jets(j)%isJet < 0 ) cycle ! not a proto-jet
          if ( (dij(i,j) < dij_min) .or. (dij_min < 0d0) ) then
            dij_min = dij(i,j)
            ij_min  = (/i,j/)
          endif
        end do
      end do

      !> we're done: nothing left to recombine
      !> for the transitions, we cannot exit yet
      if (dij_min < 0d0) qdone = .true.

      !> set the ycut transitions
      if (qytrans) then
        !> dij_min cannot become smaller than the value of the last clustering
        if (dij_min > dij_last) then
          dij_last = dij_min
        else
          dij_min = dij_last
        end if
        select case (nproto_jet - k - count_offset)
          case(4)
            list_cache(id_y45)%value   = dij_min
            list_cache(id_y45)%qcached = .true.
            ! print*, "y45", dij_min, "[",k,"/",nproto_jet,"(",count_offset,")]"
          case(3)
            list_cache(id_y34)%value   = dij_min
            list_cache(id_y34)%qcached = .true.
            ! print*, "y34", dij_min, "[",k,"/",nproto_jet,"(",count_offset,")]"
          case(2)
            list_cache(id_y23)%value   = dij_min
            list_cache(id_y23)%qcached = .true.
            ! print*, "y23", dij_min, "[",k,"/",nproto_jet,"(",count_offset,")]"
          case(1)
            list_cache(id_y12)%value   = dij_min
            list_cache(id_y12)%qcached = .true.
            ! print*, "y12", dij_min, "[",k,"/",nproto_jet,"(",count_offset,")]"
          case(0)
            list_cache(id_y01)%value   = dij_min
            list_cache(id_y01)%qcached = .true.
            ! print*, "y01", dij_min, "[",k,"/",nproto_jet,"(",count_offset,")]"
        end select
        if (qdone) cycle
      end if

      !> if we're doing a "normal" clustering, check for termination
      if (qdone) exit
      if ( .not.qytrans .and. (dij_min > ycut_val) ) then
        exit  ! the recombination do-loop
      end if

      !> merge smallest dij
      i = ij_min(1)
      j = ij_min(2)
      call merge(jets(i),jets(j))
      ! print*, "merge", j, "into", i
      ! print*, ">", i, jets(i)%isJet, jets(i)%p(:)
      !> store clustering 
      clusterHist_jet(j) = i  ! 2nd argument merged into 1st
      do j=1,nproto_jet
        if ( .not. jets(j)%isJet < 0 ) cycle ! not a proto-jet
        if (j==i) cycle
        if (i<j) dij(i,j) = del_ij(jets(i),jets(j))
        if (i>j) dij(j,i) = del_ij(jets(i),jets(j))
        ! print*, "$ dij", min(i,j), max(i,j), dij(min(i,j),max(i,j))
      end do
    end do  ! k=1,nproto_jet

    !> the clustering is complete => all proto-jets are jets
    njets_jet = 0
    do i=1,nproto_jet
      if ( .not. jets(i)%isJet < 0 ) cycle  ! not a proto-jet
      jets(i)%isJet = 1
      njets_jet = njets_jet + 1
    end do

    ! !> DEBUG consistency check
    ! if ( qytrans ) then
    !   if ( (list_cache(id_y01)%value < list_cache(id_y12)%value) .or.  &
    !      & (list_cache(id_y12)%value < list_cache(id_y23)%value) .or.  &
    !      & (list_cache(id_y23)%value < list_cache(id_y34)%value) .or.  &
    !      & (list_cache(id_y34)%value < list_cache(id_y45)%value) ) then
    !     print*, "cluster_jade_jet: hierarchy wrong:"
    !     print*, "> y01", list_cache(id_y01)%value
    !     print*, "> y12", list_cache(id_y12)%value
    !     print*, "> y23", list_cache(id_y23)%value
    !     print*, "> y34", list_cache(id_y34)%value
    !     print*, "> y45", list_cache(id_y45)%value
    !     stop
    !   end if
    ! end if

    !> DIS with JADE algorithm:
    !> we should remove the jet associated with the proton remnant
    !> this is the "+1" in the DIS "(n+1) jet" event which we don't count
    if (nonpis == 1) then
      ! find out where the entry `nproto_jet` (slot for the remnant) was clustered 
      j = nproto_jet
      do while ( jets(j)%isJet == 0 )
        j = clusterHist_jet(j)
      end do
      ! eliminate the beam jet
      jets(j)%isJet = 0
      njets_jet = njets_jet -1
    end if

    return

  contains ! functions local to `cluster_jet`

    pure double precision function del_ij(j1,j2)
      use Constants_mod
      type(Jet_t), intent(in) :: j1,j2
      double precision :: Ei,Ej,cos_ij
      Ei = j1%p(4)
      Ej = j2%p(4)
      cos_ij = dot_eucl(j1%p(1:3), j2%p(1:3))  &
      &        / sqrt( dot_eucl(j1%p(1:3), j1%p(1:3)) )  &
      &        / sqrt( dot_eucl(j2%p(1:3), j2%p(1:3)) )
      del_ij = 2d0*Ei*Ej*(1d0-cos_ij) / M2_norm
    end function del_ij

    pure subroutine merge(j1,j2)
      type(Jet_t), intent(inout) :: j1,j2
      if(.not.Et_recom) then
        j1%p   = j1%p + j2%p ! add four-momenta
      else
        j1%p   = Etrecom(j1,j2) ! Use Et-weighted recombination
      endif
      j1%pt2 = v1_pt2(j1%p)
      j1%y   = v1_y  (j1%p)
      j1%phi = v1_phi(j1%p)
      ! second jet eliminated
      j2%isJet = 0
    end subroutine merge

    pure function Etrecom(j1,j2)
      use Constants_mod
      type(Jet_t), intent(in) :: j1,j2
      double precision, dimension (4) :: Etrecom
      double precision :: Et1,Et2,Etsum,eta,phi,paux1,paux2

      Et1=v1_et(j1%p)
      Et2=v1_et(j2%p)
      Etsum= Et1+Et2
      paux1=v1_phi(j1%p)
      paux2=v1_phi(j2%p)
      eta=(Et1*v1_eta(j1%p)+Et2*v1_eta(j2%p))/Etsum
      if(abs(paux1-paux2).gt.pi) then
        if(paux1.lt.paux2) then
          paux1=2d0*pi+paux1
        else
          paux2=2d0*pi+paux2
        endif
      endif

      phi=(Et1*paux1+Et2*paux2)/Etsum
      if(phi.gt.2d0*pi) phi=phi-2d0*pi
      
      !> combination such that final state is massless
      Etrecom(1)=Etsum*dcos(phi)
      Etrecom(2)=Etsum*dsin(phi)
      Etrecom(3)=Etsum*dsinh(eta)
      Etrecom(4)=Etsum*dcosh(eta)
    end function Etrecom

  end subroutine cluster_jade_jet


#ifdef USEFASTJET
  !-----------------------------------------------------------------------------
  !> @brief
  !> Apply the jet recombination and store the jets in descending order of
  !> transverse momentum in kin%pjets. Wrapper using FastJet as the backend
  !
  !> @param[in] npar kinematics specified by the number of particles
  !-----------------------------------------------------------------------------
  subroutine cluster_fastjet(npar)
    integer, intent(in) :: npar
    type(Kin_t), pointer :: kd
    double precision :: pproto(4,nproto_jet), pfjets(4,nproto_jet)
    integer :: i, nfjets
    double precision :: obs_val
    ! !>----- DEBUG vs. NNLOJET implementation
    ! logical :: qmismatch
    ! integer :: njets_nnlojet
    ! double precision, parameter :: sum_abs_threshold = 1d-11
    ! double precision :: sum_abs_delta, p_tmp(4)
    ! double precision :: pjets_nnlojet(4,nproto_jet), delta(4,nproto_jet)
    ! !>----- /DEBUG

    ! !>----- DEBUG: first run the NNLOJET clustering and store results
    ! call cluster_jet()
    ! call sort_jet()
    ! njets_nnlojet = njets_jet
    ! do i=1,njets_jet
    !   pjets_nnlojet(:,i) = jets(pt_sort_jet(i))%p(:)
    ! end do
    ! !>----- /DEBUG

    !> prettify code
    kd => kin(npar)

    !----- copy over the momenta
    nproto_jet = 0
    do i = kd%ijets_lower, kd%ijets_upper
      nproto_jet = nproto_jet + 1
      pproto(:,nproto_jet) = kd%p(:,i)
    end do

    if (nproto_jet .ne. kd%nproto) then
      print*, "cluster_fastjet: nproto don't match"
      stop
    end if

    !----- run the algorithm --> sorted in descending order wrt pt(jet)
    call fastjetppgenkt(pproto,nproto_jet,Rcone,DBLE(ip_jet),pfjets,nfjets)
    !call fastjetppgenkt(pproto,nproto,Rcone,-DBLE(ip_jet),pfjets,nfjets)

    !> copy over result into local bufer and set up information arrays
    npar_jet   = npar
    njets_jet  = nfjets
    do i=1,nfjets
      jets(i)%isJet = 1
      jets(i)%p(:)  =         pfjets(:,i)
      jets(i)%pt2   = v1_pt2( pfjets(:,i) )
      jets(i)%y     = v1_y  ( pfjets(:,i) )
      jets(i)%phi   = v1_phi( pfjets(:,i) )
      !> already sorted (note: jet cuts can still kill entries!)
      pt_sort_jet(i) = i
    end do
    !> no clustering history available
    clusterHist_jet = -1

    ! !>----- DEBUG: now check fastjet vs NNLOJET clustering results
    ! qmismatch = .false.
    ! ! test # jets
    ! if ( nfjets /= njets_nnlojet ) then
    !   print*, "DEBUG: mismatch in number of reconstructed jets:", nfjets, njets_nnlojet
    !   qmismatch = .true.
    ! end if
    ! ! test momenta of jets
    ! delta = 0d0
    ! do i=1,njets_jet
    !   delta(:,i) = pjets_nnlojet(:,i) - pfjets(:,i)
    ! end do
    ! sum_abs_delta = sum( abs(delta) )
    ! if ( .not. qmismatch .and. sum_abs_delta > sum_abs_threshold ) then
    !   !print*, "DEBUG: mismatch in momenta of jets:", sum_abs_delta
    !   qmismatch = .true.
    !   ! try fixing it by swapping j1 & j2
    !   if (njets_jet == 2) then
    !     ! swap
    !     p_tmp(:) = pjets_nnlojet(:,1)
    !     pjets_nnlojet(:,1) = pjets_nnlojet(:,2)
    !     pjets_nnlojet(:,2) = p_tmp(:)
    !     ! recompute delta
    !     delta = 0d0
    !     do i=1,njets_jet
    !       delta(:,i) = pjets_nnlojet(:,i) - pfjets(:,i)
    !     end do
    !     sum_abs_delta = sum( abs(delta) ) 
    !     if (sum_abs_delta < sum_abs_threshold) then
    !       !print*, "DEBUG: fixed by a swap j1 <-> j2"
    !       qmismatch = .false.
    !     end if
    !   end if
    ! end if
    ! ! output is there is a mismatch
    ! if (qmismatch) then
    !   print*, "njets", nfjets, njets_nnlojet
    !   do i=1,max(nfjets, njets_nnlojet)
    !     if (i<=njets_nnlojet) print*, i, "NNLOJET:", pjets_nnlojet(:,i) 
    !     if (i<=nfjets)        print*, i, "fastjet:", pfjets(:,i)
    !   end do
    !   read*,
    ! end if
    ! !>----- /DEBUG

  end subroutine cluster_fastjet
#endif

  !-----------------------------------------------------------------------------
  !> @brief
  !> Apply the jet recombination
  !
  !> @param[in] RSq optional argument to override the Rcone value (squared!)
  !> @param[out] njets optional argument to get the # of reconstructed jets
  !> @param[out] ptjets optional argument to get the pt's of the jets
  !-----------------------------------------------------------------------------
  subroutine applyCuts_jet()
    integer :: i, k
    double precision :: obs_val

    do i=1,nproto_jet
      if ( .not. jets(i)%isJet > 0 ) cycle ! not a jet
      do k=1,n_obsj
        if (.not.list_obsj(k)%qactive) cycle
        obs_val = eval_obsj(k, npar_jet, jets(i)%p)
        if ( .not. list_obsj(k)%qinvert .and.  &
        & (  &
        &   (list_obsj(k)%qlower .and. obs_val<list_obsj(k)%lower) .or.  &
        &   (list_obsj(k)%qupper .and. obs_val>list_obsj(k)%upper)  &
        & )  &
        & .or.  &
        & list_obsj(k)%qinvert .and.  &
        & (  &
        &   (list_obsj(k)%qlower .and. obs_val>=list_obsj(k)%lower .and.  &
        &    list_obsj(k)%qupper .and. obs_val<=list_obsj(k)%upper) .or.  &
        &   (list_obsj(k)%qlower .and. obs_val>=list_obsj(k)%lower .and.  &
        &    .not. list_obsj(k)%qupper ) .or.  &
        &   (list_obsj(k)%qupper .and. obs_val<=list_obsj(k)%upper .and.  &
        &    .not. list_obsj(k)%qlower )  &
        & ) ) then 
          jets(i)%isJet = 0
          clusterHist_jet(i) = 0
          njets_jet = njets_jet - 1
          exit
        end if
      end do
    end do

  end subroutine applyCuts_jet


  !-----------------------------------------------------------------------------
  !> @brief
  !> Apply the jet recombination
  !
  !> @param[in] RSq optional argument to override the Rcone value (squared!)
  !> @param[out] njets optional argument to get the # of reconstructed jets
  !> @param[out] ptjets optional argument to get the pt's of the jets
  !-----------------------------------------------------------------------------
  subroutine sort_jet()
    integer :: i,k
    double precision :: pt2_max
    integer :: i_max

    pt_sort_jet = 0
    do k=1,njets_jet
      ! find highest pt
      pt2_max = -1d0
      i_max  = -1
      do i=1,nproto_jet
        ! not a jet or already written out
        if ( .not. jets(i)%isJet > 0 .or. jets(i)%isJet == 2 ) cycle
        if ( jets(i)%pt2 > pt2_max ) then
          pt2_max = jets(i)%pt2
          i_max   = i
        end if
      end do
      ! write out
      !print*, "sort:", k, i_max, pt2_max
      jets(i_max)%isJet = 2  ! flag with "2" for sorted jets
      pt_sort_jet(k) = i_max
    end do

  end subroutine sort_jet


  !----- lepton 1 (next-to-last entry)

  double precision function eval_pt_l1(npar)
    integer, intent(in) :: npar
    eval_pt_l1 = v1_pt(kin(npar)%p(:,npar-1))
  end function eval_pt_l1

  double precision function eval_y_l1(npar) ! same as eta
    integer, intent(in) :: npar
    eval_y_l1 = v1_y(kin(npar)%p(:,npar-1))
  end function eval_y_l1


  !----- lepton 2 (last entry)

  double precision function eval_pt_l2(npar)
    integer, intent(in) :: npar
    eval_pt_l2 = v1_pt(kin(npar)%p(:,npar))
  end function eval_pt_l2

  double precision function eval_y_l2(npar) ! same as eta
    integer, intent(in) :: npar
    eval_y_l2 = v1_y(kin(npar)%p(:,npar))
  end function eval_y_l2


  !----- lepton 3 (next-to-next-to-next-to-last entry)

  double precision function eval_pt_l3(npar)
    integer, intent(in) :: npar
    eval_pt_l3 = v1_pt(kin(npar)%p(:,npar-3))
  end function eval_pt_l3

  double precision function eval_y_l3(npar) ! same as eta
    integer, intent(in) :: npar
    eval_y_l3 = v1_y(kin(npar)%p(:,npar-3))
  end function eval_y_l3


  !----- lepton 4 (next-to-next-to-last entry)

  double precision function eval_pt_l4(npar)
    integer, intent(in) :: npar
    eval_pt_l4 = v1_pt(kin(npar)%p(:,npar-2))
  end function eval_pt_l4

  double precision function eval_y_l4(npar) ! same as eta
    integer, intent(in) :: npar
    eval_y_l4 = v1_y(kin(npar)%p(:,npar-2))
  end function eval_y_l4

  double precision function eval_epem_costh_ej1(npar)
    integer, intent(in) :: npar
    double precision :: pi(1:4),pj(1:4),cos_ij,Emax,E
    integer i,imax

    !> Find the highest energy jet
    Emax=-1d0
    do i=1,kin(npar)%njets
       E=kin(npar)%pjets(4,i)
       if(E.gt.Emax) then
          Emax = E
          imax = i
       endif
    enddo

    !> leading jet
    pi = kin(npar)%pjets(:,imax)
    !> electron momentum
    pj = kin(npar)%p(:,2)
    !> the angle
    cos_ij = dot_eucl(pi(1: 3), pj(1: 3))  &
         &        / sqrt( dot_eucl(pi(1: 3), pi(1: 3)) )  &
         &        / sqrt( dot_eucl(pj(1: 3), pj(1: 3)) )
    !> only positive values 
    eval_epem_costh_ej1 = abs(cos_ij)
  end function eval_epem_costh_ej1

  double precision function eval_epem_costh_en3j(npar)
    integer, intent(in) :: npar
    double precision :: pj1(1:4),pj2(1:4),pi(1:3),pj(1:4),cos_ij
    !> leading jet
    pj1 = kin(npar)%pjets(:,1)
    !> subleading jet
    pj2 = kin(npar)%pjets(:,2)
    !> cross product j1xj2
    pi(1)=pj1(2)*pj2(3)-pj1(3)*pj2(2)
    pi(2)=pj1(3)*pj2(1)-pj1(1)*pj2(3)
    pi(3)=pj1(1)*pj2(2)-pj1(2)*pj2(1)
    !> electron momentum
    pj = kin(npar)%p(:,2)
    !> the angle
    cos_ij = dot_eucl(pi(1: 3), pj(1: 3))  &
         &        / sqrt( dot_eucl(pi(1: 3), pi(1: 3)) )  &
         &        / sqrt( dot_eucl(pj(1: 3), pj(1: 3)) )
    eval_epem_costh_en3j = abs(cos_ij)
  end function eval_epem_costh_en3j


  double precision function eval_epem_chi(npar)
    use Constants_mod
    integer, intent(in) :: npar
    double precision :: pj1(1:4),pj2(1:4),pe(1:4),pin1(1:3),pj(1:3),cos_ij
    double precision :: val,Emax,E
    integer i,imax,j

    !> Find the highest energy jet
    Emax=-1d0
    do i=1,kin(npar)%njets
       E=kin(npar)%pjets(4,i)
       if(E.gt.Emax) then
          Emax = E
          imax = i
       endif
    enddo

    if(imax.eq.1) then
       j=2
    else
       j=1
    endif

    !> leading jet
    pj1 = kin(npar)%pjets(:,imax)
    !> subleading jet
    pj2 = kin(npar)%pjets(:,j)
    !> cross product j1xj2
    pin1(1)=pj1(2)*pj2(3)-pj1(3)*pj2(2)
    pin1(2)=pj1(3)*pj2(1)-pj1(1)*pj2(3)
    pin1(3)=pj1(1)*pj2(2)-pj1(2)*pj2(1)
    !> electron momentum
    pe = kin(npar)%p(:,2)
    !> cross product j1xe
    pj(1)=pj1(2)*pe(3)-pj1(3)*pe(2)
    pj(2)=pj1(3)*pe(1)-pj1(1)*pe(3)
    pj(3)=pj1(1)*pe(2)-pj1(2)*pe(1)

    !> the angle
    cos_ij = dot_eucl(pin1(1: 3), pj(1: 3))  &
         &        / sqrt( dot_eucl(pin1(1: 3), pin1(1: 3)) )  &
         &        / sqrt( dot_eucl(pj(1: 3), pj(1: 3)) )

    val = acos(cos_ij)
    !> range 0<val<pi/2
    if(val.gt.pi/2d0) val = pi - val
    if (val.gt.pi/2d0 .or. val.lt.0d0) then
       write(*,*) val 
       stop
    endif
    eval_epem_chi = val
  end function eval_epem_chi


  !----- C parameter in ep em
  double precision function eval_epem_C(npar)
    integer, intent(in) :: npar
    double precision :: t11,t12,t13,t22,t23,t33
    double precision :: Cpar
    t11=theta(1,1,npar)
    t12=theta(1,2,npar)
    t13=theta(1,3,npar)
    t22=theta(2,2,npar)
    t23=theta(2,3,npar)
    t33=theta(3,3,npar)
    Cpar=-t12**2-t13**2-t23**2+t11*t22+t11*t33+t22*t33
    eval_epem_C = 3d0*Cpar
  end function eval_epem_C

  !----- D parameter in ep em
  double precision function eval_epem_D(npar)
    integer, intent(in) :: npar
    double precision :: t11,t12,t13,t22,t23,t33
    double precision :: Dpar
    t11=theta(1,1,npar)
    t12=theta(1,2,npar)
    t13=theta(1,3,npar)
    t22=theta(2,2,npar)
    t23=theta(2,3,npar)
    t33=theta(3,3,npar)
    Dpar=t11*t22*t33-t11*t23**2-t22*t13**2-t33*t12**2+2d0*t12*t13*t23
    Dpar=27d0*Dpar
  end function eval_epem_D
  

  !----- ordering the four lepton system in H --> 4l according to pt of leptons

  double precision function eval_pt_4l(npar,i)
    integer, intent(in) :: npar
    integer :: j,i
    integer, dimension(4) :: order
    double precision, dimension(4) :: input_pt_4l,order_pt_4l

       input_pt_4l(1) = eval_obs(id_pt_l1, npar)
       input_pt_4l(2) = eval_obs(id_pt_l2, npar)
       input_pt_4l(3) = eval_obs(id_pt_l3, npar)
       input_pt_4l(4) = eval_obs(id_pt_l4, npar)

       order(1)=0
       order(2)=0
       order(3)=0
       order(4)=0

       order_pt_4l(1)=0d0
       order_pt_4l(2)=0d0
       order_pt_4l(3)=0d0
       order_pt_4l(4)=0d0


    do j=1,4
       if (input_pt_4l(j) > order_pt_4l(1)) then
          order(1)=j
          order_pt_4l(1)=input_pt_4l(j)
       endif
    enddo
 
    do j=1,4
       if (j.ne.order(1)) then
          if (input_pt_4l(j) > order_pt_4l(2)) then
             order_pt_4l(2)=input_pt_4l(j)
             order(2)=j
          endif
       endif
    enddo

    do j=1,4
       if ((j.ne.order(1)).and.(j.ne.order(2))) then
          if (input_pt_4l(j) > order_pt_4l(3)) then
             order_pt_4l(3)=input_pt_4l(j)
             order(3)=j
          endif
       endif
    enddo

    do j=1,4
       if ((j.ne.order(1)).and.(j.ne.order(2)).and.(j.ne.order(3))) then
             order_pt_4l(4)=input_pt_4l(j)
             order(4)=j
       endif
    enddo

    eval_pt_4l=order_pt_4l(i)

  end function eval_pt_4l


  double precision function eval_y_4l(npar,i)
    integer, intent(in) :: npar
    integer :: j,i
    integer, dimension(4) :: order
    double precision, dimension(4) :: input_pt_4l,order_pt_4l,y_4l

       input_pt_4l(1) = eval_obs(id_pt_l1, npar)
       input_pt_4l(2) = eval_obs(id_pt_l2, npar)
       input_pt_4l(3) = eval_obs(id_pt_l3, npar)
       input_pt_4l(4) = eval_obs(id_pt_l4, npar)

       y_4l(1) = eval_obs(id_y_l1, npar)
       y_4l(2) = eval_obs(id_y_l2, npar)
       y_4l(3) = eval_obs(id_y_l3, npar)
       y_4l(4) = eval_obs(id_y_l4, npar)

       order(1)=0
       order(2)=0
       order(3)=0
       order(4)=0

       order_pt_4l(1)=0d0
       order_pt_4l(2)=0d0
       order_pt_4l(3)=0d0
       order_pt_4l(4)=0d0


    do j=1,4
       if (input_pt_4l(j) > order_pt_4l(1)) then
          order(1)=j
          order_pt_4l(1)=input_pt_4l(j)
       endif
    enddo
 
    do j=1,4
       if (j.ne.order(1)) then
          if (input_pt_4l(j) > order_pt_4l(2)) then
             order_pt_4l(2)=input_pt_4l(j)
             order(2)=j
          endif
       endif
    enddo

    do j=1,4
       if ((j.ne.order(1)).and.(j.ne.order(2))) then
          if (input_pt_4l(j) > order_pt_4l(3)) then
             order_pt_4l(3)=input_pt_4l(j)
             order(3)=j
          endif
       endif
    enddo

    do j=1,4
       if ((j.ne.order(1)).and.(j.ne.order(2)).and.(j.ne.order(3))) then
             order_pt_4l(4)=input_pt_4l(j)
             order(4)=j
       endif
    enddo

    eval_y_4l=y_4l(order(i))

  end function eval_y_4l


  !----- leading & sub-leading photon / leptons in Z (last two entries)

  !> this routine sets the indices of which is the leading/sub-leading one
  recursive subroutine set_idx_g12(npar)
    integer, intent(in) :: npar
    double precision :: pt_l1, pt_l2
    double precision :: rnd
    integer :: npar_from
    !integer :: idx1_old, idx2_old, idx1_new, idx2_new
    !> get pt's
    pt_l1 = eval_obs(id_pt_l1, npar)
    pt_l2 = eval_obs(id_pt_l2, npar)
    !> get mapping situation
    npar_from = kin(npar)%npar_from

    ! print*, ">>>>> set_idx_g12", npar, ",", npar_from
    ! print*, "pt_l1", pt_l1
    ! print*, "pt_l2", pt_l2
    ! print*, "y_l1 ", eval_obs(id_y_l1, npar)
    ! print*, "y_l2 ", eval_obs(id_y_l2, npar)

    if ( abs((pt_l1-pt_l2)/(pt_l1+pt_l2)) > 1d-6 .or. npar > 4 ) then
      !> if they differ in the ~6th digit or more we use the normal ordering
      !> if I'm in a 2->n configuration with n>2 we never randomise
      if (pt_l1 > pt_l2) then
        list_cache(id_idx_g1)%value = DBLE(1)
        list_cache(id_idx_g2)%value = DBLE(2)
      else
        list_cache(id_idx_g1)%value = DBLE(2)
        list_cache(id_idx_g2)%value = DBLE(1)
      end if
    else
      !> if the pt's are the "same" the assignment becomes ambiguous
      if ( npar < npar_from ) then
        !> we're in a subtraction term
        !>  => determine assignment from real-emission kinemtics
        ! print*, "JUMP UP", npar_from
        call clearCache_obs()
        call set_idx_g12(npar_from)
        call clearCache_obs()
      else
        !> out of options
        !>  => randomly choose an assignemnt  
        call random_number(rnd)
        ! print*, "rnd", rnd
        if (rnd < 0.5d0) then
          list_cache(id_idx_g1)%value = DBLE(1)
          list_cache(id_idx_g2)%value = DBLE(2)
        else
          list_cache(id_idx_g1)%value = DBLE(2)
          list_cache(id_idx_g2)%value = DBLE(1)
        end if
      end if
    end if
    list_cache(id_idx_g1)%qcached = .true.
    list_cache(id_idx_g2)%qcached = .true.
    ! print*, "idx1", INT( list_cache(id_idx_g1)%value )
    ! print*, "idx2", INT( list_cache(id_idx_g2)%value )
    ! read*,

    !> debugging
    ! if (npar_max-npar .eq. 1) then
    !   idx1_old = INT( list_cache(id_idx_g1)%value )
    !   idx2_old = INT( list_cache(id_idx_g2)%value )
    !   call clearCache_obs()
    !   call set_idx_g12(npar_max)
    !   call clearCache_obs()
    !   idx1_new = INT( list_cache(id_idx_g1)%value )
    !   idx2_new = INT( list_cache(id_idx_g2)%value )
    !   if (idx1_new /= idx1_old) then
    !     print*, "-----------------------------------------"
    !     print*, "old", idx1_old, idx2_old
    !     print*, "new", idx1_new, idx2_new
    !     print*, "-----------------------------------------"
    !     call print_momenta(npar_max)
    !     call print_momenta(npar)
    !     read*,
    !   end if
    ! end if

  end subroutine set_idx_g12

  double precision function eval_pt_g1(npar)
    integer, intent(in) :: npar
    integer :: idx
    idx = INT( eval_obs(id_idx_g1, npar) )
    if (idx == 1) then
      eval_pt_g1 = eval_obs(id_pt_l1, npar)
    else
      eval_pt_g1 = eval_obs(id_pt_l2, npar)
    end if
  end function eval_pt_g1

  double precision function eval_pt_g2(npar)
    integer, intent(in) :: npar
    integer :: idx
    idx = INT( eval_obs(id_idx_g2, npar) )
    if (idx == 1) then
      eval_pt_g2 = eval_obs(id_pt_l1, npar)
    else
      eval_pt_g2 = eval_obs(id_pt_l2, npar)
    end if
  end function eval_pt_g2

  double precision function eval_y_g1(npar) ! same as eta
    integer, intent(in) :: npar
    integer :: idx
    idx = INT( eval_obs(id_idx_g1, npar) )
    if (idx == 1) then
      eval_y_g1 = eval_obs(id_y_l1, npar)
    else
      eval_y_g1 = eval_obs(id_y_l2, npar)
    end if
  end function eval_y_g1

  double precision function eval_y_g2(npar) ! same as eta
    integer, intent(in) :: npar
    integer :: idx
    idx = INT( eval_obs(id_idx_g2, npar) )
    if (idx == 1) then
      eval_y_g2 = eval_obs(id_y_l1, npar)
    else
      eval_y_g2 = eval_obs(id_y_l2, npar)
    end if
  end function eval_y_g2

  double precision function eval_deltaphi_g1g2(npar)
    integer, intent(in) :: npar
    eval_deltaphi_g1g2 = v2_delta_phi( kin(npar)%p(:,npar) , & 
                                     & kin(npar)%p(:,npar-1) )
  end function eval_deltaphi_g1g2

  double precision function eval_deltay_g1g2(npar)
    integer, intent(in) :: npar
    double precision :: y_g1, y_g2
    y_g1 = eval_obs(id_y_g1, npar)  ! leading
    y_g2 = eval_obs(id_y_g2, npar)  ! sub-leading
    eval_deltay_g1g2 = y_g1 - y_g2  ! no abs(...) !
  end function eval_deltay_g1g2

  double precision function eval_ptt_g1g2(npar)
    integer, intent(in) :: npar
    double precision :: pt_g1, pt_g2
    pt_g1=v1_pt(kin(npar)%p(:,npar-1))
    pt_g2=v1_pt(kin(npar)%p(:,npar))
    if (pt_g1.ge.pt_g2) then
      eval_ptt_g1g2=(kin(npar)%p(1,npar-1)*kin(npar)%p(2,npar) &
                 &  -kin(npar)%p(1,npar)*kin(npar)%p(2,npar-1))/2/(pt_g1-pt_g2)
    else
      eval_ptt_g1g2=(kin(npar)%p(1,npar)*kin(npar)%p(2,npar-1) &
                 &  -kin(npar)%p(1,npar-1)*kin(npar)%p(2,npar))/2/(pt_g2-pt_g1)
    endif
  end function eval_ptt_g1g2


  double precision function eval_deltaphi_Vj1(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_deltaphi_Vj1 = v2_delta_phi(pV , kin(npar)%pjets(:,1) )
  end function eval_deltaphi_Vj1

  double precision function eval_deltay_Vj1(npar)
    integer, intent(in) :: npar
    double precision :: yV, yj1
    yV  = eval_obs(id_y_V, npar)
    yj1 = eval_obs(id_y_j1, npar)
    eval_deltay_Vj1 = abs(yV-yj1)
  end function eval_deltay_Vj1

  double precision function eval_ydiff_Vj1(npar)
    integer, intent(in) :: npar
    double precision :: deltay
    deltay = eval_obs(id_deltay_Vj1, npar)
    eval_ydiff_Vj1 = deltay / 2d0
  end function eval_ydiff_Vj1

  double precision function eval_ysum_Vj1(npar)
    integer, intent(in) :: npar
    double precision :: yV, yj1
    yV  = eval_obs(id_y_V, npar)
    yj1 = eval_obs(id_y_j1, npar)
    eval_ysum_Vj1 = abs(yV+yj1) / 2d0
  end function eval_ysum_Vj1


  subroutine set_Collins_Soper(npar)
    use Constants_mod
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    double precision :: QV,QV2, PTV,PTV2, XTV
    double precision :: pl_x,pl_y,pl_z
    integer :: i
    ! !> naive implementation for cross checks
    ! double precision, dimension(4) :: bp,bm,pl
    ! double precision, dimension(3) :: ep,em, ex,ey,ez
    ! double precision :: pl_x_alt,pl_y_alt,pl_z_alt, dummy

    !> momentum of intermediate gauge boson
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do

    !> Eq.(23) of hep-ph/9406381 
    !> with corrected signs and generalised to pV(2) /= 0
    !> 'n-1':  l1 = lm (for Z)
    QV2  = dot(pV,pV)
    PTV2 = pV(1)**2 + pV(2)**2
    QV   = sqrt(QV2)
    PTV  = sqrt(PTV2)
    XTV  = sqrt(QV2 + PTV2)
    pl_x = +PTV/QV/XTV * (  &
         &    -pV(4)*kin(npar)%p(4,npar-1)  &
         &    +pV(3)*kin(npar)%p(3,npar-1) )  &
         & +XTV/QV/PTV * (  &
         &    +pV(1)*kin(npar)%p(1,npar-1)  &
         &    +pV(2)*kin(npar)%p(2,npar-1) ) 
    pl_y = sign(1d0,pV(3))/PTV * (  &
         &    -pV(2)*kin(npar)%p(1,npar-1)  &
         &    +pV(1)*kin(npar)%p(2,npar-1) ) 
    pl_z = sign(1d0,pV(3))/XTV * ( &
         &    -pV(3)*kin(npar)%p(4,npar-1)  &
         &    +pV(4)*kin(npar)%p(3,npar-1) )

    ! !----- BEGIN: naive implementation
    ! !> beam directions
    ! ! plus (+)
    ! bp(1) =  0d0
    ! bp(2) =  0d0
    ! bp(3) =  1d0
    ! bp(4) =  1d0
    ! ! minus (-)
    ! bm(1) =  0d0
    ! bm(2) =  0d0
    ! bm(3) = -1d0
    ! bm(4) =  1d0
    ! !> lepton momentum 
    ! pl(1:4) = kin(npar)%p(1:4,npar-1)  ! = l1 = lm (for Z)
    ! !> boost into the gauge-boson rest frame
    ! call boost_rest(bp,pV)
    ! call boost_rest(bm,pV)
    ! call boost_rest(pl,pV)
    ! !> only interested in directions
    ! ep(1:3) = bp(1:3)
    ! em(1:3) = bm(1:3)
    ! ! normalise 
    ! ep = ep / sqrt( dot_eucl(ep,ep) )
    ! em = em / sqrt( dot_eucl(em,em) )
    ! !> set up the Collins-Spoer coordinate system
    ! ! the z-axis is the bisector of (bp,-bm) (direction of pV in lab frame)
    ! if (pV(3) > 0d0) then
    !   ez = ep-em
    ! else
    !   ez = em-ep
    ! end if
    ! ez = ez / sqrt( dot_eucl(ez,ez) )
    ! ! orthogonal to z & in the ep-em-plane (opposite direction to bp+bm)
    ! ex = -ep-em
    ! ex = ex / sqrt( dot_eucl(ex,ex) )
    ! ! y-axis to get a right-handed coordinate system: ey = ez x ex
    ! ey(1) = ez(2)*ex(3) - ez(3)*ex(2)
    ! ey(2) = ez(3)*ex(1) - ez(1)*ex(3)
    ! ey(3) = ez(1)*ex(2) - ez(2)*ex(1)
    ! !> project the lepton momentum onto the axes
    ! pl_x_alt = dot_eucl(pl(1:3),ex)
    ! pl_y_alt = dot_eucl(pl(1:3),ey)
    ! pl_z_alt = dot_eucl(pl(1:3),ez)
    ! dummy = atan2(pl_y_alt, pl_x_alt)
    ! if (dummy < 0d0) dummy = dummy + 2d0*pi    
    ! !read*,
    ! if ( abs(pl_x/pl_x_alt-1d0)>1d-6 .or. &
    !    & abs(pl_y/pl_y_alt-1d0)>1d-6 .or. &
    !    & abs(pl_z/pl_z_alt-1d0)>1d-6 ) then
    !   print*, "pl_x", pl_x, pl_x_alt, pl_x/pl_x_alt
    !   print*, "pl_y", pl_y, pl_y_alt, pl_y/pl_y_alt
    !   print*, "pl_z", pl_z, pl_z_alt, pl_z/pl_z_alt
    !   read*,
    ! end if
    ! !----- END: naive implementation

    list_cache(id_cos_theta_CS)%value = pl_z / sqrt(pl_x**2 + pl_y**2 + pl_z**2)
    
    list_cache(id_theta_CS)%value = atan2( sqrt(pl_x**2+pl_y**2) , pl_z )

    list_cache(id_phi_CS)%value = atan2(pl_y, pl_x)
    if (list_cache(id_phi_CS)%value < 0d0) then
      list_cache(id_phi_CS)%value = list_cache(id_phi_CS)%value+2d0*pi
    end if

    list_cache(id_cos_theta_CS)%qcached = .true.
    list_cache(id_theta_CS)%qcached     = .true.
    list_cache(id_phi_CS)%qcached       = .true.
  end subroutine set_Collins_Soper


  !> routine to find the thrust axis for ep em collisions

  subroutine set_Taxis(npar)
    integer, intent(in) :: npar
    integer :: i,j,k,n,ithrust
    double precision, dimension(4) :: Taxis
    double precision :: pp(5,15),Energy,smoms
    double precision :: bplus,bminus,bot,bsum,bmax,bmin
    double precision :: em2h,em2l,t,tmp,tmp1,tmp2,Tpar,absp
    
    !> vector of possible thrust directions
    pp=0d0       
    n=0
    !> sum of final state energy
    Energy=0d0
    !> sum of final state scalar momenta
    smoms=0d0

    do j=3,npar
       n=n+1
       Energy = Energy+kin(npar)%p(4,j)
       smoms  = smoms+sqrt(kin(npar)%p(1,j)**2+kin(npar)%p(2,j)**2+kin(npar)%p(3,j)**2)
       do i=1,4
          pp(i,n)=kin(npar)%p(i,j)
       enddo
       pp(5,n)=pp(4,n)
    end do

    if(npar.eq.6)then
       do k=4,npar 
          n=n+1
          do i=1,4
             pp(i,n)=kin(npar)%p(i,3)+kin(npar)%p(i,k)
          enddo
          pp(5,n)=sqrt(pp(1,n)**2+pp(2,n)**2+pp(3,n)**2)
       enddo
    elseif(npar.eq.7)then
       do j=3,npar-1
          do k=j+1,npar 
             n=n+1
             do i=1,4
                pp(i,n)=kin(npar)%p(i,j)+kin(npar)%p(i,k)
             enddo
             pp(5,n)=sqrt(pp(1,n)**2+pp(2,n)**2+pp(3,n)**2)
          enddo
       enddo
    endif
    
    !>
    !> Find thrust axis and calculate thrust
    !> ithrust gives index for thrust axis
    ithrust=0
    Tpar=0d0
    do i=1,15
       t=2d0*pp(5,i)
       if(t.gt.Tpar) then
          Tpar=t
          ithrust=i
       endif
    enddo
    !> thrust axis (unit vector)
    Taxis(4) = Tpar/smoms
    do j=1,3
       Taxis(j) = 2d0*pp(j,ithrust)/Tpar
    enddo
    
    !> angle of thrust axis
    list_cache(id_epem_cosT)%value=dsqrt(Taxis(1)**2+Taxis(2)**2)
    !> writing results to cache
    list_cache(id_epem_T)%value   = Taxis(4)
    list_cache(id_epem_1mT)%value = 1d0-Taxis(4)
   
    list_cache(id_epem_T)%qcached    = .true.
    list_cache(id_epem_1mT)%qcached  = .true.
    list_cache(id_epem_cosT)%qcached = .true.

    !> evaluate all related results
    !> HJM (heavy jet mass)
    em2h=pp(4,ithrust)**2-pp(5,ithrust)**2
    em2l=(Energy-pp(4,ithrust))**2-pp(5,ithrust)**2
    if(em2l.gt.em2h) then
       tmp=em2l
       em2l=em2h
       em2h=tmp
    endif
    list_cache(id_epem_HJM)%value   = em2h/energy**2
    list_cache(id_epem_HJM)%qcached = .true.
    list_cache(id_epem_SHM)%value   = (em2h+em2l)/energy**2
    list_cache(id_epem_SHM)%qcached = .true.

    !>
    !> Hemisphere broadening
    !>
    bplus=0d0
    bminus=0d0
    bot=0d0
    absp = 0d0
    !> 1 to npar-2 are the normal momenta
    do i=1,npar-2
       absp = sqrt(pp(1,i)**2 + pp(2,i)**2 + pp(3,i)**2)
       bot=bot+absp
       if(i.ne.ithrust)then
          !> |p|cos(theta_(iT))
          tmp1 = pp(1,i)*Taxis(1)+pp(2,i)*Taxis(2)+pp(3,i)*Taxis(3)
          !> |p|sin(theta_(iT))
          tmp2 = sqrt(absp**2-tmp1**2)
          if(tmp1.gt.0d0)then
             bplus=bplus+tmp2 
          elseif(tmp1.lt.0d0)then
             bminus=bminus+tmp2 
          endif
       endif
    enddo
    bplus=bplus/bot/2d0
    bminus=bminus/bot/2d0
    bmax=max(bplus,bminus)
    bmin=min(bplus,bminus)
    bsum =bmax+bmin 

    list_cache(id_epem_TJB)%value = bsum
    list_cache(id_epem_WJB)%value = bmax

    list_cache(id_epem_TJB)%qcached = .true.
    list_cache(id_epem_WJB)%qcached = .true.
    
  end subroutine set_Taxis

  !> Calculate entries in theta matrix to determine C parameter
  pure double precision function theta(ii,ij,npar)
    integer, intent(in) :: ii,ij,npar
    double precision :: top,bot,absp
    integer :: i

    top  = 0d0
    bot  = 0d0
    !> magnitude of impulse
    absp = 0d0

    do i=3,npar
       absp = sqrt(kin(npar)%p(1,i)**2 + &
            & kin(npar)%p(2,i)**2 + kin(npar)%p(3,i)**2 )
   
       top=top+kin(npar)%p(ii,i)*kin(npar)%p(ij,i)/absp
       bot=bot+absp
    enddo
    theta=top/bot
    return
  end function theta

  !> Calculate entries in theta matrix to determine C parameter for DIS
  pure double precision function thetadis(ii,ij,npar)
    integer, intent(in) :: ii,ij,npar
    double precision :: top,bot,absp,ehtot
    integer :: i,ncount

    top  = 0d0
    bot  = 0d0
    !> magnitude of impulse
    absp = 0d0
    ehtot = 0d0
    ncount=0

    do i=4,npar
        ! Select only particles in the current hemisphere (pz > 0)        
        if (kin(npar)%p(3,i) > 0d0) then       
           ncount=ncount+1
           absp = sqrt(kin(npar)%p(1,i)**2 + &
                 & kin(npar)%p(2,i)**2 + kin(npar)%p(3,i)**2 )
            top=top+kin(npar)%p(ii,i)*kin(npar)%p(ij,i)/absp
            bot=bot+absp
        end if
    enddo
    if(ncount.eq.0) then
       thetadis=-1d0
    else
       thetadis=top/bot
    endif
    return
  end function thetadis

  !> projectors from [1606.00689]
  subroutine set_proj_Ai(npar)
    integer, intent(in) :: npar
    !> read in cached observables
    double precision :: cos_th_cs, th_cs, phi_cs
    !> locally cache more
    double precision :: cos2_th_cs, sin2_th_cs, sin_th_cs, sin_2th_cs
    double precision :: cos_phi_cs, sin_phi_cs, cos_2phi_cs, sin_2phi_cs

    cos_th_cs   = eval_obs(id_cos_theta_CS, npar)
    th_cs       = eval_obs(id_theta_CS, npar)
    phi_cs      = eval_obs(id_phi_CS, npar)

    cos2_th_cs  = cos_th_cs**2
    sin2_th_cs  = (1d0-cos_th_cs)*(1d0+cos_th_cs)
    sin_th_cs   = sqrt(sin2_th_cs)
    sin_2th_cs  = 2d0*cos_th_cs*sin_th_cs
    cos_phi_cs  = cos(phi_cs)
    sin_phi_cs  = sin(phi_cs)
    cos_2phi_cs = (cos_phi_cs-sin_phi_cs)*(cos_phi_cs+sin_phi_cs)
    sin_2phi_cs = 2d0*cos_phi_cs*sin_phi_cs

    ! !> sanity checks:
    ! print*, "cos_th_cs ", cos_th_cs     / cos(th_cs)
    ! print*, "cos2_th_cs", cos2_th_cs    / cos(th_cs)**2
    ! print*, "sin2_th_cs", sin2_th_cs    / sin(th_cs)**2
    ! print*, "sin_th_cs ", sin_th_cs     / sin(th_cs)
    ! print*, "sin_2th_cs", sin_2th_cs    / sin(2d0*th_cs)
    ! print*, "sin_2phi_cs", sin_2phi_cs  / sin(2d0*phi_cs)
    ! print*, "cos_2phi_cs", cos_2phi_cs  / cos(2d0*phi_cs)
    ! read*,

    list_cache(id_proj_A0)%value       = 4d0-10d0*cos2_th_cs
    list_cache(id_proj_A1)%value       = 5d0*sin_2th_cs*cos_phi_cs
    list_cache(id_proj_A2)%value       = 10d0*sin2_th_cs*cos_2phi_cs
    list_cache(id_proj_A3)%value       = 4d0*sin_th_cs*cos_phi_cs
    list_cache(id_proj_A4)%value       = 4d0*cos_th_cs
    list_cache(id_proj_A5)%value       = 5d0*sin2_th_cs*sin_2phi_cs
    list_cache(id_proj_A6)%value       = 5d0*sin_2th_cs*sin_phi_cs
    list_cache(id_proj_A7)%value       = 4d0*sin_th_cs*sin_phi_cs
    list_cache(id_proj_A0mA2)%value    =  list_cache(id_proj_A0)%value  &
                                       & -list_cache(id_proj_A2)%value
    list_cache(id_proj_A0mA2alt)%value =  -6d0 +20d0*sin2_th_cs*sin_phi_cs**2

    list_cache(id_proj_A0)%qcached       = .true.
    list_cache(id_proj_A1)%qcached       = .true.
    list_cache(id_proj_A2)%qcached       = .true.
    list_cache(id_proj_A3)%qcached       = .true.
    list_cache(id_proj_A4)%qcached       = .true.
    list_cache(id_proj_A5)%qcached       = .true.
    list_cache(id_proj_A6)%qcached       = .true.
    list_cache(id_proj_A7)%qcached       = .true.
    list_cache(id_proj_A0mA2)%qcached    = .true.
    list_cache(id_proj_A0mA2alt)%qcached = .true.

  end subroutine set_proj_Ai

  double precision function eval_cos_theta_CS_0PT(npar)
    !> Definition of cos(theta*) that is also valid at ptZ = 0
    !> so can be used for Z+0j production/ cases where ptZ = 0
    !> Tested for ZJ(J)LO against eval_cos_theta_CS and agrees perfectly
    integer, intent(in) :: npar
    integer :: i
    double precision, dimension(4) :: pV
    double precision :: mll, pzll, ptll, p1p, p1m, p2p, p2m
    !> 'n-1':  l1 = lm (for Z)
    !> 3 = z component
    !> 4 = E component
    do i=1,4 ! Only need first 3 components
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    mll = v1_minv(pV)
    pzll = pV(3)
    ptll = sqrt(pV(1)**2+pV(2)**2)
    p1p = kin(npar)%p(4,npar-1)+kin(npar)%p(3,npar-1)
    p1m = kin(npar)%p(4,npar-1)-kin(npar)%p(3,npar-1)
    p2p = kin(npar)%p(4,npar)+kin(npar)%p(3,npar)
    p2m = kin(npar)%p(4,npar)-kin(npar)%p(3,npar)
    eval_cos_theta_CS_0PT = sign(1d0,pzll)*((p1p*p2m)-(p1m*p2p))/mll/sqrt(mll**2+ptll**2)
  end function eval_cos_theta_CS_0PT
  

  double precision function eval_costhetastar(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pg1,pg2,pH
    double precision :: ptg1,ptg2,ptH
    integer :: i
    do i=1,4
      pg1(i) = kin(npar)%p(i,npar) 
      pg2(i) = kin(npar)%p(i,npar-1) 
      pH(i)  = kin(npar)%p(i,npar)+kin(npar)%p(i,npar-1) 
    end do
      ptg1=v1_pt(pg1)
      ptg2=v1_pt(pg2)
      ptH =v1_pt(pH)
      eval_costhetastar= abs(dsinh(v1_y(pg1)-v1_y(pg2)))*2d0*ptg1*ptg2 &
                       & /v1_minv2(pH)/dsqrt(1d0+(ptH/v1_minv(pH))**2d0)
  end function eval_costhetastar

  !----- gauge boson (sum of the related leptons)

  double precision function eval_pt_V(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_pt_V = v1_pt(pV)
    !> this is slow because it induces an allocation
    !eval_pt_V = v1_pt( kin(npar)%p(:,npar-1) + kin(npar)%p(:,npar) )
  end function eval_pt_V

  double precision function eval_pt_V2(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2)
    end do
    eval_pt_V2 = v1_pt(pV)
    !> this is slow because it induces an allocation
    !eval_pt_V2 = v1_pt( kin(npar)%p(:,npar-1) + kin(npar)%p(:,npar) )
  end function eval_pt_V2

  double precision function eval_pt_V4l(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) =   kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2) &
            & + kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_pt_V4l = v1_pt(pV)
  end function eval_pt_V4l

  double precision function eval_minv_V(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_minv_V = v1_minv(pV)
    !> this is slow because it induces an allocation
    !eval_minv_V = v1_minv( kin(npar)%p(:,npar-1) + kin(npar)%p(:,npar) )
  end function eval_minv_V

  double precision function eval_minv_V2(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2)
    end do
    eval_minv_V2 = v1_minv(pV)
    !> this is slow because it induces an allocation
    !eval_minv_V2 = v1_minv( kin(npar)%p(:,npar-1) + kin(npar)%p(:,npar) )
  end function eval_minv_V2

  double precision function eval_minv_V4l(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) =   kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2) &
            & + kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_minv_V4l = v1_minv(pV)
  end function eval_minv_V4l

  double precision function eval_minv_V23(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-2) + kin(npar)%p(i,npar-1)
    end do
    eval_minv_V23 = v1_minv(pV)
  end function eval_minv_V23

  double precision function eval_minv_V14(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar) + kin(npar)%p(i,npar-3)
    end do
    eval_minv_V14 = v1_minv(pV)
  end function eval_minv_V14

  double precision function eval_minv_mZ1(npar)
    integer, intent(in) :: npar
    double precision :: mV1, mV2
    double precision emz, ezwidth
    common /Zmass/ emz, ezwidth
    mV1 = eval_obs(id_minv_V,  npar)
    mV2 = eval_obs(id_minv_V2, npar)
    if ( abs(mV1-emz) < abs(mV2-emz) ) then
      eval_minv_mZ1 = mV1
    else
      eval_minv_mZ1 = mV2
    end if
  end function eval_minv_mZ1

  double precision function eval_minvSFOS_vs_mZ(npar)
    integer, intent(in) :: npar
    double precision :: m2e, m2mu
    double precision emz, ezwidth
    common /Zmass/ emz, ezwidth
    m2mu = eval_obs(id_minv_V,  npar)
    m2e = eval_obs(id_minv_V2, npar)
    if ( abs(m2mu-emz) < abs(m2e-emz) ) then
      eval_minvSFOS_vs_mZ = 2d0
    else
      eval_minvSFOS_vs_mZ = 1d0
    end if
  end function eval_minvSFOS_vs_mZ

  double precision function eval_minv_mZSFOS(npar,i)
    integer, intent(in) :: npar
    integer :: j,i
    integer, dimension(4) :: order
    double precision, dimension(4) :: input_minv_4l,order_minv_4l,diff_minv_4l
    double precision emz, ezwidth
    common /Zmass/ emz, ezwidth

       input_minv_4l(1) = eval_obs(id_minv_V, npar)
       input_minv_4l(2) = eval_obs(id_minv_V2, npar)
       input_minv_4l(3) = eval_obs(id_minv_V23, npar)
       input_minv_4l(4) = eval_obs(id_minv_V14, npar)

       diff_minv_4l(1) = abs(input_minv_4l(1)-emz)
       diff_minv_4l(2) = abs(input_minv_4l(2)-emz)
       diff_minv_4l(3) = abs(input_minv_4l(3)-emz)
       diff_minv_4l(4) = abs(input_minv_4l(4)-emz)

       order(1)=0
       order(2)=0
       order(3)=0
       order(4)=0

       order_minv_4l(1)=0d0
       order_minv_4l(2)=0d0
       order_minv_4l(3)=0d0
       order_minv_4l(4)=0d0


    do j=1,4
       if (diff_minv_4l(j) > order_minv_4l(1)) then
          order(1)=j
          order_minv_4l(1)=input_minv_4l(j)
       endif
    enddo
 
    do j=1,4
       if (j.ne.order(1)) then
          if (diff_minv_4l(j) > order_minv_4l(2)) then
             order_minv_4l(2)=input_minv_4l(j)
             order(2)=j
          endif
       endif
    enddo

    do j=1,4
       if ((j.ne.order(1)).and.(j.ne.order(2))) then
          if (diff_minv_4l(j) > order_minv_4l(3)) then
             order_minv_4l(3)=input_minv_4l(j)
             order(3)=j
          endif
       endif
    enddo

    do j=1,4
       if ((j.ne.order(1)).and.(j.ne.order(2)).and.(j.ne.order(3))) then
             order_minv_4l(4)=input_minv_4l(j)
             order(4)=j
       endif
    enddo

    eval_minv_mZSFOS=order_minv_4l(i)

  end function eval_minv_mZSFOS

  double precision function eval_minv_mZ2(npar)
    integer, intent(in) :: npar
    double precision :: mV1, mV2
    double precision emz, ezwidth
    common /Zmass/ emz, ezwidth
    mV1 = eval_obs(id_minv_V,  npar)
    mV2 = eval_obs(id_minv_V2, npar)
    if ( abs(mV1-emz) > abs(mV2-emz) ) then
      eval_minv_mZ2 = mV1
    else
      eval_minv_mZ2 = mV2
    end if
  end function eval_minv_mZ2

  double precision function eval_y_V(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_y_V = v1_y(pV)
    !> this is slow because it induces an allocation
    !eval_y_V = v1_y( kin(npar)%p(:,npar-1) + kin(npar)%p(:,npar) )
  end function eval_y_V

  double precision function eval_y_V2(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2)
    end do
    eval_y_V2 = v1_y(pV)
    !> this is slow because it induces an allocation
    !eval_y_V2 = v1_y( kin(npar)%p(:,npar-1) + kin(npar)%p(:,npar) )
  end function eval_y_V2

  double precision function eval_y_V4l(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) =   kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2) &
            & + kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_y_V4l = v1_y(pV)
  end function eval_y_V4l

  double precision function eval_eta_V(npar) ! not same as y_eta
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_eta_V = v1_eta(pV)
    !> this is slow because it induces an allocation
    !eval_eta_V = v1_eta( kin(npar)%p(:,npar-1) + kin(npar)%p(:,npar) )
  end function eval_eta_V

  double precision function eval_eta_V2(npar) ! not same as y_eta
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2)
    end do
    eval_eta_V2 = v1_eta(pV)
    !> this is slow because it induces an allocation
    !eval_eta_V2 = v1_eta( kin(npar)%p(:,npar-1) + kin(npar)%p(:,npar) )
  end function eval_eta_V2

  double precision function eval_eta_V4l(npar) ! not same as y_eta
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) =   kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2) &
            & + kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_eta_V4l = v1_eta(pV)
  end function eval_eta_V4l

  double precision function eval_mt_V(npar)
    integer, intent(in) :: npar
    !> we assume massless daugthers (lepton, neutrino) here
    eval_mt_V = sqrt( 2d0*v1_pt(kin(npar)%p(:,npar-1))*v1_pt(kin(npar)%p(:,npar)) &
              & *(1d0-cos(v2_delta_phi(kin(npar)%p(:,npar-1),kin(npar)%p(:,npar)))) )
  end function eval_mt_V

  double precision function eval_mt_V2(npar)
    integer, intent(in) :: npar
    !> we assume massless daugthers (lepton, neutrino) here
    eval_mt_V2 = sqrt( 2d0*v1_pt(kin(npar)%p(:,npar-3))*v1_pt(kin(npar)%p(:,npar-2)) &
              & *(1d0-cos(v2_delta_phi(kin(npar)%p(:,npar-3),kin(npar)%p(:,npar-2)))) )
  end function eval_mt_V2

  ! ET = sqrt(mll^2 + ptz^2)
  double precision function eval_et_V(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_et_V = sqrt( v1_minv2(pV) + v1_pt2(pV) )
  end function eval_et_V

  ! ET = sqrt(mll^2 + ptz^2)
  double precision function eval_et_V2(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2)
    end do
    eval_et_V2 = sqrt( v1_minv2(pV) + v1_pt2(pV) )
  end function eval_et_V2

  double precision function eval_pt_Vjj(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pVjj
    integer :: i
    do i=1,4
      pVjj(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar) &
            & + kin(npar)%pjets(i,1) + kin(npar)%pjets(i,2)
    end do
    eval_pt_Vjj = v1_pt(pVjj)
  end function eval_pt_Vjj

  double precision function eval_pt_Vj(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pVj
    integer :: i
    do i=1,4
      pVj(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar) &
            & + kin(npar)%pjets(i,1)
    end do
    eval_pt_Vj = v1_pt(pVj)
  end function eval_pt_Vj

  ! Added by KR
  double precision function eval_ystar_Vj(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_ystar_Vj = (v1_y(pV) - v1_y(kin(npar)%pjets(:,1))) / 2.d0 
  end function eval_ystar_Vj

  double precision function eval_yboost_Vj(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_yboost_Vj = (v1_y(pV) + v1_y(kin(npar)%pjets(:,1))) / 2.d0 
  end function eval_yboost_Vj
  ! End: Added by KR


  !----- jets: 1,2,...  <-->  leading, sub-leading, ...

  ! after recombination eta_j[i] in general not the same as y_j[i]
  double precision function eval_pt_j1(npar)
    integer, intent(in) :: npar
    eval_pt_j1 = v1_pt(kin(npar)%pjets(:,1))
  end function eval_pt_j1

  double precision function eval_y_j1(npar)
    integer, intent(in) :: npar
    eval_y_j1 = v1_y(kin(npar)%pjets(:,1))
  end function eval_y_j1

  double precision function eval_eta_j1(npar)
    integer, intent(in) :: npar
    eval_eta_j1 = v1_eta(kin(npar)%pjets(:,1))
  end function eval_eta_j1


  double precision function eval_pt_j2(npar)
    integer, intent(in) :: npar
    eval_pt_j2 = v1_pt(kin(npar)%pjets(:,2))
  end function eval_pt_j2

  double precision function eval_y_j2(npar)
    integer, intent(in) :: npar
    eval_y_j2 = v1_y(kin(npar)%pjets(:,2))
  end function eval_y_j2

  double precision function eval_eta_j2(npar)
    integer, intent(in) :: npar
    eval_eta_j2 = v1_eta(kin(npar)%pjets(:,2))
  end function eval_eta_j2


  double precision function eval_pt_j3(npar)
    integer, intent(in) :: npar
    eval_pt_j3 = v1_pt(kin(npar)%pjets(:,3))
  end function eval_pt_j3

  double precision function eval_y_j3(npar)
    integer, intent(in) :: npar
    eval_y_j3 = v1_y(kin(npar)%pjets(:,3))
  end function eval_y_j3

  double precision function eval_eta_j3(npar)
    integer, intent(in) :: npar
    eval_eta_j3 = v1_eta(kin(npar)%pjets(:,3))
  end function eval_eta_j3


  double precision function eval_pt_j4(npar)
    integer, intent(in) :: npar
    eval_pt_j4 = v1_pt(kin(npar)%pjets(:,4))
  end function eval_pt_j4

  double precision function eval_y_j4(npar)
    integer, intent(in) :: npar
    eval_y_j4 = v1_y(kin(npar)%pjets(:,4))
  end function eval_y_j4

  double precision function eval_eta_j4(npar)
    integer, intent(in) :: npar
    eval_eta_j4 = v1_eta(kin(npar)%pjets(:,4))
  end function eval_eta_j4

  double precision function eval_sum_pt_jets(npar)
    integer, intent(in) :: npar
    integer :: i
    double precision :: sum_pt
    sum_pt=0d0
    do i=0,kin(npar)%njets
       if(i.eq.0) cycle
       sum_pt=sum_pt+v1_pt(kin(npar)%pjets(:,i))
    enddo
    eval_sum_pt_jets = sum_pt
  end function eval_sum_pt_jets

  double precision function eval_deltapt_j12(npar)
    integer, intent(in) :: npar
    double precision :: ptj1, ptj2
    ptj1 = eval_obs(id_pt_j1, npar)
    ptj2 = eval_obs(id_pt_j2, npar)
    eval_deltapt_j12 = abs(ptj1 - ptj2)
  end function eval_deltapt_j12

  double precision function eval_deltaeta_j12(npar)
    integer, intent(in) :: npar
    double precision :: etaj1, etaj2
    etaj1 = eval_obs(id_eta_j1, npar)
    etaj2 = eval_obs(id_eta_j2, npar)
    eval_deltaeta_j12 = abs(etaj1 - etaj2)
  end function eval_deltaeta_j12

  double precision function eval_deltay_j12(npar)
    integer, intent(in) :: npar
    double precision :: yj1, yj2
    yj1 = eval_obs(id_y_j1, npar)
    yj2 = eval_obs(id_y_j2, npar)
    eval_deltay_j12 = abs(yj1 - yj2)
  end function eval_deltay_j12

  double precision function eval_deltaphi_j12(npar)
    integer, intent(in) :: npar
    eval_deltaphi_j12 = v2_delta_phi( kin(npar)%pjets(:,1),  &
    &                                 kin(npar)%pjets(:,2) )
  end function eval_deltaphi_j12

  double precision function eval_chi_j12(npar)
    integer, intent(in) :: npar
    double precision :: abs_del_y_j12
    abs_del_y_j12 = eval_obs(id_deltay_j12, npar)
    eval_chi_j12 = exp( abs_del_y_j12 )
  end function eval_chi_j12

  double precision function eval_ystar_j12(npar)
    integer, intent(in) :: npar
    double precision :: abs_del_y_j12
    abs_del_y_j12 = eval_obs(id_deltay_j12, npar)
    eval_ystar_j12 = abs_del_y_j12 / 2d0
  end function eval_ystar_j12

  double precision function eval_yboost_j12(npar)
    integer, intent(in) :: npar
    double precision :: yj1, yj2
    yj1 = eval_obs(id_y_j1, npar)
    yj2 = eval_obs(id_y_j2, npar)
    eval_yboost_j12 = abs(yj1 + yj2) / 2d0
  end function eval_yboost_j12

  double precision function eval_max_y_j12(npar)
    integer, intent(in) :: npar
    double precision :: yj1, yj2
    yj1 = eval_obs(id_y_j1, npar)
    yj2 = eval_obs(id_y_j2, npar)
    eval_max_y_j12 = max(yj1, yj2)
  end function eval_max_y_j12

  double precision function eval_minv_j12(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pj12
    pj12(1:4) = kin(npar)%pjets(1:4,1) + kin(npar)%pjets(1:4,2)
    eval_minv_j12 = v1_minv(pj12)
    !> this is slow because it induces an allocation
    !eval_minv_j12 = v1_minv( kin(npar)%pjets(:,1) + kin(npar)%pjets(:,2) )
  end function eval_minv_j12

  double precision function eval_pt_j1expYstar(npar)
    integer, intent(in) :: npar
    eval_pt_j1expYstar = v1_pt(kin(npar)%pjets(:,1))*exp(0.3d0*eval_ystar_j12(npar))
  end function eval_pt_j1expYstar

  double precision function eval_minv_j123(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pj123
    pj123(1:4) = kin(npar)%pjets(1:4,1)  &
    &          + kin(npar)%pjets(1:4,2)  &
    &          + kin(npar)%pjets(1:4,3)
    eval_minv_j123 = v1_minv(pj123)
  end function eval_minv_j123

  double precision function eval_ptavg_j12(npar)
    integer, intent(in) :: npar
    eval_ptavg_j12 = ( v1_pt(kin(npar)%pjets(:,1)) &
                   & + v1_pt(kin(npar)%pjets(:,2)) ) / 2d0
  end function eval_ptavg_j12

  double precision function eval_ptavg_j123(npar)
    integer, intent(in) :: npar
    eval_ptavg_j123 = ( v1_pt(kin(npar)%pjets(:,1)) &
                    & + v1_pt(kin(npar)%pjets(:,2)) &
                    & + v1_pt(kin(npar)%pjets(:,3)) ) / 3d0
  end function eval_ptavg_j123

  double precision function eval_ptavg_jall(npar)
    integer, intent(in) :: npar
    integer :: i
    eval_ptavg_jall = 0d0
    do i=1,kin(npar)%njets
      eval_ptavg_jall = eval_ptavg_jall + v1_pt(kin(npar)%pjets(:,i))
    end do
    eval_ptavg_jall = eval_ptavg_jall / DBLE(kin(npar)%njets)
  end function eval_ptavg_jall

  double precision function eval_ptavg_geom_jall(npar)
    integer, intent(in) :: npar
    integer :: i
    eval_ptavg_geom_jall = 1d0
    do i=1,kin(npar)%njets
      eval_ptavg_geom_jall = eval_ptavg_geom_jall * v1_pt(kin(npar)%pjets(:,i))
    end do
    eval_ptavg_geom_jall = eval_ptavg_geom_jall ** (1d0/DBLE(kin(npar)%njets))
  end function eval_ptavg_geom_jall

  double precision function eval_etavg_j12(npar)
    integer, intent(in) :: npar
    eval_etavg_j12 = ( v1_et(kin(npar)%pjets(:,1)) &
                   & + v1_et(kin(npar)%pjets(:,2)) ) / 2d0
  end function eval_etavg_j12

  double precision function eval_etavg_j123(npar)
    integer, intent(in) :: npar
    eval_etavg_j123 = ( v1_et(kin(npar)%pjets(:,1)) &
                    & + v1_et(kin(npar)%pjets(:,2)) &
                    & + v1_et(kin(npar)%pjets(:,3)) ) / 3d0
  end function eval_etavg_j123

  double precision function eval_etavg_jall(npar)
    integer, intent(in) :: npar
    integer :: i
    eval_etavg_jall = 0d0
    do i=1,kin(npar)%njets
      eval_etavg_jall = eval_etavg_jall + v1_et(kin(npar)%pjets(:,i))
    end do
    eval_etavg_jall = eval_etavg_jall / DBLE(kin(npar)%njets)
  end function eval_etavg_jall

  double precision function eval_etavg_geom_jall(npar)
    integer, intent(in) :: npar
    integer :: i
    eval_etavg_geom_jall = 1d0
    do i=1,kin(npar)%njets
      eval_etavg_geom_jall = eval_etavg_geom_jall * v1_et(kin(npar)%pjets(:,i))
    end do
    eval_etavg_geom_jall = eval_etavg_geom_jall ** (1d0/DBLE(kin(npar)%njets))
  end function eval_etavg_geom_jall

  double precision function eval_et2_j12(npar)
    integer, intent(in) :: npar
    eval_et2_j12 = ( v1_et(kin(npar)%pjets(:,1))**2 &
                   & + v1_et(kin(npar)%pjets(:,2))**2 )
  end function eval_et2_j12


  double precision function eval_photon_pt(npar)
    integer, intent(in) :: npar
    !> photon is always last entry in momentum set
    eval_photon_pt = v1_pt(kin(npar)%p(:,npar))
  end function eval_photon_pt

  double precision function eval_photon_y(npar)
    integer, intent(in) :: npar
    !> photon is always last entry in momentum set
    eval_photon_y = v1_y(kin(npar)%p(:,npar))
  end function eval_photon_y

  double precision function eval_tau_j1(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pj1
    double precision, dimension(4) :: ph
    pj1(1:4) = kin(npar)%pjets(1:4,1)  
    ph(1:4)  = kin(npar)%p(1:4,npar-1)+kin(npar)%p(1:4,npar) 
    eval_tau_j1 = sqrt(v1_pt2(pj1)+v1_minv2(pj1))/2/dcosh(v1_y(pj1)-v1_y(ph))
  end function eval_tau_j1

  double precision function eval_sum_tau_jets(npar)
    integer, intent(in) :: npar
    integer :: i
    double precision, dimension(4) :: pjet
    double precision, dimension(4) :: ph  
    ph(1:4)  = kin(npar)%p(1:4,npar-1)+kin(npar)%p(1:4,npar) 
    eval_sum_tau_jets = 0d0
    do i=1,kin(npar)%njets
       pjet(1:4) = kin(npar)%pjets(1:4,kin(npar)%njets)
       eval_sum_tau_jets = eval_sum_tau_jets &
         &   + sqrt(v1_pt2(pjet)+v1_minv2(pjet))/2/dcosh(v1_y(pjet)-v1_y(ph))
    end do
  end function eval_sum_tau_jets

  double precision function eval_max_tau_jet(npar)
    integer, intent(in) :: npar
    integer :: i
    double precision :: eval_tau_jet  
    double precision, dimension(4) :: pjet
    double precision, dimension(4) :: ph  
    ph(1:4)  = kin(npar)%p(1:4,npar-1)+kin(npar)%p(1:4,npar) 
    eval_max_tau_jet = 0d0
    eval_tau_jet = 0d0
    do i=1,kin(npar)%njets
       pjet(1:4) = kin(npar)%pjets(1:4,kin(npar)%njets)
       eval_tau_jet = sqrt(v1_pt2(pjet)+v1_minv2(pjet))/2/dcosh(v1_y(pjet)-v1_y(ph))
       if (eval_max_tau_jet.le.eval_tau_jet) then
           eval_max_tau_jet=eval_tau_jet
       endif
    end do
  end function eval_max_tau_jet


  !----- "perfect" jet observables
  subroutine set_jets_nocut(npar)
    integer, intent(in) :: npar
    double precision, parameter :: fixedRSq = 1d0**2
    integer :: i
    double precision :: pt_1,pt_2,pt_i, sum_pt, prod_pt

    !> initialise the recombination algorithm
    call initRecomb_jet(npar)
    !> cluster the jets
    call cluster_jet()
    !> sort in pt
    call sort_jet()

    !> if njets < 2 throw exception
    if (njets_jet < 2) then
      print*, "set_jets_nocut: less than two reconstructed jets: ", njets_jet
      stop
    end if

    !> compute sums and products
    sum_pt  = 0d0
    prod_pt = 1d0
    do i=1,njets_jet
      pt_i = sqrt( jets(pt_sort_jet(i))%pt2 )
      if (i==1) pt_1 = pt_i
      if (i==2) pt_2 = pt_i
      sum_pt  = sum_pt  + pt_i
      prod_pt = prod_pt * pt_i
    end do

    !> compute the obs & cache
    list_cache(id_pt_j1_nocut)%value             = pt_1
    list_cache(id_pt_j2_nocut)%value             = pt_2
    list_cache(id_ht_jall_nocut)%value           = sum_pt
    list_cache(id_ptavg_j12_nocut)%value         = (pt_1+pt_2)/2d0
    list_cache(id_ptavg_jall_nocut)%value        = sum_pt / DBLE(njets_jet)
    list_cache(id_ptavg_geom_j12_nocut)%value    = sqrt(pt_1*pt_2)
    list_cache(id_ptavg_geom_jall_nocut)%value   = prod_pt**(1d0/DBLE(njets_jet))

    list_cache(id_pt_j1_nocut)%qcached           = .true.
    list_cache(id_pt_j2_nocut)%qcached           = .true.
    list_cache(id_ht_jall_nocut)%qcached         = .true.
    list_cache(id_ptavg_j12_nocut)%qcached       = .true.
    list_cache(id_ptavg_jall_nocut)%qcached      = .true.
    list_cache(id_ptavg_geom_j12_nocut)%qcached  = .true.
    list_cache(id_ptavg_geom_jall_nocut)%qcached = .true.

  end subroutine set_jets_nocut



  !----- R=1 jet observables
  subroutine set_jets_R1(npar)
    integer, intent(in) :: npar
    double precision, parameter :: fixedRSq = 1d0**2
    integer :: i
    double precision :: pt_1,pt_2,pt_i, sum_pt, prod_pt

    !> initialise the recombination algorithm
    call initRecomb_jet(npar)
    !> cluster the jets
    call cluster_jet(RSq = fixedRSq)
    !> sort in pt
    call sort_jet()

    !> if njets < 2 throw exception
    if (njets_jet < 2) then
      print*, "set_jets_R1: less than two reconstructed jets: ", njets_jet
      stop
    end if

    !> compute sums and products
    sum_pt  = 0d0
    prod_pt = 1d0
    do i=1,njets_jet
      pt_i = sqrt( jets(pt_sort_jet(i))%pt2 )
      if (i==1) pt_1 = pt_i
      if (i==2) pt_2 = pt_i
      sum_pt  = sum_pt  + pt_i
      prod_pt = prod_pt * pt_i
    end do

    !> compute the obs & cache
    list_cache(id_pt_j1_R1)%value             = pt_1
    list_cache(id_pt_j2_R1)%value             = pt_2
    list_cache(id_ht_jall_R1)%value           = sum_pt
    list_cache(id_ptavg_j12_R1)%value         = (pt_1+pt_2)/2d0
    list_cache(id_ptavg_jall_R1)%value        = sum_pt / DBLE(njets_jet)
    list_cache(id_ptavg_geom_j12_R1)%value    = sqrt(pt_1*pt_2)
    list_cache(id_ptavg_geom_jall_R1)%value   = prod_pt**(1d0/DBLE(njets_jet))

    list_cache(id_pt_j1_R1)%qcached           = .true.
    list_cache(id_pt_j2_R1)%qcached           = .true.
    list_cache(id_ht_jall_R1)%qcached         = .true.
    list_cache(id_ptavg_j12_R1)%qcached       = .true.
    list_cache(id_ptavg_jall_R1)%qcached      = .true.
    list_cache(id_ptavg_geom_j12_R1)%qcached  = .true.
    list_cache(id_ptavg_geom_jall_R1)%qcached = .true.

  end subroutine set_jets_R1


  !> y_{n n+1} ycut transitions in the JADE algorithm
  !> this is *before* any jet cuts are applied!!!
  subroutine set_jade_ytrans(npar)
    integer, intent(in) :: npar

    !> initialise the recombination algorithm
    call initRecomb_jet(npar)
    !> run clustering till the end & compute y_{n n+1}
    call cluster_jade_jet(compute_ytrans = .true.)

  end subroutine set_jade_ytrans


  !----- DIS observables

  double precision function eval_dis_dscl(npar)
    integer, intent(in) :: npar
    double precision :: ptavg_j12, q2
    ptavg_j12 = eval_obs(id_ptavg_j12, npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dscl = sqrt((q2+ptavg_j12**2)/2d0)
  end function eval_dis_dscl

  double precision function eval_dis_dscl3j(npar)
    integer, intent(in) :: npar
    double precision :: ptavg_j123, q2
    ptavg_j123 = eval_obs(id_ptavg_j123, npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dscl3j = sqrt((q2+ptavg_j123**2)/2d0)
  end function eval_dis_dscl3j

  double precision function eval_dis_dsclj1(npar)
    integer, intent(in) :: npar
    double precision :: pt_jet, q2
    pt_jet = eval_obs(id_pt_j1,npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclj1 = sqrt((q2+pt_jet**2)/2d0)
  end function eval_dis_dsclj1

  double precision function eval_dis_dsclj2(npar)
    integer, intent(in) :: npar
    double precision :: pt_jet, q2
    pt_jet = eval_obs(id_pt_j2,npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclj2 = sqrt((q2+pt_jet**2)/2d0)
  end function eval_dis_dsclj2

  double precision function eval_dis_dsclj3(npar)
    integer, intent(in) :: npar
    double precision :: pt_jet, q2
    pt_jet = eval_obs(id_pt_j3,npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclj3 = sqrt((q2+pt_jet**2)/2d0)
  end function eval_dis_dsclj3

  double precision function eval_dis_dsclj4(npar)
    integer, intent(in) :: npar
    double precision :: pt_jet, q2
    pt_jet = eval_obs(id_pt_j4,npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclj4 = sqrt((q2+pt_jet**2)/2d0)
  end function eval_dis_dsclj4


  double precision function eval_dis_dsclj1l(npar)
    integer, intent(in) :: npar
    double precision :: pt_jet, q2
    pt_jet = eval_dis_hera_etj1(npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclj1l = sqrt((q2+pt_jet**2)/2d0)
  end function eval_dis_dsclj1l

  double precision function eval_dis_dsclj2l(npar)
    integer, intent(in) :: npar
    double precision :: pt_jet, q2
    pt_jet = eval_dis_hera_etj2(npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclj2l = sqrt((q2+pt_jet**2)/2d0)
  end function eval_dis_dsclj2l

  double precision function eval_dis_dsclj3l(npar)
    integer, intent(in) :: npar
    double precision :: pt_jet, q2
    pt_jet = eval_dis_hera_etj3(npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclj3l = sqrt((q2+pt_jet**2)/2d0)
  end function eval_dis_dsclj3l

  double precision function eval_dis_dsclj4l(npar)
    integer, intent(in) :: npar
    double precision :: pt_jet, q2
    pt_jet = eval_dis_hera_etj4(npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclj4l = sqrt((q2+pt_jet**2)/2d0)
  end function eval_dis_dsclj4l



  double precision function eval_dis_dsclZEUS(npar)
    integer, intent(in) :: npar
    double precision :: etavg_j12, q2
    etavg_j12 = eval_obs(id_etavg_j12,npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclZEUS = sqrt(q2+etavg_j12**2)
  end function eval_dis_dsclZEUS

  double precision function eval_dis_dsclZEUS2(npar)
    integer, intent(in) :: npar
    double precision :: etavg_jall, q2
    etavg_jall = eval_obs(id_etavg_jall,npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_dsclZEUS2 = sqrt(q2+etavg_jall**2)
  end function eval_dis_dsclZEUS2

  double precision function eval_dis_hera_deltaeta_j12()
    use DIS_mod
    double precision :: eta1,eta2
    eta1=v1_eta(dis%HERAjets(:,1))
    eta2=v1_eta(dis%HERAjets(:,2))
    eval_dis_hera_deltaeta_j12=abs(eta1-eta2)
  end function eval_dis_hera_deltaeta_j12

  double precision function eval_dis_hera_etaavg_j12()
    use DIS_mod, only : dis
    double precision :: eta1,eta2    
    eta1=v1_eta(dis%HERAjets(:,1))
    eta2=v1_eta(dis%HERAjets(:,2))
    eval_dis_hera_etaavg_j12=0.5d0*(eta1+eta2)
  end function eval_dis_hera_etaavg_j12

  double precision function eval_dis_hera_etaj1(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_etaj1=v1_eta(dis%HERAjets(1:4,1))
  end function eval_dis_hera_etaj1

  double precision function eval_dis_hera_etaj2(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_etaj2=v1_eta(dis%HERAjets(1:4,2))
  end function eval_dis_hera_etaj2

  double precision function eval_dis_hera_etaj3(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_etaj3=v1_eta(dis%HERAjets(1:4,3))
  end function eval_dis_hera_etaj3

  double precision function eval_dis_hera_etaj4(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_etaj4=v1_eta(dis%HERAjets(1:4,4))
  end function eval_dis_hera_etaj4

  double precision function eval_dis_hera_etj1(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_etj1=v1_et(dis%HERAjets(1:4,1))
  end function eval_dis_hera_etj1

  double precision function eval_dis_hera_etj2(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_etj2=v1_et(dis%HERAjets(1:4,2))
  end function eval_dis_hera_etj2

  double precision function eval_dis_hera_etj3(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_etj3=v1_et(dis%HERAjets(1:4,3))
  end function eval_dis_hera_etj3

  double precision function eval_dis_hera_etj4(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_etj4=v1_et(dis%HERAjets(1:4,4))
  end function eval_dis_hera_etj4

  double precision function eval_dis_hera_et2j1(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_et2j1=v1_et2(dis%HERAjets(1:4,1))
  end function eval_dis_hera_et2j1

  double precision function eval_dis_hera_et2j2(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_et2j2=v1_et2(dis%HERAjets(1:4,2))
  end function eval_dis_hera_et2j2

  double precision function eval_dis_hera_et2j3(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_et2j3=v1_et2(dis%HERAjets(1:4,3))
  end function eval_dis_hera_et2j3

  double precision function eval_dis_hera_et2j4(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_et2j4=v1_et2(dis%HERAjets(1:4,4))
  end function eval_dis_hera_et2j4

  double precision function eval_dis_hera_ET2byQ2j1(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    double precision q2
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_hera_ET2byQ2j1=v1_et2(dis%HERAjets(1:4,1))/q2
  end function eval_dis_hera_ET2byQ2j1

  double precision function eval_dis_hera_ET2byQ2j2(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    double precision q2
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_hera_ET2byQ2j2=v1_et2(dis%HERAjets(1:4,2))/q2
  end function eval_dis_hera_ET2byQ2j2

  double precision function eval_dis_hera_ET2byQ2j3(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    double precision q2
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_hera_ET2byQ2j3=v1_et2(dis%HERAjets(1:4,3))/q2
  end function eval_dis_hera_ET2byQ2j3

  double precision function eval_dis_hera_ET2byQ2j4(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    double precision q2
    q2  = eval_obs(id_DIS_Q2, npar)
    eval_dis_hera_ET2byQ2j4=v1_et2(dis%HERAjets(1:4,4))/q2
  end function eval_dis_hera_ET2byQ2j4

  double precision function eval_dis_hera_eta_j1_0()
    use DIS_mod, only : dis
    eval_dis_hera_eta_j1_0=v1_eta(dis%raw_HERAjets(:,1))
  end function eval_dis_hera_eta_j1_0

  double precision function eval_dis_hera_xgamma(npar)
  use DIS_mod
  integer, intent(in) :: npar
  double precision :: Etj1,Etj2,Pzj1,Pzj2,y,Ee  
  Etj1=v1_et(dis%HERAjets(:,1))
  Etj2=v1_et(dis%HERAjets(:,2))
  Pzj1=dis%HERAjets(3,1)
  Pzj2=dis%HERAjets(3,2)
  y=eval_obs(id_DIS_y,npar)
  Ee=getEe()
  eval_dis_hera_xgamma=(Etj1-Pzj1+Etj2-Pzj2)/(2d0*Ee*y)  
  end function eval_dis_hera_xgamma
  
  double precision function eval_dis_hera_eta_j2_0()
    use DIS_mod, only : dis
    eval_dis_hera_eta_j2_0=v1_eta(dis%raw_HERAjets(:,2))
  end function eval_dis_hera_eta_j2_0

  double precision function eval_gammap_eta_j1()
    use DIS_mod, only : dis
    eval_gammap_eta_j1=v1_eta(dis%HERApgammajets(:,1))
  end function eval_gammap_eta_j1

  double precision function eval_gammap_eta_j2()
    use DIS_mod, only : dis
    eval_gammap_eta_j2=v1_eta(dis%HERApgammajets(:,2))
  end function eval_gammap_eta_j2

  double precision function eval_DIS_hera_FJ_reweight(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    double precision :: et2sum,q2
    integer :: i
    et2sum=0d0
    q2  = eval_obs(id_DIS_Q2, npar)
    if(kin(npar)%njets.ge.1) then
       do i=1,kin(npar)%njets
          et2sum=et2sum+v1_et(dis%HERAjets(1:4,i))**2
       enddo
    endif
    eval_DIS_hera_FJ_reweight=q2+et2sum
  end function eval_DIS_hera_FJ_reweight

  double precision function eval_dis_hera_deta1(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_deta1=v1_eta(dis%HERAplus2jets(1:4,1))&
         & -v1_eta(dis%HERAplus2jets(1:4,2))
  end function eval_dis_hera_deta1

  double precision function eval_dis_hera_deta2(npar)
    use DIS_mod, only : dis
    integer, intent(in) :: npar
    eval_dis_hera_deta2=v1_eta(dis%HERAjets(1:4,1))&
         & -v1_eta(dis%HERAplus2jets(1:4,1))
  end function eval_dis_hera_deta2

  double precision function eval_dis_delta_phis()
    use DIS_mod, only : dis
    use Constants_mod
    double precision :: phi1,phi2,phi
    phi1 = v1_phi(dis%HERApgammajets(:,1))
    phi2 = v1_phi(dis%HERApgammajets(:,2))
    phi = abs(phi1-phi2)
    if(phi.gt.pi) phi = 2d0*pi-phi
    eval_dis_delta_phis= phi
  end function eval_dis_delta_phis
  
  double precision function eval_dis_xi2(npar)
    integer, intent(in) :: npar
    double precision :: m12, q2, x
    m12 = eval_obs(id_minv_j12, npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    x   = eval_obs(id_DIS_x, npar)
    eval_dis_xi2 = x*(1d0+m12**2/q2)
  end function eval_dis_xi2

  double precision function eval_dis_thrust(npar)
     integer, intent(in) :: npar
     double precision :: top,bot,absp
     integer :: i,ncount
    
     top = 0d0
     bot = 0d0
     absp = 0d0
     ncount = 0
     do i=4, npar
         ! Select only particles in the current hemisphere (pz > 0)        
         if (kin(npar)%p(3,i) > 0d0) then
             absp = sqrt(kin(npar)%p(1,i)**2 + &
                     & kin(npar)%p(2,i)**2 + kin(npar)%p(3,i)**2 )
             top = top + sqrt(kin(npar)%p(3,i)**2)
             bot = bot + absp
             ncount=ncount+1
         end if
     end do
     if(ncount.eq.0) then
        eval_dis_thrust = -1d0
     else
        eval_dis_thrust = 1d0 - top / bot
     endif
  end function eval_dis_thrust
 
  double precision function eval_dis_thrust_c(npar)
     integer, intent(in) :: npar
     double precision :: top,bot
     double precision, dimension(3) :: dis_Taxis
     integer :: i,j,ncount
    
     top = 0d0
     bot = 0d0
     dis_Taxis = (/0d0, 0d0, 0d0/)
     ncount=0
     ! thrust axis is just in the direction of the sum of all particles momenta    
     do j=4,npar
             ! Select only particles in the current hemisphere (pz > 0)    
        if (kin(npar)%p(3,j) > 0d0) then
           ncount=ncount+1
           do i=1,3
              dis_Taxis(i) = dis_Taxis(i) + kin(npar)%p(i,j)
           enddo
        end if
     enddo

     if(ncount.eq.0) then
        eval_dis_thrust_c = -1d0
        return
     endif

     dis_Taxis = dis_Taxis / sqrt(dis_Taxis(1)**2+dis_Taxis(2)**2+dis_Taxis(3)**2)
    
     do i=4,npar
         if (kin(npar)%p(3,i) > 0d0) then
             top = top + sqrt(dot_eucl(kin(npar)%p(1:3,i),dis_Taxis)**2)
             bot = bot + sqrt(kin(npar)%p(1,i)**2+kin(npar)%p(2,i)**2+kin(npar)%p(3,i)**2)
         end if
     end do
    
     eval_dis_thrust_c = 1d0 - top / bot
  end function eval_dis_thrust_c
 
  double precision function eval_dis_JB(npar)
     integer, intent(in) :: npar
     double precision :: top,bot,absp
     integer :: i,ncount
    
     top = 0d0
     bot = 0d0
     absp = 0d0
     ncount = 0
     do i=4, npar
         ! Select only particles in the current hemisphere (pz > 0)        
         if (kin(npar)%p(3,i) > 0d0) then
            ncount=ncount+1
             absp = sqrt(kin(npar)%p(1,i)**2 + &
                     & kin(npar)%p(2,i)**2 + kin(npar)%p(3,i)**2 )
             top = top + sqrt(kin(npar)%p(1,i)**2 + kin(npar)%p(2,i)**2)
             bot = bot + absp
         end if
     end do
    if(ncount.eq.0) then
       eval_dis_JB = -1d0
    else
       eval_dis_JB = top / (2d0*bot)
    endif
  end function eval_dis_JB
 
  double precision function eval_dis_JM2(npar)
     integer, intent(in) :: npar
     double precision :: absptot,absp!,ehtot
     double precision :: pxtot,pytot,pztot
     integer :: i,ncount
 
     pxtot = 0d0
     pytot = 0d0
     pztot = 0d0
     !ehtot = 0d0 !rho_0: replace eh with modulus of 3-momentum
     absp = 0d0
     absptot = 0d0
     ncount=0

     do i=4, npar
         ! Select only particles in the current hemisphere (pz > 0)        
         if (kin(npar)%p(3,i) > 0d0) then
            ncount=ncount+1
            pxtot = pxtot + kin(npar)%p(1,i)
            pytot = pytot + kin(npar)%p(2,i)
            pztot = pztot + kin(npar)%p(3,i)
            absp = sqrt(kin(npar)%p(1,i)**2 + &
                 & kin(npar)%p(2,i)**2 + kin(npar)%p(3,i)**2 )
            !ehtot = ehtot + kin(npar)%p(4,i)
            absptot = absptot + absp
         end if
     end do    
     if(ncount.eq.0) then
        eval_dis_JM2 = -1d0
     else
        eval_dis_JM2 = (absptot**2 - pxtot**2 - pytot**2 - pztot**2) &
             &/ (2d0 * absptot)**2
     endif
  end function eval_dis_JM2
 
  double precision function eval_dis_C(npar)
    integer, intent(in) :: npar
    double precision :: t11,t12,t13,t22,t23,t33
    double precision :: Cpar
    t11=thetadis(1,1,npar)
    t12=thetadis(1,2,npar)
    t13=thetadis(1,3,npar)
    t22=thetadis(2,2,npar)
    t23=thetadis(2,3,npar)
    t33=thetadis(3,3,npar)
    Cpar=-t12**2-t13**2-t23**2+t11*t22+t11*t33+t22*t33
    eval_dis_C = 3d0*Cpar
  end function eval_dis_C
 
  double precision function eval_dis_Lgxi2(npar)
    integer, intent(in) :: npar
    double precision :: m12, q2, x
    m12 = eval_obs(id_minv_j12, npar)
    q2  = eval_obs(id_DIS_Q2, npar)
    x   = eval_obs(id_DIS_x, npar)
    eval_dis_Lgxi2 = log10(x*(1d0+m12**2/q2))
  end function eval_dis_Lgxi2

  double precision function eval_dis_etas(npar)
    integer, intent(in) :: npar
    double precision :: deta    
    deta=eval_obs(id_deltaeta_j12,npar)
    eval_dis_etas=abs(deta)/2d0
  end function eval_dis_etas

  double precision function eval_dis_visbyW(npar)
    integer, intent(in) :: npar
    double precision :: Esum,W
    integer i
    W    = eval_obs(id_dis_W,npar)
    Esum = 0d0
    do i=1,kin(npar)%njets
       Esum=Esum+kin(npar)%p(4,i)
    enddo
    eval_dis_visbyW=Esum/W
  end function eval_dis_visbyW

  double precision function eval_dis_xi3(npar)
    integer, intent(in) :: npar
    double precision :: m123, q2, x
    m123 = eval_obs(id_minv_j123, npar)
    q2   = eval_obs(id_DIS_Q2, npar)
    x    = eval_obs(id_DIS_x, npar)
    eval_dis_xi3 = x*(1d0+m123**2/q2)
  end function eval_dis_xi3

  double precision function eval_dis_ptl(npar)
    integer, intent(in) :: npar
    eval_dis_ptl = v1_pt(kin(npar)%p(:,3))
  end function eval_dis_ptl


  !----- more complicated guys

  ! scalar pt sum over reconstructed jets
  double precision function eval_ht_jets(npar)
    integer, intent(in) :: npar
    integer :: i
    eval_ht_jets = 0d0
    do i=1,kin(npar)%njets
      eval_ht_jets = eval_ht_jets + v1_pt(kin(npar)%pjets(:,i))
    end do
  end function eval_ht_jets

  ! scalar pt sum over the partons
  double precision function eval_ht_part(npar)
    integer, intent(in) :: npar
    integer :: i
    eval_ht_part = 0d0
    do i=kin(npar)%ijets_lower,kin(npar)%ijets_upper
      eval_ht_part = eval_ht_part + v1_pt(kin(npar)%p(:,i))
    end do
  end function eval_ht_part

  ! scalar pt sum over leptons and reconstructed jets
  double precision function eval_ht_full(npar)
    integer, intent(in) :: npar
    integer :: i
    eval_ht_full = v1_pt(kin(npar)%p(:,npar-1)) + v1_pt(kin(npar)%p(:,npar))
    do i=1,kin(npar)%njets
      eval_ht_full = eval_ht_full + v1_pt(kin(npar)%pjets(:,i))
    end do
  end function eval_ht_full


  ! smallest angular separation between reconstructed jets and l1
  double precision function eval_min_dR_l1j(npar)
    integer, intent(in) :: npar
    integer :: i
    double precision :: dR_l1i
    eval_min_dR_l1j = -1d0
    do i=1,kin(npar)%njets
      dR_l1i = v2_delta_R( kin(npar)%p(:,npar-1) , kin(npar)%pjets(:,i) )
      if ( (dR_l1i < eval_min_dR_l1j) .or. (eval_min_dR_l1j < 0d0) ) then
        eval_min_dR_l1j = dR_l1i
      end if
    end do
  end function eval_min_dR_l1j

  ! smallest angular separation between reconstructed jets and l2
  double precision function eval_min_dR_l2j(npar)
    integer, intent(in) :: npar
    integer :: i
    double precision :: dR_l2i
    eval_min_dR_l2j = -1d0
    do i=1,kin(npar)%njets
      dR_l2i = v2_delta_R( kin(npar)%p(:,npar) , kin(npar)%pjets(:,i) )
      if ( (dR_l2i < eval_min_dR_l2j) .or. (eval_min_dR_l2j < 0d0) ) then
        eval_min_dR_l2j = dR_l2i
      end if
    end do
  end function eval_min_dR_l2j

  ! smallest angular separation between reconstructed jets and l3
  double precision function eval_min_dR_l3j(npar)
    integer, intent(in) :: npar
    integer :: i
    double precision :: dR_l3i
    eval_min_dR_l3j = -1d0
    do i=1,kin(npar)%njets
      dR_l3i = v2_delta_R( kin(npar)%p(:,npar-3) , kin(npar)%pjets(:,i) )
      if ( (dR_l3i < eval_min_dR_l3j) .or. (eval_min_dR_l3j < 0d0) ) then
        eval_min_dR_l3j = dR_l3i
      end if
    end do
  end function eval_min_dR_l3j

  ! smallest angular separation between reconstructed jets and l4
  double precision function eval_min_dR_l4j(npar)
    integer, intent(in) :: npar
    integer :: i
    double precision :: dR_l4i
    eval_min_dR_l4j = -1d0
    do i=1,kin(npar)%njets
      dR_l4i = v2_delta_R( kin(npar)%p(:,npar-2) , kin(npar)%pjets(:,i) )
      if ( (dR_l4i < eval_min_dR_l4j) .or. (eval_min_dR_l4j < 0d0) ) then
        eval_min_dR_l4j = dR_l4i
      end if
    end do
  end function eval_min_dR_l4j

  ! smallest angular separation between reconstructed jets and l1/l2
  double precision function eval_min_dR_l12j(npar)
    integer, intent(in) :: npar
    eval_min_dR_l12j = min( eval_min_dR_l1j(npar) , eval_min_dR_l2j(npar) )
  end function eval_min_dR_l12j

  ! smallest angular separation between reconstructed jets and l3/l4
  double precision function eval_min_dR_l34j(npar)
    integer, intent(in) :: npar
    eval_min_dR_l34j = min( eval_min_dR_l3j(npar) , eval_min_dR_l4j(npar) )
  end function eval_min_dR_l34j

  ! smallest angular separation between reconstructed jets and l1/l2/l3/l4
  double precision function eval_min_dR_l1234j(npar)
    integer, intent(in) :: npar
    eval_min_dR_l1234j = min( eval_min_dR_l1j(npar) , eval_min_dR_l2j(npar) &
                           &, eval_min_dR_l3j(npar) , eval_min_dR_l4j(npar) )
  end function eval_min_dR_l1234j

  ! angular separation between l1 & l2
  double precision function eval_dR_l12(npar)
    integer, intent(in) :: npar
    eval_dR_l12 = v2_delta_R( kin(npar)%p(:,npar) , kin(npar)%p(:,npar-1) )
  end function eval_dR_l12

  ! angular separation between l3 & l4
  double precision function eval_dR_l34(npar)
    integer, intent(in) :: npar
    eval_dR_l34 = v2_delta_R( kin(npar)%p(:,npar-2) , kin(npar)%p(:,npar-3) )
  end function eval_dR_l34

  ! smallest angular separation between same-flavour leptons (l1,l2), (l3,l4)
  double precision function eval_min_dR_l_sf(npar)
    integer, intent(in) :: npar
    double precision :: dR_l12, dR_l34
    dR_l12 = eval_obs(id_dR_l12, npar)
    dR_l34 = eval_obs(id_dR_l34, npar)
    eval_min_dR_l_sf = min(dR_l12, dR_l34)
  end function eval_min_dR_l_sf

  ! smallest angular separation between different-flavour leptons (l1,l3), (l1,l4), (l2,l3), (l2,l4)
  double precision function eval_min_dR_l_df(npar)
    integer, intent(in) :: npar
    integer :: i,j
    double precision :: dR_lij
    eval_min_dR_l_df = -1d0
    do i=0,1
    do j=2,3
      dR_lij = v2_delta_R( kin(npar)%p(:,npar-i) , kin(npar)%p(:,npar-j) )
      if ( (dR_lij < eval_min_dR_l_df) .or. (eval_min_dR_l_df < 0d0) ) then
        eval_min_dR_l_df = dR_lij
      end if
    end do
    end do
  end function eval_min_dR_l_df

  ! angular separation between j1 and the i_th lepton ordered by their pT
  double precision function eval_dR_lj1(npar,i)
    integer, intent(in) :: npar
    integer :: j,i
    integer, dimension(4) :: order
    double precision, dimension(4) :: input_pt_4l,order_pt_4l,dR_lj1

       input_pt_4l(1) = eval_obs(id_pt_l1, npar)
       input_pt_4l(2) = eval_obs(id_pt_l2, npar)
       input_pt_4l(3) = eval_obs(id_pt_l3, npar)
       input_pt_4l(4) = eval_obs(id_pt_l4, npar)

       dR_lj1(1) = v2_delta_R( kin(npar)%p(:,npar-1) , kin(npar)%pjets(:,1) )
       dR_lj1(2) = v2_delta_R( kin(npar)%p(:,npar) , kin(npar)%pjets(:,1) )
       dR_lj1(3) = v2_delta_R( kin(npar)%p(:,npar-3) , kin(npar)%pjets(:,1) )
       dR_lj1(4) = v2_delta_R( kin(npar)%p(:,npar-2) , kin(npar)%pjets(:,1) )

       order(1)=0
       order(2)=0
       order(3)=0
       order(4)=0

       order_pt_4l(1)=0d0
       order_pt_4l(2)=0d0
       order_pt_4l(3)=0d0
       order_pt_4l(4)=0d0


    do j=1,4
       if (input_pt_4l(j) > order_pt_4l(1)) then
          order(1)=j
          order_pt_4l(1)=input_pt_4l(j)
       endif
    enddo
 
    do j=1,4
       if (j.ne.order(1)) then
          if (input_pt_4l(j) > order_pt_4l(2)) then
             order_pt_4l(2)=input_pt_4l(j)
             order(2)=j
          endif
       endif
    enddo

    do j=1,4
       if ((j.ne.order(1)).and.(j.ne.order(2))) then
          if (input_pt_4l(j) > order_pt_4l(3)) then
             order_pt_4l(3)=input_pt_4l(j)
             order(3)=j
          endif
       endif
    enddo

    do j=1,4
       if ((j.ne.order(1)).and.(j.ne.order(2)).and.(j.ne.order(3))) then
             order_pt_4l(4)=input_pt_4l(j)
             order(4)=j
       endif
    enddo

    eval_dR_lj1=dR_lj1(order(i))

  end function eval_dR_lj1

  ! \phi^*_\eta
  ! defined in Sect. 3 of arXiv:1009.1580 [hep-ex]
  double precision function eval_phi_star_eta(npar)
    use Constants_mod
    integer, intent(in) :: npar
    double precision :: del_phi, del_eta, cos_th,sin_th
    del_phi = v2_delta_phi( kin(npar)%p(:,npar-1), kin(npar)%p(:,npar) )
    del_eta = v2_delta_y  ( kin(npar)%p(:,npar-1), kin(npar)%p(:,npar) ) ! massless: eta = y
    cos_th = tanh(del_eta/2d0)
    sin_th = sqrt( (1d0+cos_th)*(1d0-cos_th) )
    eval_phi_star_eta = tan( (pi-del_phi)/2d0 ) * sin_th
  end function eval_phi_star_eta

  ! dynamical scale for ZJ production
  ! sqrt(mll^2 + sum ptj{i}^2) [arXiv:1512.01291]
  double precision function eval_ZJ_dscale1(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_ZJ_dscale1 = v1_minv2(pV)
    do i=1,kin(npar)%njets
      eval_ZJ_dscale1 = eval_ZJ_dscale1 + v1_pt2(kin(npar)%pjets(:,i))
    end do
    eval_ZJ_dscale1 = sqrt(eval_ZJ_dscale1)
  end function eval_ZJ_dscale1

  ! dynamical scale for ZJ production
  ! sqrt(mV^2 + ptV^2) + sum pt{partons}
  double precision function eval_ZJ_HT(npar)
    integer, intent(in) :: npar
    integer :: i
    eval_ZJ_HT = eval_obs(id_et_V, npar)
    do i=kin(npar)%ijets_lower,kin(npar)%ijets_upper
      eval_ZJ_HT = eval_ZJ_HT + v1_pt(kin(npar)%p(:,i))
    end do
  end function eval_ZJ_HT

  ! dynamical scale for ZJ production & shape uncertainties
  ! eps * mV + 1/eps * sum pt{V,partons}
  double precision function eval_ZJ_HT_shape(npar, eps)
    integer, intent(in) :: npar
    double precision, intent(in) :: eps
    double precision :: sum_pt
    integer :: i
    sum_pt = eval_obs(id_pt_V, npar)
    do i=kin(npar)%ijets_lower,kin(npar)%ijets_upper
      sum_pt = sum_pt + v1_pt(kin(npar)%p(:,i))
    end do
    eval_ZJ_HT_shape = eps * eval_obs(id_minv_V, npar) + sum_pt / eps 
  end function eval_ZJ_HT_shape


  ! dynamical scale for ZJ production
  ! sqrt(m4l^2 + sum ptj{i}^2) [arXiv:1512.01291]
  double precision function eval_H4lJ_dscale1(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) =   kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2) &
            & + kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_H4lJ_dscale1 = v1_minv2(pV)
    do i=1,kin(npar)%njets
      eval_H4lJ_dscale1 = eval_H4lJ_dscale1 + v1_pt2(kin(npar)%pjets(:,i))
    end do
    eval_H4lJ_dscale1 = sqrt(eval_H4lJ_dscale1)
  end function eval_H4lJ_dscale1

  ! dynamical scale for inclusive Z production
  ! sqrt(m4l^2 + ptz^2)
  double precision function eval_H4lJ_dscale2(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) =   kin(npar)%p(i,npar-3) + kin(npar)%p(i,npar-2) &
            & + kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_H4lJ_dscale2 = sqrt( v1_minv2(pV) + v1_pt2(pV) )
  end function eval_H4lJ_dscale2

  ! dynamical scale for Z(H)J production
  ! sqrt(mll^2 + ptZ^2) + sum ptj{i} [for Les Houches 2017]
  double precision function eval_ZJ_dscale3(npar)
    integer, intent(in) :: npar
    double precision, dimension(4) :: pV
    integer :: i
    do i=1,4
      pV(i) = kin(npar)%p(i,npar-1) + kin(npar)%p(i,npar)
    end do
    eval_ZJ_dscale3 = sqrt( v1_minv2(pV) + v1_pt2(pV) )
    do i=1,kin(npar)%njets
      eval_ZJ_dscale3 = eval_ZJ_dscale3 + v1_pt(kin(npar)%pjets(:,i))
    end do
  end function eval_ZJ_dscale3

  ! VFH selectors for VBF cuts
  double precision function eval_VFH_deltay(npar)
    integer, intent(in) :: npar
    eval_VFH_deltay = abs(v2_delta_y( kin(npar)%pjets(:,1), kin(npar)%pjets(:,2) ))
  end function eval_VFH_deltay

  double precision function eval_VFH_y1xy2(npar)
    integer, intent(in) :: npar
    eval_VFH_y1xy2 = v1_y( kin(npar)%pjets(:,1) ) * v1_y( kin(npar)%pjets(:,2) )
  end function eval_VFH_y1xy2

  ! dynamical scale for VFH production
  ! mu0^2 = mH/2 * sqrt(mH^2/4 + ptH^2)
  double precision function eval_VFH_dscale(npar)
    integer, intent(in)             :: npar
    double precision                :: mHo2, ptH2
    double precision, dimension (4) :: pV
    pV   = kin(npar)%p(:,npar-1) + kin(npar)%p(:,npar)
    mHo2 = v1_minv(pV)/2d0
    ptH2 = v1_pt2(pV)
    eval_VFH_dscale = dsqrt(mHo2*dsqrt(mHo2**2 + ptH2))
  end function eval_VFH_dscale

  ! normalised central rapidity
  double precision function eval_VFH_z3(npar)
     integer, intent(in) :: npar
     double precision :: y1, y2, y3
     y1 = v1_y(kin(npar)%pjets(:,1))
     y2 = v1_y(kin(npar)%pjets(:,2))
     y3 = v1_y(kin(npar)%pjets(:,3))
     eval_VFH_z3 = (y3 - (y1+y2)/2d0)/(y1-y2)
  end function eval_VFH_z3

  !>----- Spike plots
  double precision function eval_sij(npar, i, j)
    integer, intent(in) :: npar, i, j
    integer             :: n
    double precision    :: shat
    !> get the highest momentum set by traversing npar_from
    n = kin(npar)%npar_from
    if (n /= kin(n)%npar_from) n = kin(npar)%npar_from
    shat = dabs(kin(n)%s(1,2))
    eval_sij = dabs(kin(n)%s(i,j)) / shat
  end function eval_sij


  !--------------------------------!
  !  some useful helper functions  !
  !--------------------------------!


  pure double precision function v1_pt2(v)
    double precision, intent(in), dimension(4) :: v
    v1_pt2 = v(1)**2 + v(2)**2
  end function v1_pt2

  pure double precision function v1_pt(v)
    double precision, intent(in), dimension(4) :: v
    v1_pt = sqrt(v1_pt2(v))
  end function v1_pt

  pure double precision function v1_p2(v)
    double precision, intent(in), dimension(4) :: v
    v1_p2 = v(1)**2 + v(2)**2 + v(3)**2
  end function v1_p2

  pure double precision function v1_p(v)
    double precision, intent(in), dimension(4) :: v
    v1_p = sqrt(v1_p2(v))
  end function v1_p

  pure double precision function v1_minv2(v)
    use FPA_mod
    double precision, intent(in), dimension(4) :: v
    double precision :: err, acc_i,err_i
    integer :: i
    v1_minv2 = v(4)**2
    err      = 0d0
    do i=1,3
      call twoSum(v1_minv2,-v(i)**2, acc_i,err_i)
      v1_minv2 = acc_i
      err      = err + err_i
    end do
    v1_minv2 = v1_minv2 + err
    !v1_minv2 = v(4)**2 - v(1)**2 - v(2)**2 - v(3)**2
  end function v1_minv2

  pure double precision function v1_minv(v)
    double precision, intent(in), dimension(4) :: v
    v1_minv = sqrt(v1_minv2(v))
  end function v1_minv

  pure double precision function v1_et2(v)
    double precision, intent(in), dimension(4) :: v
    !> MCFM seems to define ET = E * sqrt(px2+py2) / sqrt(px2+py2+pz2)
    v1_et2 = v1_pt2(v) + v1_minv2(v)
  end function v1_et2

  pure double precision function v1_et(v)
    double precision, intent(in), dimension(4) :: v
    v1_et = sqrt(v1_et2(v))
  end function v1_et

  pure double precision function v1_y(v)
    double precision, intent(in), dimension(4) :: v
    v1_y = log((v(4)+v(3))/(v(4)-v(3)))/2d0
  end function v1_y

  pure double precision function v1_eta(v)
    double precision, intent(in), dimension(4) :: v
    double precision ::p
    p = v1_p(v)
    v1_eta = log((p+v(3))/(p-v(3)))/2d0
  end function v1_eta

  pure double precision function v1_phi(v) ! in [0,2pi]
    use Constants_mod
    double precision, intent(in), dimension(4) :: v
    v1_phi = atan2(v(2),v(1))
    if (v1_phi.lt.0d0) v1_phi = v1_phi+2d0*pi
  end function v1_phi

  pure double precision function v2_delta_y(v1,v2)
    double precision, intent(in), dimension(4) :: v1,v2
    v2_delta_y = v1_y(v1) - v1_y(v2)
  end function v2_delta_y

  pure double precision function v2_delta_eta(v1,v2)
    double precision, intent(in), dimension(4) :: v1,v2
    v2_delta_eta = v1_eta(v1) - v1_eta(v2)
  end function v2_delta_eta

  pure double precision function v2_delta_phi(v1,v2) ! in [0,pi]
    use Constants_mod
    double precision, intent(in), dimension(4) :: v1,v2
    v2_delta_phi = abs(v1_phi(v1) - v1_phi(v2))
    if(v2_delta_phi.gt.pi) v2_delta_phi = 2d0*pi-v2_delta_phi
  end function v2_delta_phi

  pure double precision function v2_delta_R(v1,v2)
    double precision, intent(in), dimension(4) :: v1,v2
    v2_delta_R = sqrt( v2_delta_phi(v1,v2)**2 + v2_delta_y(v1,v2)**2 )
  end function v2_delta_R

end module EvalObs_mod


!-------------------------------------------------------------------------------
module Observables_mod
!-------------------------------------------------------------------------------
  use EvalObs_mod
  use KinData_mod
  use Spikes_mod, only: initialiseSpike
  implicit none
  private


  !---------------!
  !  Observables  !
  !---------------!

  public :: init_obs, destroy_obs
  public :: digestActive_obs
  public :: initBuffer_obs, destroyBuffer_obs
  public :: nUsed_obs, setUsed_obs, purgeUsed_obs, isUsed_obs
  public :: getNameFromId_obs, getIdFromName_obs
  public :: evalAll_obs, getValue_obs, getValid_obs
  public :: printAll_obs, compCache_obs
  public :: recombine_jet
  public :: npar_max, npar_current


  type Observable_t
    logical :: qactive = .false.
    logical :: qused   = .false.
    character(len=max_obs_name_length) :: name
    character(len=max_obs_desc_length) :: desc
  end type Observable_t

  type(Observable_t), dimension(n_obs)        :: list_obs
  integer,            dimension(n_obs)        :: mapAll2Used_obs
  integer                                     :: nUsed_obs = 0
  integer,          allocatable, dimension(:) :: mapUsed2All_obs
  double precision, allocatable, dimension(:) :: valueUsed_obs
  logical,          allocatable, dimension(:) :: validUsed_obs
!$OMP THREADPRIVATE(valueUsed_obs,validUsed_obs)


  !> store some information on the kinematics:
  !>  * npar_max: highest momentum set in the current calcualtion
  integer :: npar_max
  integer :: npar_current
  !>  * npar_current: currently active npar, can be different for the threads!
  !>    (i.e. associated with the cached obs)
!$OMP THREADPRIVATE(npar_current)


contains


  !-----------------------------------------------------------------------------
  !> @brief
  !> Initialization routine for Observables. Register all allowed Observables
  !> with their names and associate it to the corresponding evaluation
  !> function. Position in list_obs corresponds to the unique identifyer of the
  !> observable.
  !
  !> @param[in] process the process under consideration
  !-----------------------------------------------------------------------------
  subroutine init_obs(process)
    character(len=*), intent(in) :: process
    !> note: depending on the process we can select the subset of allowed
    !> observables and change their respective string identifyer

    !> init to invalid states
    npar_max     = -1
    npar_current = -1

    !------------------------
    ! set boson observables !
    !------------------------
    !> `npar-1` : anti-lepton
    !> `npar`   : lepton
    !> @todo check that the lepton assignments are correct


    select case (process)
      case ("Z","ZJ","ZJJ","ZJJJ")
        ! leptons from the Z -> lp lm decay
        call set_active(id_pt_l1,       "ptlm", &
                       & "transverse momentum of l-")
        call set_active(id_y_l1,        "ylm", &
                       & "(pseudo-)rapidity of l-")
        call set_active(id_abs_y_l1,    "abs_ylm", &
                       & "|ylm|")
        call set_active(id_pt_l2,       "ptlp", &
                       & "transverse momentum of l+")
        call set_active(id_y_l2,        "ylp", &
                       & "(pseudo-)rapidity of l+")
        call set_active(id_abs_y_l2,    "abs_ylp", &
                       & "|ylp|")
        ! leading & sub-leading lepton
        call set_active(id_pt_g1,       "ptl1", &
                       & "transverse momentum of leading lepton")
        call set_active(id_y_g1,        "yl1", &
                       & "(pseudo-)rapidity of leading lepton")
        call set_active(id_abs_y_g1,    "abs_yl1", &
                       & "|yl1|")
        call set_active(id_pt_g2,       "ptl2", &
                       & "transverse momentum of sub-leading lepton")
        call set_active(id_y_g2,        "yl2", &
                       & "(pseudo-)rapidity of sub-leading lepton")
        call set_active(id_abs_y_g2,    "abs_yl2", &
                       & "|yl2|")
        call set_active(id_deltaphi_g1g2,     "dphi_l1l2", &
                       & "Azimuthal angle difference of di-lepton")
        ! Z-boson = lp + lm
        call set_active(id_pt_V,        "ptz", &
                       & "transverse momentum of Z")
        call set_active(id_minv_V,      "mll", &
                       & "invariant mass of Z")
        call set_active(id_et_V,     "etz", &
                       & "transverse energy: sqrt(mll^2 + ptz^2)")
        call set_active(id_y_V,         "yz", &
                       & "rapidity of Z")
        call set_active(id_abs_y_V,     "abs_yz", &
                       & "|yz|")
        call set_active(id_eta_V,       "etaz", &
                       & "pseudo-rapidity of Z")
        call set_active(id_abs_eta_V,   "abs_etaz", &
                       & "|etaz|")
        call set_active(id_dR_l12, "dr_l12", &
                       & "Delta R (lp, lm)")
        ! Z : two charged leptons
        call set_active(id_min_dR_l12j, "min_dr_lj", &
                       & "minimum of Delta R (jets, leptons)")
        ! including non-coloured final states
        call set_active(id_ht_full,     "ht_full", &
                       & "full HT")

        ! Added by KR
        call set_active(id_ystar_Vj, "ystar_Zj", &
                   & "ystar_Zj = |yj1-yz|/2")
        call set_active(id_yboost_Vj, "yboost_Zj", &
                   & "yboost_Zj = |yj1+yz|/2")
        ! End: Added by KR

        call set_active(id_ZJ_HT,     "z_ht", &
                       & "etz + sum pt{partons}")
        call set_active(id_ZJ_HTshape05,     "z_ht_shape05", &
                       & "mll/2 + sum pt{Z,partons} * 2")        
        call set_active(id_ZJ_HTshape20,     "z_ht_shape20", &
                       & "mll*2 + sum pt{Z,partons} / 2")

        call set_active(id_phi_star_eta,     "phi_star", &
                       & "\phi^*_\eta [arXiv:1009.1580]")
        ! Collins-Soper frame (angles wrt l1   = 47
        call set_active(id_cos_theta_CS, "cos_theta_CS", &
                      & "cos(theta) in Collins-Soper frame")
        call set_active(id_theta_CS,     "theta_CS", &
                      & "theta in Collins-Soper frame")
        call set_active(id_phi_CS,       "phi_CS", &
                      & "phi in Collins-Soper frame")
        ! projectors for the angular coefficients A_i [1606.00689]
        call set_active(id_proj_A0, "proj_A0", &
                      & "projector for the angular coefficient A0")
        call set_active(id_proj_A1, "proj_A1", &
                      & "projector for the angular coefficient A1")
        call set_active(id_proj_A2, "proj_A2", &
                      & "projector for the angular coefficient A2")
        call set_active(id_proj_A3, "proj_A3", &
                      & "projector for the angular coefficient A3")
        call set_active(id_proj_A4, "proj_A4", &
                      & "projector for the angular coefficient A4")
        call set_active(id_proj_A5, "proj_A5", &
                      & "projector for the angular coefficient A5")
        call set_active(id_proj_A6, "proj_A6", &
                      & "projector for the angular coefficient A6")
        call set_active(id_proj_A7, "proj_A7", &
                      & "projector for the angular coefficient A7")
        call set_active(id_proj_A0mA2, "proj_A0mA2", &
                      & "projector for the difference A0-A2")
        call set_active(id_cos_theta_CS_0PT, "cos_theta_CS_0PT", &
                      & "cos(theta) in Collins-Soper frame, usable at PTZ=0")
        ! dynamical scales
        call set_active(id_ZJ_dscale1,     "dscale1", &
                       & "sqrt(mll^2 + sum ptj{i}^2) [arXiv:1512.01291]")
        call set_active(id_ZJ_dscale3,     "dscale3", &
                       & "sqrt(mll^2 + ptZ^2) + sum ptj{i} [Les Houches 2017]")
        ! MiNLO 
        call set_active(id_minlo,     "minlo", &
                       & "MiNLO reweighting factor")
        call set_active(id_sudakov,     "sudakov", &
                       & "sudakov factor")


        ! ! debugging
        ! call set_active(id_idx_g1,     "idx_g1", &
        !                & "idx_g1")
        ! call set_active(id_idx_g2,     "idx_g2", &
        !                & "idx_g2")


      case ("Wm","WmJ","WmJJ","WmJJJ")
        ! leptons from the Wm -> nu-bar lm decay
        call set_active(id_pt_l1,       "ptlm", &
                       & "transverse momentum of l-")
        call set_active(id_y_l1,        "ylm", &
                       & "(pseudo-)rapidity of l-")
        call set_active(id_abs_y_l1,    "abs_ylm", &
                       & "|ylm|")
        call set_active(id_pt_l2,       "etmiss", &
                       & "missing transverse energy/momentum")
        call set_active(id_y_l2,        "ynu", &
                       & "(pseudo-)rapidity of nu-bar")
        call set_active(id_abs_y_l2,    "abs_ynu", &
                       & "|ynu|")
        ! Wm-boson = nu-bar + lm
        call set_active(id_pt_V,        "ptw", &
                       & "transverse momentum of W")
        call set_active(id_mt_V,        "mt", &
                       & "transverse mass of W")
        ! Wm : only one charged lepton at `npar`
        !      note, same name as in Z but different eval function
        call set_active(id_min_dR_l1j,  "min_dr_lj", &
                       & "minimum of Delta R (jets, l-)")
        ! including non-coloured final states
        call set_active(id_ht_full,     "ht_full", &
                       & "full HT")        
        call set_active(id_y_V,         "yw", &
                       & "rapidity of W")
        call set_active(id_abs_y_V,     "abs_yw", &
                       & "|yw|")
        call set_active(id_eta_V,       "etaw", &
                       & "pseudo-rapidity of W")
        call set_active(id_minv_V,      "mln", &
                       & "invariant mass of W")
        call set_active(id_et_V,     "etw", &
                       & "transverse energy: sqrt(mln^2 + ptw^2)")

      case ("Wp","WpJ","WpJJ","WpJJJ")
        ! leptons from the Wp -> lp nu decay
        call set_active(id_pt_l1,       "etmiss", &
                       & "missing transverse energy/momentum")
        call set_active(id_pt_l2,       "ptlp", &
                       & "transverse momentum of l+")
        call set_active(id_y_l2,        "ylp", &
                       & "(pseudo-)rapidity of l+")
        call set_active(id_abs_y_l2,    "abs_ylp", &
                       & "|ylp|")
        call set_active(id_y_l1,        "ynu", &
                       & "(pseudo-)rapidity of nu")
        call set_active(id_abs_y_l1,    "abs_ynu", &
                       & "|ynu|")
        ! Wp-boson = lp + nu
        call set_active(id_pt_V,        "ptw", &
                       & "transverse momentum of W")
        call set_active(id_mt_V,        "mt", &
                       & "transverse mass of W")
        ! Wm : only one charged lepton at `npar-1`
        !      note, same name as in Z but different eval function
        call set_active(id_min_dR_l2j,  "min_dr_lj", &
                       & "minimum of Delta R (jets, l+)")
        ! including non-coloured final states
        call set_active(id_ht_full,     "ht_full", &
                       & "full HT")
        call set_active(id_deltay_Vj1,     "deltay_zj1", &
                       & "|yz - yj1|")
        call set_active(id_y_V,         "yw", &
                       & "rapidity of W")
        call set_active(id_abs_y_V,     "abs_yw", &
                       & "|yw|")
        call set_active(id_eta_V,       "etaw", &
                       & "pseudo-rapidity of W")
        call set_active(id_minv_V,      "mln", &
                       & "invariant mass of W")
        call set_active(id_et_V,     "etw", &
                       & "transverse energy: sqrt(mln^2 + ptw^2)")

      case ("H","HJ","HJJ","HJJJ","Hm","HJm")
        ! H-boson
        call set_active(id_pt_V,        "pth", &
                       & "transverse momentum of H")
        call set_active(id_y_V,         "yh", &
                       & "rapidity of H")
        call set_active(id_abs_y_V,     "abs_yh", &
                       & "|yh|")
        call set_active(id_eta_V,       "etah", &
                       & "pseudo-rapidity of H")
        call set_active(id_abs_eta_V,   "abs_etah", &
                       & "|etah|")
        ! including non-coloured final states
        call set_active(id_ht_full,     "ht_full", &
                       & "full HT")
        call set_active(id_pt_Vjj,     "pthjj", &
                       & "transverse momentum of H plus two leading jets")
        call set_active(id_pt_Vj,     "pthj", &
                       & "transverse momentum of H plus the leading jet")
        call set_active(id_deltay_Vj1,     "deltay_hj1", &
                       & "|yh - yj1|")
        call set_active(id_deltaphi_Vj1,     "dphi_hj1", &
                       & "Azimuthal angle difference of H and 1st jet")
        call set_active(id_ZJ_dscale1,     "dscale1", &
                       & "sqrt(mgg(H)^2 + sum ptj{i}^2) [arXiv:1605.04692]")
        call set_active(id_et_V,     "dscale2", &
                       & "sqrt(mgg(H)^2 + ptz^2)")
        call set_active(id_ZJ_dscale3,     "dscale3", &
                       & "sqrt(mll^2 + ptZ^2) + sum ptj{i} [Les Houches 2017]")
        call set_active(id_phi_star_eta,     "phi_star", &
                       & "\phi^*_\eta [arXiv:1009.1580]")

      case ("Hto2pm","Hto2pJm","Hto2p","Hto2pJ","Hto2pJJ","Hto2pJJJ")
        ! H-boson
        call set_active(id_pt_V,        "pth", &
                       & "transverse momentum of H")
        call set_active(id_y_V,         "yh", &
                       & "rapidity of H")
        call set_active(id_abs_y_V,     "abs_yh", &
                       & "|yh|")
        call set_active(id_eta_V,       "etah", &
                       & "pseudo-rapidity of H")
        call set_active(id_abs_eta_V,   "abs_etah", &
                       & "|etah|")
        ! photon from the H -> gamma gamma decay
        call set_active(id_pt_g1,       "ptg1", &
                       & "transverse momentum of leading photon")
        call set_active(id_y_g1,        "yg1", &
                       & "(pseudo-)rapidity of leading photon")
        call set_active(id_abs_y_g1,    "abs_yg1", &
                       & "|yg1|")
        call set_active(id_pt_g2,       "ptg2", &
                       & "transverse momentum of sub-leading photon")
        call set_active(id_y_g2,        "yg2", &
                       & "(pseudo-)rapidity of sub-leading photon")
        call set_active(id_abs_y_g2,    "abs_yg2", &
                       & "|yg2|")
        call set_active(id_min_dR_l12j, "min_dr_gj", &
                       & "minimum of Delta R (jets, photons)")
        call set_active(id_deltaphi_g1g2,     "dphi_g1g2", &
                       & "Azimuthal angle difference of di-photon")
        call set_active(id_deltay_g1g2,     "dy_g1g2", &
                       & "(pseudo-)rapidity difference of di-photon")
        call set_active(id_abs_deltay_g1g2,     "abs_dy_g1g2", &
                       & "|dy_g1g2|")
        call set_active(id_costhetastar,      "costar_gg", &
                       & "photon decay angles in the Collins-Soper frame")
        call set_active(id_phi_star_eta,     "phi_star", &
                       & "\phi^*_\eta [arXiv:1009.1580]")
        call set_active(id_ptt_g1g2,     "ptt_g1g2", &
                       & "pTt of di-photon")
        ! including non-coloured final states
        call set_active(id_ht_full,     "ht_full", &
                       & "full HT")
        call set_active(id_pt_Vjj,     "pthjj", &
                       & "transverse momentum of H plus two leading jets")
        call set_active(id_pt_Vj,     "pthj", &
                       & "transverse momentum of H plus the leading jet")
        call set_active(id_deltay_Vj1,     "deltay_hj1", &
                       & "|yh - yj1|")
        call set_active(id_deltaphi_Vj1,     "dphi_hj1", &
                       & "Azimuthal angle difference of H and 1st jet")
        call set_active(id_ZJ_dscale1,     "dscale1", &
                       & "sqrt(mgg(H)^2 + sum ptj{i}^2) [arXiv:1605.04692]")
        call set_active(id_et_V,     "dscale2", &
                       & "sqrt(mgg(H)^2 + ptz^2)")
        call set_active(id_tau_j1,     "tau_j1", &
                       & "tau-jet of the leading jet")
        call set_active(id_max_tau_jet,     "max_tau_jet", &
                       & "the maximum of all tau-jet among all jets")
        call set_active(id_sum_tau_jets,     "sum_tau_jets", &
                       & "sum of all jet's tau-jet")

      case ("Hto2taum","Hto2tauJm","Hto2tau","Hto2tauJ","Hto2tauJJ","Hto2tauJJJ")
        ! H-boson
        call set_active(id_pt_V,        "pth", &
                       & "transverse momentum of H")
        call set_active(id_y_V,         "yh", &
                       & "rapidity of H")
        call set_active(id_abs_y_V,     "abs_yh", &
                       & "|yh|")
        call set_active(id_eta_V,       "etah", &
                       & "pseudo-rapidity of H")
        call set_active(id_abs_eta_V,   "abs_etah", &
                       & "|etah|")
        ! tau from the H -> tau+ tau- decay
        call set_active(id_pt_g1,       "pttau1", &
                       & "transverse momentum of leading tau")
        call set_active(id_y_g1,        "ytau1", &
                       & "(pseudo-)rapidity of leading tau")
        call set_active(id_abs_y_g1,    "abs_ytau1", &
                       & "|ytau1|")
        call set_active(id_pt_g2,       "pttau2", &
                       & "transverse momentum of sub-leading tau")
        call set_active(id_y_g2,        "ytau2", &
                       & "(pseudo-)rapidity of sub-leading tau")
        call set_active(id_abs_y_g2,    "abs_ytau2", &
                       & "|ytau2|")
        call set_active(id_min_dR_l12j, "min_dr_tauj", &
                       & "minimum of Delta R (jets, taus)")
        call set_active(id_deltaphi_g1g2,     "dphi_tau1tau2", &
                       & "Azimuthal angle difference of tau+ and tau-")
        call set_active(id_deltay_g1g2,     "dy_tau1tau2", &
                       & "(pseudo-)rapidity difference of tau+ and tau-")
        call set_active(id_abs_deltay_g1g2,     "abs_dy_tau1tau2", &
                       & "|dy_tau1tau2|")
        call set_active(id_costhetastar,      "costar_tautau", &
                       & "tau decay angles in the Collins-Soper frame")
        call set_active(id_phi_star_eta,     "phi_star", &
                       & "\phi^*_\eta [arXiv:1009.1580]")
        call set_active(id_ptt_g1g2,     "ptt_tau1tau2", &
                       & "pTt of tau+ and tau-")
        ! including non-coloured final states
        call set_active(id_ht_full,     "ht_full", &
                       & "full HT")
        call set_active(id_pt_Vjj,     "pthjj", &
                       & "transverse momentum of H plus two leading jets")
        call set_active(id_pt_Vj,     "pthj", &
                       & "transverse momentum of H plus the leading jet")
        call set_active(id_deltay_Vj1,     "deltay_hj1", &
                       & "|yh - yj1|")
        call set_active(id_deltaphi_Vj1,     "dphi_hj1", &
                       & "Azimuthal angle difference of H and 1st jet")
        call set_active(id_ZJ_dscale1,     "dscale1", &
                       & "sqrt(mgg(H)^2 + sum ptj{i}^2) [arXiv:1605.04692]")
        call set_active(id_et_V,     "dscale2", &
                       & "sqrt(mgg(H)^2 + ptz^2)")
        call set_active(id_tau_j1,     "tau_j1", &
                       & "tau-jet of the leading jet")
        call set_active(id_max_tau_jet,     "max_tau_jet", &
                       & "the maximum of all tau-jet among all jets")
        call set_active(id_sum_tau_jets,     "sum_tau_jets", &
                       & "sum of all jet's tau-jet")

      case ("Hto2l1pm","Hto2l1pJm","Hto2l1p","Hto2l1pJ","Hto2l1pJJ","Hto2l1pJJJ")
        ! H-boson
        call set_active(id_pt_V,        "pth", &
                       & "transverse momentum of H")
        call set_active(id_y_V,         "yh", &
                       & "rapidity of H")
        call set_active(id_abs_y_V,     "abs_yh", &
                       & "|yh|")
        call set_active(id_eta_V,       "etah", &
                       & "pseudo-rapidity of H")
        call set_active(id_abs_eta_V,   "abs_etah", &
                       & "|etah|")
        call set_active(id_minv_V23,      "m23", &
                       & "invariant mass of di-lepton system")
        call set_active(id_pt_l4,       "ptgamma", &
                       & "transverse momentum of photon")

      case ("Hto4em","Hto4eJm","Hto4e","Hto4eJ","Hto4eJJ","Hto4eJJJ")
        ! leptons from the Z -> lp lm decay
        call set_active(id_pt_l1,       "ptem1", &
                       & "transverse momentum of one e-")
        call set_active(id_y_l1,        "yem1", &
                       & "(pseudo-)rapidity of one e-")
        call set_active(id_abs_y_l1,    "abs_yem1", &
                       & "|ye1m|")
        call set_active(id_pt_l2,       "ptep1", &
                       & "transverse momentum of one e+")
        call set_active(id_y_l2,        "yep1", &
                       & "(pseudo-)rapidity of one e+")
        call set_active(id_abs_y_l2,    "abs_yep1", &
                       & "|ye1p|")
        call set_active(id_pt_l3,       "ptem2", &
                       & "transverse momentum of another e-")
        call set_active(id_y_l3,        "yem2", &
                       & "(pseudo-)rapidity of another e-")
        call set_active(id_abs_y_l3,    "abs_yem2", &
                       & "|ye2m|")
        call set_active(id_pt_l4,       "ptep2", &
                       & "transverse momentum of another e+")
        call set_active(id_y_l4,        "yep2", &
                       & "(pseudo-)rapidity of another e+")
        call set_active(id_abs_y_l4,    "abs_yep2", &
                       & "|ye2p|")
        ! ordered leptons from the two Z -> lp lm decay
        call set_active(id_pt_l1st,       "pte1st", &
                       & "transverse momentum of the leading lepton in the 4e system")
        call set_active(id_y_l1st,        "ye1st", &
                       & "(pseudo-)rapidity of the leading lepton in the 4e system")
        call set_active(id_abs_y_l1st,    "abs_ye1st", &
                       & "|ye1st|")
        call set_active(id_pt_l2nd,       "pte2nd", &
                       & "transverse momentum of  the sub-leading lepton in the 4e system")
        call set_active(id_y_l2nd,        "ye2nd", &
                       & "(pseudo-)rapidity of  the sub-leading lepton in the 4e system")
        call set_active(id_abs_y_l2nd,    "abs_ye2nd", &
                       & "|ye2nd|")
        call set_active(id_pt_l3rd,       "pte3rd", &
                       & "transverse momentum of the thrid-leading lepton in the 4e system")
        call set_active(id_y_l3rd,        "ye3rd", &
                       & "(pseudo-)rapidity of  the thrid-leading lepton in the 4e system")
        call set_active(id_abs_y_l3rd,    "abs_ye3rd", &
                       & "|ye3rd|")
        call set_active(id_pt_l4th,       "pte4th", &
                       & "transverse momentum of the least-leading lepton in the 4e system")
        call set_active(id_y_l4th,        "ye4th", &
                       & "(pseudo-)rapidity of the least-leading lepton in the 4e system")
        call set_active(id_abs_y_l4th,    "abs_ye4th", &
                       & "|ye4th|")
        ! Z-boson = one e+ + e- system
        call set_active(id_pt_V,        "ptepem1", &
                       & "transverse momentum of one Z-->e+e- ")
        call set_active(id_y_V,         "yepem1", &
                       & "rapidity of one Z-->e+e-")
        call set_active(id_abs_y_V,     "abs_yepem1", &
                       & "|yepem1|")
        call set_active(id_eta_V,       "etaepem1", &
                       & "pseudo-rapidity of one Z-->e+e-")
        call set_active(id_abs_eta_V,   "abs_etaepem1", &
                       & "|etaepem1|")
        call set_active(id_minv_V,      "mepem1", &
                       & "invariant mass of one Z-->e+e-")
        call set_active(id_et_V,     "etepem1", &
                       & "transverse energy of one e+e- pair: sqrt(mepem1^2 + ptz^2)")
        ! Z-boson = another e+ + e- system
        call set_active(id_pt_V2,        "ptepem2", &
                       & "transverse momentum of another Z-->e+e- ")
        call set_active(id_y_V2,         "yepem2", &
                       & "rapidity of another Z-->e+e-")
        call set_active(id_abs_y_V2,     "abs_yepem2", &
                       & "|ymupmum2|")
        call set_active(id_eta_V2,       "etaepem2", &
                       & "pseudo-rapidity of another Z-->e+e-")
        call set_active(id_abs_eta_V2,   "abs_etaepem2", &
                       & "|etamupmum2|")
        call set_active(id_minv_V2,      "mepem2", &
                       & "invariant mass of another Z-->e+e-")
        call set_active(id_et_V2,     "etmupmum2", &
                       & "transverse energy: sqrt(mepem2^2 + ptz^2)")
        ! angular separation
        call set_active(id_min_dR_l_sf,     "min_dr_epem", &
                       & "minimum Delta R_e+e- of the two e+e- systems")
        call set_active(id_min_dR_l_df,     "min_dr_e1e2", &
                       & "minimum Delta R of ep1_ep2, ep1_em2, em1_ep2, em1_em2")
        call set_active(id_min_dR_l1234j, "min_dr_lj", &
                       & "minimum of Delta R (jets, leptons)")
        call set_active(id_R_l1st_j1, "dr_e1st_j1", &
                       & "Delta R distance between (leading jet, leading lepton)")
        ! Z-boson = 4e ordering system
        call set_active(id_minv_mZSFOS1st,      "m2lfarestmZ", &
                       & "invariant mass of SFOS e+e- pair farest to the Z-mass")
        call set_active(id_minv_mZSFOS2nd,      "m2lfarmZ", &
                       & "invariant mass of SFOS e+e- pair third close to the Z-mass")
        call set_active(id_minv_mZSFOS3rd,      "m2lclosemZ", &
                       & "invariant mass of SFOS e+e- pair second close to the Z-mass")
        call set_active(id_minv_mZSFOS4th,      "m2lclosestmZ", &
                       & "invariant mass of SFOS e+e- pair closest to the Z-mass")
        ! H-boson = l1 + l2 +l3 +l4
        call set_active(id_pt_V4l,        "pth", &
                       & "transverse momentum of H")
        call set_active(id_minv_V4l,      "m4l", &
                       & "invariant mass of H")
        call set_active(id_y_V4l,         "yh", &
                       & "rapidity of H")
        call set_active(id_abs_y_V4l,     "abs_yh", &
                       & "|yh|")
        call set_active(id_eta_V4l,       "etah", &
                       & "pseudo-rapidity of H")
        call set_active(id_abs_eta_V4l,   "abs_etah", &
                       & "|etah|")
        ! dynamical scales
        call set_active(id_H4lJ_dscale1,     "dscale1", &
                       & "sqrt(m4l^2 + sum ptj{i}^2) [arXiv:1512.01291]")
        call set_active(id_H4lJ_dscale2,     "dscale2", &
                       & "sqrt(m4l^2 + ptz^2)")


      case ("Hto2e2mum","Hto2e2muJm","Hto2e2mu","Hto2e2muJ","Hto2e2muJJ","Hto2e2muJJJ")
        ! leptons from the Z -> lp lm decay
        call set_active(id_pt_l1,       "ptmum", &
                       & "transverse momentum of mu-")
        call set_active(id_y_l1,        "ymum", &
                       & "(pseudo-)rapidity of mu-")
        call set_active(id_abs_y_l1,    "abs_ymum", &
                       & "|ymum|")
        call set_active(id_pt_l2,       "ptmup", &
                       & "transverse momentum of mu+")
        call set_active(id_y_l2,        "ymup", &
                       & "(pseudo-)rapidity of mu+")
        call set_active(id_abs_y_l2,    "abs_ymup", &
                       & "|ymup|")
        call set_active(id_pt_l3,       "ptem", &
                       & "transverse momentum of e-")
        call set_active(id_y_l3,        "yem", &
                       & "(pseudo-)rapidity of e-")
        call set_active(id_abs_y_l3,    "abs_yem", &
                       & "|yem|")
        call set_active(id_pt_l4,       "ptep", &
                       & "transverse momentum of e+")
        call set_active(id_y_l4,        "yep", &
                       & "(pseudo-)rapidity of e+")
        call set_active(id_abs_y_l4,    "abs_yep", &
                       & "|yep|")
        ! ordered leptons from the two Z -> lp lm decay
        call set_active(id_pt_l1st,       "ptl1st", &
                       & "transverse momentum of the leading lepton in the 2e2mu system")
        call set_active(id_y_l1st,        "yl1st", &
                       & "(pseudo-)rapidity of the leading lepton in the 2e2mu system")
        call set_active(id_abs_y_l1st,    "abs_yl1st", &
                       & "|yl1st|")
        call set_active(id_pt_l2nd,       "ptl2nd", &
                       & "transverse momentum of  the sub-leading lepton in the 2e2mu system")
        call set_active(id_y_l2nd,        "yl2nd", &
                       & "(pseudo-)rapidity of  the sub-leading lepton in the 2e2mu system")
        call set_active(id_abs_y_l2nd,    "abs_yl2nd", &
                       & "|yl2nd|")
        call set_active(id_pt_l3rd,       "ptl3rd", &
                       & "transverse momentum of the thrid-leading lepton in the 2e2mu system")
        call set_active(id_y_l3rd,        "yl3rd", &
                       & "(pseudo-)rapidity of  the thrid-leading lepton in the 2e2mu system")
        call set_active(id_abs_y_l3rd,    "abs_yl3rd", &
                       & "|yl3rd|")
        call set_active(id_pt_l4th,       "ptl4th", &
                       & "transverse momentum of the least-leading lepton in the 2e2mu system")
        call set_active(id_y_l4th,        "yl4th", &
                       & "(pseudo-)rapidity of the least-leading lepton in the 2e2mu system")
        call set_active(id_abs_y_l4th,    "abs_yl4th", &
                       & "|yl4th|")
        ! Z-boson = mu+ + mu- system
        call set_active(id_pt_V,        "ptmupmum", &
                       & "transverse momentum of Z-->mu+mu- ")
        call set_active(id_y_V,         "ymupmum", &
                       & "rapidity of Z-->mu+mu-")
        call set_active(id_abs_y_V,     "abs_ymupmum", &
                       & "|ymupmum|")
        call set_active(id_eta_V,       "etaepem", &
                       & "pseudo-rapidity of Z-->e+e-")
        call set_active(id_abs_eta_V,   "abs_etamupmum", &
                       & "|etamupmum|")
        call set_active(id_minv_V,      "mmupmum", &
                       & "invariant mass of Z-->mu+mu-")
        call set_active(id_et_V,     "etmupmum", &
                       & "transverse energy: sqrt(mmupmum^2 + ptz^2)")
        call set_active(id_dR_l12, "dr_mupmum", &
                       & "Delta R (mu+, mu-)")
        ! Z-boson = e+ + e- system
        call set_active(id_pt_V2,        "ptepem", &
                       & "transverse momentum of Z-->e+e- ")
        call set_active(id_y_V2,         "yepem", &
                       & "rapidity of Z-->e+e-")
        call set_active(id_abs_y_V2,     "abs_yepem", &
                       & "|yepem|")
        call set_active(id_eta_V2,       "etaepem", &
                       & "pseudo-rapidity of Z-->e+e-")
        call set_active(id_abs_eta_V2,   "abs_etaepem", &
                       & "|etaepem|")
        call set_active(id_minv_V2,      "mepem", &
                       & "invariant mass of Z-->e+e-")
        call set_active(id_et_V2,     "etepem", &
                       & "transverse energy: sqrt(mepem^2 + ptz^2)")
        call set_active(id_dR_l34, "dr_epem", &
                       & "Delta R (e+, e-)")
        ! Z-boson alternative l+l- system
        call set_active(id_minv_V23,      "mepmum", &
                       & "invariant mass of Z-->e+mu-")
        call set_active(id_minv_V14,      "memmup", &
                       & "invariant mass of Z-->e-mu+")
        ! angular separation
        call set_active(id_min_dR_l_sf,     "min_dr_sf", &
                       & "minimum of Delta R (same flavour leptons)")
        call set_active(id_min_dR_l_df,     "min_dr_df", &
                       & "minimum of Delta R (different flavour leptons)")
        call set_active(id_min_dR_l1234j, "min_dr_lj", &
                       & "minimum of Delta R (jets, leptons)")
        call set_active(id_min_dR_l12j, "min_dr_muj", &
                       & "minimum of Delta R (jets, muons)")
        call set_active(id_min_dR_l34j, "min_dr_ej", &
                       & "minimum of Delta R (jets, electrons)")
        call set_active(id_R_l1st_j1, "dr_l1st_j1", &
                       & "Delta R distance between (leading jet, leading lepton)")
        ! Z-boson = four lepton ordering system
        call set_active(id_minv_mZ1,      "m2lnearmZ", &
                       & "invariant mass of lepton pair closest to the Z-mass")
        call set_active(id_minv_mZ2,      "m2lfarmZ", &
                       & "invariant mass of lepton pair fartest to the Z-mass")
        call set_active(id_minvSFOS_vs_mZ,    "m2em2muvsmZ", &
                       & "identify the SFOS minv closest to mZ and label 1 (2) for m2e (m2mu)")
        ! H-boson = l1 + l2 +l3 +l4
        call set_active(id_pt_V4l,        "pth", &
                       & "transverse momentum of H")
        call set_active(id_minv_V4l,      "m4l", &
                       & "invariant mass of H")
        call set_active(id_y_V4l,         "yh", &
                       & "rapidity of H")
        call set_active(id_abs_y_V4l,     "abs_yh", &
                       & "|yh|")
        call set_active(id_eta_V4l,       "etah", &
                       & "pseudo-rapidity of H")
        call set_active(id_abs_eta_V4l,   "abs_etah", &
                       & "|etah|")
        ! dynamical scales
        call set_active(id_H4lJ_dscale1,     "dscale1", &
                       & "sqrt(m4l^2 + sum ptj{i}^2) [arXiv:1512.01291]")
        call set_active(id_H4lJ_dscale2,     "dscale2", &
                       & "sqrt(m4l^2 + ptz^2)")

      case ("Hto2l2nm","Hto2l2nJm","Hto2l2n","Hto2l2nJ","Hto2l2nJJ","Hto2l2nJJJ")
        ! H-boson = l1 + l2 +l3 +l4
        call set_active(id_pt_V4l,        "pth", &
                       & "transverse momentum of H")
        call set_active(id_minv_V4l,      "m2l2n", &
                       & "invariant mass of H")
        call set_active(id_y_V4l,         "yh", &
                       & "rapidity of H")
        call set_active(id_abs_y_V4l,     "abs_yh", &
                       & "|yh|")
        call set_active(id_eta_V4l,       "etah", &
                       & "pseudo-rapidity of H")
        call set_active(id_abs_eta_V4l,   "abs_etah", &
                       & "|etah|")
        ! dynamical scales
        call set_active(id_H4lJ_dscale1,     "dscale1", &
                       & "sqrt(m2l2n^2 + sum ptj{i}^2) [arXiv:1512.01291]")
        call set_active(id_H4lJ_dscale2,     "dscale2", &
                       & "sqrt(m2l2n^2 + ptz^2)")

      case ("VFH", "VFHJ")
        ! H-boson
        call set_active(id_pt_V,        "pth", &
                       & "transverse momentum of H")
        call set_active(id_y_V,         "yh", &
                       & "rapidity of H")
        call set_active(id_abs_y_V,     "abs_yh", &
                       & "|yh|")
        call set_active(id_eta_V,       "etah", &
                       & "pseudo-rapidity of H")
        call set_active(id_abs_eta_V,   "abs_etah", &
                       & "|etah|")
        ! including non-coloured final states
        call set_active(id_ht_full,     "ht_full", &
                       & "full HT")
        ! VFH cuts
        call set_active(id_VFH_deltay,  "deltay", &
                       & "|y1 - y2|")
        call set_active(id_VFH_y1xy2,  "y1xy2", &
                       & "y1*y2")
        call set_active(id_VFH_dscale,  "dscale", &
                       & "sqrt(mH/2 * sqrt(mH^2/4 + ptH^2))")
        call set_active(id_VFH_z3, "z3", &  ! normalised central rapidity distribution of 3rd jet:
                       & " z3 = (y3 - (y1+y2)/2)/(y1-y2)")

      case("GJ")
         call set_active(id_photon_pt,  "pt_gam", &
                       & "pt of photon")
         call set_active(id_photon_y,  "y_gam", &
                       & "rapidity of photon")

      case ("DIS1j","SFDIS","RESDIS")
        ! DIS electron cuts
        call set_active(id_DIS_Q2,        "q2", &
                       & "momentum transfer of the lepton")
        call set_active(id_DIS_sqrtQ2,     "q",  &
                       & "sqrt of q2")
        call set_active(id_DIS_x,         "x", &
                       & "Bjorken x")
        call set_active(id_DIS_y,         "y", &
                       & "fraction of electron energy transferred to proton")
        call set_active(id_y01, "y01", &
                       & "0 to 1 transition rate (JADE) *before* jet cuts")
        call set_active(id_y12, "y12", &
                       & "1 to 2 transition rate (JADE) *before* jet cuts")
        call set_active(id_y23, "y23", &
                       & "2 to 3 transition rate (JADE) *before* jet cuts")
        call set_active(id_y34, "y34", &
                       & "3 to 4 transition rate (JADE) *before* jet cuts")
        call set_active(id_DIS_resrho, "resrho", &
                       & "resummed jet mass (rho)")

      case ("DIS")
        ! DIS electron cuts
        call set_active(id_DIS_dscl,        "dscl", &
                       & "sqrt((Q2+ptavg_j12**2)/2)")
        call set_active(id_DIS_dscl3j,        "dscl3j", &
                       & "sqrt((Q2+ptavg_j123**2)/2)")
        call set_active(id_DIS_dsclj1,        "dsclj1", &
                       & "sqrt((Q2+pt_jet1**2)/2)")
        call set_active(id_DIS_dsclj2,        "dsclj2", &
                       & "sqrt((Q2+pt_jet2**2)/2)")
        call set_active(id_DIS_dsclj3,        "dsclj3", &
                       & "sqrt((Q2+pt_jet3**2)/2)")
        call set_active(id_DIS_dsclj4,        "dsclj4", &
                       & "sqrt((Q2+pt_jet4**2)/2)")
        call set_active(id_DIS_dsclj1l,        "dsclj1l", &
                       & "sqrt((Q2+pt_jet1_lab**2)/2)")
        call set_active(id_DIS_dsclj2l,        "dsclj2l", &
                       & "sqrt((Q2+pt_jet2_lab**2)/2)")
        call set_active(id_DIS_dsclj3l,        "dsclj3l", &
                       & "sqrt((Q2+pt_jet3_lab**2)/2)")
        call set_active(id_DIS_dsclj4l,        "dsclj4l", &
                       & "sqrt((Q2+pt_jet4_lab**2)/2)")
        call set_active(id_DIS_dsclZEUS,        "dsclZEUS", &
                       & "sqrt(Q2+etavg_j12**2)")
        call set_active(id_DIS_dsclZEUS2,        "dsclZEUS2", &
                       & "sqrt(Q2+etavg_jall**2)")
        call set_active(id_DIS_Q2,        "q2", &
                       & "momentum transfer of the lepton")
        call set_active(id_DIS_W2,        "W2", &
                       & "(photon-proton energy)^2")
        call set_active(id_DIS_W,        "W", &
                       & "photon-proton energy")
        call set_active(id_DIS_nu,        "disnu", &
                       & "energy of photon in proton restframe")
        call set_active(id_DIS_visbyW,    "visbyW", &
                       & "E_visible/W")
        call set_active(id_DIS_cosgammah,        "cosgammah", &
                       & "cos(angle) of scattered quark ")
        call set_active(id_DIS_sqrtQ2,     "q",  &
                       & "sqrt of q2")
        call set_active(id_DIS_x,         "x", &
                       & "Bjorken x")
        call set_active(id_DIS_y,         "y", &
                       & "fraction of electron energy transferred to proton")
        call set_active(id_DIS_xi2,         "xi2", &
                       & "xbj*(1d0+m12**2/q2)")
        call set_active(id_DIS_Thrust,         "dis_thrust", &
                       & "Thrust in DIS (tau=1-T)")
        call set_active(id_DIS_Thrust_c,         "dis_thrust_c", &
                       & "Thrust tau w.r.t. thrust axis in DIS")
        call set_active(id_DIS_JB,         "dis_JB", &
                       & "Jet Broadening in DIS")
        call set_active(id_DIS_JM2,         "dis_JM2", &
                       & "Squared Jet Mass in DIS")
        call set_active(id_DIS_C,         "dis_C", &
                       & "C-parameter in DIS")                       
        call set_active(id_DIS_hera_eta_j1_0,         "eta1_0", &
                       & "HERA eta of leading jet in Breit frame before cuts")
        call set_active(id_DIS_hera_eta_j2_0,         "eta2_0", &
                       & "HERA eta of subleading jet in Breit frame before cuts")
        call set_active(id_DIS_Lgxi2,         "lgxi2", &
                       & "log(xi2)")
        call set_active(id_DIS_etas,         "etas", &
                       & "log(xi2)")
        call set_active(id_DIS_xi3,         "xi3", &
                       & "abs(eta_j1-eta_j2)/2")
        call set_active(id_ptavg_j12, "ptavg_12", &
                       & "(ptj1+ptj2)/2")
        call set_active(id_ptavg_j123, "ptavg_123", &
                       & "(ptj1+ptj2+ptj3)/3")
        call set_active(id_ptavg_jall, "ptavg_all", &
                       & "(ptj1+ptj2+...)/njets")
        call set_active(id_etavg_j12, "etavg_12", &
                       & "(etj1+etj2)/2")
        call set_active(id_et2_j12, "et2_12", &
                       & "(etj1+etj2)/2")
        call set_active(id_etavg_j123, "etavg_123", &
                       & "(etj1+etj2+etj3)/3")
        call set_active(id_etavg_jall, "etavg_all", &
                       & "(etj1+etj2+...)/njets")

        call set_active(id_hera_ptj, "hera_jets_pt", &
                       & "jet pt in HERA frame")
        call set_active(id_hera_yj, "hera_jets_y", &
                       & "jet rapidity in HERA frame")
        call set_active(id_hera_abs_yj, "hera_jets_abs_y", &
                       & "|hera_jets_y|")
        call set_active(id_hera_etaj, "hera_jets_eta", &
                       & "jet pseudo-rapidity in HERA frame")
        call set_active(id_hera_abs_etaj, "hera_jets_abs_eta", &
                       & "|etaj_hera|")
        call set_active(id_hera_etj, "hera_jets_et", &
                       & "jet Et in HERA frame")
        call set_active(id_hera_xjets, "hera_xjets", &
                       & "x jet in HERA frame: p_{z,jet}/E_P")
        call set_active(id_hera_H1xjets, "hera_H1xjets", &
                       & "x jet in HERA frame: E_{jet}/E_P")
        call set_active(id_dis_hera_etj1, "hera_etj1", &
                       & "Et of leading jet in hera frame")
        call set_active(id_dis_hera_etj2, "hera_etj2", &
                       & "Et of subleading jet in hera frame")
        call set_active(id_dis_hera_etj3, "hera_etj3", &
                       & "Et of 3rd jet in hera frame")
        call set_active(id_dis_hera_etj4, "hera_etj4", &
                       & "Et of 4th jet in hera frame")
        call set_active(id_dis_hera_et2j1, "hera_et2j1", &
                       & "Et2 of leading jet in hera frame")
        call set_active(id_dis_hera_et2j2, "hera_et2j2", &
                       & "Et2 of subleading jet in hera frame")
        call set_active(id_dis_hera_et2j3, "hera_et2j3", &
                       & "Et2 of 3rd jet in hera frame")
        call set_active(id_dis_hera_et2j4, "hera_et2j4", &
                       & "Et2 of 4th jet in hera frame")
        call set_active(id_dis_hera_etaj1, "hera_etaj1", &
                       & "Eta of leading jet in hera frame")
        call set_active(id_dis_hera_etaj2, "hera_etaj2", &
                       & "Eta of subleading jet in hera frame")
        call set_active(id_dis_hera_etaj3, "hera_etaj3", &
                       & "Eta of 3rd jet in hera frame")
        call set_active(id_dis_hera_etaj4, "hera_etaj4", &
                       & "Eta of 4th jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j1, "hera_et2byQ2j1", &
                       & "Et^2/Q2 of leading jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j2, "hera_et2byQ2j2", &
                       & "Et^2/Q2 of subleading jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j3, "hera_et2byQ2j3", &
                       & "Et^2/Q2 of 3rd jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j4, "hera_et2byQ2j4", &
                       & "Et^2/Q2 of 4th jet in hera frame")
        call set_active(id_DIS_hera_FJ_reweight, "hera_FJ_reweight", &
                       & "Q2+sumET2")
        call set_active(id_dis_hera_deta1, "hera_deta1", &
                       & "etaj1-etaj2 for +2 jets")
        call set_active(id_dis_hera_deta2, "hera_deta2", &
                       & "etafj-etaj1 for +2 jets")
        call set_active(id_DIS_hera_etaavg_j12, "hera_etaavg_12", &
                       & "0.5*(etaj1+etaj2)")
        call set_active(id_DIS_hera_deltaeta_j12, "hera_deltaeta_12", &
                       & "|etaj1-etaj2|")
        call set_active(id_gammap_eta_j1, "gammap_eta1",&
                       & "eta1 in gamma proton frame")
        call set_active(id_gammap_eta_j2, "gammap_eta2",&
                       & "eta2 in gamma proton frame")
        call set_active(id_dis_delta_phis, "delta_phis",&
                       & "phi1-phi2 in proton gamma frame")
        call set_active(id_DIS_hera_xgamma, "xgamma",&
                       & "(E-pZ)_LabJets / (2*Ee*y)")

        call set_active(id_y01, "y01", &
                       & "0 to 1 transition rate (JADE) *before* jet cuts")
        call set_active(id_y12, "y12", &
                       & "1 to 2 transition rate (JADE) *before* jet cuts")
        call set_active(id_y23, "y23", &
                       & "2 to 3 transition rate (JADE) *before* jet cuts")
        call set_active(id_y34, "y34", &
                       & "3 to 4 transition rate (JADE) *before* jet cuts")
        call set_active(id_DIS_ptl, "ptl", &
                       & "transverse momentum of lepton")

      case ("DISWm")
        ! DIS electron cuts
        call set_active(id_DIS_dscl,        "dscl", &
                       & "sqrt((Q2+ptavg_j12**2)/2)")
        call set_active(id_DIS_dscl3j,        "dscl3j", &
                       & "sqrt((Q2+ptavg_j123**2)/2)")
        call set_active(id_DIS_dsclj1,        "dsclj1", &
                       & "sqrt((Q2+pt_jet1**2)/2)")
        call set_active(id_DIS_dsclj2,        "dsclj2", &
                       & "sqrt((Q2+pt_jet2**2)/2)")
        call set_active(id_DIS_dsclj3,        "dsclj3", &
                       & "sqrt((Q2+pt_jet3**2)/2)")
        call set_active(id_DIS_dsclj4,        "dsclj4", &
                       & "sqrt((Q2+pt_jet4**2)/2)")
        call set_active(id_DIS_dsclZEUS,        "dsclZEUS", &
                       & "sqrt(Q2+etavg_j12**2)")
        call set_active(id_DIS_dsclZEUS2,        "dsclZEUS2", &
                       & "sqrt(Q2+etavg_jall**2)")
        call set_active(id_DIS_Q2,        "q2", &
                       & "momentum transfer of the lepton")
        call set_active(id_DIS_W2,        "W2", &
                       & "(photon-proton energy)^2")
        call set_active(id_DIS_W,        "W", &
                       & "photon-proton energy")
        call set_active(id_DIS_nu,        "disnu", &
                       & "energy of photon in proton restframe")
        call set_active(id_DIS_visbyW,    "visbyW", &
                       & "E_visible/W")
        call set_active(id_DIS_cosgammah,        "cosgammah", &
                       & "cos(angle) of scattered quark ")
        call set_active(id_DIS_sqrtQ2,     "q",  &
                       & "sqrt of q2")
        call set_active(id_DIS_x,         "x", &
                       & "Bjorken x")
        call set_active(id_DIS_y,         "y", &
                       & "fraction of electron energy transferred to proton")
        call set_active(id_DIS_xi2,         "xi2", &
                       & "xbj*(1d0+m12**2/q2)")
        call set_active(id_DIS_Thrust,         "dis_thrust", &
                       & "Thrust in DIS (tau=1-T)")
        call set_active(id_DIS_Thrust_c,         "dis_thrust_c", &
                       & "Thrust tau w.r.t. thrust axis in DIS")
        call set_active(id_DIS_JB,         "dis_JB", &
                       & "Jet Broadening in DIS")
        call set_active(id_DIS_JM2,         "dis_JM2", &
                       & "Squared Jet Mass in DIS")
        call set_active(id_DIS_C,         "dis_C", &
                       & "C-parameter in DIS")                       
        call set_active(id_DIS_hera_eta_j1_0,         "eta1_0", &
                       & "HERA eta of leading jet in Breit frame before cuts")
        call set_active(id_DIS_hera_eta_j2_0,         "eta2_0", &
                       & "HERA eta of subleading jet in Breit frame before cuts")
        call set_active(id_DIS_Lgxi2,         "lgxi2", &
                       & "log(xi2)")
        call set_active(id_DIS_etas,         "etas", &
                       & "log(xi2)")
        call set_active(id_DIS_xi3,         "xi3", &
                       & "abs(eta_j1-eta_j2)/2")
        call set_active(id_ptavg_j12, "ptavg_12", &
                       & "(ptj1+ptj2)/2")
        call set_active(id_ptavg_j123, "ptavg_123", &
                       & "(ptj1+ptj2+ptj3)/3")
        call set_active(id_ptavg_jall, "ptavg_all", &
                       & "(ptj1+ptj2+...)/njets")
        call set_active(id_etavg_j12, "etavg_12", &
                       & "(etj1+etj2)/2")
        call set_active(id_et2_j12, "et2_12", &
                       & "(etj1+etj2)/2")
        call set_active(id_etavg_j123, "etavg_123", &
                       & "(etj1+etj2+etj3)/3")
        call set_active(id_etavg_jall, "etavg_all", &
                       & "(etj1+etj2+...)/njets")

        call set_active(id_hera_ptj, "hera_jets_pt", &
                       & "jet pt in HERA frame")
        call set_active(id_hera_yj, "hera_jets_y", &
                       & "jet rapidity in HERA frame")
        call set_active(id_hera_abs_yj, "hera_jets_abs_y", &
                       & "|hera_jets_y|")
        call set_active(id_hera_etaj, "hera_jets_eta", &
                       & "jet pseudo-rapidity in HERA frame")
        call set_active(id_hera_abs_etaj, "hera_jets_abs_eta", &
                       & "|etaj_hera|")
        call set_active(id_hera_etj, "hera_jets_et", &
                       & "jet Et in HERA frame")
        call set_active(id_hera_xjets, "hera_xjets", &
                       & "x jet in HERA frame: p_{z,jet}/E_P")
        call set_active(id_hera_H1xjets, "hera_H1xjets", &
                       & "x jet in HERA frame: E_{jet}/E_P")
        call set_active(id_dis_hera_etj1, "hera_etj1", &
                       & "Et of leading jet in hera frame")
        call set_active(id_dis_hera_etj2, "hera_etj2", &
                       & "Et of subleading jet in hera frame")
        call set_active(id_dis_hera_etj3, "hera_etj3", &
                       & "Et of 3rd jet in hera frame")
        call set_active(id_dis_hera_etj4, "hera_etj4", &
                       & "Et of 4th jet in hera frame")
        call set_active(id_dis_hera_et2j1, "hera_et2j1", &
                       & "Et2 of leading jet in hera frame")
        call set_active(id_dis_hera_et2j2, "hera_et2j2", &
                       & "Et2 of subleading jet in hera frame")
        call set_active(id_dis_hera_et2j3, "hera_et2j3", &
                       & "Et2 of 3rd jet in hera frame")
        call set_active(id_dis_hera_et2j4, "hera_et2j4", &
                       & "Et2 of 4th jet in hera frame")
        call set_active(id_dis_hera_etaj1, "hera_etaj1", &
                       & "Eta of leading jet in hera frame")
        call set_active(id_dis_hera_etaj2, "hera_etaj2", &
                       & "Eta of subleading jet in hera frame")
        call set_active(id_dis_hera_etaj3, "hera_etaj3", &
                       & "Eta of 3rd jet in hera frame")
        call set_active(id_dis_hera_etaj4, "hera_etaj4", &
                       & "Eta of 4th jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j1, "hera_et2byQ2j1", &
                       & "Et^2/Q2 of leading jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j2, "hera_et2byQ2j2", &
                       & "Et^2/Q2 of subleading jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j3, "hera_et2byQ2j3", &
                       & "Et^2/Q2 of 3rd jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j4, "hera_et2byQ2j4", &
                       & "Et^2/Q2 of 4th jet in hera frame")
        call set_active(id_DIS_hera_FJ_reweight, "hera_FJ_reweight", &
                       & "Q2+sumET2")
        call set_active(id_dis_hera_deta1, "hera_deta1", &
                       & "etaj1-etaj2 for +2 jets")
        call set_active(id_dis_hera_deta2, "hera_deta2", &
                       & "etafj-etaj1 for +2 jets")
        call set_active(id_DIS_hera_etaavg_j12, "hera_etaavg_12", &
                       & "0.5*(etaj1+etaj2)")
        call set_active(id_DIS_hera_deltaeta_j12, "hera_deltaeta_12", &
                       & "|etaj1-etaj2|")
        call set_active(id_gammap_eta_j1, "gammap_eta1",&
                       & "eta1 in gamma proton frame")
        call set_active(id_gammap_eta_j2, "gammap_eta2",&
                       & "eta2 in gamma proton frame")
        call set_active(id_dis_delta_phis, "delta_phis",&
                       & "phi1-phi2 in proton gamma frame")
        call set_active(id_DIS_hera_xgamma, "xgamma",&
                       & "(E-pZ)_LabJets / (2*Ee*y)")
        call set_active(id_y01, "y01", &
                       & "0 to 1 transition rate (JADE)")
        call set_active(id_y12, "y12", &
                       & "1 to 2 transition rate (JADE)")
        call set_active(id_y23, "y23", &
                       & "2 to 3 transition rate (JADE)")
        call set_active(id_y34, "y34", &
                       & "3 to 4 transition rate (JADE)")
        call set_active(id_DIS_ptl, "ptl", &
                       & "transverse momentum of lepton")

      case ("DISWp")
        ! DIS electron cuts
        call set_active(id_DIS_dscl,        "dscl", &
                       & "sqrt((Q2+ptavg_j12**2)/2)")
        call set_active(id_DIS_dscl3j,        "dscl3j", &
                       & "sqrt((Q2+ptavg_j123**2)/2)")
        call set_active(id_DIS_dsclj1,        "dsclj1", &
                       & "sqrt((Q2+pt_jet1**2)/2)")
        call set_active(id_DIS_dsclj2,        "dsclj2", &
                       & "sqrt((Q2+pt_jet2**2)/2)")
        call set_active(id_DIS_dsclj3,        "dsclj3", &
                       & "sqrt((Q2+pt_jet3**2)/2)")
        call set_active(id_DIS_dsclj4,        "dsclj4", &
                       & "sqrt((Q2+pt_jet4**2)/2)")
        call set_active(id_DIS_dsclZEUS,        "dsclZEUS", &
                       & "sqrt(Q2+etavg_j12**2)")
        call set_active(id_DIS_dsclZEUS2,        "dsclZEUS2", &
                       & "sqrt(Q2+etavg_jall**2)")
        call set_active(id_DIS_Q2,        "q2", &
                       & "momentum transfer of the lepton")
        call set_active(id_DIS_W2,        "W2", &
                       & "(photon-proton energy)^2")
        call set_active(id_DIS_W,        "W", &
                       & "photon-proton energy")
        call set_active(id_DIS_nu,        "disnu", &
                       & "energy of photon in proton restframe")
        call set_active(id_DIS_visbyW,    "visbyW", &
                       & "E_visible/W")
        call set_active(id_DIS_cosgammah,        "cosgammah", &
                       & "cos(angle) of scattered quark ")
        call set_active(id_DIS_sqrtQ2,     "q",  &
                       & "sqrt of q2")
        call set_active(id_DIS_x,         "x", &
                       & "Bjorken x")
        call set_active(id_DIS_y,         "y", &
                       & "fraction of electron energy transferred to proton")
        call set_active(id_DIS_xi2,         "xi2", &
                       & "xbj*(1d0+m12**2/q2)")
        call set_active(id_DIS_Thrust,         "dis_thrust", &
                       & "Thrust in DIS (tau=1-T)")
        call set_active(id_DIS_Thrust_c,         "dis_thrust_c", &
                       & "Thrust tau w.r.t. thrust axis in DIS")
        call set_active(id_DIS_JB,         "dis_JB", &
                       & "Jet Broadening in DIS")
        call set_active(id_DIS_JM2,         "dis_JM2", &
                       & "Squared Jet Mass in DIS")
        call set_active(id_DIS_C,         "dis_C", &
                       & "C-parameter in DIS")                       
        call set_active(id_DIS_hera_eta_j1_0,         "eta1_0", &
                       & "HERA eta of leading jet in Breit frame before cuts")
        call set_active(id_DIS_hera_eta_j2_0,         "eta2_0", &
                       & "HERA eta of subleading jet in Breit frame before cuts")
        call set_active(id_DIS_Lgxi2,         "lgxi2", &
                       & "log(xi2)")
        call set_active(id_DIS_etas,         "etas", &
                       & "log(xi2)")
        call set_active(id_DIS_xi3,         "xi3", &
                       & "abs(eta_j1-eta_j2)/2")
        call set_active(id_ptavg_j12, "ptavg_12", &
                       & "(ptj1+ptj2)/2")
        call set_active(id_ptavg_j123, "ptavg_123", &
                       & "(ptj1+ptj2+ptj3)/3")
        call set_active(id_ptavg_jall, "ptavg_all", &
                       & "(ptj1+ptj2+...)/njets")
        call set_active(id_etavg_j12, "etavg_12", &
                       & "(etj1+etj2)/2")
        call set_active(id_et2_j12, "et2_12", &
                       & "(etj1+etj2)/2")
        call set_active(id_etavg_j123, "etavg_123", &
                       & "(etj1+etj2+etj3)/3")
        call set_active(id_etavg_jall, "etavg_all", &
                       & "(etj1+etj2+...)/njets")
        call set_active(id_hera_ptj, "hera_jets_pt", &
                       & "jet pt in HERA frame")
        call set_active(id_hera_yj, "hera_jets_y", &
                       & "jet rapidity in HERA frame")
        call set_active(id_hera_abs_yj, "hera_jets_abs_y", &
                       & "|hera_jets_y|")
        call set_active(id_hera_etaj, "hera_jets_eta", &
                       & "jet pseudo-rapidity in HERA frame")
        call set_active(id_hera_abs_etaj, "hera_jets_abs_eta", &
                       & "|etaj_hera|")
        call set_active(id_hera_etj, "hera_jets_et", &
                       & "jet Et in HERA frame")
        call set_active(id_hera_xjets, "hera_xjets", &
                       & "x jet in HERA frame: p_{z,jet}/E_P")
        call set_active(id_hera_H1xjets, "hera_H1xjets", &
                       & "x jet in HERA frame: E_{jet}/E_P")
        call set_active(id_dis_hera_etj1, "hera_etj1", &
                       & "Et of leading jet in hera frame")
        call set_active(id_dis_hera_etj2, "hera_etj2", &
                       & "Et of subleading jet in hera frame")
        call set_active(id_dis_hera_etj3, "hera_etj3", &
                       & "Et of 3rd jet in hera frame")
        call set_active(id_dis_hera_etj4, "hera_etj4", &
                       & "Et of 4th jet in hera frame")
        call set_active(id_dis_hera_et2j1, "hera_et2j1", &
                       & "Et2 of leading jet in hera frame")
        call set_active(id_dis_hera_et2j2, "hera_et2j2", &
                       & "Et2 of subleading jet in hera frame")
        call set_active(id_dis_hera_et2j3, "hera_et2j3", &
                       & "Et2 of 3rd jet in hera frame")
        call set_active(id_dis_hera_et2j4, "hera_et2j4", &
                       & "Et2 of 4th jet in hera frame")
        call set_active(id_dis_hera_etaj1, "hera_etaj1", &
                       & "Eta of leading jet in hera frame")
        call set_active(id_dis_hera_etaj2, "hera_etaj2", &
                       & "Eta of subleading jet in hera frame")
        call set_active(id_dis_hera_etaj3, "hera_etaj3", &
                       & "Eta of 3rd jet in hera frame")
        call set_active(id_dis_hera_etaj4, "hera_etaj4", &
                       & "Eta of 4th jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j1, "hera_et2byQ2j1", &
                       & "Et^2/Q2 of leading jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j2, "hera_et2byQ2j2", &
                       & "Et^2/Q2 of subleading jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j3, "hera_et2byQ2j3", &
                       & "Et^2/Q2 of 3rd jet in hera frame")
        call set_active(id_dis_hera_et2byQ2j4, "hera_et2byQ2j4", &
                       & "Et^2/Q2 of 4th jet in hera frame")
        call set_active(id_DIS_hera_FJ_reweight, "hera_FJ_reweight", &
                       & "Q2+sumET2")
        call set_active(id_dis_hera_deta1, "hera_deta1", &
                       & "etaj1-etaj2 for +2 jets")
        call set_active(id_dis_hera_deta2, "hera_deta2", &
                       & "etafj-etaj1 for +2 jets")
        call set_active(id_DIS_hera_etaavg_j12, "hera_etaavg_12", &
                       & "0.5*(etaj1+etaj2)")
        call set_active(id_DIS_hera_deltaeta_j12, "hera_deltaeta_12", &
                       & "|etaj1-etaj2|")
        call set_active(id_gammap_eta_j1, "gammap_eta1",&
                       & "eta1 in gamma proton frame")
        call set_active(id_gammap_eta_j2, "gammap_eta2",&
                       & "eta2 in gamma proton frame")
        call set_active(id_dis_delta_phis, "delta_phis",&
                       & "phi1-phi2 in proton gamma frame")
        call set_active(id_DIS_hera_xgamma, "xgamma",&
                       & "(E-pZ)_LabJets / (2*Ee*y)")
        call set_active(id_y01, "y01", &
                       & "0 to 1 transition rate (JADE)")
        call set_active(id_y12, "y12", &
                       & "1 to 2 transition rate (JADE)")
        call set_active(id_y23, "y23", &
                       & "2 to 3 transition rate (JADE)")
        call set_active(id_y34, "y34", &
                       & "3 to 4 transition rate (JADE)")


        
        !> ep em observables
     case ("EPEM")
        call set_active(id_epem_C, "C", &
             & "C parameter")
        call set_active(id_epem_D, "D", &
             & "D parameter")
        call set_active(id_epem_WJB, "WJB", &
             & "Wide jet broadening")
        call set_active(id_epem_TJB, "TJB", &
             & "total jet broadening")
        call set_active(id_epem_T, "T", &
             & "Thrust")
        call set_active(id_epem_1mT, "1mT", &
             & "1-T")
        call set_active(id_epem_cosT, "cosT", &
             & "angle between thrust and z axis")
        call set_active(id_epem_HJM, "HJM", &
             & "Heavy jet mass")
        call set_active(id_epem_SHM, "SHM", &
             & "Sum hemisphere masses")
        call set_active(id_y45 , "y45", &
             & "5 to 4 transition rate (JADE)")
        call set_active(id_y34, "y34", &
             & "4 to 3 transition rate (JADE)")
        call set_active(id_y23, "y23", &
             & "2 to 3 transition rate (JADE)")
        call set_active(id_epem_costh_ej1, "cost_ej1", &
             & "costh between leading jet and electron direction")
        call set_active(id_epem_costh_en3j, "costh_en3j", &
             & "costh bewteen 3 jet plane and electron direction")
        call set_active(id_epem_chi, "orchi", &
             & "angle between 3 jet and e^- - leading jet plane")

      case ("1jet", "2jet", "3jet", "4jet")
        !> multi jet observables
        call set_active(id_ptavg_j12, "ptavg_j12", &
                       & "(ptj1+ptj2)/2")
        call set_active(id_ptavg_jall, "ptavg_all", &
                       & "(ptj1+ptj2+...)/njets")
        call set_active(id_ptavg_geom_jall, "ptavg_geom_all", &
                       & "(ptj1*ptj2*...)^{1/njets}")
        !> "perfect" jet observables
        call set_active(id_pt_j1_nocut,    "ptj1_nocut", &
                       & "ptj1 with no jet cuts")
        call set_active(id_pt_j2_nocut,    "ptj2_nocut", &
                       & "ptj2 with no jet cuts")
        call set_active(id_ht_jall_nocut,    "ht_jets_nocut", &
                       & "ht_jets with no jet cuts")
        call set_active(id_ptavg_j12_nocut,    "ptavg_j12_nocut", &
                       & "(ptj1+ptj2)/2 with no jet cuts")
        call set_active(id_ptavg_jall_nocut,    "ptavg_jall_nocut", &
                       & "(ptj1+ptj2+...)/njets with no jet cuts")
        call set_active(id_ptavg_geom_j12_nocut,    "ptavg_geom_j12_nocut", &
                       & "sqrt(ptj1*ptj2) with no jet cuts")
        call set_active(id_ptavg_geom_jall_nocut,    "ptavg_geom_jall_nocut", &
                       & "(ptj1*ptj2*...)^{1/njets} with no jet cuts")
        !> R=1 jet observables (are "perfect")
        call set_active(id_pt_j1_R1,    "ptj1_R1", &
                       & "ptj1 using Rcone=1")
        call set_active(id_pt_j2_R1,    "ptj2_R1", &
                       & "ptj2 using Rcone=1")
        call set_active(id_ht_jall_R1,    "ht_jets_R1", &
                       & "ht_jets using Rcone=1")
        call set_active(id_ptavg_j12_R1,    "ptavg_j12_R1", &
                       & "(ptj1+ptj2)/2 using Rcone=1")
        call set_active(id_ptavg_jall_R1,    "ptavg_jall_R1", &
                       & "(ptj1+ptj2+...)/njets using Rcone=1")
        call set_active(id_ptavg_geom_j12_R1,    "ptavg_geom_j12_R1", &
                       & "sqrt(ptj1*ptj2) using Rcone=1")
        call set_active(id_ptavg_geom_jall_R1,    "ptavg_geom_jall_R1", &
                       & "(ptj1*ptj2*...)^{1/njets} using Rcone=1")

    end select


    ! Spikeplots
    call set_active(id_s13,     "s13", &
                   & "s13 invariant / shat")
    call set_active(id_s14,     "s14", &
                   & "s14 invariant / shat")
    call set_active(id_s15,     "s15", &
                   & "s15 invariant / shat")
    call set_active(id_s16,     "s16", &
                   & "s16 invariant / shat")
    call set_active(id_s23,     "s23", &
                   & "s23 invariant / shat")
    call set_active(id_s24,     "s24", &
                   & "s24 invariant / shat")
    call set_active(id_s25,     "s25", &
                   & "s25 invariant / shat")
    call set_active(id_s26,     "s26", &
                   & "s26 invariant / shat")
    call set_active(id_s34,     "s34", &
                   & "s34 invariant / shat")
    call set_active(id_s35,     "s35", &
                   & "s35 invariant / shat")
    call set_active(id_s36,     "s36", &
                   & "s36 invariant / shat")
    call set_active(id_s45,     "s45", &
                   & "s45 invariant / shat")
    call set_active(id_s46,     "s46", &
                   & "s46 invariant / shat")
    call set_active(id_s56,     "s56", &
                   & "s56 invariant / shat")

    !----------------------
    ! set jet observables !
    !----------------------
    ! number of reconstructed jets
    call set_active(id_njets,       "njets", &
                   &"number of reconstructed jets")
    call set_active(id_npartons,       "npartons", &
                   &"number of final-state partons (for debugging!)")
    ! leading jet
    call set_active(id_pt_j1,       "ptj1", &
                   &"transverse momentum of leading jet")
    call set_active(id_y_j1,        "yj1", &
                   &"rapidity of leading jet")
    call set_active(id_abs_y_j1,    "abs_yj1", &
                   &"|yj1|")
    call set_active(id_eta_j1,      "etaj1", &
                   &"pseudo-rapidity of leading jet")
    call set_active(id_abs_eta_j1,  "abs_etaj1", &
                   &"|etaj1|")
    ! sub-leading jet
    call set_active(id_pt_j2,       "ptj2", &
                   &"transverse momentum of sub-leading jet")
    call set_active(id_y_j2,        "yj2", &
                   &"rapidity of sub-leading jet")
    call set_active(id_abs_y_j2,    "abs_yj2", &
                   &"|yj2|")
    call set_active(id_eta_j2,      "etaj2", &
                   &"pseudo-rapidity of sub-leading jet")
    call set_active(id_abs_eta_j2,  "abs_etaj2", &
                   &"|etaj2|")
    ! 3rd jet
    call set_active(id_pt_j3,       "ptj3", &
                   &"transverse momentum of 3rd jet")
    call set_active(id_y_j3,        "yj3", &
                   &"rapidity of sub-leading jet")
    call set_active(id_abs_y_j3,    "abs_yj3", &
                   &"|yj3|")
    call set_active(id_eta_j3,      "etaj3", &
                   &"pseudo-rapidity of 3rd jet")
    call set_active(id_abs_eta_j3,  "abs_etaj3", &
                   &"|etaj3|")
    ! 4th jet
    call set_active(id_pt_j4,       "ptj4", &
                   &"transverse momentum of 4th jet")
    call set_active(id_y_j4,        "yj4", &
                   &"rapidity of sub-leading jet")
    call set_active(id_abs_y_j4,    "abs_yj4", &
                   &"|yj4|")
    call set_active(id_eta_j4,      "etaj4", &
                   &"pseudo-rapidity of 4th jet")
    call set_active(id_abs_eta_j4,  "abs_etaj4", &
                   &"|etaj4|")
    call set_active(id_sum_pt_jets, "ptsum_jets", &
                   &"sum of jet pt's")
    ! dijet invariant mass
    call set_active(id_minv_j12,    "m12", &
                   &"dijet invariant mass")
    ! more di-jet observables
    call set_active(id_deltapt_j12, "deltapt_j12", &
                   & "|ptj1-ptj2|")
    call set_active(id_deltaeta_j12, "deltaeta_j12", &
                   & "|etaj1-etaj2|")
    call set_active(id_deltay_j12, "deltay_j12", &
                   & "|yj1-yj2|")
    call set_active(id_deltaphi_j12, "deltaphi_j12", &
                   & "azimuthal angle between j1 & j2")
    call set_active(id_chi_j12, "chi", &
                   & "chi = exp(|yj1-yj2|)")
    call set_active(id_ystar_j12, "ystar", &
                   & "ystar = |yj1-yj2|/2")
    call set_active(id_yboost_j12, "yboost", &
                   & "yboost = |yj1+yj2|/2")
    call set_active(id_max_y_j12, "ymax_j12", &
                   & "max(yj1,yj2)")
    call set_active(id_pt_j1expYstar,   "ptj1EXP_03ystar", &
                   &"ptj1*dexp(0.3*ystar)")


    ! trijet invariant mass
    call set_active(id_minv_j123,    "m123", &
                   &"trijet invariant mass")
    ! HT only summing over coloured final-state particles
    call set_active(id_ht_jets,     "ht_jets", &
                   &"HT of jets")
    call set_active(id_ht_part,     "ht_part", &
                   &"HT of partons")


    !> set the intrinsic Fortran random number generator
    call random_seed()


    !> @todo set of pseudo-observables like `log(sij)` etc... for debugging

    contains

      subroutine set_active(id,name,desc)
        integer, intent(in) :: id
        character(len=*), intent(in) :: name, desc
        list_obs(id) = Observable_t(.true., .false., name, desc)
      end subroutine set_active

  end subroutine init_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> routine to deallocate global members
  !-----------------------------------------------------------------------------
  subroutine destroy_obs()
    deallocate(mapUsed2All_obs)
  end subroutine destroy_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> mark the observable (id) as used
  !
  !> @param[in] id unique observable identifier
  !-----------------------------------------------------------------------------
  recursive subroutine setUsed_obs(id)
    integer, intent(in) :: id
    integer :: i, j
    if ( .not.list_obs(id)%qactive ) then
      stop "setUsed_obs: can't set an inactive observables to a used state"
    end if
    ! if ( isExtra_obs(id) ) then
    !   select case (id)
    !     !> MiNLO needs ptz and mll
    !     case (id_minlo)
    !       call setUsed_obs(id_pt_V)
    !       call setUsed_obs(id_minv_V)
    !       !> also register sudakov observables
    !       !> to make them also "cachable"!!!
    !       call setUsed_obs(id_sudakov)
    !   end select
    !   return
    ! end if
    if (.not.list_obs(id)%qused) then
      list_obs(id)%qused = .true.
      nUsed_obs = nUsed_obs + 1
      select case (id) ! Spike plots
      case(id_s13, id_s14, id_s15, id_s16, id_s23, id_s24, id_s25, id_s26, &
           & id_s34, id_s35, id_s36, id_s45, id_s46, id_s56)
         select case (id)
         case (id_s13) 
            i = 1 ; j = 3
         case (id_s14)
            i = 1 ; j = 4
         case (id_s15)
            i = 1 ; j = 5
         case (id_s16)
            i = 1 ; j = 6
         case (id_s23)
            i = 2 ; j = 3
         case (id_s24)
            i = 2 ; j = 4
         case (id_s25)
            i = 2 ; j = 5
         case (id_s26)
            i = 2 ; j = 6
         case (id_s34)
            i = 3 ; j = 4
         case (id_s35)
            i = 3 ; j = 5
         case (id_s36)
            i = 3 ; j = 6
         case (id_s45)
            i = 4 ; j = 5
         case (id_s46)
            i = 4 ; j = 6
         case (id_s56)
            i = 5 ; j = 6
         end select
         call initialiseSpike(i,j)
      end select
    end if
  end subroutine setUsed_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> purge all un-used observables and set up the mapping arrays
  !-----------------------------------------------------------------------------
  subroutine purgeUsed_obs()
    integer :: i,j

    ! character(len=*), parameter :: fmt_act = '(1x,"*",3x,A10,"---",3x,A50,"*")'
    ! print*, "used observables"
    ! do i=1,n_obs
    !   if (.not.list_obs(i)%qused) cycle
    !     write(*,fmt_act) list_obs(i)%name, list_obs(i)%desc
    ! end do

    allocate(mapUsed2All_obs(nUsed_obs))
    mapUsed2All_obs = 0
    mapAll2Used_obs = 0
    j = 0
    do i = 1,n_obs
      if (list_obs(i)%qused) then
        j = j+1
        mapUsed2All_obs(j) = i
        mapAll2Used_obs(i) = j
      end if
    end do
    if ( j.ne.nUsed_obs ) stop "purgeUsed_obs: number mismatch"
  end subroutine purgeUsed_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> inqure if an observable is used in a run
  !
  !> @param[in] id unique observable identifier
  !> @return .true. if used, .false. otherwise
  !-----------------------------------------------------------------------------
  logical function isUsed_obs(id)
    integer, intent(in) :: id
    isUsed_obs = list_obs(id)%qused
  end function isUsed_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> print out a list of active observables
  !-----------------------------------------------------------------------------
  subroutine digestActive_obs()
    integer :: i
    character(len=*), parameter :: fmt_act = '(1x,"*",3x,A20,"---",3x,A60,"*")'
    character(len=*), parameter :: fmt_one = '(1x,"*",3x,A86,"*")'

    write(*,fmt_one) " "
    write(*,fmt_one) "  > jet veto specifiers <  "
    write(*,fmt_one) " "
    do i=1,n_obsj
      if (.not.list_obsj(i)%qactive) cycle
      write(*,fmt_one) list_obsj(i)%name
    end do

    write(*,fmt_one) " "
    write(*,fmt_one) "  > process-specific observables <  "
    write(*,fmt_one) " "
    do i=1,n_obs
      if (.not.list_obs(i)%qactive) cycle
        write(*,fmt_act) list_obs(i)%name, list_obs(i)%desc
    end do

  end subroutine digestActive_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> routine to initialize the (threadprivate) observables cache
  !-----------------------------------------------------------------------------
  subroutine initBuffer_obs()
    allocate(valueUsed_obs(nUsed_obs))
    allocate(validUsed_obs(nUsed_obs))
    npar_current = -1
  end subroutine initBuffer_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> routine to deallocate the (threadprivate) observables cache
  !-----------------------------------------------------------------------------
  subroutine destroyBuffer_obs()
    deallocate(valueUsed_obs)
    deallocate(validUsed_obs)
  end subroutine destroyBuffer_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> Convert the observable identifier to its corresponding name string
  !
  !> @param[in] id unique observable identifier
  !> @return name string
  !-----------------------------------------------------------------------------
  function getNameFromId_obs(id)
    character(len=max_obs_name_length) :: getNameFromId_obs
    integer, intent(in) :: id
    if (id < 0) then ! negative id's = jet observables in the reconstruction
      if ( -id > n_obsj .or. .not.list_obsj(-id)%qactive)  &
        &  stop "invalid jet observable id"
      getNameFromId_obs = list_obsj(-id)%name
    else
      if (id < 1 .or. id > n_obs) stop "invalid observable id"
      getNameFromId_obs = list_obs(id)%name
    end if
  end function getNameFromId_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> Convert a name string to the observable identifier
  !
  !> @param[in] name string containing the name of an observable
  !> @return observable id (0 if unsuccessful)
  !-----------------------------------------------------------------------------
  function getIdFromName_obs(name)
    use StrHelper_mod, only : to_upper_str
    integer :: getIdFromName_obs
    character(len=*), intent(in) :: name
    integer :: i
    getIdFromName_obs = 0
    do i=1,n_obs
      if (.not.list_obs(i)%qactive) cycle
      if ( trim(to_upper_str(list_obs(i)%name)) == &
         & trim(to_upper_str(name)) ) then
        getIdFromName_obs = i
        return
      endif
    end do
    ! didn't find it in the set of Observables, check jet-observables
    do i=1,n_obsj
      if (.not.list_obsj(i)%qactive) cycle
      if ( trim(to_upper_str(list_obsj(i)%name)) == &
         & trim(to_upper_str(name)) ) then
        getIdFromName_obs = -i ! negative id's for jet observables
        return
      endif
    end do
  end function getIdFromName_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> Evaluate all registered observables
  !
  !> @param[in] npar kinematics specified by the number of external partlicles
  !> @param[in] skipJetAlgo optional flag to skip jet clustering
  !> @param[in] skipCache optional flag to skip reading & writing to the cache
  !-----------------------------------------------------------------------------
  subroutine evalAll_obs(npar, skipJetAlgo, skipCache)
    use KinData_mod
    integer, intent(in) :: npar
    logical, intent(in), optional :: skipJetAlgo, skipCache
    integer :: i
    logical :: skipJ, skipC

    !> save the current npar
    npar_current =  npar

    skipC = .false.  ! default
    if (present(skipCache)) skipC = skipCache

    !> clean up the local cache before calling `eval_obs`
    call clearCache_obs()

    if ( kin(npar)%qcachedObs .and. .not.skipC ) then
      !> we have the observables cached, read in
      call readCacheObs_kin(npar, valueUsed_obs)
      do i=1,nUsed_obs
        validUsed_obs(i) = isValid_obs(mapUsed2All_obs(i) , npar)
      end do

    else
      !> need to actually compute the observables
      skipJ = .false.
      if (present(skipJetAlgo)) skipJ = skipJetAlgo
      if (.not.skipJ) call recombine_jet(npar)

      !> set manual observables here...
      ! call evalAll_manual(npar)

      !> compute all used observables
      do i=1,nUsed_obs
        !> validUsed must be set *AFTER* the jet algorithm
        validUsed_obs(i) = isValid_obs(mapUsed2All_obs(i) , npar)
        if ( validUsed_obs(i) ) then
          valueUsed_obs(i) = eval_obs( mapUsed2All_obs(i) , npar)
        else
          valueUsed_obs(i) = 0d0
        end if
      end do
      !> store in the kin cache
      if (.not.skipC) call writeCacheObs_kin(npar, valueUsed_obs)
    end if

  end subroutine evalAll_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> Print all registered observables
  !
  !> @param[in] npar kinematics specified by the number of external partlicles
  !> @param[in] skipCache optional flag to skip reading from the cache
  !> @param[in] unit optional logical unit connected for formatted output (file, screen)
  !-----------------------------------------------------------------------------
  subroutine printAll_obs(npar, skipCache, unit)
    use KinData_mod
    integer, intent(in) :: npar
    logical, intent(in), optional :: skipCache
    integer, intent(in), optional :: unit
    double precision, dimension(nUsed_obs) :: valueUsed_npar
    logical,          dimension(nUsed_obs) :: validUsed_npar
    integer :: i, iout
    logical :: skipC
    character(len=*), parameter :: fmt_int     = '(1x,"*",3x,A,T30,I10,T51,"*")'
    character(len=*), parameter :: fmt_dble    = '(1x,"*",3x,A,T30,1e18.12,T51,"*")'
    character(len=*), parameter :: fmt_invalid = '(1x,"*",3x,A,T30,"<invalid>",T51,"*")'

    skipC = .false.  ! default
    if (present(skipCache)) skipC = skipCache

    if (present(unit)) then
      iout = unit
    else
      iout = 6 ! output to screen
    endif

    !> get the cached values or recompute and store in local arrays
    if ( kin(npar)%qcachedObs .and. .not.skipC ) then
      call readCacheObs_kin(npar, valueUsed_npar)
      do i=1,nUsed_obs
        validUsed_npar(i) = isValid_obs(mapUsed2All_obs(i) , npar)
      end do
    else
      write(iout,*) ">>> re-evaluate observables <<<"
      !> clean up the local cache before calling `eval_obs`
      call clearCache_obs()
      !> then re-compute them
      do i=1,nUsed_obs
        !> depending on where the write statement is, the global varibale might not been set yet
        !> will then be printed as invalid!
        if(isGlobal_obs(mapUsed2All_obs(i)).and. .not.isCached_obs(mapUsed2All_obs(i)))  cycle
        validUsed_npar(i) = isValid_obs(mapUsed2All_obs(i) , npar)
        if ( validUsed_npar(i)) then
          valueUsed_npar(i) = eval_obs( mapUsed2All_obs(i) , npar)
        else
          valueUsed_npar(i) = 0d0
        end if
      end do
    endif
    
    !> pretty printout
    write(iout,'(/,1x,50("*"),T15," printAll_obs(",I1,") ")') npar
    write(iout,fmt_int) "cacheID", kin(npar)%cacheID
    write(iout,fmt_int) "npar_from", kin(npar)%npar_from
    do i=1,nUsed_obs
        if ( validUsed_npar(i) ) then
          if ( isInteger_obs(mapUsed2All_obs(i)) ) then
            write(iout,fmt_int) trim(getNameFromId_obs(mapUsed2All_obs(i))), INT(valueUsed_npar(i))
          else
            write(iout,fmt_dble) trim(getNameFromId_obs(mapUsed2All_obs(i))), valueUsed_npar(i)
          end if
        else
          write(iout,fmt_invalid) trim(getNameFromId_obs(mapUsed2All_obs(i)))
        end if
    end do
    write(iout,'(1x,50("*"))')

  end subroutine printAll_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> Compare observables between currently cached ones and cache in KinData
  !
  !> @param[in] npar kinematics specified by the number of external partlicles
  !> @param[in] threshold optional value for the threshold
  !> @param[in] verbose optional integer flag to control verbosity
  !> @param[out] fail optional flag to check if obs agree within threshold
  !-----------------------------------------------------------------------------
  subroutine compCache_obs(npar, threshold, verbose, fail)
    use KinData_mod
    integer, intent(in) :: npar
    double precision, intent(in), optional :: threshold
    integer, intent(in), optional :: verbose
    logical, intent(out), optional :: fail
    double precision :: thresh
    integer :: vrb
    integer :: i
    double precision :: diff,sum,reldiff, maxreldiff
    double precision, dimension(nUsed_obs) :: valueUsed_npar
    logical,          dimension(nUsed_obs) :: validUsed_npar
    character(len=*), parameter :: fmt_int      = '(1x,"*",3x,A,T30,I10,T95,"*")'
    character(len=*), parameter :: fmt_dble     = '(1x,"*",3x,A,T30,1e18.12,T95,"*")'
    character(len=*), parameter :: fmt_one_dble = '(1e18.12)'
    character(len=*), parameter :: fmt_comp     = '(1x,"*",3x,A,T20,A,T45,A,T70,A,T95,"*")'
    character(len=*), parameter :: fmt_invalid  = '(1x,"*",3x,A,T30,"<invalid>",T95,"*")'
    character(len=30) :: v1, v2, v3


    !> check if cached
    if ( .not. kin(npar)%qcachedObs ) then
      print*, "compCache_obs: no cached observables to compare to:", npar
      return
    end if

    !> threshold
    thresh = 1d-3  ! default value: ~ 3-4 digit agreement
    if (present(threshold)) thresh = threshold

    !> verbosity
    vrb = 1
    if (present(verbose)) vrb = verbose

    !> get the npar-cache values
    call readCacheObs_kin(npar, valueUsed_npar)
    do i=1,nUsed_obs
      validUsed_npar(i) = isValid_obs(mapUsed2All_obs(i) , npar)
    end do

    !> compute maxreldiff
    maxreldiff = 0d0
    do i=1,nUsed_obs
      if ( .not.validUsed_obs(i) .and. .not.validUsed_npar(i) ) cycle
      if ( .not.validUsed_obs(i) .or.  .not.validUsed_npar(i) ) then
        reldiff = thresh
      else
        diff = valueUsed_obs(i) - valueUsed_npar(i)
        sum  = valueUsed_obs(i) + valueUsed_npar(i)
        reldiff = 0d0
        if ( sum .ne. 0d0 ) reldiff = diff / sum
      end if
      if ( abs(reldiff) > maxreldiff ) maxreldiff = abs(reldiff)
    end do

    !> flag compatibility
    if (present(fail)) fail = ( maxreldiff >= thresh )

    !> optional output
    if (vrb > 0) then
      !> pretty printout
      write(*,'(/,1x,94("*"),T35," compCache_obs(",I1,") ")') npar
      write(*,fmt_int) "cacheID", kin(npar)%cacheID
      write(*,fmt_int) "npar_from", kin(npar)%npar_from
      write(*,fmt_dble) "threshold", thresh
      write(*,fmt_dble) "maxreldiff", maxreldiff
      if ( maxreldiff >= thresh ) then
        write(*,fmt_comp) "----------", "Observables", "KinData", "rel. diff"
        do i=1,nUsed_obs
          if ( .not.validUsed_obs(i) .and. .not.validUsed_npar(i) ) cycle
          write(v1,fmt_one_dble) valueUsed_obs(i)
          write(v2,fmt_one_dble) valueUsed_npar(i)
          if ( .not.validUsed_obs(i) .or.  .not.validUsed_npar(i) ) then
            v3 = '<invalid>'
            if ( .not.validUsed_obs(i) )  v1 = '<invalid>'
            if ( .not.validUsed_npar(i) ) v2 = '<invalid>'
          else
            diff = valueUsed_obs(i) - valueUsed_npar(i)
            sum  = valueUsed_obs(i) + valueUsed_npar(i)
            reldiff = 0d0
            if ( sum .ne. 0d0 ) reldiff = diff / sum
            write(v3,fmt_one_dble) reldiff
          end if
          write(*,fmt_comp) trim(getNameFromId_obs(mapUsed2All_obs(i))), &
            & v1, v2, v3
        end do
      end if
    end if
    
  end subroutine compCache_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> retrieve the value for the observable (id)
  !
  !> @param[in] id unique observable identifier
  !> @param[out] valid flag
  !> @return value of the observable
  !-----------------------------------------------------------------------------
  function getValue_obs(id, valid)
    double precision :: getValue_obs
    integer, intent(in) :: id
    logical, intent(out) :: valid
    
    if (id < 1) stop "getValue_obs: invalid input"

    if (.not.list_obs(id)%qused) stop "getValue_obs: inactive observable"

    if ( isGlobal_obs(id) ) then
      !> global observables default to true
      !> getManual_obs will error out if value not set
      !> we should *not* retrieve this value from the array
      !> reason is that once isCached is true, it will never
      !> be reset. Since this function can be called before
      !> evalAll_obs which fills the array it has to be treated
      !> separately!
      !> manual observables are determined on a kinematics basis 
      !> and should therefore be read in from the array, the "local"
      !> value in the EvalObs_mod is not necessarily the up-to-date 
      !> one if we read in from the kinematics cache
      valid        = .true.
      getValue_obs = getManual_obs(id)
    else if ( isExtraj_obs(id) ) then
      !> extras should never be accessed this way
      valid        = .false.
      getValue_obs = 0d0
    else
      valid        = validUsed_obs(mapAll2Used_obs(id))
      getValue_obs = valueUsed_obs(mapAll2Used_obs(id))
    end if 
  end function getValue_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> retrieve the value for the observable (id)
  !
  !> @param[in] id unique observable identifier
  !> @return value of the observable
  !-----------------------------------------------------------------------------
  function getValid_obs(id)
    logical :: getValid_obs
    integer, intent(in) :: id
    if (id < 1) stop "getValid_obs: invalid input"
    if ( isGlobal_obs(id) ) then
      !> global observables default to true
      !> getManual_obs will error out if value not set
      getValid_obs = .true.
    else if ( isExtraj_obs(id) ) then
      !> extras should never be accessed this way
      getValid_obs = .false.
    else
      getValid_obs = validUsed_obs(mapAll2Used_obs(id))
    end if 
  end function getValid_obs


  !-----------------------------------------------------------------------------
  !> @brief
  !> Apply the jet recombination and store the jets in descending order of
  !> transverse momentum in kd%pjets.
  !
  !> @param[in] npar kinematics specified by the number of particles
  !-----------------------------------------------------------------------------
  recursive subroutine recombine_jet(npar)
    integer, intent(in) :: npar
    type(Kin_t), pointer :: kd
    integer :: i
    !>----- for special 2->2 config in dijets
    double precision :: pt2j1, pt2j2, del_pt2_abs, del_pt2_rel
    double precision, parameter :: threshold_del_pt2_abs = 0.001d0**2  ! GeV**2
    double precision, parameter :: threshold_del_pt2_rel = 1d-8
    integer :: npar_copy

    !> prettify code
    kd => kin(npar)

    !> recombination has been performed already in the past
    if ( kd%qcachedObs .and. kd%njets >= 0 ) then
      return
    end if

    !> no jet algo
    if (jetalg == 0) then
      kd%njets = 0
      return
    end if

    call initRecomb_jet(npar)
#ifdef USEFASTJET
    call cluster_fastjet(npar)
#else
    if (jetalg == 4) then
      !> JADE
      call cluster_jade_jet()
    else
      call cluster_jet()
    end if
#endif
    call applyCuts_jet()
    call sort_jet()

    !------ WRITE OUT INTO KINDATA (need to do this before any matching)
    kd%RconeSq = RconeSq
    kd%njets = njets_jet
    do i=1,njets_jet
      kd%pjets(:,i) = jets(pt_sort_jet(i))%p(:)
    end do

    !------ treat ambiguous 2->2(+soft) configuration where pt2j1~pt2j2
    if ( njets_jet >= 2 ) then
      ! note: these are pt^2 [GeV^2]
      pt2j1 = jets(pt_sort_jet(1))%pt2
      pt2j2 = jets(pt_sort_jet(2))%pt2
      del_pt2_abs = abs( pt2j1-pt2j2 )
      del_pt2_rel = del_pt2_abs / abs( pt2j1+pt2j2 )
      if ( del_pt2_abs > threshold_del_pt2_abs .and.  &
      &    del_pt2_rel > threshold_del_pt2_rel ) then
        return
      end if

      ! !_______ DEBUG
      ! call randJJ(kd)
      ! return
      ! !_______ DEBUG

      npar_copy = kd%npar_from

      if ( kd%npar < npar_copy ) then
        !> we're in a subtraction term
        !>  => determine assignment from real-emission kinemtics
        !> make sure clustering was performed in the higher momentum set
        call recombine_jet(npar_copy)
        if ( kin(npar_copy)%njets == 0 ) then
          !> there is no reconstructed jet... randomise!
          if (del_pt2_abs == 0d0) call randJJ
        else if (   deltaJJ(kd%pjets(:,1),kin(npar_copy)%pjets(:,1)) &
        &         > deltaJJ(kd%pjets(:,2),kin(npar_copy)%pjets(:,1)) ) then
          !> align the leading jets
          call flipJJ
        end if
      else
        !> out of options
        !>  => randomly choose an assignemnt  
        if (del_pt2_abs == 0d0) call randJJ
      end if
      
    end if

    return

  contains ! functions local to `recombine_jet`

    !> funtion to measure "closeness" of jet momenta (for matching)
    pure double precision function deltaJJ(p1, p2)
      double precision, dimension(4), intent(in) :: p1, p2
      integer :: i

      deltaJJ = 0d0
      do i=1,4
        deltaJJ = deltaJJ + (p1(i)-p2(i))**2
      end do

      !> attention: only angles are no good!
      ! deltaJJ = v2_delta_R(p1,p2)

      ! deltaJJ = 0d0
      ! do i=1,3
      !   deltaJJ = deltaJJ + (p1(i)-p2(i))**2
      ! end do

      !deltaJJ = abs( dot(p1,p2) )

    end function deltaJJ

    !> routine to flip j1,j2
    subroutine flipJJ
      double precision :: pi
      integer :: i
      !print*, "FLIP"
      do i=1,4
        pi            = kd%pjets(i,1)
        kd%pjets(i,1) = kd%pjets(i,2)
        kd%pjets(i,2) = pi
      end do
    end subroutine flipJJ

    !> routine to randomise j1,j2
    subroutine randJJ
      double precision :: rnd
      call random_number(rnd)
      !print*, "RAND"
      !> default assignment is deterministic: j1=3, j2=4
      !> flip the two with a 50% probability
      if (rnd < 0.5d0) call flipJJ
    end subroutine randJJ

  end subroutine recombine_jet

end module Observables_mod
