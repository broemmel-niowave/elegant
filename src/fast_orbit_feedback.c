/*************************************************************************\
* Copyright (c) 2026 The University of Chicago, as Operator of Argonne
* National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution.
\*************************************************************************/

/* file: fast_orbit_feedback.c
 * purpose: fast orbit feedback (FOFB) simulation command.
 *
 * The command combines realistic beam dynamics (turn-by-turn tracking, with
 * synchrotron radiation) with a digital orbit-feedback loop:
 *
 *   for each FOFB iteration (run_control n_steps):
 *     track the (single, persistent) beam for run_control n_passes turns.
 *       - each turn every BPM's reading is pushed through its own digital IIR
 *         filter to form a running "tick" readout;
 *       - each turn every steering actuator's physically-applied kick is
 *         advanced one step of its own z-transform (IIR step response) toward
 *         the currently-held setpoint u[] -- so a setpoint that changes at the
 *         tick is NOT seen instantly by the beam, it ramps in turn-by-turn.
 *     at the end of the iteration the BPM ticks are collected into an orbit
 *     vector Q, projected into actuator space through the inverse response
 *     matrix T=-C^-1, and passed through a per-actuator-class PID controller
 *     (with anti-windup) to update the held setpoints u[].
 *
 * The steering-corrector core is complemented by an optional RF-frequency
 * actuator (include_rf_frequency, rf_Kp/Ki/Kd): the selected RFCA cavities form
 * one joint frequency knob appended to the x-plane response matrix as the
 * analytic dispersive column dx_i/df = -eta_i/(alpha_c*f0).  It closes the
 * steering<->energy loop that pure horizontal steering cannot, removing the slow
 * residual horizontal drift.  Which cavities are driven is selected through
 * &steering_element target=fofb, item=FREQ (else every RFCA is auto-discovered).
 *
 * A third, independent actuator is the RF-PHASE energy loop (include_rf_phase,
 * energy_Kp/Ki/Kd): a STANDALONE scalar loop (not a response-matrix column) that deduces
 * the common (mode-0) energy offset from the horizontal BPMs by dispersion projection
 * delta_est=sum(x_i*eta_i)/sum(eta_i^2) and nulls it by modulating the main-RF phase.
 * Run fast (step << synchrotron period) it damps the energy oscillation directly -- the
 * fast counterpart to the slow frequency radial loop.  subtract_dispersion feeds the x
 * steering solve the betatron-only residual so it can drop the long synchrotron boxcar.
 * Cavities are selected through &steering_element target=fofb, item=PHASE.
 *
 * See fast_orbit_feedback.nl for the namelist and the plan file for the full
 * design rationale.
 */

#include "mdb.h"
#include "track.h"
#include "correctDefs.h"
#include "fast_orbit_feedback.h"

/* nonzero while a fast_orbit_feedback command is running; read by do_tracking
   to enable the per-turn BPM store, the actuator z-transform hook, and by
   fofbStoreBpmTick/fofbUpdateActuators as a guard. */
long fofbActive = 0;

/* Total passes tracked across all FOFB steps (control->n_passes*control->n_steps).
   i_pass is continuous across steps (passOffset accumulates), so a WATCH in
   centroid/parameter mode accumulates one row per pass over the whole feedback run,
   not just n_passes.  dump_watch_parameters sizes its SDDS table from this when
   nonzero; it is 0 whenever FOFB is not running. */
long fofbTotalPasses = 0;

#define FOFB_MAX_FILTERS 100

double noise_value(double xamplitude, double xcutoff, long xerror_type);
long add_steer_type_to_lists(STEERING_LIST *SL, long plane, long type, char *item, double tweek, double limit,
                             LINE_LIST *beamline, RUN *run, long forceQuads);
long add_steer_elem_to_lists(STEERING_LIST *SL, long plane, char *name, char *item,
                             char *element_type, double tweek, double limit,
                             long start_occurence, long end_occurence, long occurence_step,
                             double s_start, double s_end,
                             LINE_LIST *beamline, RUN *run, long forceQuadsBends, long verbose);

/* FOFB-owned steering lists, populated by fofbAddSteerElem when a &steering_element
   command sets target="fast_orbit_feedback".  These are independent of the global
   &correct SLx/SLy.  fofbSL[0]=x correctors, fofbSL[1]=y correctors; fofbRFsl holds
   the RF cavities (item=FREQ).  setupPlane consumes them (auto-adding the standard
   dedicated correctors / all RFCAs when a class is not declared). */
static STEERING_LIST fofbSL[2];
static STEERING_LIST fofbRFsl;
static STEERING_LIST fofbRFPhasesl;
static long fofbSLDeclared[2] = {0, 0};
static long fofbRFDeclared = 0;
static long fofbRFPhaseDeclared = 0;

/* Copies of the closed-orbit-start namelist flags, set in setupFastOrbitFeedback and
   read by elegant.c (through the accessors below) to center/offset the beam on the
   computed closed orbit before tracking begins, mirroring the &track command. */
static long fofbCenterOnOrbitFlag = 0;
static long fofbCenterMomentumAlsoFlag = 1;
static long fofbOffsetByOrbitFlag = 0;
static long fofbOffsetMomentumAlsoFlag = 1;

long fofbCenterOnOrbit(void) { return fofbCenterOnOrbitFlag; }
long fofbCenterMomentumAlso(void) { return fofbCenterMomentumAlsoFlag; }
long fofbOffsetByOrbit(void) { return fofbOffsetByOrbitFlag; }
long fofbOffsetMomentumAlso(void) { return fofbOffsetMomentumAlsoFlag; }

/* ------------------------------------------------------------------ */
/* FOFB-owned steering-element intake: routed here from add_steering_element when a
   &steering_element command sets target="fast_orbit_feedback".  item=FREQ selects
   RF cavities for the joint RF-frequency knob; item=PHASE selects RF cavities for the
   RF-phase energy actuator (a separate scalar loop); anything else is a plane corrector. */

long fofbAddSteerElem(long plane, char *name, char *item, char *element_type, double tweek, double limit,
                      long start_occurence, long end_occurence, long occurence_step,
                      double s_start, double s_end, LINE_LIST *beamline, RUN *run, long verbose) {
  long found;
  if (item && strcmp(item, "FREQ") == 0) {
    found = add_steer_elem_to_lists(&fofbRFsl, plane, name, item, element_type, tweek, limit,
                                    start_occurence, end_occurence, occurence_step, s_start, s_end,
                                    beamline, run, 0, verbose);
    if (found)
      fofbRFDeclared = 1;
  } else if (item && strcmp(item, "PHASE") == 0) {
    found = add_steer_elem_to_lists(&fofbRFPhasesl, plane, name, item, element_type, tweek, limit,
                                    start_occurence, end_occurence, occurence_step, s_start, s_end,
                                    beamline, run, 0, verbose);
    if (found)
      fofbRFPhaseDeclared = 1;
  } else {
    long idx = (plane == 2) ? 1 : 0;
    found = add_steer_elem_to_lists(&fofbSL[idx], plane, name, item ? item : (idx ? "VKICK" : "HKICK"),
                                    element_type, tweek, limit,
                                    start_occurence, end_occurence, occurence_step, s_start, s_end,
                                    beamline, run, 0, verbose);
    if (found)
      fofbSLDeclared[idx] = 1;
  }
  return found;
}

/* release a STEERING_LIST's arrays (mirrors the reset branch of add_steer_elem_to_lists) */
static void freeSteeringList(STEERING_LIST *SL) {
  long i;
  if (SL->corr_param) {
    for (i = 0; i < SL->n_corr_types; i++)
      if (SL->corr_param[i])
        free(SL->corr_param[i]);
    free(SL->corr_param);
  }
  if (SL->elem)
    free(SL->elem);
  if (SL->corr_tweek)
    free(SL->corr_tweek);
  if (SL->corr_limit)
    free(SL->corr_limit);
  if (SL->param_offset)
    free(SL->param_offset);
  if (SL->param_index)
    free(SL->param_index);
  memset(SL, 0, sizeof(*SL));
}

/* ------------------------------------------------------------------ */
/* per-monitor runtime state, attached through ELEMENT_LIST->p_elem->fofbData */

typedef struct {
  IIRFILTER *xFilter, *yFilter; /* independent IIR banks (own state) for this BPM */
  long nxFilter, nyFilter;
  double xTick, yTick; /* running filtered readout (value at the last turn = the tick) */
  /* separate x-filter for the RF-frequency loop (x/dispersive plane only): the frequency
     integrator must average over the synchrotron period (notch f_s) while steering keeps a
     short bpm_filter_file; populated only when rf_bpm_filter_file is set. */
  IIRFILTER *xFilterRF;
  long nxFilterRF;
  double xTickRF; /* running RF-loop filtered readout; == raw reading when no RF bank */
} FOFB_BPM_DATA;

/* per-plane control state, all arrays indexed by corrector 0..CM.ncor-1 */

typedef struct {
  long coord;   /* 0 = x, 2 = y */
  long active;  /* nonzero if this plane has both monitors and correctors */
  CORMON_DATA CM;
  STEERING_LIST SL;
  double *base;        /* actuator base setpoint (parameter units, p_elem) */
  double *u;           /* currently-held actuator setpoint (parameter units) */
  double *Iacc;        /* PID integral accumulator (parameter units) */
  double *ePrev;       /* previous PID error (parameter units) */
  double *lastApplied; /* last physically-applied parameter value (for change test) */
  short *actPegged;    /* nonzero while actuator setpoint is pegged at the limit */
  long *kickOffset;    /* SL.param_offset[sl_index] for each corrector */
  IIRFILTER **actFilter; /* per-corrector actuator IIR bank (own state) */
  long *nActFilter;
  double Kp, Ki, Kd;   /* steering-class PID gains */
  /* RF-frequency actuator (x plane only): the selected RFCAs form one joint knob,
     registered as one extra actuator column at index rfIndex (-1 if inactive). */
  long rfIndex;        /* actuator index of the joint RF-frequency knob, or -1 */
  double rfKp, rfKi, rfKd; /* RF-class PID gains */
  double s_rf;         /* SVD conditioning scale applied to the RF response column */
  ELEMENT_LIST **rfElem; /* selected RFCA elements driven with a common detuning */
  double *rfBaseFreq;  /* each selected cavity's base frequency (Hz) */
  long *rfFreqOffset;  /* offset of FREQ within each cavity's p_elem */
  long nRF;            /* number of selected RF cavities */
  /* RF-phase energy actuator (x/dispersive plane only): a standalone scalar loop,
     NOT a response-matrix column.  The common (mode-0) energy offset is deduced from
     the BPMs by dispersion projection delta_est=sum(x_i*eta_i)/sum(eta_i^2) and nulled
     by modulating the main-RF phase (energy kick), without touching phase_fiducial. */
  long rfPhaseActive;         /* nonzero when include_rf_phase drives this plane */
  ELEMENT_LIST **rfPhaseElem; /* selected RFCA cavities whose PHASE is modulated */
  double *rfBasePhase;        /* each cavity's base PHASE (deg) */
  long *rfPhaseOffset;        /* offset of PHASE within each cavity's p_elem */
  long nRFPhase;              /* number of selected cavities */
  double energyKp, energyKi, energyKd; /* scalar PID gains on delta_est */
  double phaseScale;          /* deg of RF phase per unit delta_est (signed) */
  double phaseCmd;            /* current commanded phase offset (deg), held over a step */
  double lastPhaseApplied;    /* last phase offset physically written (deg) */
  double ePhaseAcc;           /* PID integral accumulator */
  double ePhasePrev;          /* previous delta_est (derivative term) */
  double deltaEst;            /* last dispersion-projected energy offset (for output) */
  short phasePegged;          /* nonzero while phaseCmd is clamped at rf_phase_limit */
} FOFB_PLANE;

static FOFB_PLANE planeData[2]; /* [0]=x, [1]=y */
static long fofbSetupDone = 0;

/* diagnostic output.  One SDDS page is written per feedback step and flushed to
   disk as the run proceeds (mirroring do_transport_analysis in analyze.c), so the
   file can be examined live while the code is still running rather than only at the
   end.  fofbHasDfrfColumn records whether the RF-frequency detuning column exists. */
static SDDS_DATASET SDDS_fofb;
static long fofbOutputActive = 0;
static long fofbHasDfrfColumn = 0;
static long fofbHasPhaseColumns = 0;

/* master-only guard for output under MPI */
#if USE_MPI
#  define FOFB_IS_MASTER (myid == 0)
#else
#  define FOFB_IS_MASTER 1
#endif

/* ------------------------------------------------------------------ */

static void **fofbDataPtr(ELEMENT_LIST *eptr) {
  switch (eptr->type) {
  case T_HMON:
    return &(((HMON *)eptr->p_elem)->fofbData);
  case T_VMON:
    return &(((VMON *)eptr->p_elem)->fofbData);
  case T_MONI:
    return &(((MONI *)eptr->p_elem)->fofbData);
  default:
    return NULL;
  }
}

/* allocate a fresh IIR bank and load it from file (own state); returns count */
static IIRFILTER *loadFilterBank(char *file, long *nFilter) {
  IIRFILTER *bank;
  *nFilter = 0;
  if (!file)
    return NULL;
  bank = tmalloc(sizeof(*bank) * FOFB_MAX_FILTERS);
  memset(bank, 0, sizeof(*bank) * FOFB_MAX_FILTERS);
  *nFilter = readIIRFilter(bank, FOFB_MAX_FILTERS, file);
  if (*nFilter <= 0) {
    free(bank);
    return NULL;
  }
  return bank;
}

/* zero the running state of an IIR bank without discarding coefficients */
static void resetFilterBank(IIRFILTER *bank, long nFilter) {
  long i, j;
  for (i = 0; i < nFilter; i++) {
    bank[i].iBuffer = 0;
    for (j = 0; j < bank[i].nTerms; j++)
      bank[i].xn[j] = bank[i].yn[j] = 0;
  }
}

/* ------------------------------------------------------------------ */
/* attach FOFB_BPM_DATA to every monitor used by either plane */

static void attachMonitorData(FOFB_PLANE *plane) {
  long i;
  for (i = 0; i < plane->CM.nmon; i++) {
    ELEMENT_LIST *moni = plane->CM.umoni[i];
    void **pp = fofbDataPtr(moni);
    FOFB_BPM_DATA *bd;
    if (!pp)
      continue;
    if (!(bd = *pp)) {
      bd = tmalloc(sizeof(*bd));
      memset(bd, 0, sizeof(*bd));
      *pp = bd;
    }
    if (plane->coord == 0 && !bd->xFilter)
      bd->xFilter = loadFilterBank(bpm_filter_file, &bd->nxFilter);
    if (plane->coord == 2 && !bd->yFilter)
      bd->yFilter = loadFilterBank(bpm_filter_file, &bd->nyFilter);
    /* dedicated BPM filter for the RF-frequency loop (x plane only) */
    if (plane->coord == 0 && include_rf_frequency && rf_bpm_filter_file && !bd->xFilterRF)
      bd->xFilterRF = loadFilterBank(rf_bpm_filter_file, &bd->nxFilterRF);
  }
}

/* ------------------------------------------------------------------ */
/* invert CM->C into CM->T = -C^-1, mirroring the COMPUTE_RESPONSE_INVERT block of
   compute_orbcor_matrices (correct.c).  Used after the RF column is appended. */

static void fofbInvertCM(CORMON_DATA *CM) {
  double conditionNumber;
  if (CM->T) {
    matrix_free(CM->T);
    CM->T = NULL;
  }
  if (CM->auto_limit_SVs && (CM->C->m < CM->C->n) && CM->remove_smallest_SVs < (long)(CM->C->n - CM->C->m))
    CM->remove_smallest_SVs = CM->C->n - CM->C->m;
  CM->T = matrix_invert(CM->C, CM->equalW ? NULL : CM->weight,
                        (int32_t)CM->keep_largest_SVs, (int32_t)CM->remove_smallest_SVs,
                        CM->minimum_SV_ratio, CM->Tikhonov_relative_alpha, CM->Tikhonov_n,
                        0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &conditionNumber);
  matrix_scmul(CM->T, -1);
  if (verbosity)
    printf("fast_orbit_feedback: RF-augmented correction matrix condition number %e\n", conditionNumber);
}

/* Append the joint RF-frequency actuator to the x-plane response matrix.  Builds the
   analytic dispersive column dx_i/df = -eta_i/(alpha_c*f0), scales it by s_rf for SVD
   conditioning, appends it to CM->C, and extends the CM per-corrector arrays by one
   (kick_coef=1/s_rf, sl_index=-1 as the RF sentinel).  Populates plane->rf* state. */

static void appendRFActuator(FOFB_PLANE *plane, RUN *run, LINE_LIST *beamline) {
  CORMON_DATA *CM = &plane->CM;
  long nSteer = CM->ncor, nmon = CM->nmon, i, j, k;
  double alpha_c, f0, *rfcol, sC = 0, sR = 0, rmsC, rmsR, s_rf;
  MAT *Caug;

  /* resolve the cavity list: declared entries, else auto-add every RFCA.
     add_steer_elem_to_lists uppercases its name/item/element_type in place, so it
     must receive mutable copies (not string literals in read-only memory). */
  if (!fofbRFDeclared || fofbRFsl.n_corr_types == 0) {
    char *anyName, *freqItem, *rfcaType;
    cp_str(&anyName, "*");
    cp_str(&freqItem, "FREQ");
    cp_str(&rfcaType, "RFCA");
    if (!add_steer_elem_to_lists(&fofbRFsl, 0, anyName, freqItem, rfcaType, 1e-6, 0,
                                 0, 0, 1, -1, -1, beamline, run, 0, verbosity > 1))
      bombElegant("fast_orbit_feedback: include_rf_frequency set but no RFCA cavities found", NULL);
    free(anyName);
    free(freqItem);
    free(rfcaType);
  }
  plane->nRF = fofbRFsl.n_corr_types;

  alpha_c = beamline->alpha[0];
  if (alpha_c == 0)
    bombElegant("fast_orbit_feedback: momentum compaction alpha_c is zero; cannot form RF-frequency response (need &twiss_output)", NULL);
  f0 = ((RFCA *)(fofbRFsl.elem[0]->p_elem))->freq;
  if (f0 <= 0)
    bombElegant("fast_orbit_feedback: RF cavity has non-positive FREQ", NULL);

  /* analytic per-Hz dispersive column at the monitors */
  rfcol = tmalloc(sizeof(*rfcol) * nmon);
  for (i = 0; i < nmon; i++)
    rfcol[i] = -CM->umoni[i]->twiss->etax / (alpha_c * f0);

  /* conditioning scale so the RF column carries weight comparable to the correctors */
  for (j = 0; j < nSteer; j++)
    for (i = 0; i < nmon; i++)
      sC += Mij(CM->C, i, j) * Mij(CM->C, i, j);
  for (i = 0; i < nmon; i++)
    sR += rfcol[i] * rfcol[i];
  rmsC = nSteer ? sqrt(sC / ((double)nmon * nSteer)) : 1.0;
  rmsR = sqrt(sR / nmon);
  if (rf_response_scale > 0)
    s_rf = rf_response_scale;
  else
    s_rf = (rmsR > 0 && rmsC > 0) ? rmsC / rmsR : 1.0;
  plane->s_rf = s_rf;

  /* augment C with the scaled RF column as the last column */
  Caug = matrix_get(nmon, nSteer + 1);
  for (j = 0; j < nSteer; j++)
    memcpy(Caug->me[j], CM->C->me[j], sizeof(double) * nmon);
  for (i = 0; i < nmon; i++)
    Mij(Caug, i, nSteer) = rfcol[i] * s_rf;
  free(rfcol);
  matrix_free(CM->C);
  CM->C = Caug;

  /* extend the per-corrector CM arrays by one for the RF entry */
  CM->kick_coef = trealloc(CM->kick_coef, sizeof(*CM->kick_coef) * (nSteer + 1));
  CM->sl_index = trealloc(CM->sl_index, sizeof(*CM->sl_index) * (nSteer + 1));
  CM->ucorr = trealloc(CM->ucorr, sizeof(*CM->ucorr) * (nSteer + 1));
  CM->pegged = trealloc(CM->pegged, sizeof(*CM->pegged) * (nSteer + 1));
  CM->kick_coef[nSteer] = 1.0 / s_rf; /* e = Mij(dK,rf)/kick_coef is then Hz */
  CM->sl_index[nSteer] = -1;          /* RF sentinel */
  CM->ucorr[nSteer] = fofbRFsl.elem[0];
  CM->pegged[nSteer] = 0;
  CM->ncor = nSteer + 1;
  plane->rfIndex = nSteer;

  /* cache the cavity elements / base frequencies for the turn-by-turn apply */
  plane->rfElem = tmalloc(sizeof(*plane->rfElem) * plane->nRF);
  plane->rfBaseFreq = tmalloc(sizeof(*plane->rfBaseFreq) * plane->nRF);
  plane->rfFreqOffset = tmalloc(sizeof(*plane->rfFreqOffset) * plane->nRF);
  for (k = 0; k < plane->nRF; k++) {
    plane->rfElem[k] = fofbRFsl.elem[k];
    plane->rfFreqOffset[k] = fofbRFsl.param_offset[k];
    plane->rfBaseFreq[k] = *((double *)(fofbRFsl.elem[k]->p_elem + fofbRFsl.param_offset[k]));
  }

  if (verbosity)
    printf("fast_orbit_feedback: RF-frequency actuator active over %ld cavity(ies), f0=%.6e Hz, alpha_c=%.6e, s_rf=%.3e\n",
           plane->nRF, f0, alpha_c, s_rf);
  fflush(stdout);
}

/* ------------------------------------------------------------------ */
/* Resolve and cache the cavity list driven by the RF-phase energy actuator, and
   compute the phase-per-delta scale.  This is a STANDALONE scalar loop -- it does NOT
   append a response-matrix column (unlike appendRFActuator).  Called only for the x
   (dispersive) plane when include_rf_phase is set.  Mirrors appendRFActuator's cavity
   discovery (declared entries, else auto-add every RFCA). */

static void setupPhaseActuator(FOFB_PLANE *plane, RUN *run, VARY *control, LINE_LIST *beamline) {
  long k;
  double voltScaled, cosPhiS, phiS;
  RFCA *rfca0;

  /* resolve the cavity list: declared entries, else auto-add every RFCA (mutable
     copies required -- add_steer_elem_to_lists uppercases its args in place). */
  if (!fofbRFPhaseDeclared || fofbRFPhasesl.n_corr_types == 0) {
    char *anyName, *phaseItem, *rfcaType;
    cp_str(&anyName, "*");
    cp_str(&phaseItem, "PHASE");
    cp_str(&rfcaType, "RFCA");
    if (!add_steer_elem_to_lists(&fofbRFPhasesl, 0, anyName, phaseItem, rfcaType, 1e-6, 0,
                                 0, 0, 1, -1, -1, beamline, run, 0, verbosity > 1))
      bombElegant("fast_orbit_feedback: include_rf_phase set but no RFCA cavities found", NULL);
    free(anyName);
    free(phaseItem);
    free(rfcaType);
  }
  plane->nRFPhase = fofbRFPhasesl.n_corr_types;

  plane->rfPhaseElem = tmalloc(sizeof(*plane->rfPhaseElem) * plane->nRFPhase);
  plane->rfBasePhase = tmalloc(sizeof(*plane->rfBasePhase) * plane->nRFPhase);
  plane->rfPhaseOffset = tmalloc(sizeof(*plane->rfPhaseOffset) * plane->nRFPhase);
  for (k = 0; k < plane->nRFPhase; k++) {
    plane->rfPhaseElem[k] = fofbRFPhasesl.elem[k];
    plane->rfPhaseOffset[k] = fofbRFPhasesl.param_offset[k];
    plane->rfBasePhase[k] = *((double *)(fofbRFPhasesl.elem[k]->p_elem + fofbRFPhasesl.param_offset[k]));
  }

  plane->energyKp = energy_Kp;
  plane->energyKi = energy_Ki;
  plane->energyKd = energy_Kd;
  plane->phaseCmd = 0;
  plane->lastPhaseApplied = 0; /* offset relative to base; base is written at pass 0 */
  plane->ePhaseAcc = plane->ePhasePrev = plane->deltaEst = 0;
  plane->phasePegged = 0;

  /* phaseScale: deg of RF phase per unit energy offset delta.  energy_response_scale,
     if given, is authoritative (signed, user-calibrated).  Otherwise form the analytic
     estimate -(180/pi)*pCentral/(n_passes*voltScaled*cos phi_s), where voltScaled matches
     the dgamma=volt*sin(phase) convention of simple_rfca.c (volt in gamma units) and
     pCentral is gamma0.  The analytic SIGN can be wrong (it depends on sign(cos phi_s) and
     above/below transition, and the true phi_s = phase_fiducial + rfca->phase, with
     phase_fiducial not yet known at setup), so it is a convenience default only. */
  rfca0 = (RFCA *)(plane->rfPhaseElem[0]->p_elem);
  if (energy_response_scale != 0) {
    plane->phaseScale = energy_response_scale;
  } else {
    voltScaled = rfca0->volt / (1e6 * particleMassMV * particleRelSign);
    phiS = rfca0->phase * PI / 180.0;
    cosPhiS = cos(phiS);
    if (voltScaled != 0 && cosPhiS != 0 && control->n_passes > 0)
      plane->phaseScale = -(180.0 / PI) * run->p_central /
                          ((double)control->n_passes * voltScaled * cosPhiS);
    else
      plane->phaseScale = 0;
    printWarning("fast_orbit_feedback: using an analytic estimate for the RF-phase energy "
                 "response scale.",
                 "Its sign may be wrong (it anti-damps if so).  Set energy_response_scale "
                 "explicitly (deg per unit delta) after calibrating with a short run.");
  }

  /* efficacy warning: a phase offset held across a step loses its energy-damping action
     once the beam begins to re-phase, i.e. once the step spans a sizable fraction of the
     synchrotron period.  Estimate nu_s from the standard formula and warn if the step is
     long compared with it.  Best-effort only (skipped silently if the estimate is not a
     sane positive number). */
  {
    double h, Eev, nus2, Tsynch;
    h = rfca0->freq * beamline->revolution_length / c_mks; /* beta ~ 1 */
    Eev = run->p_central * particleMassMV * 1e6;
    nus2 = (Eev > 0) ? beamline->alpha[0] * h * rfca0->volt * fabs(cos(rfca0->phase * PI / 180.0)) /
                         (PIx2 * Eev)
                     : 0;
    if (nus2 > 0) {
      Tsynch = 1.0 / sqrt(nus2);
      if (control->n_passes > Tsynch / 5.0)
        printWarning("fast_orbit_feedback: the RF-phase actuator step is long compared with "
                     "the synchrotron period.",
                     "n_passes exceeds ~T_synch/5; a held phase offset loses energy-damping "
                     "efficacy as the beam re-phases.  Use a faster loop (smaller n_passes).");
    }
  }

  plane->rfPhaseActive = 1;

  if (verbosity)
    printf("fast_orbit_feedback: RF-phase energy actuator active over %ld cavity(ies), "
           "phaseScale=%.6e deg/delta%s\n",
           plane->nRFPhase, plane->phaseScale,
           energy_response_scale != 0 ? " (user)" : " (analytic)");
  fflush(stdout);
}

/* ------------------------------------------------------------------ */

static void setupPlane(FOFB_PLANE *plane, long coord, RUN *run, VARY *control, LINE_LIST *beamline) {
  CORMON_DATA *CM;
  STEERING_LIST *SL;
  long i, found = 0, idx = (coord == 0) ? 0 : 1, rfActive;
  unsigned long flags;
  char *item;

  memset(plane, 0, sizeof(*plane));
  plane->coord = coord;
  plane->rfIndex = -1;
  CM = &plane->CM;

  /* select the corrector steering list: if this plane's correctors were declared
     through &steering_element target=fofb, use exactly those; otherwise auto-add the
     standard dedicated-corrector element types (mirrors correction_setup defaults,
     minus quad/bend steering). */
  if (fofbSLDeclared[idx]) {
    SL = &fofbSL[idx];
    found = SL->n_corr_types;
    if (verbosity)
      printf("fast_orbit_feedback: %s plane uses %ld user-declared corrector(s)\n",
             coord == 0 ? "x" : "y", found);
  } else {
    SL = &plane->SL;
    if (coord == 0) {
      cp_str(&item, "KICK");
      found += add_steer_type_to_lists(SL, 0, T_HCOR, item, 1e-6, corrector_limit, beamline, run, 0);
      found += add_steer_type_to_lists(SL, 0, T_EHCOR, item, 1e-6, corrector_limit, beamline, run, 0);
      cp_str(&item, "HKICK");
      found += add_steer_type_to_lists(SL, 0, T_HVCOR, item, 1e-6, corrector_limit, beamline, run, 0);
      found += add_steer_type_to_lists(SL, 0, T_EHVCOR, item, 1e-6, corrector_limit, beamline, run, 0);
    } else {
      cp_str(&item, "KICK");
      found += add_steer_type_to_lists(SL, 2, T_VCOR, item, 1e-6, corrector_limit, beamline, run, 0);
      found += add_steer_type_to_lists(SL, 2, T_EVCOR, item, 1e-6, corrector_limit, beamline, run, 0);
      cp_str(&item, "VKICK");
      found += add_steer_type_to_lists(SL, 2, T_HVCOR, item, 1e-6, corrector_limit, beamline, run, 0);
      found += add_steer_type_to_lists(SL, 2, T_EHVCOR, item, 1e-6, corrector_limit, beamline, run, 0);
    }
  }

  /* the RF-frequency actuator is a horizontal (dispersive) knob only */
  rfActive = (coord == 0 && include_rf_frequency);

  if (!found && !rfActive) {
    printWarning("No steering correctors found for a plane.",
                 coord == 0 ? "fast_orbit_feedback x plane inactive." : "fast_orbit_feedback y plane inactive.");
    plane->active = 0;
    return;
  }

  /* SV / inversion controls consumed by compute_orbcor_matrices */
  CM->nmon = CM->ncor = 0;
  CM->C = CM->T = NULL;
  /* fixed_length is the correction-time orbit model (unused: FOFB tracks the beam and
     never runs a closed-orbit correction step). fixed_length_matrix selects the orbit
     model used to build the response matrix (0 vary RF frequency, 1 vary energy,
     2 full 6-D); honored by both the analytic and computed-orbit builders below. */
  CM->fixed_length = 0;
  CM->fixed_length_matrix = fixed_length_matrix;
  CM->auto_limit_SVs = 1;
  CM->keep_largest_SVs = keep_largest_SVs;
  CM->remove_smallest_SVs = remove_smallest_SVs;
  CM->minimum_SV_ratio = minimum_SV_ratio;
  CM->Tikhonov_relative_alpha = Tikhonov_relative_alpha;
  CM->Tikhonov_n = Tikhonov_n;

  /* When the RF actuator is active we append its column before inverting, so defer
     the inversion; otherwise invert here exactly as increment 1. */
  flags = (rfActive ? 0 : COMPUTE_RESPONSE_INVERT) | (verbosity < 2 ? COMPUTE_RESPONSE_SILENT : 0);
  if (use_response_from_computed_orbits) {
    /* Build the response matrix by tweaking each corrector and differencing perturbed
       closed orbits (mirrors &correct use_response_from_computed_orbits). Reuse the
       serial per-plane builder compute_orbcor_matrices1 on all ranks: it is compiled
       unconditionally, matrix-based (single test particle) for closed_orbit_tracking_turns=0,
       and MPI-safe run redundantly -- the batched compute_orbcor_matrices1p does both
       planes in one collective call, incompatible with FOFB's per-plane + x-only RF-append
       flow. A local CORRECTION carries only the closed-orbit knobs and this plane's CM/SL. */
    CORRECTION corrTmp;
    memset(&corrTmp, 0, sizeof(corrTmp));
    corrTmp.clorb_accuracy = closed_orbit_accuracy;
    corrTmp.clorb_accuracy_requirement = closed_orbit_accuracy_requirement;
    corrTmp.clorb_iterations = closed_orbit_iterations;               /* struct field is double */
    corrTmp.clorb_iter_fraction = closed_orbit_iteration_fraction;
    corrTmp.clorb_fraction_multiplier = closed_orbit_fraction_multiplier;
    corrTmp.clorb_multiplier_interval = closed_orbit_multiplier_interval; /* struct field is short */
    corrTmp.clorb_track_for_orbit = closed_orbit_tracking_turns;      /* struct field is short */
    corrTmp.rpn_store_response_matrix = 0;
    /* shallow struct copy: SL's pointer members are read-only during response
       computation, and corrTmp is never freed. */
    if (coord == 0) {
      corrTmp.CMFx = CM;
      corrTmp.SLx = *SL;
    } else {
      corrTmp.CMFy = CM;
      corrTmp.SLy = *SL;
    }
    compute_orbcor_matrices1(&corrTmp, coord, run, beamline, flags, 0, NULL, NULL);
  } else {
    compute_orbcor_matrices(CM, SL, coord, run, beamline, flags, 0);
  }

  if (CM->nmon == 0) {
    plane->active = 0;
    return;
  }

  if (rfActive) {
    appendRFActuator(plane, run, beamline);
    fofbInvertCM(CM);
  }

  if (CM->ncor == 0 || CM->T == NULL) {
    plane->active = 0;
    return;
  }
  plane->active = 1;

  /* copy BPM noise controls into CM (used by fofbStoreBpmTick) */
  CM->bpm_noise = bpm_noise;
  CM->bpm_noise_cutoff = bpm_noise_cutoff;
  CM->bpm_noise_distribution = 1; /* gaussian; uniform via bpm_noise_distribution string below */
  if (bpm_noise_distribution && strncmp(bpm_noise_distribution, "uniform", 7) == 0)
    CM->bpm_noise_distribution = 2;

  /* per-actuator state */
  plane->base = tmalloc(sizeof(*plane->base) * CM->ncor);
  plane->u = tmalloc(sizeof(*plane->u) * CM->ncor);
  plane->Iacc = tmalloc(sizeof(*plane->Iacc) * CM->ncor);
  plane->ePrev = tmalloc(sizeof(*plane->ePrev) * CM->ncor);
  plane->lastApplied = tmalloc(sizeof(*plane->lastApplied) * CM->ncor);
  plane->actPegged = tmalloc(sizeof(*plane->actPegged) * CM->ncor);
  plane->kickOffset = tmalloc(sizeof(*plane->kickOffset) * CM->ncor);
  plane->actFilter = tmalloc(sizeof(*plane->actFilter) * CM->ncor);
  plane->nActFilter = tmalloc(sizeof(*plane->nActFilter) * CM->ncor);
  for (i = 0; i < CM->ncor; i++) {
    long sl_index = CM->sl_index[i];
    if (sl_index < 0) {
      /* joint RF-frequency knob: base setpoint is the reference cavity frequency */
      plane->kickOffset[i] = plane->rfFreqOffset[0];
      plane->base[i] = plane->rfBaseFreq[0];
      plane->actFilter[i] = loadFilterBank(rf_filter_file, &plane->nActFilter[i]);
    } else {
      ELEMENT_LIST *corr = CM->ucorr[i];
      plane->kickOffset[i] = SL->param_offset[sl_index];
      plane->base[i] = *((double *)(corr->p_elem + plane->kickOffset[i]));
      plane->actFilter[i] = loadFilterBank(steering_filter_file, &plane->nActFilter[i]);
    }
    plane->u[i] = plane->base[i];
    plane->lastApplied[i] = plane->base[i];
    plane->Iacc[i] = plane->ePrev[i] = 0;
    plane->actPegged[i] = 0;
  }

  plane->Kp = steering_Kp;
  plane->Ki = steering_Ki;
  plane->Kd = steering_Kd;
  plane->rfKp = rf_Kp;
  plane->rfKi = rf_Ki;
  plane->rfKd = rf_Kd;

  /* RF-phase energy actuator rides on the x (dispersive) plane's active BPM set */
  if (coord == 0 && include_rf_phase)
    setupPhaseActuator(plane, run, control, beamline);

  attachMonitorData(plane);

  if (verbosity)
    printf("fast_orbit_feedback: %s plane has %ld monitors, %ld correctors\n",
           coord == 0 ? "x" : "y", (long)CM->nmon, (long)CM->ncor);
  fflush(stdout);
}

/* ------------------------------------------------------------------ */

void setupFastOrbitFeedback(NAMELIST_TEXT *nltext, RUN *run, VARY *control, LINE_LIST *beamline) {
  /* process namelist */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&fast_orbit_feedback, nltext) == NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists)
    print_namelist(stdout, &fast_orbit_feedback);

  if (output_interval < 1)
    output_interval = 1;

  /* response-matrix orbit model (mirrors &correct fixed_length_matrix guards) */
  if (fixed_length_matrix < 0 || fixed_length_matrix > 2)
    bombElegant("fixed_length_matrix must be 0 (vary RF frequency), 1 (vary energy, momentum secant), or 2 (full 6-D closed orbit).", NULL);
  if (fixed_length_matrix == 1 && checkChangeT(beamline))
    bombElegant("change_t is nonzero on one or more RF cavities. This is incompatible with fixed_length_matrix=1 orbit computations.", NULL);
  if (fixed_length_matrix && include_rf_frequency)
    printWarning("fast_orbit_feedback: fixed_length_matrix is set together with include_rf_frequency.",
                 "The RF-frequency actuator moves the orbit by varying the RF frequency, which is inconsistent with a fixed-length response matrix. Use fixed_length_matrix=0 (variable-length, vary-RF-frequency model) or disable include_rf_frequency.");
  /* The RF-frequency column folds the dispersive signal eta*delta INTO the x solve;
     subtract_dispersion REMOVES exactly that signal from the x solve.  The two are the
     complementary halves of the energy channel -- use one or the other, never both. */
  if (subtract_dispersion && include_rf_frequency)
    bombElegant("fast_orbit_feedback: subtract_dispersion and include_rf_frequency are mutually exclusive "
                "(the frequency column needs the eta*delta signal that subtract_dispersion removes). "
                "Use the RF-phase energy loop (include_rf_phase) with subtract_dispersion for a fast energy "
                "channel, or include_rf_frequency alone for the slow radial loop.", NULL);
  /* rf_bpm_filter_file gives the RF-frequency loop its own BPM filter (to notch the
     synchrotron line that destabilizes the frequency integrator at any gain) -- it only
     has meaning for that loop. */
  if (rf_bpm_filter_file && !include_rf_frequency)
    bombElegant("fast_orbit_feedback: rf_bpm_filter_file is set but include_rf_frequency is 0. "
                "The dedicated RF-loop BPM filter only applies to the RF-frequency actuator; "
                "enable include_rf_frequency or remove rf_bpm_filter_file.", NULL);
  if (subtract_dispersion && !include_rf_phase)
    printWarning("fast_orbit_feedback: subtract_dispersion is set without include_rf_phase.",
                 "The dispersive energy offset is subtracted from the x steering solve but no actuator "
                 "corrects it; the common energy error will persist.  Enable include_rf_phase (or "
                 "include_rf_frequency without subtract_dispersion) to act on it.");
  if (use_response_from_computed_orbits) {
    if (closed_orbit_accuracy <= 0)
      bombElegant("closed_orbit_accuracy must be > 0", NULL);
    if (closed_orbit_accuracy_requirement <= 0)
      bombElegant("closed_orbit_accuracy_requirement must be > 0", NULL);
    if (closed_orbit_iteration_fraction <= 0 || closed_orbit_iteration_fraction > 1)
      bombElegant("closed_orbit_iteration_fraction must be on (0, 1]", NULL);
    if (closed_orbit_iterations <= 0)
      bombElegant("closed_orbit_iterations <= 0", NULL);
  }

  /* stash the closed-orbit-start flags for elegant.c to act on before tracking */
  fofbCenterOnOrbitFlag = center_on_orbit;
  fofbCenterMomentumAlsoFlag = center_momentum_also;
  fofbOffsetByOrbitFlag = offset_by_orbit;
  fofbOffsetMomentumAlsoFlag = offset_momentum_also;

#if USE_MPI
  printWarning("fast_orbit_feedback is validated for serial (elegant) running.",
               "Under Pelegant with a distributed beam, BPM readings are not reduced across ranks.");
#endif

  /* Warn if the deck defined &steering_element correctors but none target the feedback.
     Such commands default to target="correct" and populate the &correct SLx/SLy lists,
     which fast_orbit_feedback does not read -- so the feedback silently falls back to
     auto-discovering dedicated correctors (and ignores DQCOR-type or other user-chosen
     correctors).  This is the common footgun when adapting an orbit-correction deck. */
  if (!fofbSLDeclared[0] && !fofbSLDeclared[1] && !fofbRFDeclared &&
      nNonFOFBSteeringElementsSeen() > 0)
    printWarning("fast_orbit_feedback: steering_element commands were given but none target the feedback.",
                 "Add target=\"fofb\" to the &steering_element commands so their correctors drive the feedback; "
                 "otherwise those definitions are ignored and the feedback auto-discovers dedicated correctors.");

  setupPlane(&planeData[0], 0, run, control, beamline);
  setupPlane(&planeData[1], 2, run, control, beamline);

  if (!planeData[0].active && !planeData[1].active)
    bombElegant("fast_orbit_feedback: no active correction plane (no monitors and/or correctors found)", NULL);

  if (include_rf_phase && !planeData[0].rfPhaseActive)
    bombElegant("fast_orbit_feedback: include_rf_phase is set but the x plane is inactive (no monitors). "
                "The RF-phase energy actuator rides on the horizontal BPM set; declare H correctors "
                "(or include_rf_frequency) so the x plane is built.", NULL);

  /* diagnostic output setup (master only) */
  fofbOutputActive = 0;
  if (output && FOFB_IS_MASTER) {
    output = compose_filename(output, run->rootname);
    /* SDDS_DefineParameter/SDDS_DefineColumn return the new element's index (>=0)
       and -1 on error -- unlike the DefineSimple* forms (1/0), so the first
       element's index of 0 would make a "!" test spuriously fail; check "< 0". */
    if (!SDDS_InitializeOutput(&SDDS_fofb, SDDS_BINARY, 1, NULL, "fast orbit feedback", output) ||
        SDDS_DefineParameter(&SDDS_fofb, "nMonitorsX", NULL, NULL,
                             "Number of horizontal beam-position monitors read by the feedback",
                             NULL, SDDS_LONG, NULL) < 0 ||
        SDDS_DefineParameter(&SDDS_fofb, "nCorrectorsX", NULL, NULL,
                             "Number of horizontal actuators driven by the feedback",
                             NULL, SDDS_LONG, NULL) < 0 ||
        SDDS_DefineParameter(&SDDS_fofb, "nMonitorsY", NULL, NULL,
                             "Number of vertical beam-position monitors read by the feedback",
                             NULL, SDDS_LONG, NULL) < 0 ||
        SDDS_DefineParameter(&SDDS_fofb, "nCorrectorsY", NULL, NULL,
                             "Number of vertical steering correctors driven by the feedback",
                             NULL, SDDS_LONG, NULL) < 0 ||
        SDDS_DefineColumn(&SDDS_fofb, "Step", NULL, NULL,
                          "Feedback iteration index (0-based); one row per iteration",
                          NULL, SDDS_LONG, 0) < 0 ||
        SDDS_DefineColumn(&SDDS_fofb, "Pass", NULL, NULL,
                          "Cumulative tracking turns completed through the end of this iteration",
                          NULL, SDDS_LONG, 0) < 0 ||
        SDDS_DefineColumn(&SDDS_fofb, "xBpmRms", NULL, "m",
                          "RMS over the horizontal BPM errors of the filtered orbit reading measured at this iteration",
                          NULL, SDDS_DOUBLE, 0) < 0 ||
        SDDS_DefineColumn(&SDDS_fofb, "yBpmRms", NULL, "m",
                          "RMS over the vertical BPM errors of the filtered orbit reading measured at this iteration",
                          NULL, SDDS_DOUBLE, 0) < 0 ||
        SDDS_DefineColumn(&SDDS_fofb, "MaxCorrX", NULL, "rad",
                          "Largest horizontal steering-corrector kick excursion from its base value "
                          "at this iteration (the RF-frequency actuator is excluded)",
                          NULL, SDDS_DOUBLE, 0) < 0 ||
        SDDS_DefineColumn(&SDDS_fofb, "MaxCorrY", NULL, "rad",
                          "Largest vertical steering-corrector kick excursion from its base value "
                          "at this iteration",
                          NULL, SDDS_DOUBLE, 0) < 0 ||
        SDDS_DefineColumn(&SDDS_fofb, "PeggedX", NULL, NULL,
                          "Number of horizontal actuators saturated at their limit at this iteration",
                          NULL, SDDS_LONG, 0) < 0 ||
        SDDS_DefineColumn(&SDDS_fofb, "PeggedY", NULL, NULL,
                          "Number of vertical actuators saturated at their limit at this iteration",
                          NULL, SDDS_LONG, 0) < 0) {
      SDDS_SetError("Unable to set up fast_orbit_feedback output file");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    /* the RF-frequency detuning column exists only when the RF actuator is active,
       so include_rf_frequency=0 runs keep the increment-1 output layout unchanged */
    if (include_rf_frequency &&
        SDDS_DefineColumn(&SDDS_fofb, "DeltaFrf", NULL, "Hz",
                          "Applied RF-frequency detuning from the base frequency at this iteration "
                          "(one common fractional detuning shared by all selected cavities)",
                          NULL, SDDS_DOUBLE, 0) < 0) {
      SDDS_SetError("Unable to set up fast_orbit_feedback DeltaFrf column");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    /* RF-phase energy-actuator diagnostics exist only when include_rf_phase is set, so
       include_rf_phase=0 runs keep the output layout unchanged */
    if (include_rf_phase &&
        (SDDS_DefineColumn(&SDDS_fofb, "DeltaPhase", NULL, "deg",
                           "Applied RF-phase offset from the base phase at this iteration "
                           "(one common offset shared by all selected cavities), the energy actuator command",
                           NULL, SDDS_DOUBLE, 0) < 0 ||
         SDDS_DefineColumn(&SDDS_fofb, "deltaEnergyEst", NULL, NULL,
                           "Dispersion-projected common energy offset sum(x_i*eta_i)/sum(eta_i^2) "
                           "estimated from the horizontal BPMs at this iteration",
                           NULL, SDDS_DOUBLE, 0) < 0)) {
      SDDS_SetError("Unable to set up fast_orbit_feedback RF-phase columns");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WriteLayout(&SDDS_fofb)) {
      SDDS_SetError("Unable to write fast_orbit_feedback output layout");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    fofbOutputActive = 1;
    fofbHasDfrfColumn = include_rf_frequency ? 1 : 0;
    fofbHasPhaseColumns = include_rf_phase ? 1 : 0;
  }

  fofbSetupDone = 1;
}

/* ------------------------------------------------------------------ */
/* per-turn hook: advance each actuator's applied kick one z-transform step
   toward its held setpoint u[].  Called from do_tracking's pass loop.
   Returns nonzero if any element matrix was changed. */

long fofbUpdateActuators(LINE_LIST *beamline, RUN *run, long i_pass) {
  long ip, i, changed = 0;
  for (ip = 0; ip < 2; ip++) {
    FOFB_PLANE *plane = &planeData[ip];
    if (!plane->active)
      continue;
    for (i = 0; i < plane->CM.ncor; i++) {
      double delta, applied, *pval;
      delta = plane->u[i] - plane->base[i];
      if (plane->nActFilter[i] > 0)
        applied = plane->base[i] + applyIIRFilter(plane->actFilter[i], plane->nActFilter[i], delta);
      else
        applied = plane->u[i]; /* no step response -> setpoint applied immediately */
      if (applied == plane->lastApplied[i])
        continue;

      if (i == plane->rfIndex) {
        /* joint RF-frequency knob: every selected cavity gets the same fractional
           detuning applied/base.  The RFCA reads freq live from p_elem during
           tracking, so no matrix recompute is needed; but a bare freq write would
           leave phase_fiducial (= -omega_old*t0) stale, so rescale it by the
           frequency ratio to hold the fiducial time t0 fixed (subtlety (a)). */
        double fidRatio = (plane->lastApplied[i] != 0) ? applied / plane->lastApplied[i] : 1.0;
        double frac = (plane->base[i] != 0) ? applied / plane->base[i] : 1.0;
        long k;
        for (k = 0; k < plane->nRF; k++) {
          RFCA *rfca = (RFCA *)(plane->rfElem[k]->p_elem);
          *((double *)(plane->rfElem[k]->p_elem + plane->rfFreqOffset[k])) = plane->rfBaseFreq[k] * frac;
          if (rfca->fiducial_seen) {
            rfca->phase_fiducial *= fidRatio; /* preserves t0 = -phase_fiducial/omega */
            set_phase_reference(rfca->phase_reference, rfca->phase_fiducial);
          }
        }
        plane->lastApplied[i] = applied;
        continue; /* no compute_matrix / assert_element_links for the RF entry */
      }

      {
        ELEMENT_LIST *corr = plane->CM.ucorr[i];
        pval = (double *)(corr->p_elem + plane->kickOffset[i]);
        *pval = applied;
        plane->lastApplied[i] = applied;
        /* refresh the element matrix so matrix tracking sees the new kick */
        if (corr->matrix) {
          free_matrices(corr->matrix);
          tfree(corr->matrix);
          corr->matrix = NULL;
        }
        compute_matrix(corr, run, NULL);
        changed = 1;
      }
    }

    /* RF-phase energy actuator: a SEPARATE scalar knob (NOT in CM->ncor).  Write
       rfca->phase = rfBasePhase[k] + phaseCmd on every selected cavity.  Unlike the
       RF-frequency block this does NOT touch phase_fiducial (t0 = -phase_fiducial/omega
       is independent of rfca->phase) and does NOT set changed=1: the RFCA live-reads
       rfca->phase each pass (simple_rfca.c), so no compute_matrix is needed.  Pass 0
       writes the base phase (phaseCmd=0) so fiducialization happens at the base. */
    if (plane->rfPhaseActive && plane->phaseCmd != plane->lastPhaseApplied) {
      long k;
      for (k = 0; k < plane->nRFPhase; k++)
        *((double *)(plane->rfPhaseElem[k]->p_elem + plane->rfPhaseOffset[k])) =
          plane->rfBasePhase[k] + plane->phaseCmd;
      plane->lastPhaseApplied = plane->phaseCmd;
    }
  }
  if (changed && beamline->links)
    assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
  return changed;
}

/* ------------------------------------------------------------------ */
/* per-turn hook: push a BPM reading through its IIR filter to update the tick */

void fofbStoreBpmTick(ELEMENT_LIST *eptr, double xReading, double yReading) {
  void **pp = fofbDataPtr(eptr);
  FOFB_BPM_DATA *bd;
  if (!pp || !(bd = *pp))
    return;
  if (eptr->type == T_HMON || eptr->type == T_MONI) {
    if (bpm_noise)
      xReading += noise_value(bpm_noise, bpm_noise_cutoff,
                              (bpm_noise_distribution && strncmp(bpm_noise_distribution, "uniform", 7) == 0) ? 2 : 1);
    bd->xTick = bd->nxFilter > 0 ? applyIIRFilter(bd->xFilter, bd->nxFilter, xReading) : xReading;
    /* same reading (incl. noise) through the RF-loop's own filter, when configured */
    bd->xTickRF = bd->nxFilterRF > 0 ? applyIIRFilter(bd->xFilterRF, bd->nxFilterRF, xReading) : xReading;
  }
  if (eptr->type == T_VMON || eptr->type == T_MONI) {
    if (bpm_noise)
      yReading += noise_value(bpm_noise, bpm_noise_cutoff,
                              (bpm_noise_distribution && strncmp(bpm_noise_distribution, "uniform", 7) == 0) ? 2 : 1);
    bd->yTick = bd->nyFilter > 0 ? applyIIRFilter(bd->yFilter, bd->nyFilter, yReading) : yReading;
  }
}

/* ------------------------------------------------------------------ */
/* PID update in actuator space, run once per FOFB iteration */

static void updateSetpoints(FOFB_PLANE *plane, double *rmsOrbit, double *maxCorr, long *nPegged) {
  CORMON_DATA *CM = &plane->CM;
  MAT *Q, *dK;
  MAT *Q_rf = NULL, *dK_rf = NULL; /* separate solve feeding ONLY the RF-frequency row */
  long i;
  double sum2 = 0, maxc = 0;
  *nPegged = 0;

  /* Dispersion projection delta_est = sum(x_i*eta_i)/sum(eta_i^2) from the RAW x ticks,
     for the RF-phase energy loop and/or subtract_dispersion (x plane only).  sums->centroid
     (hence xTick) is already MPI_Allreduce'd, so delta_est is bit-identical on every rank. */
  if (plane->coord == 0 && (plane->rfPhaseActive || subtract_dispersion)) {
    double sxe = 0, see = 0;
    for (i = 0; i < CM->nmon; i++) {
      void **pp = fofbDataPtr(CM->umoni[i]);
      FOFB_BPM_DATA *bd = pp ? *pp : NULL;
      double eta = CM->umoni[i]->twiss ? CM->umoni[i]->twiss->etax : 0;
      double reading = bd ? bd->xTick : 0;
      sxe += reading * eta;
      see += eta * eta;
    }
    plane->deltaEst = (see > 0) ? sxe / see : 0;
  }

  /* collect BPM ticks into the orbit vector.  xBpmRms always reports the RAW orbit; when
     subtract_dispersion is set the CONTROL vector Q uses the betatron-only residual
     (x_i - eta_i*delta_est) so the x steering loop need not average out the dispersive
     synchrotron signal -- that channel is handled by the RF-phase energy actuator. */
  Q = matrix_get(CM->nmon, 1);
  for (i = 0; i < CM->nmon; i++) {
    void **pp = fofbDataPtr(CM->umoni[i]);
    FOFB_BPM_DATA *bd = pp ? *pp : NULL;
    double reading = 0, control;
    if (bd)
      reading = plane->coord == 0 ? bd->xTick : bd->yTick;
    control = reading;
    if (plane->coord == 0 && subtract_dispersion) {
      double eta = CM->umoni[i]->twiss ? CM->umoni[i]->twiss->etax : 0;
      control -= eta * plane->deltaEst;
    }
    Mij(Q, i, 0) = control;
    sum2 += reading * reading;
  }
  *rmsOrbit = CM->nmon ? sqrt(sum2 / CM->nmon) : 0;

  /* dK = T * Q = -C^-1 * Q : deadbeat correction demand (kick units) */
  dK = matrix_mult(CM->T, Q);

  /* Separate BPM filter for the RF-frequency loop (x plane only): build a second control
     vector from the RF-filtered ticks and solve with the SAME T, so the RF-frequency row's
     error comes from an f_s-averaged orbit while the steering rows keep the fast bpm_filter
     (see dK below).  Collapses to dK exactly when rf_bpm_filter_file is unset. */
  if (plane->coord == 0 && plane->rfIndex >= 0 && rf_bpm_filter_file) {
    Q_rf = matrix_get(CM->nmon, 1);
    for (i = 0; i < CM->nmon; i++) {
      void **pp = fofbDataPtr(CM->umoni[i]);
      FOFB_BPM_DATA *bd = pp ? *pp : NULL;
      Mij(Q_rf, i, 0) = bd ? bd->xTickRF : 0;
    }
    dK_rf = matrix_mult(CM->T, Q_rf);
  }

  for (i = 0; i < CM->ncor; i++) {
    double e, d, uProp, uNew, applied, Kp, Ki, Kd, limitVal, limitRef;
    long isRF = (i == plane->rfIndex);
    MAT *dKsrc = (isRF && dK_rf) ? dK_rf : dK;
    if (isRF) {
      /* RF-frequency class: gains rf_*, and the limit is a detuning about the base
         frequency (|f-f0| <= rf_frequency_limit), not an absolute value. */
      Kp = plane->rfKp;
      Ki = plane->rfKi;
      Kd = plane->rfKd;
      limitVal = rf_frequency_limit;
      limitRef = plane->base[i];
    } else {
      Kp = plane->Kp;
      Ki = plane->Ki;
      Kd = plane->Kd;
      limitVal = corrector_limit;
      limitRef = 0; /* correctors clamp on the absolute kick */
    }
    /* the RF-frequency row draws its error from the RF-filtered solve when a dedicated
       rf_bpm_filter_file is configured; steering rows always use the fast-filtered dK */
    e = Mij(dKsrc, i, 0) / CM->kick_coef[i]; /* param-space error (Hz for the RF entry) */
    d = e - plane->ePrev[i];
    if (!(anti_windup && plane->actPegged[i]))
      plane->Iacc[i] += Ki * e;
    uProp = plane->base[i] + Kp * e + Kd * d;
    uNew = uProp + plane->Iacc[i];
    if (limitVal > 0 && fabs(uNew - limitRef) > limitVal) {
      double lim = limitRef + (uNew > limitRef ? 1.0 : -1.0) * limitVal;
      if (anti_windup)
        plane->Iacc[i] = lim - uProp; /* back-calculate: freeze integrator at saturation */
      uNew = lim;
      plane->actPegged[i] = 1;
      (*nPegged)++;
    } else
      plane->actPegged[i] = 0;
    plane->u[i] = uNew;
    plane->ePrev[i] = e;
    if (!isRF) {
      /* MaxCorr* reports the steering excursion (rad); the RF detuning is a
         different quantity and is excluded here. */
      applied = fabs(plane->u[i] - plane->base[i]);
      if (applied > maxc)
        maxc = applied;
    }
  }
  *maxCorr = maxc;

  matrix_free(dK);
  matrix_free(Q);
  if (dK_rf)
    matrix_free(dK_rf);
  if (Q_rf)
    matrix_free(Q_rf);
}

/* ------------------------------------------------------------------ */
/* RF-phase energy loop: a scalar PID on the dispersion-projected energy offset
   delta_est (populated by updateSetpoints), producing the held phase command phaseCmd
   (deg).  Run once per FOFB iteration, on ALL ranks (the slave holding the particle must
   apply the phase), right after updateSetpoints for the x plane.  deltaEst and the gains
   are identical across ranks, so phaseCmd is bit-identical. */

static void updateEnergyLoop(FOFB_PLANE *plane) {
  double delta, c, cmd, prop;
  if (!plane->active || !plane->rfPhaseActive)
    return;
  delta = plane->deltaEst;
  /* PID in delta; energy_Ki defaults 0 (DC delta stays the frequency loop's job). */
  prop = plane->energyKp * delta + plane->energyKd * (delta - plane->ePhasePrev);
  if (!(anti_windup && plane->phasePegged))
    plane->ePhaseAcc += plane->energyKi * delta;
  c = prop + plane->ePhaseAcc;
  cmd = plane->phaseScale * c; /* sign of phaseScale sets damping vs anti-damping */
  if (rf_phase_limit > 0 && fabs(cmd) > rf_phase_limit) {
    double lim = (cmd > 0 ? 1.0 : -1.0) * rf_phase_limit;
    if (anti_windup && plane->phaseScale != 0)
      /* back-calculate the integrator so it holds at the clamp (phaseScale maps c->deg) */
      plane->ePhaseAcc = lim / plane->phaseScale - prop;
    cmd = lim;
    plane->phasePegged = 1;
  } else
    plane->phasePegged = 0;
  plane->phaseCmd = cmd;
  plane->ePhasePrev = delta;
}

/* ------------------------------------------------------------------ */

long doFastOrbitFeedback(RUN *run, VARY *control, LINE_LIST *beamline, BEAM *beam, OUTPUT_FILES *output_files) {
  long i_step, ip, passOffset = 0;
  double pCentral, finalCharge = 0;
  double clockBase = 0.0;
  unsigned long baseFlags, savedBeamlineFidFlag;

  if (!fofbSetupDone)
    bombElegant("fast_orbit_feedback: setup was not performed", NULL);

  pCentral = run->p_central;
  fofbActive = 1;
  fofbTotalPasses = control->n_passes * control->n_steps;

  /* Unlike &track (one do_tracking call per step), FOFB tracks ONE persistent beam
     across all n_steps, calling do_tracking once per step so tracking is continuous.
     do_tracking OR's beamline->fiducial_flag into its flags (do_tracking.c), and that
     flag carries RESET_RF_FOR_EACH_STEP by default (run_control's reset_rf_for_each_step,
     copied to the beamline in elegant.c). That would make do_tracking delete the phase
     references and re-fiducialize the RF at the start of EVERY step -- recomputing each
     cavity's fiducial time t0 from the mid-flight beam and shifting the RF energy kick at
     the step boundary (a spurious centroid-energy jump). The RF must be fiducialized once,
     on the step-0 beam, and then held for the life of the continuous beam. We therefore
     save the flag, let step 0 fiducialize exactly as a standalone first step would, and
     clear RESET_RF_FOR_EACH_STEP before the remaining steps; the flag is restored on exit. */
  savedBeamlineFidFlag = beamline->fiducial_flag;

  /* SUPPRESS_TWISS_UPDATE: the actuators change corrector kicks and rf frequency
     every turn, which invalidates the cached twiss and would otherwise force
     do_tracking to re-derive the twiss parameters on every pass (see
     do_tracking.c).  Those thin kicks move only the closed orbit, not the linear
     optics, so the recomputation is wasted work; suppress it automatically. */
  baseFlags = FINAL_SUMS_ONLY | ALLOW_MPI_ABORT_TRACKING | SUPPRESS_TWISS_UPDATE;

  for (i_step = 0; i_step < control->n_steps; i_step++) {
    unsigned long flags;
    double xrms = 0, yrms = 0, xcor = 0, ycor = 0;
    long xpeg = 0, ypeg = 0;

    /* optional reset of the BPM filter running state (never the actuator/Iacc state) */
    if (reset_filters_each_step) {
      for (ip = 0; ip < 2; ip++) {
        FOFB_PLANE *plane = &planeData[ip];
        long i;
        if (!plane->active)
          continue;
        for (i = 0; i < plane->CM.nmon; i++) {
          void **pp = fofbDataPtr(plane->CM.umoni[i]);
          FOFB_BPM_DATA *bd = pp ? *pp : NULL;
          if (!bd)
            continue;
          if (bd->nxFilter > 0)
            resetFilterBank(bd->xFilter, bd->nxFilter);
          if (bd->nyFilter > 0)
            resetFilterBank(bd->yFilter, bd->nyFilter);
          if (bd->nxFilterRF > 0)
            resetFilterBank(bd->xFilterRF, bd->nxFilterRF);
        }
      }
    }

    flags = baseFlags |
            (control->fiducial_flag &
             (FIRST_BEAM_IS_FIDUCIAL + FIDUCIAL_BEAM_SEEN + RESTRICT_FIDUCIALIZATION +
              LINEAR_CHROMATIC_MATRIX + LONGITUDINAL_RING_ONLY));

    /* FOFB tracks one persistent beam across n_steps do_tracking calls, but each
       do_tracking call resets the tracking clock accumulator (do_tracking.c
       resetTrackingClock()).  Standard &track makes a single do_tracking call, so its
       CHANGE_T/step_frequency macro-time offset accumulates continuously over all passes.
       To reproduce that continuity here, seed each step's base with the total macro time
       elapsed through the previous step; the upcoming do_tracking zeroes the accumulator
       before any consumer (e.g. applyElementModulations) reads trackingClockOffset(), so
       the offset stays continuous across step boundaries instead of sawtoothing to 0. */
    setTrackingClockBase(clockBase);

    do_tracking(beam, NULL, 0, NULL, beamline, &pCentral, beam->accepted, NULL, NULL, NULL,
                run, i_step, flags, control->n_passes, passOffset, NULL, NULL, &finalCharge, NULL, NULL);

    /* capture this step's accumulated macro time (base + per-step accum) so the next
       step continues from it */
    clockBase = trackingClockOffset();

    /* hold the RF fiducial across the remaining iterations */
    if (control->fiducial_flag & FIRST_BEAM_IS_FIDUCIAL) {
      control->fiducial_flag |= FIDUCIAL_BEAM_SEEN;
      beamline->fiducial_flag |= FIDUCIAL_BEAM_SEEN;
    }
    /* After step 0 has established the RF fiducial on the continuous beam, stop
       do_tracking from deleting phase references / re-fiducializing the RF at each
       subsequent step boundary (see the header comment above). */
    beamline->fiducial_flag &= ~RESET_RF_FOR_EACH_STEP;

    passOffset += control->n_passes;

    /* run the controller for each active plane */
    if (planeData[0].active) {
      updateSetpoints(&planeData[0], &xrms, &xcor, &xpeg);
      updateEnergyLoop(&planeData[0]); /* RF-phase energy loop (sets phaseCmd for next step) */
    }
    if (planeData[1].active)
      updateSetpoints(&planeData[1], &yrms, &ycor, &ypeg);

    if (verbosity)
      printf("fast_orbit_feedback step %ld: xBpmRms=%.6e m, yBpmRms=%.6e m (pegged x=%ld, y=%ld)\n",
             i_step, xrms, yrms, xpeg, ypeg);

    /* Write (and flush) one output page for this step so the file can be examined
       while the run is still in progress, rather than only at the end. */
    if (fofbOutputActive && FOFB_IS_MASTER && (i_step % output_interval == 0)) {
      if (!SDDS_StartPage(&SDDS_fofb, 1) ||
          !SDDS_SetParameters(&SDDS_fofb, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                              "nMonitorsX", (int32_t)planeData[0].CM.nmon, "nCorrectorsX", (int32_t)planeData[0].CM.ncor,
                              "nMonitorsY", (int32_t)planeData[1].CM.nmon, "nCorrectorsY", (int32_t)planeData[1].CM.ncor, NULL) ||
          !SDDS_SetRowValues(&SDDS_fofb, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, 0,
                             "Step", (int32_t)i_step, "Pass", (int32_t)passOffset,
                             "xBpmRms", xrms, "yBpmRms", yrms, "MaxCorrX", xcor, "MaxCorrY", ycor,
                             "PeggedX", (int32_t)xpeg, "PeggedY", (int32_t)ypeg, NULL)) {
        SDDS_SetError("Unable to write fast_orbit_feedback output page");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
      }
      if (fofbHasDfrfColumn) {
        double dfrf = (planeData[0].active && planeData[0].rfIndex >= 0)
                        ? planeData[0].u[planeData[0].rfIndex] - planeData[0].base[planeData[0].rfIndex]
                        : 0.0;
        if (!SDDS_SetRowValues(&SDDS_fofb, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, 0,
                               "DeltaFrf", dfrf, NULL)) {
          SDDS_SetError("Unable to set fast_orbit_feedback DeltaFrf value");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
        }
      }
      if (fofbHasPhaseColumns) {
        double dphase = (planeData[0].active && planeData[0].rfPhaseActive) ? planeData[0].phaseCmd : 0.0;
        double dest = (planeData[0].active && planeData[0].rfPhaseActive) ? planeData[0].deltaEst : 0.0;
        if (!SDDS_SetRowValues(&SDDS_fofb, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, 0,
                               "DeltaPhase", dphase, "deltaEnergyEst", dest, NULL)) {
          SDDS_SetError("Unable to set fast_orbit_feedback RF-phase values");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
        }
      }
      if (!SDDS_WritePage(&SDDS_fofb)) {
        SDDS_SetError("Unable to write fast_orbit_feedback output page");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
      }
      if (!inhibitFileSync)
        SDDS_DoFSync(&SDDS_fofb);
    }
  }

  fofbActive = 0;
  fofbTotalPasses = 0;
  beamline->fiducial_flag = savedBeamlineFidFlag;
  return 1;
}

/* ------------------------------------------------------------------ */

static void freePlane(FOFB_PLANE *plane) {
  long i;
  if (!plane->active) {
    if (plane->CM.C)
      matrix_free(plane->CM.C);
    if (plane->CM.T)
      matrix_free(plane->CM.T);
    plane->CM.C = plane->CM.T = NULL;
    return;
  }
  /* detach & free per-monitor data */
  for (i = 0; i < plane->CM.nmon; i++) {
    void **pp = fofbDataPtr(plane->CM.umoni[i]);
    FOFB_BPM_DATA *bd = pp ? *pp : NULL;
    if (!bd)
      continue;
    if (bd->xFilter) {
      freeIIRFilterMemory(bd->xFilter, bd->nxFilter);
      free(bd->xFilter);
    }
    if (bd->yFilter) {
      freeIIRFilterMemory(bd->yFilter, bd->nyFilter);
      free(bd->yFilter);
    }
    if (bd->xFilterRF) {
      freeIIRFilterMemory(bd->xFilterRF, bd->nxFilterRF);
      free(bd->xFilterRF);
    }
    free(bd);
    *pp = NULL;
  }
  for (i = 0; i < plane->CM.ncor; i++) {
    if (plane->actFilter[i]) {
      freeIIRFilterMemory(plane->actFilter[i], plane->nActFilter[i]);
      free(plane->actFilter[i]);
    }
  }
  free(plane->base);
  free(plane->u);
  free(plane->Iacc);
  free(plane->ePrev);
  free(plane->lastApplied);
  free(plane->actPegged);
  free(plane->kickOffset);
  free(plane->actFilter);
  free(plane->nActFilter);
  if (plane->rfElem)
    free(plane->rfElem);
  if (plane->rfBaseFreq)
    free(plane->rfBaseFreq);
  if (plane->rfFreqOffset)
    free(plane->rfFreqOffset);
  plane->rfElem = NULL;
  plane->rfBaseFreq = NULL;
  plane->rfFreqOffset = NULL;
  if (plane->rfPhaseElem)
    free(plane->rfPhaseElem);
  if (plane->rfBasePhase)
    free(plane->rfBasePhase);
  if (plane->rfPhaseOffset)
    free(plane->rfPhaseOffset);
  plane->rfPhaseElem = NULL;
  plane->rfBasePhase = NULL;
  plane->rfPhaseOffset = NULL;
  plane->rfPhaseActive = 0;
  if (plane->CM.C)
    matrix_free(plane->CM.C);
  if (plane->CM.T)
    matrix_free(plane->CM.T);
  plane->CM.C = plane->CM.T = NULL;
  plane->active = 0;
}

void finishFastOrbitFeedback(void) {
  fofbActive = 0;
  fofbTotalPasses = 0;
  if (!fofbSetupDone)
    return;

  if (fofbOutputActive && FOFB_IS_MASTER) {
    /* every step already wrote and flushed its own page during tracking; just
       close the file out here. */
    if (!SDDS_Terminate(&SDDS_fofb)) {
      SDDS_SetError("Unable to finish fast_orbit_feedback output file");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
  }
  fofbOutputActive = 0;
  fofbHasDfrfColumn = 0;
  fofbHasPhaseColumns = 0;

  /* restore each phase-actuator cavity to its base phase for a clean end state (the
     RFCA otherwise retains the last commanded offset); do this before freePlane frees
     the cached arrays.  The RFCA live-reads rfca->phase, so a plain write suffices. */
  {
    long ip, k;
    for (ip = 0; ip < 2; ip++) {
      FOFB_PLANE *plane = &planeData[ip];
      if (plane->active && plane->rfPhaseActive && plane->rfPhaseElem)
        for (k = 0; k < plane->nRFPhase; k++)
          *((double *)(plane->rfPhaseElem[k]->p_elem + plane->rfPhaseOffset[k])) = plane->rfBasePhase[k];
    }
  }

  freePlane(&planeData[0]);
  freePlane(&planeData[1]);

  /* release the FOFB-owned steering intake so a later &fast_orbit_feedback
     re-declares its actuators cleanly */
  freeSteeringList(&fofbSL[0]);
  freeSteeringList(&fofbSL[1]);
  freeSteeringList(&fofbRFsl);
  freeSteeringList(&fofbRFPhasesl);
  fofbSLDeclared[0] = fofbSLDeclared[1] = 0;
  fofbRFDeclared = 0;
  fofbRFPhaseDeclared = 0;

  fofbSetupDone = 0;
}
