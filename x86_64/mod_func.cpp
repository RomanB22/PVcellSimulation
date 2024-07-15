#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _ar_traub_reg(void);
extern void _cadad_reg(void);
extern void _cadyn_reg(void);
extern void _cagk_reg(void);
extern void _cal_mh_reg(void);
extern void _cal_mig_reg(void);
extern void _canin_reg(void);
extern void _can_mig_reg(void);
extern void _catcb_reg(void);
extern void _cat_mig_reg(void);
extern void _cat_traub_reg(void);
extern void _gabab_reg(void);
extern void _h_BS_reg(void);
extern void _HCN1_reg(void);
extern void _h_harnett_reg(void);
extern void _hin_reg(void);
extern void _h_kole_reg(void);
extern void _h_migliore_reg(void);
extern void _ican_sidi_reg(void);
extern void _IC_reg(void);
extern void _IKsin_reg(void);
extern void _kap_BS_reg(void);
extern void _kapcb_reg(void);
extern void _kapin_reg(void);
extern void _kBK_reg(void);
extern void _kctin_reg(void);
extern void _kdmc_BS_reg(void);
extern void _kdr_BS_reg(void);
extern void _kdrin_reg(void);
extern void _MyExp2SynBB_reg(void);
extern void _MyExp2SynNMDABB_reg(void);
extern void _nafx_reg(void);
extern void _nap_sidi_reg(void);
extern void _nax_BS_reg(void);
extern void _savedist_reg(void);
extern void _vecevent_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"mod/ar_traub.mod\"");
    fprintf(stderr, " \"mod/cadad.mod\"");
    fprintf(stderr, " \"mod/cadyn.mod\"");
    fprintf(stderr, " \"mod/cagk.mod\"");
    fprintf(stderr, " \"mod/cal_mh.mod\"");
    fprintf(stderr, " \"mod/cal_mig.mod\"");
    fprintf(stderr, " \"mod/canin.mod\"");
    fprintf(stderr, " \"mod/can_mig.mod\"");
    fprintf(stderr, " \"mod/catcb.mod\"");
    fprintf(stderr, " \"mod/cat_mig.mod\"");
    fprintf(stderr, " \"mod/cat_traub.mod\"");
    fprintf(stderr, " \"mod/gabab.mod\"");
    fprintf(stderr, " \"mod/h_BS.mod\"");
    fprintf(stderr, " \"mod/HCN1.mod\"");
    fprintf(stderr, " \"mod/h_harnett.mod\"");
    fprintf(stderr, " \"mod/hin.mod\"");
    fprintf(stderr, " \"mod/h_kole.mod\"");
    fprintf(stderr, " \"mod/h_migliore.mod\"");
    fprintf(stderr, " \"mod/ican_sidi.mod\"");
    fprintf(stderr, " \"mod/IC.mod\"");
    fprintf(stderr, " \"mod/IKsin.mod\"");
    fprintf(stderr, " \"mod/kap_BS.mod\"");
    fprintf(stderr, " \"mod/kapcb.mod\"");
    fprintf(stderr, " \"mod/kapin.mod\"");
    fprintf(stderr, " \"mod/kBK.mod\"");
    fprintf(stderr, " \"mod/kctin.mod\"");
    fprintf(stderr, " \"mod/kdmc_BS.mod\"");
    fprintf(stderr, " \"mod/kdr_BS.mod\"");
    fprintf(stderr, " \"mod/kdrin.mod\"");
    fprintf(stderr, " \"mod/MyExp2SynBB.mod\"");
    fprintf(stderr, " \"mod/MyExp2SynNMDABB.mod\"");
    fprintf(stderr, " \"mod/nafx.mod\"");
    fprintf(stderr, " \"mod/nap_sidi.mod\"");
    fprintf(stderr, " \"mod/nax_BS.mod\"");
    fprintf(stderr, " \"mod/savedist.mod\"");
    fprintf(stderr, " \"mod/vecevent.mod\"");
    fprintf(stderr, "\n");
  }
  _ar_traub_reg();
  _cadad_reg();
  _cadyn_reg();
  _cagk_reg();
  _cal_mh_reg();
  _cal_mig_reg();
  _canin_reg();
  _can_mig_reg();
  _catcb_reg();
  _cat_mig_reg();
  _cat_traub_reg();
  _gabab_reg();
  _h_BS_reg();
  _HCN1_reg();
  _h_harnett_reg();
  _hin_reg();
  _h_kole_reg();
  _h_migliore_reg();
  _ican_sidi_reg();
  _IC_reg();
  _IKsin_reg();
  _kap_BS_reg();
  _kapcb_reg();
  _kapin_reg();
  _kBK_reg();
  _kctin_reg();
  _kdmc_BS_reg();
  _kdr_BS_reg();
  _kdrin_reg();
  _MyExp2SynBB_reg();
  _MyExp2SynNMDABB_reg();
  _nafx_reg();
  _nap_sidi_reg();
  _nax_BS_reg();
  _savedist_reg();
  _vecevent_reg();
}

#if defined(__cplusplus)
}
#endif
