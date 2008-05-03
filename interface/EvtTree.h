#ifndef TagAndProbe_EvtTree_h
#define TagAndProbe_EvtTree_h
// -*- C++ -*-
//
// Package:     MuonAnalysis
// Class  :     EvtTree
// 
/**\class EvtTree EvtTree.h MuonAnalysis/TagAndProbe/interface/EvtTree.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author: Nadia Adam
//         Created:  Tue Oct 31 10:23:32 CST 2006
// $Id: EvtTree.h,v 1.2 2008/04/14 18:22:39 neadam Exp $
//

// system include files

// user include files

// forward declarations

class EvtTree
{

   public:

      enum { MAXNTRG = 500 };
      enum { MAXNJET = 100 };
      enum { MAXNMET = 100 };
      enum { MAXNMU  = 100 };
      enum { MAXNTRK = 1000 };
      enum { MAXNE  = 100 };
      enum { MAXNMCP = 20 };
      enum { MAXNMCT = 20 };
      enum { MAXNPAIR = 100 };
      enum { MAXNPHO = 100 };

      EvtTree();
      virtual ~EvtTree();

      // Public member data
      int run;
      int event;
      float num;
      float xs;

      // Trigger Info
      int nl1;
      int nhlt;

      // Information about the MC particles
      int nmcpart;
      int nmc[MAXNMCP];
      int mc_bc[MAXNMCP][MAXNMCT];
      int mc_pid[MAXNMCP][MAXNMCT];
      int mc_pbc[MAXNMCP][MAXNMCT];
      int mc_ppid[MAXNMCP][MAXNMCT];
      float mc_mass[MAXNMCP][MAXNMCT];
      float mc_p[MAXNMCP][MAXNMCT];
      float mc_px[MAXNMCP][MAXNMCT];
      float mc_py[MAXNMCP][MAXNMCT];
      float mc_pz[MAXNMCP][MAXNMCT];
      float mc_pt[MAXNMCP][MAXNMCT];
      float mc_e[MAXNMCP][MAXNMCT];
      float mc_eta[MAXNMCP][MAXNMCT];
      float mc_phi[MAXNMCP][MAXNMCT];


      // Information about the reconstructed tracks
      int nrtrk;
      int   trk_type[MAXNTRK];
      int   trk_id[MAXNTRK];
      float trk_chi2[MAXNTRK];
      float trk_ndof[MAXNTRK];
      float trk_p[MAXNTRK];
      float trk_px[MAXNTRK];
      float trk_py[MAXNTRK];
      float trk_pz[MAXNTRK];
      float trk_pt[MAXNTRK];
      float trk_e[MAXNTRK];
      float trk_q[MAXNTRK];
      float trk_eta[MAXNTRK];
      float trk_phi[MAXNTRK];
      float trk_dxy[MAXNTRK];
      float trk_d0[MAXNTRK];
      float trk_dsz[MAXNTRK];
      float trk_dz[MAXNTRK];
      float trk_vx[MAXNTRK];
      float trk_vy[MAXNTRK];
      float trk_vz[MAXNTRK];

      // Information about the tag-probe muon pairs
      int nrtp;                   /* Total number of TP pairs stored in this event. */
      int   tp_type[MAXNPAIR];    /* TP type, i.e. type of efficiency. */
      int   tp_true[MAXNPAIR];    /* Is the TP pair a true Z? */
      int   tp_ppass[MAXNPAIR];   /* Did the probe pass eff criteria? (if relevant) */
      int   tp_l1[MAXNPAIR];      /* Did the tag cause a L1 trigger? */
      int   tp_hlt[MAXNPAIR];     /* Did the tag cause a HLT trigger? */ 
      float tp_mass[MAXNPAIR];    /* Invariant mass of the TP pair. */
      float tp_p[MAXNPAIR];       /* Momentum of the TP pair. */
      float tp_pt[MAXNPAIR];      /* Transverse momentum of the TP pair. */
      float tp_px[MAXNPAIR];      /* X momentum of the TP pair. */
      float tp_py[MAXNPAIR];      /* Y momentum of the TP pair. */
      float tp_pz[MAXNPAIR];      /* Z momentum of the TP pair. */
      float tp_e[MAXNPAIR];       /* Energy of the TP pair. */
      float tp_et[MAXNPAIR];      /* Transverse energy of the TP pair. */
      float tp_dp[MAXNPAIR][2];   /* Daughter muon momentum: [0]=Tag, [1]=Probe */
      float tp_dpx[MAXNPAIR][2];  /* Daughter muon X momentum: [0]=Tag, [1]=Probe */
      float tp_dpy[MAXNPAIR][2];  /* Daughter muon Y momentum: [0]=Tag, [1]=Probe */
      float tp_dpz[MAXNPAIR][2];  /* Daughter muon Z momentum: [0]=Tag, [1]=Probe */
      float tp_dpt[MAXNPAIR][2];  /* Daughter muon transverse momentum: [0]=Tag, [1]=Probe */
      float tp_de[MAXNPAIR][2];   /* Daughter muon energy: [0]=Tag, [1]=Probe */
      float tp_det[MAXNPAIR][2];  /* Daughter muon transverse energy: [0]=Tag, [1]=Probe */
      float tp_dq[MAXNPAIR][2];   /* Daughter muon charge: [0]=Tag, [1]=Probe */
      float tp_deta[MAXNPAIR][2]; /* Daughter muon eta: [0]=Tag, [1]=Probe */
      float tp_dphi[MAXNPAIR][2]; /* Daughter muon phi: [0]=Tag, [1]=Probe */

      // Information about the true muon/electron efficiency
      // From all generator level muons in the event
      int ncnd;
      int   cnd_type[MAXNMU];
      int   cnd_tag[MAXNMU];
      int   cnd_aprobe[MAXNMU];
      int   cnd_pprobe[MAXNMU];
      int   cnd_moid[MAXNMU];
      int   cnd_gmid[MAXNMU];
      float cnd_p[MAXNMU];
      float cnd_px[MAXNMU];
      float cnd_py[MAXNMU];
      float cnd_pz[MAXNMU];
      float cnd_pt[MAXNMU];
      float cnd_e[MAXNMU];
      float cnd_et[MAXNMU];
      float cnd_q[MAXNMU];
      float cnd_eta[MAXNMU];
      float cnd_phi[MAXNMU];
      float cnd_rp[MAXNMU];
      float cnd_rpx[MAXNMU];
      float cnd_rpy[MAXNMU];
      float cnd_rpz[MAXNMU];
      float cnd_rpt[MAXNMU];
      float cnd_re[MAXNMU];
      float cnd_ret[MAXNMU];
      float cnd_rq[MAXNMU];
      float cnd_reta[MAXNMU];
      float cnd_rphi[MAXNMU];


      // Information about primary "Vtx's"
      int nvtx;
      int   vtx_IndexPrimaryVx;
      int   vtx_PvnTracks;
      float vtx_Chi2;
      float vtx_Ndof;
      float vtx_NormalizedChi2;
      float vtx_PVx;
      float vtx_PVy;
      float vtx_PVz;
      float vtx_PVdx;
      float vtx_PVdy;
      float vtx_PVdz;
      float vtx_PvPtsum;


      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void init()
      {
	 run   = -999;
	 event = -999;
	 xs = -999;
	 num = -999;

	 nl1 = -999;
	 nhlt = -999;

	 nmcpart = 0;
	 for( int imc=0; imc<MAXNMCP; imc++ )
	 {
	    nmc[imc] = 0;

	    for( int jmc=0; jmc<MAXNMCT; jmc++ )
	    {
	       mc_bc[imc][jmc] = -999;
	       mc_pid[imc][jmc] = -999;
	       mc_pbc[imc][jmc] = -999;
	       mc_ppid[imc][jmc] = -999;
	       mc_bc[imc][jmc] = -999;
	       mc_mass[imc][jmc] = -999;
	       mc_p[imc][jmc] = -999;
	       mc_px[imc][jmc] = -999;
	       mc_py[imc][jmc] = -999;
	       mc_pz[imc][jmc] = -999;
	       mc_pt[imc][jmc] = -999;
	       mc_e[imc][jmc] = -999;
	       mc_eta[imc][jmc] = -999;
	       mc_phi[imc][jmc] = -999;
	    }
	 }


	 // Tracks
	 nrtrk = 0;
	 for( int ij = 0; ij<MAXNTRK; ij++ )
	 {
	    trk_type[ij] = -999;
	    trk_id[ij] = -999;
	    trk_chi2[ij] = -999;
	    trk_ndof[ij] = 1;
	    trk_p[ij] = -999;
	    trk_px[ij] = -999;
	    trk_py[ij] = -999;
	    trk_pz[ij] = -999;
	    trk_pt[ij] = -999;
	    trk_e[ij] = -999;
	    trk_q[ij] = -999;
	    trk_eta[ij] = -999;
	    trk_phi[ij] = -999;
	    trk_dxy[ij] = -999;
	    trk_d0[ij] = -999;
	    trk_dsz[ij] = -999;
	    trk_dz[ij] = -999;
	    trk_vx[ij] = -999;
	    trk_vy[ij] = -999;
	    trk_vz[ij] = -999;

	 }

	 // Tag-probe muon pairs
	 nrtp = 0;
	 for( int ij=0; ij<MAXNPAIR; ++ij )
	 {
	    tp_type[ij] = -999;
	    tp_true[ij] = -999;
	    tp_ppass[ij] = -999;
	    tp_l1[ij] = -999;
	    tp_hlt[ij] = -999;
	    tp_mass[ij] = -999;
	    tp_p[ij] = -999;
	    tp_pt[ij] = -999;
	    tp_px[ij] = -999;
	    tp_py[ij] = -999;
	    tp_pz[ij] = -999;
	    tp_e[ij] = -999;
	    tp_et[ij] = -999;

	    tp_dp[ij][0] = -999;
	    tp_dpx[ij][0] = -999;
	    tp_dpy[ij][0] = -999;
	    tp_dpz[ij][0] = -999;
	    tp_dpt[ij][0] = -999;
	    tp_de[ij][0] = -999;
	    tp_det[ij][0] = -999;
	    tp_dq[ij][0] = -999;
	    tp_deta[ij][0] = -999;
	    tp_dphi[ij][0] = -999;

	    tp_dp[ij][1] = -999;
	    tp_dpx[ij][1] = -999;
	    tp_dpy[ij][1] = -999;
	    tp_dpz[ij][1] = -999;
	    tp_dpt[ij][1] = -999;
	    tp_de[ij][1] = -999;
	    tp_det[ij][1] = -999;
	    tp_dq[ij][1] = -999;
	    tp_deta[ij][1] = -999;
	    tp_dphi[ij][1] = -999;
	 }

 	 // Muon Truth and Matching
	 ncnd = 0;
	 for( int ij = 0; ij<MAXNMU; ij++ )
	 {
	    cnd_type[ij] = -999;
	    cnd_tag[ij] = -999;
	    cnd_aprobe[ij] = -999;
	    cnd_pprobe[ij] = -999;
	    cnd_moid[ij] = -999;
	    cnd_gmid[ij] = -999;
	    cnd_p[ij] = -999;
	    cnd_px[ij] = -999;
	    cnd_py[ij] = -999;
	    cnd_pz[ij] = -999;
	    cnd_pt[ij] = -999;
	    cnd_e[ij] = -999;
	    cnd_et[ij] = -999;
	    cnd_q[ij] = -999;
	    cnd_eta[ij] = -999;
	    cnd_phi[ij] = -999;
	    cnd_rp[ij] = -999;
	    cnd_rpx[ij] = -999;
	    cnd_rpy[ij] = -999;
	    cnd_rpz[ij] = -999;
	    cnd_rpt[ij] = -999;
	    cnd_re[ij] = -999;
	    cnd_ret[ij] = -999;
	    cnd_rq[ij] = -999;
	    cnd_reta[ij] = -999;
	    cnd_rphi[ij] = -999;
	 }
	 
	 // reset "Vtx's to zero"
	 nvtx=0;
	 vtx_IndexPrimaryVx=-99;
	 vtx_PvnTracks = -99;
	 vtx_Chi2=-99.;
	 vtx_Ndof=-99.;
	 vtx_NormalizedChi2=-99.;
	 vtx_PVx=-99.;
	 vtx_PVy=-99.;
	 vtx_PVz=-99.;
	 vtx_PVdx=-99.;
	 vtx_PVdy=-99.;
	 vtx_PVdz=-99.;
	 vtx_PvPtsum=-99.;

      }

      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------

   private:
      EvtTree(const EvtTree&); // stop default

      const EvtTree& operator=(const EvtTree&); // stop default

      // ---------- member data --------------------------------

};


#endif
