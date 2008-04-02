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
// $Id$
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

      // Information about the reconstructed muons
      int nrmu;
      int   mu_id[MAXNMU];
      int   mu_gbl[MAXNMU];
      int   mu_sta[MAXNMU];
      int   mu_trk[MAXNMU];
      float mu_p[MAXNMU];
      float mu_px[MAXNMU];
      float mu_py[MAXNMU];
      float mu_pz[MAXNMU];
      float mu_pt[MAXNMU];
      float mu_e[MAXNMU];
      float mu_et[MAXNMU];
      float mu_q[MAXNMU];
      float mu_eta[MAXNMU];
      float mu_phi[MAXNMU];

      // Information about the tag-probe muon pairs
      int nrtp;
      int   tp_type[MAXNPAIR];
      int   tp_true[MAXNPAIR];
      int   tp_ppass[MAXNPAIR];
      int   tp_l1[MAXNPAIR];
      int   tp_hlt[MAXNPAIR];
      float tp_mass[MAXNPAIR];
      float tp_p[MAXNPAIR];
      float tp_pt[MAXNPAIR];
      float tp_px[MAXNPAIR];
      float tp_py[MAXNPAIR];
      float tp_pz[MAXNPAIR];
      float tp_e[MAXNPAIR];
      float tp_et[MAXNPAIR];
      float tp_dp[MAXNPAIR][2];
      float tp_dpx[MAXNPAIR][2];
      float tp_dpy[MAXNPAIR][2];
      float tp_dpz[MAXNPAIR][2];
      float tp_dpt[MAXNPAIR][2];
      float tp_de[MAXNPAIR][2];
      float tp_det[MAXNPAIR][2];
      float tp_dq[MAXNPAIR][2];
      float tp_deta[MAXNPAIR][2];
      float tp_dphi[MAXNPAIR][2];

      // Information about the true muon efficiency
      // From all generator level muons in the event
      int ngmu;
      int   gmu_tag[MAXNMU];
      int   gmu_aprobe[MAXNMU];
      int   gmu_pprobe[MAXNMU];
      int   gmu_moid[MAXNMU];
      int   gmu_gmid[MAXNMU];
      float gmu_p[MAXNMU];
      float gmu_px[MAXNMU];
      float gmu_py[MAXNMU];
      float gmu_pz[MAXNMU];
      float gmu_pt[MAXNMU];
      float gmu_e[MAXNMU];
      float gmu_et[MAXNMU];
      float gmu_q[MAXNMU];
      float gmu_eta[MAXNMU];
      float gmu_phi[MAXNMU];
      float gmu_rp[MAXNMU];
      float gmu_rpx[MAXNMU];
      float gmu_rpy[MAXNMU];
      float gmu_rpz[MAXNMU];
      float gmu_rpt[MAXNMU];
      float gmu_re[MAXNMU];
      float gmu_ret[MAXNMU];
      float gmu_rq[MAXNMU];
      float gmu_reta[MAXNMU];
      float gmu_rphi[MAXNMU];

      // Information about muon pair "Z's"
      // z_type specifies the type of Z tag-probe etc
      int nrz;
      int   z_type[MAXNPAIR];
      int   z_true[MAXNPAIR];
      int   z_tid1[MAXNPAIR];
      int   z_tid2[MAXNPAIR];
      float z_mass[MAXNPAIR];
      float z_p[MAXNPAIR];
      float z_pt[MAXNPAIR];
      float z_px[MAXNPAIR];
      float z_py[MAXNPAIR];
      float z_pz[MAXNPAIR];
      float z_e[MAXNPAIR];
      int   z_vvalid[MAXNPAIR];
      float z_vd0[MAXNPAIR];
      float z_vx[MAXNPAIR];
      float z_vy[MAXNPAIR];
      float z_vz[MAXNPAIR];
      float z_vdx[MAXNPAIR];
      float z_vdy[MAXNPAIR];
      float z_vdz[MAXNPAIR];
      float z_vchi2[MAXNPAIR];
      float z_vndof[MAXNPAIR];
      float z_dp[MAXNPAIR][2];
      float z_dpx[MAXNPAIR][2];
      float z_dpy[MAXNPAIR][2];
      float z_dpz[MAXNPAIR][2];
      float z_dpt[MAXNPAIR][2];
      float z_de[MAXNPAIR][2];
      float z_det[MAXNPAIR][2];
      float z_dq[MAXNPAIR][2];
      float z_deta[MAXNPAIR][2];
      float z_dphi[MAXNPAIR][2];


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

 	 // Muons
	 nrmu = 0;
	 for( int ij = 0; ij<MAXNMU; ij++ )
	 {
	    mu_id[ij] = -999;
	    mu_gbl[ij] = -999;
	    mu_sta[ij] = -999;
	    mu_trk[ij] = -999;
	    mu_p[ij] = -999;
	    mu_px[ij] = -999;
	    mu_py[ij] = -999;
	    mu_pz[ij] = -999;
	    mu_pt[ij] = -999;
	    mu_e[ij] = -999;
	    mu_et[ij] = -999;
	    mu_q[ij] = -999;
	    mu_eta[ij] = -999;
	    mu_phi[ij] = -999;

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
	 ngmu = 0;
	 for( int ij = 0; ij<MAXNMU; ij++ )
	 {
	    gmu_tag[ij] = -999;
	    gmu_aprobe[ij] = -999;
	    gmu_pprobe[ij] = -999;
	    gmu_moid[ij] = -999;
	    gmu_gmid[ij] = -999;
	    gmu_p[ij] = -999;
	    gmu_px[ij] = -999;
	    gmu_py[ij] = -999;
	    gmu_pz[ij] = -999;
	    gmu_pt[ij] = -999;
	    gmu_e[ij] = -999;
	    gmu_et[ij] = -999;
	    gmu_q[ij] = -999;
	    gmu_eta[ij] = -999;
	    gmu_phi[ij] = -999;
	    gmu_rp[ij] = -999;
	    gmu_rpx[ij] = -999;
	    gmu_rpy[ij] = -999;
	    gmu_rpz[ij] = -999;
	    gmu_rpt[ij] = -999;
	    gmu_re[ij] = -999;
	    gmu_ret[ij] = -999;
	    gmu_rq[ij] = -999;
	    gmu_reta[ij] = -999;
	    gmu_rphi[ij] = -999;
	 }

	 nrz = 0;
	 for( int iz = 0; iz<MAXNPAIR; ++iz )
	 {
	    z_type[iz] = -999;
	    z_true[iz] = -999;
	    z_tid1[iz] = -999;
	    z_tid2[iz] = -999;
	    z_mass[iz] = -999;
	    z_p[iz] = -999;
	    z_pt[iz] = -999;
	    z_px[iz] = -999;
	    z_py[iz] = -999;
	    z_pz[iz] = -999;
	    z_e[iz] = -999;
	    z_vvalid[iz] = -999;
	    z_vd0[iz] = -999;
	    z_vx[iz] = -999;
	    z_vy[iz] = -999;
	    z_vz[iz] = -999;
	    z_vdx[iz] = -999;
	    z_vdy[iz] = -999;
	    z_vdz[iz] = -999;
	    z_vchi2[iz] = -999;
	    z_vndof[iz] = -999;

	    z_dp[iz][0] = -999;
	    z_dpx[iz][0] = -999;
	    z_dpy[iz][0] = -999;
	    z_dpz[iz][0] = -999;
	    z_dpt[iz][0] = -999;
	    z_de[iz][0] = -999;
	    z_det[iz][0] = -999;
	    z_dq[iz][0] = -999;
	    z_deta[iz][0] = -999;
	    z_dphi[iz][0] = -999;

	    z_dp[iz][1] = -999;
	    z_dpx[iz][1] = -999;
	    z_dpy[iz][1] = -999;
	    z_dpz[iz][1] = -999;
	    z_dpt[iz][1] = -999;
	    z_de[iz][1] = -999;
	    z_det[iz][1] = -999;
	    z_dq[iz][1] = -999;
	    z_deta[iz][1] = -999;
	    z_dphi[iz][1] = -999;
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
