// -*- C++ -*-
//
// Package:    TagAndProbe
// Class  :     FillEvtTree
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author: Nadia Adam 
//         Created:  Tue Oct 31 10:23:53 CST 2006
// $Id: FillEvtTree.cc,v 1.2 2008/04/14 15:51:29 neadam Exp $
//

// system include files
#include <iostream>

// user include files
#include "MuonAnalysis/TagAndProbe/interface/FillEvtTree.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FillEvtTree::FillEvtTree()
{
}

// FillEvtTree::FillEvtTree(const FillEvtTree& rhs)
// {
//    // do actual copying here;
// }

FillEvtTree::~FillEvtTree()
{
}

//
// assignment operators
//
// const FillEvtTree& FillEvtTree::operator=(const FillEvtTree& rhs)
// {
//   //An exception safe implementation is
//   FillEvtTree temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

void 
FillEvtTree::init(std::string filename, EvtTree* treePtr) {

   m_treePtr = treePtr;

   m_file = new TFile( filename.c_str(), "RECREATE" );

   m_tree = new TTree( "evttree", "Event Info", 1 );

   m_tree->Branch("run",        &(m_treePtr->run),        "run/I");
   m_tree->Branch("event",      &(m_treePtr->event),      "event/I");
   m_tree->Branch("num",        &(m_treePtr->num),        "num/F");
   m_tree->Branch("xs",         &(m_treePtr->xs),         "xs/F");

   m_tree->Branch("nl1",        &(m_treePtr->nl1),        "nl1/I");
   m_tree->Branch("nhlt",       &(m_treePtr->nhlt),       "nhlt/I");

   m_tree->Branch("nmcpart",    &(m_treePtr->nmcpart),    "nmcpart/I");
   m_tree->Branch("nmc",        &(m_treePtr->nmc),        "nmc[nmcpart]/I");
   m_tree->Branch("mc_bc",      &(m_treePtr->mc_bc),      "mc_bc[nmcpart][20]/I");
   m_tree->Branch("mc_pid",     &(m_treePtr->mc_pid),     "mc_pid[nmcpart][20]/I");
   m_tree->Branch("mc_pbc",     &(m_treePtr->mc_pbc),     "mc_pbc[nmcpart][20]/I");
   m_tree->Branch("mc_ppid",    &(m_treePtr->mc_ppid),    "mc_ppid[nmcpart][20]/I");
   m_tree->Branch("mc_mass",    &(m_treePtr->mc_mass),    "mc_mass[nmcpart][20]/F");
   m_tree->Branch("mc_p",       &(m_treePtr->mc_p),       "mc_p[nmcpart][20]/F");
   m_tree->Branch("mc_px",      &(m_treePtr->mc_px),      "mc_px[nmcpart][20]/F");
   m_tree->Branch("mc_py",      &(m_treePtr->mc_py),      "mc_py[nmcpart][20]/F");
   m_tree->Branch("mc_pz",      &(m_treePtr->mc_pz),      "mc_pz[nmcpart][20]/F");
   m_tree->Branch("mc_pt",      &(m_treePtr->mc_pt),      "mc_pt[nmcpart][20]/F");
   m_tree->Branch("mc_e",       &(m_treePtr->mc_e),       "mc_e[nmcpart][20]/F");
   m_tree->Branch("mc_phi",     &(m_treePtr->mc_phi),     "mc_phi[nmcpart][20]/F");
   m_tree->Branch("mc_eta",     &(m_treePtr->mc_eta),     "mc_eta[nmcpart][20]/F");

   // Tracks
   m_tree->Branch("nrtrk",       &(m_treePtr->nrtrk),       "nrtrk/I");
   m_tree->Branch("trk_type",    &(m_treePtr->trk_type),    "trk_type[nrtrk]/I");
   m_tree->Branch("trk_id",      &(m_treePtr->trk_id),      "trk_id[nrtrk]/I");
   m_tree->Branch("trk_chi2",    &(m_treePtr->trk_chi2),    "trk_chi2[nrtrk]/F");
   m_tree->Branch("trk_ndof",    &(m_treePtr->trk_ndof),    "trk_ndof[nrtrk]/F");
   m_tree->Branch("trk_p",       &(m_treePtr->trk_p),       "trk_p[nrtrk]/F");
   m_tree->Branch("trk_px",      &(m_treePtr->trk_px),      "trk_px[nrtrk]/F");
   m_tree->Branch("trk_py",      &(m_treePtr->trk_py),      "trk_py[nrtrk]/F");
   m_tree->Branch("trk_pz",      &(m_treePtr->trk_pz),      "trk_pz[nrtrk]/F");
   m_tree->Branch("trk_pt",      &(m_treePtr->trk_pt),      "trk_pt[nrtrk]/F");
   m_tree->Branch("trk_e",       &(m_treePtr->trk_e),       "trk_e[nrtrk]/F");
   m_tree->Branch("trk_q",       &(m_treePtr->trk_q),       "trk_q[nrtrk]/F");
   m_tree->Branch("trk_phi",     &(m_treePtr->trk_phi),     "trk_phi[nrtrk]/F");
   m_tree->Branch("trk_eta",     &(m_treePtr->trk_eta),     "trk_eta[nrtrk]/F");
   m_tree->Branch("trk_dxy",     &(m_treePtr->trk_dxy),     "trk_dxy[nrtrk]/F");
   m_tree->Branch("trk_d0",      &(m_treePtr->trk_d0),      "trk_d0[nrtrk]/F");
   m_tree->Branch("trk_dsz",     &(m_treePtr->trk_dsz),     "trk_dsz[nrtrk]/F");
   m_tree->Branch("trk_dz",      &(m_treePtr->trk_dz),      "trk_dz[nrtrk]/F");
   m_tree->Branch("trk_vx",      &(m_treePtr->trk_vx),      "trk_vx[nrtrk]/F");
   m_tree->Branch("trk_vy",      &(m_treePtr->trk_vy),      "trk_vy[nrtrk]/F");
   m_tree->Branch("trk_vz",      &(m_treePtr->trk_vz),      "trk_vz[nrtrk]/F");

   // Tag-probe pairs
   m_tree->Branch("nrtp",       &(m_treePtr->nrtp),       "nrtp/I");
   m_tree->Branch("tp_true",    &(m_treePtr->tp_true),    "tp_true[nrtp]/I");
   m_tree->Branch("tp_ppass",   &(m_treePtr->tp_ppass),   "tp_ppass[nrtp]/I");
   m_tree->Branch("tp_l1",      &(m_treePtr->tp_l1),      "tp_l1[nrtp]/I");
   m_tree->Branch("tp_hlt",     &(m_treePtr->tp_hlt),     "tp_hlt[nrtp]/I");
   m_tree->Branch("tp_p",       &(m_treePtr->tp_p),       "tp_p[nrtp]/F");
   m_tree->Branch("tp_px",      &(m_treePtr->tp_px),      "tp_px[nrtp]/F");
   m_tree->Branch("tp_py",      &(m_treePtr->tp_py),      "tp_py[nrtp]/F");
   m_tree->Branch("tp_pz",      &(m_treePtr->tp_pz),      "tp_pz[nrtp]/F");
   m_tree->Branch("tp_pt",      &(m_treePtr->tp_pt),      "tp_pt[nrtp]/F");
   m_tree->Branch("tp_e",       &(m_treePtr->tp_e),       "tp_e[nrtp]/F");
   m_tree->Branch("tp_et",      &(m_treePtr->tp_et),      "tp_et[nrtp]/F");
   m_tree->Branch("tp_mass",    &(m_treePtr->tp_mass),    "tp_mass[nrtp]/F");

   m_tree->Branch("tp_dp",       &(m_treePtr->tp_dp),       "tp_dp[nrtp][2]/F");
   m_tree->Branch("tp_dpx",      &(m_treePtr->tp_dpx),      "tp_dpx[nrtp][2]/F");
   m_tree->Branch("tp_dpy",      &(m_treePtr->tp_dpy),      "tp_dpy[nrtp][2]/F");
   m_tree->Branch("tp_dpz",      &(m_treePtr->tp_dpz),      "tp_dpz[nrtp][2]/F");
   m_tree->Branch("tp_dpt",      &(m_treePtr->tp_dpt),      "tp_dpt[nrtp][2]/F");
   m_tree->Branch("tp_de",       &(m_treePtr->tp_de),       "tp_de[nrtp][2]/F");
   m_tree->Branch("tp_det",      &(m_treePtr->tp_det),      "tp_det[nrtp][2]/F");
   m_tree->Branch("tp_dq",       &(m_treePtr->tp_dq),       "tp_dq[nrtp][2]/F");
   m_tree->Branch("tp_dphi",     &(m_treePtr->tp_dphi),     "tp_dphi[nrtp][2]/F");
   m_tree->Branch("tp_deta",     &(m_treePtr->tp_deta),     "tp_deta[nrtp][2]/F");

   // True Candidate (usually electron/muon) Efficiency
   m_tree->Branch("ncnd",       &(m_treePtr->ncnd),       "ncnd/I");
   m_tree->Branch("cnd_tag",    &(m_treePtr->cnd_tag),    "cnd_tag[ncnd]/I");
   m_tree->Branch("cnd_aprobe", &(m_treePtr->cnd_aprobe), "cnd_aprobe[ncnd]/I");
   m_tree->Branch("cnd_pprobe", &(m_treePtr->cnd_pprobe), "cnd_pprobe[ncnd]/I");
   m_tree->Branch("cnd_moid",   &(m_treePtr->cnd_moid),   "cnd_moid[ncnd]/I");
   m_tree->Branch("cnd_gmid",   &(m_treePtr->cnd_gmid),   "cnd_gmid[ncnd]/I");
   m_tree->Branch("cnd_p",      &(m_treePtr->cnd_p),      "cnd_p[ncnd]/F");
   m_tree->Branch("cnd_px",     &(m_treePtr->cnd_px),     "cnd_px[ncnd]/F");
   m_tree->Branch("cnd_py",     &(m_treePtr->cnd_py),     "cnd_py[ncnd]/F");
   m_tree->Branch("cnd_pz",     &(m_treePtr->cnd_pz),     "cnd_pz[ncnd]/F");
   m_tree->Branch("cnd_pt",     &(m_treePtr->cnd_pt),     "cnd_pt[ncnd]/F");
   m_tree->Branch("cnd_e",      &(m_treePtr->cnd_e),      "cnd_e[ncnd]/F");
   m_tree->Branch("cnd_et",     &(m_treePtr->cnd_et),     "cnd_et[ncnd]/F");
   m_tree->Branch("cnd_q",      &(m_treePtr->cnd_q),      "cnd_q[ncnd]/F");
   m_tree->Branch("cnd_phi",    &(m_treePtr->cnd_phi),    "cnd_phi[ncnd]/F");
   m_tree->Branch("cnd_eta",    &(m_treePtr->cnd_eta),    "cnd_eta[ncnd]/F");
   m_tree->Branch("cnd_rp",     &(m_treePtr->cnd_rp),     "cnd_rp[ncnd]/F");
   m_tree->Branch("cnd_rpx",    &(m_treePtr->cnd_rpx),    "cnd_rpx[ncnd]/F");
   m_tree->Branch("cnd_rpy",    &(m_treePtr->cnd_rpy),    "cnd_rpy[ncnd]/F");
   m_tree->Branch("cnd_rpz",    &(m_treePtr->cnd_rpz),    "cnd_rpz[ncnd]/F");
   m_tree->Branch("cnd_rpt",    &(m_treePtr->cnd_rpt),    "cnd_rpt[ncnd]/F");
   m_tree->Branch("cnd_re",     &(m_treePtr->cnd_re),     "cnd_re[ncnd]/F");
   m_tree->Branch("cnd_ret",    &(m_treePtr->cnd_ret),    "cnd_ret[ncnd]/F");
   m_tree->Branch("cnd_rq",     &(m_treePtr->cnd_rq),     "cnd_rq[ncnd]/F");
   m_tree->Branch("cnd_rphi",   &(m_treePtr->cnd_rphi),   "cnd_rphi[ncnd]/F");
   m_tree->Branch("cnd_reta",   &(m_treePtr->cnd_reta),   "cnd_reta[ncnd]/F");

   // "Vtx's"
   m_tree->Branch("nvtx",       &(m_treePtr->nvtx),       "nvtx/I");
   m_tree->Branch("vtx_IndexPrimaryVx", &(m_treePtr->vtx_IndexPrimaryVx), "vtx_IndexPrimaryVx/I");
   m_tree->Branch("vtx_PvnTracks",&(m_treePtr->vtx_PvnTracks),"vtx_PvnTracks/I");
   m_tree->Branch("vtx_Chi2",     &(m_treePtr->vtx_Chi2),     "vtx_Chi2/F");
   m_tree->Branch("vtx_Ndof",     &(m_treePtr->vtx_Ndof),     "vtx_Ndof/F");
   m_tree->Branch("vtx_NormalizedChi2",     &(m_treePtr->vtx_NormalizedChi2), "vtx_NormalizedChi2/F");
   m_tree->Branch("vtx_PVx",      &(m_treePtr->vtx_PVx),     "vtx_PVx/F");
   m_tree->Branch("vtx_PVy",      &(m_treePtr->vtx_PVy),     "vtx_PVy/F");
   m_tree->Branch("vtx_PVz",      &(m_treePtr->vtx_PVz),     "vtx_PVz/F");
   m_tree->Branch("vtx_PVdx",     &(m_treePtr->vtx_PVdx),     "vtx_PVdx/F");
   m_tree->Branch("vtx_PVdy",     &(m_treePtr->vtx_PVdy),     "vtx_PVdy/F");
   m_tree->Branch("vtx_PVdz",     &(m_treePtr->vtx_PVdz),     "vtx_PVdz/F");
   m_tree->Branch("vtx_PvPtsum",  &(m_treePtr->vtx_PvPtsum),  "vtx_PvPtsum/F");



   m_tree->Print();
}

EvtTree* 
FillEvtTree::getTreePtr() {
   return m_treePtr;
}

void 
FillEvtTree::fill() {
   m_tree->Fill();
}

void 
FillEvtTree::finalize() {
   m_file->Write();
   m_file->Close();
}

//
// member functions
//

//
// const member functions
//

//
// static member functions
//
