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
// $Id: FillEvtTree.cc,v 1.1 2008/04/02 19:59:03 neadam Exp $
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

   // Muons
   m_tree->Branch("nrmu",       &(m_treePtr->nrmu),       "nrmu/I");
   m_tree->Branch("mu_id",      &(m_treePtr->mu_id),      "mu_id[nrmu]/I");
   m_tree->Branch("mu_gbl",     &(m_treePtr->mu_gbl),     "mu_gbl[nrmu]/I");
   m_tree->Branch("mu_sta",     &(m_treePtr->mu_sta),     "mu_sta[nrmu]/I");
   m_tree->Branch("mu_trk",     &(m_treePtr->mu_trk),     "mu_trk[nrmu]/I");
   m_tree->Branch("mu_p",       &(m_treePtr->mu_p),       "mu_p[nrmu]/F");
   m_tree->Branch("mu_px",      &(m_treePtr->mu_px),      "mu_px[nrmu]/F");
   m_tree->Branch("mu_py",      &(m_treePtr->mu_py),      "mu_py[nrmu]/F");
   m_tree->Branch("mu_pz",      &(m_treePtr->mu_pz),      "mu_pz[nrmu]/F");
   m_tree->Branch("mu_pt",      &(m_treePtr->mu_pt),      "mu_pt[nrmu]/F");
   m_tree->Branch("mu_e",       &(m_treePtr->mu_e),       "mu_e[nrmu]/F");
   m_tree->Branch("mu_et",      &(m_treePtr->mu_et),      "mu_et[nrmu]/F");
   m_tree->Branch("mu_q",       &(m_treePtr->mu_q),       "mu_q[nrmu]/F");
   m_tree->Branch("mu_phi",     &(m_treePtr->mu_phi),     "mu_phi[nrmu]/F");
   m_tree->Branch("mu_eta",     &(m_treePtr->mu_eta),     "mu_eta[nrmu]/F");

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

   // Muon Efficiency
   m_tree->Branch("ngmu",       &(m_treePtr->ngmu),       "ngmu/I");
   m_tree->Branch("gmu_tag",    &(m_treePtr->gmu_tag),    "gmu_tag[ngmu]/I");
   m_tree->Branch("gmu_aprobe", &(m_treePtr->gmu_aprobe), "gmu_aprobe[ngmu]/I");
   m_tree->Branch("gmu_pprobe", &(m_treePtr->gmu_pprobe), "gmu_pprobe[ngmu]/I");
   m_tree->Branch("gmu_moid",   &(m_treePtr->gmu_moid),   "gmu_moid[ngmu]/I");
   m_tree->Branch("gmu_gmid",   &(m_treePtr->gmu_gmid),   "gmu_gmid[ngmu]/I");
   m_tree->Branch("gmu_p",      &(m_treePtr->gmu_p),      "gmu_p[ngmu]/F");
   m_tree->Branch("gmu_px",     &(m_treePtr->gmu_px),     "gmu_px[ngmu]/F");
   m_tree->Branch("gmu_py",     &(m_treePtr->gmu_py),     "gmu_py[ngmu]/F");
   m_tree->Branch("gmu_pz",     &(m_treePtr->gmu_pz),     "gmu_pz[ngmu]/F");
   m_tree->Branch("gmu_pt",     &(m_treePtr->gmu_pt),     "gmu_pt[ngmu]/F");
   m_tree->Branch("gmu_e",      &(m_treePtr->gmu_e),      "gmu_e[ngmu]/F");
   m_tree->Branch("gmu_et",     &(m_treePtr->gmu_et),     "gmu_et[ngmu]/F");
   m_tree->Branch("gmu_q",      &(m_treePtr->gmu_q),      "gmu_q[ngmu]/F");
   m_tree->Branch("gmu_phi",    &(m_treePtr->gmu_phi),    "gmu_phi[ngmu]/F");
   m_tree->Branch("gmu_eta",    &(m_treePtr->gmu_eta),    "gmu_eta[ngmu]/F");
   m_tree->Branch("gmu_rp",     &(m_treePtr->gmu_rp),     "gmu_rp[ngmu]/F");
   m_tree->Branch("gmu_rpx",    &(m_treePtr->gmu_rpx),    "gmu_rpx[ngmu]/F");
   m_tree->Branch("gmu_rpy",    &(m_treePtr->gmu_rpy),    "gmu_rpy[ngmu]/F");
   m_tree->Branch("gmu_rpz",    &(m_treePtr->gmu_rpz),    "gmu_rpz[ngmu]/F");
   m_tree->Branch("gmu_rpt",    &(m_treePtr->gmu_rpt),    "gmu_rpt[ngmu]/F");
   m_tree->Branch("gmu_re",     &(m_treePtr->gmu_re),     "gmu_re[ngmu]/F");
   m_tree->Branch("gmu_ret",    &(m_treePtr->gmu_ret),    "gmu_ret[ngmu]/F");
   m_tree->Branch("gmu_rq",     &(m_treePtr->gmu_rq),     "gmu_rq[ngmu]/F");
   m_tree->Branch("gmu_rphi",   &(m_treePtr->gmu_rphi),   "gmu_rphi[ngmu]/F");
   m_tree->Branch("gmu_reta",   &(m_treePtr->gmu_reta),   "gmu_reta[ngmu]/F");

   // "Z's"
   m_tree->Branch("nrz",       &(m_treePtr->nrz),       "nrz/I");
   m_tree->Branch("z_type",    &(m_treePtr->z_type),    "z_type[nrz]/I");
   m_tree->Branch("z_true",    &(m_treePtr->z_true),    "z_true[nrz]/I");
   m_tree->Branch("z_tid1",    &(m_treePtr->z_tid1),    "z_tid1[nrz]/I");
   m_tree->Branch("z_tid2",    &(m_treePtr->z_tid2),    "z_tid2[nrz]/I");
   m_tree->Branch("z_p",       &(m_treePtr->z_p),       "z_p[nrz]/F");
   m_tree->Branch("z_px",      &(m_treePtr->z_px),      "z_px[nrz]/F");
   m_tree->Branch("z_py",      &(m_treePtr->z_py),      "z_py[nrz]/F");
   m_tree->Branch("z_pz",      &(m_treePtr->z_pz),      "z_pz[nrz]/F");
   m_tree->Branch("z_pt",      &(m_treePtr->z_pt),      "z_pt[nrz]/F");
   m_tree->Branch("z_e",       &(m_treePtr->z_e),       "z_e[nrz]/F");
   m_tree->Branch("z_mass",    &(m_treePtr->z_mass),    "z_mass[nrz]/F");
   m_tree->Branch("z_vvalid",  &(m_treePtr->z_vvalid),  "z_vvalid[nrz]/I");
   m_tree->Branch("z_vchi2",   &(m_treePtr->z_vchi2),   "z_vchi2[nrz]/F");
   m_tree->Branch("z_vndof",   &(m_treePtr->z_vndof),   "z_vndof[nrz]/F");
   m_tree->Branch("z_vd0",     &(m_treePtr->z_vd0),     "z_vd0[nrz]/F");
   m_tree->Branch("z_vx",      &(m_treePtr->z_vx),      "z_vx[nrz]/F");
   m_tree->Branch("z_vy",      &(m_treePtr->z_vy),      "z_vy[nrz]/F");
   m_tree->Branch("z_vz",      &(m_treePtr->z_vz),      "z_vz[nrz]/F");
   m_tree->Branch("z_vdx",     &(m_treePtr->z_vdx),     "z_vdx[nrz]/F");
   m_tree->Branch("z_vdy",     &(m_treePtr->z_vdy),     "z_vdy[nrz]/F");
   m_tree->Branch("z_vdz",     &(m_treePtr->z_vdz),     "z_vdz[nrz]/F");

   m_tree->Branch("z_dp",       &(m_treePtr->z_dp),       "z_dp[nrz][2]/F");
   m_tree->Branch("z_dpx",      &(m_treePtr->z_dpx),      "z_dpx[nrz][2]/F");
   m_tree->Branch("z_dpy",      &(m_treePtr->z_dpy),      "z_dpy[nrz][2]/F");
   m_tree->Branch("z_dpz",      &(m_treePtr->z_dpz),      "z_dpz[nrz][2]/F");
   m_tree->Branch("z_dpt",      &(m_treePtr->z_dpt),      "z_dpt[nrz][2]/F");
   m_tree->Branch("z_de",       &(m_treePtr->z_de),       "z_de[nrz][2]/F");
   m_tree->Branch("z_det",      &(m_treePtr->z_det),      "z_det[nrz][2]/F");
   m_tree->Branch("z_dq",       &(m_treePtr->z_dq),       "z_dq[nrz][2]/F");
   m_tree->Branch("z_dphi",     &(m_treePtr->z_dphi),     "z_dphi[nrz][2]/F");
   m_tree->Branch("z_deta",     &(m_treePtr->z_deta),     "z_deta[nrz][2]/F");

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
