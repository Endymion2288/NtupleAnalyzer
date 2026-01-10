#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
J/psi + J/psi + Phi window study for DPS: focus on 0.3 < DeltaPhi(Jpsi2, Phi) < 0.7
and 0 < DeltaY(Jpsi2, Phi) < 0.4. Keeps all baseline JJP selections identical to
existing ntuple analysis, then records kinematics of Jpsi1, Jpsi2, Phi and their
constituent muons/kaons inside that window.
"""
import argparse
import glob
import math
import os
import time
import tempfile
import multiprocessing

import ROOT
from ROOT import TFile, TChain, TH1F, TH2F, TLorentzVector

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
JJP_DATA_PATH_DEFAULT = "/eos/user/x/xcheng/learn_MC/JJP_DPS_MC_output/Ntuple/"
TREE_NAME = "mkcands/X_data"

JPSI1_MASS_MIN, JPSI1_MASS_MAX = 2.9, 3.3
JPSI2_MASS_MIN, JPSI2_MASS_MAX = 2.9, 3.3
PHI_MASS_MIN, PHI_MASS_MAX = 0.99, 1.10

JPSI_PT_MIN = 3.0
JPSI_VTXPROB_MIN = 0.05
PHI_PT_MIN = 2.0
PHI_VTXPROB_MIN = 0.05
PHI_K_PT_MIN = 2.0

JPSI_MUON_ID = 'soft'
OUTPUT_DIR = "/eos/user/x/xcheng/CMSSW_14_0_18/src/NtupleAnalyzer/output/"
PLOT_DIR_DEFAULT = None
DEFAULT_WORKERS = max(1, min(8, multiprocessing.cpu_count()))

MUON_MASS = 0.105658
KAON_MASS = 0.493677

# Window definition for Jpsi2 vs Phi
DPHI_MIN, DPHI_MAX = 0.3, 0.7
DY_MIN, DY_MAX = 0.0, 0.4


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def setup_root():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)


def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2 * math.pi
    while dphi < -math.pi:
        dphi += 2 * math.pi
    return dphi


def check_muon_id(chain, mu_idx, id_type):
    """Muon ID check, mirroring analyze_ntuple_JJP.py."""
    if id_type is None:
        return True
    try:
        idx = int(mu_idx)
        if idx < 0:
            return False
        if id_type == 'loose':
            return bool(chain.muIsPatLooseMuon.at(idx))
        if id_type == 'medium':
            return bool(chain.muIsPatMediumMuon.at(idx))
        if id_type == 'tight':
            return bool(chain.muIsPatTightMuon.at(idx))
        if id_type == 'soft':
            return bool(chain.muIsPatSoftMuon.at(idx))
        return True
    except Exception:
        return False


def build_vec_from_pxpypz(px, py, pz, mass):
    e = math.sqrt(px * px + py * py + pz * pz + mass * mass)
    vec = TLorentzVector()
    vec.SetPxPyPzE(px, py, pz, e)
    return vec


def create_histograms():
    h = {}
    # Jpsi/Phi kinematics in window
    h['h_jpsi1_pt'] = TH1F('h_jpsi1_pt', 'J/#psi_{1} p_{T};p_{T} [GeV];Events', 18, 0, 40)
    h['h_jpsi1_eta'] = TH1F('h_jpsi1_eta', 'J/#psi_{1} #eta;#eta;Events', 18, -3, 3)
    h['h_jpsi1_phi'] = TH1F('h_jpsi1_phi', 'J/#psi_{1} #phi;#phi;Events', 18, -math.pi, math.pi)

    h['h_jpsi2_pt'] = TH1F('h_jpsi2_pt', 'J/#psi_{2} p_{T};p_{T} [GeV];Events', 18, 0, 40)
    h['h_jpsi2_eta'] = TH1F('h_jpsi2_eta', 'J/#psi_{2} #eta;#eta;Events', 18, -3, 3)
    h['h_jpsi2_phi'] = TH1F('h_jpsi2_phi', 'J/#psi_{2} #phi;#phi;Events', 18, -math.pi, math.pi)

    h['h_phi_pt'] = TH1F('h_phi_pt', '#phi p_{T};p_{T} [GeV];Events', 18, 0, 40)
    h['h_phi_eta'] = TH1F('h_phi_eta', '#phi #eta;#eta;Events', 18, -3, 3)
    h['h_phi_phi'] = TH1F('h_phi_phi', '#phi #phi;#phi;Events', 18, -math.pi, math.pi)

    # Muon kinematics (Jpsi1 mu1/mu2, Jpsi2 mu1/mu2)
    h['h_mu_jpsi1_mu1_pt'] = TH1F('h_mu_jpsi1_mu1_pt', 'Muon (J/#psi_{1} #mu_{1}) p_{T};p_{T} [GeV];Events', 18, 0, 40)
    h['h_mu_jpsi1_mu1_eta'] = TH1F('h_mu_jpsi1_mu1_eta', 'Muon (J/#psi_{1} #mu_{1}) #eta;#eta;Events', 18, -3, 3)
    h['h_mu_jpsi1_mu1_phi'] = TH1F('h_mu_jpsi1_mu1_phi', 'Muon (J/#psi_{1} #mu_{1}) #phi;#phi;Events', 18, -math.pi, math.pi)

    h['h_mu_jpsi1_mu2_pt'] = TH1F('h_mu_jpsi1_mu2_pt', 'Muon (J/#psi_{1} #mu_{2}) p_{T};p_{T} [GeV];Events', 18, 0, 40)
    h['h_mu_jpsi1_mu2_eta'] = TH1F('h_mu_jpsi1_mu2_eta', 'Muon (J/#psi_{1} #mu_{2}) #eta;#eta;Events', 18, -3, 3)
    h['h_mu_jpsi1_mu2_phi'] = TH1F('h_mu_jpsi1_mu2_phi', 'Muon (J/#psi_{1} #mu_{2}) #phi;#phi;Events', 18, -math.pi, math.pi)

    h['h_mu_jpsi2_mu1_pt'] = TH1F('h_mu_jpsi2_mu1_pt', 'Muon (J/#psi_{2} #mu_{1}) p_{T};p_{T} [GeV];Events', 18, 0, 40)
    h['h_mu_jpsi2_mu1_eta'] = TH1F('h_mu_jpsi2_mu1_eta', 'Muon (J/#psi_{2} #mu_{1}) #eta;#eta;Events', 18, -3, 3)
    h['h_mu_jpsi2_mu1_phi'] = TH1F('h_mu_jpsi2_mu1_phi', 'Muon (J/#psi_{2} #mu_{1}) #phi;#phi;Events', 18, -math.pi, math.pi)

    h['h_mu_jpsi2_mu2_pt'] = TH1F('h_mu_jpsi2_mu2_pt', 'Muon (J/#psi_{2} #mu_{2}) p_{T};p_{T} [GeV];Events', 18, 0, 40)
    h['h_mu_jpsi2_mu2_eta'] = TH1F('h_mu_jpsi2_mu2_eta', 'Muon (J/#psi_{2} #mu_{2}) #eta;#eta;Events', 18, -3, 3)
    h['h_mu_jpsi2_mu2_phi'] = TH1F('h_mu_jpsi2_mu2_phi', 'Muon (J/#psi_{2} #mu_{2}) #phi;#phi;Events', 18, -math.pi, math.pi)

    # Kaon kinematics
    h['h_k1_pt'] = TH1F('h_k1_pt', 'K_{1} p_{T};p_{T} [GeV];Events', 18, 0, 40)
    h['h_k1_eta'] = TH1F('h_k1_eta', 'K_{1} #eta;#eta;Events', 18, -3, 3)
    h['h_k1_phi'] = TH1F('h_k1_phi', 'K_{1} #phi;#phi;Events', 18, -math.pi, math.pi)

    h['h_k2_pt'] = TH1F('h_k2_pt', 'K_{2} p_{T};p_{T} [GeV];Events', 18, 0, 40)
    h['h_k2_eta'] = TH1F('h_k2_eta', 'K_{2} #eta;#eta;Events', 18, -3, 3)
    h['h_k2_phi'] = TH1F('h_k2_phi', 'K_{2} #phi;#phi;Events', 18, -math.pi, math.pi)

    # Delta y vs Delta phi (Jpsi2 vs Phi) inside window, zoom to window range
    # x: Delta y (6 bins), y: Delta phi (8 bins), both restricted to window range
    h['h2_dy_dphi_jpsi2_phi'] = TH2F('h2_dy_dphi_jpsi2_phi',
        'J/#psi_{2} - #phi: #Delta y vs #Delta#phi;#Delta y;#Delta#phi',
        6, DY_MIN, DY_MAX, 8, DPHI_MIN, DPHI_MAX)
    return h


def save_plots(histos, plot_dir):
    os.makedirs(plot_dir, exist_ok=True)
    canvas = ROOT.TCanvas("c", "c", 800, 600)
    for name, hist in histos.items():
        canvas.Clear()
        if isinstance(hist, TH2F):
            hist.Draw("COLZ")
        else:
            hist.Draw("HIST")
        canvas.SaveAs(os.path.join(plot_dir, f"{name}.png"))


def merge_histograms(dest, src):
    for name, hdest in dest.items():
        hsrc = src.Get(name)
        if hsrc:
            hdest.Add(hsrc)


def process_file_batch(file_list, max_events, muon_id, tree_name):
    """Worker: process a batch of files, return temp ROOT path and stats."""
    chain = TChain(tree_name)
    for f in file_list:
        chain.Add(f)

    histos = create_histograms()
    n_total = chain.GetEntries()
    n_to_process = n_total if max_events < 0 else min(max_events, n_total)
    n_window = 0
    n_has_cand = 0
    n_pass_baseline = 0
    n_fail_mass = 0
    n_fail_pt = 0
    n_fail_vtx = 0
    n_fail_kpt = 0
    n_fail_mu = 0
    n_mu_fill_ok = 0
    n_mu_fill_fail = 0
    n_mu_size_zero = 0
    n_mu_idx_neg = 0
    n_mu_idx_ge_size = 0
    mu_size_min = None
    mu_size_max = None
    mu_idx_max_seen = None
    debug_fail_prints = 0

    for i_evt in range(n_to_process):
        chain.GetEntry(i_evt)
        try:
            n_cand = chain.Jpsi_1_mass.size()
        except Exception:
            continue
        if n_cand == 0:
            continue
        n_has_cand += 1

        best_cand = None
        best_score = -1

        for i_cand in range(n_cand):
            try:
                jpsi1_mass = chain.Jpsi_1_mass.at(i_cand)
                jpsi2_mass = chain.Jpsi_2_mass.at(i_cand)
                phi_mass = chain.Phi_mass.at(i_cand)
                if not (JPSI1_MASS_MIN < jpsi1_mass < JPSI1_MASS_MAX):
                    n_fail_mass += 1
                    continue
                if not (JPSI2_MASS_MIN < jpsi2_mass < JPSI2_MASS_MAX):
                    n_fail_mass += 1
                    continue
                if not (PHI_MASS_MIN < phi_mass < PHI_MASS_MAX):
                    n_fail_mass += 1
                    continue

                jpsi1_pt = chain.Jpsi_1_pt.at(i_cand)
                jpsi2_pt = chain.Jpsi_2_pt.at(i_cand)
                phi_pt = chain.Phi_pt.at(i_cand)
                if jpsi1_pt < JPSI_PT_MIN or jpsi2_pt < JPSI_PT_MIN or phi_pt < PHI_PT_MIN:
                    n_fail_pt += 1
                    continue

                jpsi1_vtxprob = chain.Jpsi_1_VtxProb.at(i_cand)
                jpsi2_vtxprob = chain.Jpsi_2_VtxProb.at(i_cand)
                phi_vtxprob = chain.Phi_VtxProb.at(i_cand)
                if jpsi1_vtxprob < JPSI_VTXPROB_MIN or jpsi2_vtxprob < JPSI_VTXPROB_MIN or phi_vtxprob < PHI_VTXPROB_MIN:
                    n_fail_vtx += 1
                    continue

                phi_k1_pt = chain.Phi_K_1_pt.at(i_cand)
                phi_k2_pt = chain.Phi_K_2_pt.at(i_cand)
                if phi_k1_pt < PHI_K_PT_MIN or phi_k2_pt < PHI_K_PT_MIN:
                    n_fail_kpt += 1
                    continue

                jpsi1_mu1_idx = int(chain.Jpsi_1_mu_1_Idx.at(i_cand))
                jpsi1_mu2_idx = int(chain.Jpsi_1_mu_2_Idx.at(i_cand))
                jpsi2_mu1_idx = int(chain.Jpsi_2_mu_1_Idx.at(i_cand))
                jpsi2_mu2_idx = int(chain.Jpsi_2_mu_2_Idx.at(i_cand))

                if muon_id and muon_id.lower() != 'none':
                    if not (check_muon_id(chain, jpsi1_mu1_idx, muon_id) and
                            check_muon_id(chain, jpsi1_mu2_idx, muon_id) and
                            check_muon_id(chain, jpsi2_mu1_idx, muon_id) and
                            check_muon_id(chain, jpsi2_mu2_idx, muon_id)):
                        n_fail_mu += 1
                        continue

                jpsi1_eta = chain.Jpsi_1_eta.at(i_cand)
                jpsi1_phi = chain.Jpsi_1_phi.at(i_cand)
                jpsi2_eta = chain.Jpsi_2_eta.at(i_cand)
                jpsi2_phi = chain.Jpsi_2_phi.at(i_cand)
                phi_eta = chain.Phi_eta.at(i_cand)
                phi_phi = chain.Phi_phi.at(i_cand)

                score = math.sqrt(jpsi1_pt * jpsi1_pt + jpsi2_pt * jpsi2_pt + phi_pt * phi_pt)
                if score > best_score:
                    best_score = score
                    best_cand = {
                        'jpsi1_pt': jpsi1_pt, 'jpsi1_eta': jpsi1_eta, 'jpsi1_phi': jpsi1_phi, 'jpsi1_mass': jpsi1_mass,
                        'jpsi2_pt': jpsi2_pt, 'jpsi2_eta': jpsi2_eta, 'jpsi2_phi': jpsi2_phi, 'jpsi2_mass': jpsi2_mass,
                        'phi_pt': phi_pt, 'phi_eta': phi_eta, 'phi_phi': phi_phi, 'phi_mass': phi_mass,
                        'jpsi1_mu1_idx': jpsi1_mu1_idx, 'jpsi1_mu2_idx': jpsi1_mu2_idx,
                        'jpsi2_mu1_idx': jpsi2_mu1_idx, 'jpsi2_mu2_idx': jpsi2_mu2_idx,
                        'phi_k1_pt': phi_k1_pt, 'phi_k1_eta': chain.Phi_K_1_eta.at(i_cand), 'phi_k1_phi': chain.Phi_K_1_phi.at(i_cand),
                        'phi_k2_pt': phi_k2_pt, 'phi_k2_eta': chain.Phi_K_2_eta.at(i_cand), 'phi_k2_phi': chain.Phi_K_2_phi.at(i_cand),
                        'phi_k1_px': chain.Phi_K_1_px.at(i_cand), 'phi_k1_py': chain.Phi_K_1_py.at(i_cand), 'phi_k1_pz': chain.Phi_K_1_pz.at(i_cand),
                        'phi_k2_px': chain.Phi_K_2_px.at(i_cand), 'phi_k2_py': chain.Phi_K_2_py.at(i_cand), 'phi_k2_pz': chain.Phi_K_2_pz.at(i_cand),
                    }
            except Exception:
                continue

        if best_cand is None:
            continue

        n_pass_baseline += 1

        jpsi1_4vec = TLorentzVector()
        jpsi2_4vec = TLorentzVector()
        phi_4vec = TLorentzVector()
        jpsi1_4vec.SetPtEtaPhiM(best_cand['jpsi1_pt'], best_cand['jpsi1_eta'], best_cand['jpsi1_phi'], best_cand['jpsi1_mass'])
        jpsi2_4vec.SetPtEtaPhiM(best_cand['jpsi2_pt'], best_cand['jpsi2_eta'], best_cand['jpsi2_phi'], best_cand['jpsi2_mass'])
        phi_4vec.SetPtEtaPhiM(best_cand['phi_pt'], best_cand['phi_eta'], best_cand['phi_phi'], best_cand['phi_mass'])

        dy_jpsi2_phi = abs(jpsi2_4vec.Rapidity() - phi_4vec.Rapidity())
        dphi_jpsi2_phi = abs(delta_phi(jpsi2_4vec.Phi(), phi_4vec.Phi()))

        if not (DPHI_MIN < dphi_jpsi2_phi < DPHI_MAX and DY_MIN < dy_jpsi2_phi < DY_MAX):
            continue

        n_window += 1

        histos['h2_dy_dphi_jpsi2_phi'].Fill(dy_jpsi2_phi, dphi_jpsi2_phi)
        histos['h_jpsi1_pt'].Fill(jpsi1_4vec.Pt())
        histos['h_jpsi1_eta'].Fill(jpsi1_4vec.Eta())
        histos['h_jpsi1_phi'].Fill(jpsi1_4vec.Phi())
        histos['h_jpsi2_pt'].Fill(jpsi2_4vec.Pt())
        histos['h_jpsi2_eta'].Fill(jpsi2_4vec.Eta())
        histos['h_jpsi2_phi'].Fill(jpsi2_4vec.Phi())
        histos['h_phi_pt'].Fill(phi_4vec.Pt())
        histos['h_phi_eta'].Fill(phi_4vec.Eta())
        histos['h_phi_phi'].Fill(phi_4vec.Phi())

        try:
            mu_size = chain.muPx.size()
            if mu_size_min is None or mu_size < mu_size_min:
                mu_size_min = mu_size
            if mu_size_max is None or mu_size > mu_size_max:
                mu_size_max = mu_size
            idxs = [int(best_cand['jpsi1_mu1_idx']), int(best_cand['jpsi1_mu2_idx']),
                    int(best_cand['jpsi2_mu1_idx']), int(best_cand['jpsi2_mu2_idx'])]
            max_idx = max(idxs)
            if mu_idx_max_seen is None or max_idx > mu_idx_max_seen:
                mu_idx_max_seen = max_idx
            invalid_neg = any(idx < 0 for idx in idxs)
            invalid_ge = any(idx >= mu_size for idx in idxs)
            if invalid_neg or invalid_ge:
                n_mu_fill_fail += 1
                if mu_size == 0:
                    n_mu_size_zero += 1
                if invalid_neg:
                    n_mu_idx_neg += 1
                if invalid_ge:
                    n_mu_idx_ge_size += 1
                if debug_fail_prints < 5:
                    print(f"[DEBUG] Mu fill fail: idxs={idxs}, mu_size={mu_size}")
                    debug_fail_prints += 1
            else:
                mu1_vec = build_vec_from_pxpypz(chain.muPx.at(best_cand['jpsi1_mu1_idx']),
                                                 chain.muPy.at(best_cand['jpsi1_mu1_idx']),
                                                 chain.muPz.at(best_cand['jpsi1_mu1_idx']), MUON_MASS)
                mu2_vec = build_vec_from_pxpypz(chain.muPx.at(best_cand['jpsi1_mu2_idx']),
                                                 chain.muPy.at(best_cand['jpsi1_mu2_idx']),
                                                 chain.muPz.at(best_cand['jpsi1_mu2_idx']), MUON_MASS)
                mu3_vec = build_vec_from_pxpypz(chain.muPx.at(best_cand['jpsi2_mu1_idx']),
                                                 chain.muPy.at(best_cand['jpsi2_mu1_idx']),
                                                 chain.muPz.at(best_cand['jpsi2_mu1_idx']), MUON_MASS)
                mu4_vec = build_vec_from_pxpypz(chain.muPx.at(best_cand['jpsi2_mu2_idx']),
                                                 chain.muPy.at(best_cand['jpsi2_mu2_idx']),
                                                 chain.muPz.at(best_cand['jpsi2_mu2_idx']), MUON_MASS)

                histos['h_mu_jpsi1_mu1_pt'].Fill(mu1_vec.Pt())
                histos['h_mu_jpsi1_mu1_eta'].Fill(mu1_vec.Eta())
                histos['h_mu_jpsi1_mu1_phi'].Fill(mu1_vec.Phi())
                histos['h_mu_jpsi1_mu2_pt'].Fill(mu2_vec.Pt())
                histos['h_mu_jpsi1_mu2_eta'].Fill(mu2_vec.Eta())
                histos['h_mu_jpsi1_mu2_phi'].Fill(mu2_vec.Phi())
                histos['h_mu_jpsi2_mu1_pt'].Fill(mu3_vec.Pt())
                histos['h_mu_jpsi2_mu1_eta'].Fill(mu3_vec.Eta())
                histos['h_mu_jpsi2_mu1_phi'].Fill(mu3_vec.Phi())
                histos['h_mu_jpsi2_mu2_pt'].Fill(mu4_vec.Pt())
                histos['h_mu_jpsi2_mu2_eta'].Fill(mu4_vec.Eta())
                histos['h_mu_jpsi2_mu2_phi'].Fill(mu4_vec.Phi())
                n_mu_fill_ok += 1
        except Exception:
            n_mu_fill_fail += 1

        k1_vec = build_vec_from_pxpypz(best_cand['phi_k1_px'], best_cand['phi_k1_py'], best_cand['phi_k1_pz'], KAON_MASS)
        k2_vec = build_vec_from_pxpypz(best_cand['phi_k2_px'], best_cand['phi_k2_py'], best_cand['phi_k2_pz'], KAON_MASS)
        histos['h_k1_pt'].Fill(k1_vec.Pt())
        histos['h_k1_eta'].Fill(k1_vec.Eta())
        histos['h_k1_phi'].Fill(k1_vec.Phi())
        histos['h_k2_pt'].Fill(k2_vec.Pt())
        histos['h_k2_eta'].Fill(k2_vec.Eta())
        histos['h_k2_phi'].Fill(k2_vec.Phi())

    fd, tmp_path = tempfile.mkstemp(suffix=".root", prefix="jjp_window_tmp_")
    os.close(fd)
    fout = TFile(tmp_path, "RECREATE")
    for h in histos.values():
        h.Write()
    fout.Close()

    return (tmp_path, n_to_process, n_window, n_has_cand, n_pass_baseline,
            n_fail_mass, n_fail_pt, n_fail_vtx, n_fail_kpt, n_fail_mu,
            n_mu_fill_ok, n_mu_fill_fail, n_mu_size_zero,
            mu_size_min if mu_size_min is not None else 0,
            mu_size_max if mu_size_max is not None else 0,
            mu_idx_max_seen if mu_idx_max_seen is not None else -1,
            n_mu_idx_neg, n_mu_idx_ge_size)


# -----------------------------------------------------------------------------
# Core analysis
# -----------------------------------------------------------------------------
def analyze_jjp_window(max_events=-1, muon_id='soft', input_dir=None, output_file=None, plot_dir=None, n_workers=1):
    global JPSI_MUON_ID
    JPSI_MUON_ID = muon_id

    setup_root()

    data_path = input_dir if input_dir else JJP_DATA_PATH_DEFAULT
    data_files = [os.path.join(data_path, os.path.basename(f)) for f in glob.glob(os.path.join(data_path, '*.root')) if os.path.exists(os.path.join(data_path, os.path.basename(f)))]
    n_files = len(data_files)
    print(f"[INFO] 输入目录: {data_path}")
    print(f"[INFO] 加载 {n_files} 个文件")

    if n_files == 0:
        print("[ERROR] 未找到输入文件")
        return None

    if n_workers <= 1:
        batch_files = [data_files]
    else:
        batch_files = [[] for _ in range(n_workers)]
        for idx, f in enumerate(data_files):
            batch_files[idx % n_workers].append(f)
        batch_files = [b for b in batch_files if b]
        n_workers = len(batch_files)

    if max_events < 0:
        per_batch_events = -1
    else:
        per_batch_events = math.ceil(max_events / len(batch_files))

    # Run workers
    if len(batch_files) == 1:
        result = process_file_batch(batch_files[0], per_batch_events, muon_id, TREE_NAME)
        temp_files = [result[0]]
        total_processed = result[1]
        total_in_window = result[2]
        total_has_cand = result[3]
        total_pass_baseline = result[4]
        total_fail_mass = result[5]
        total_fail_pt = result[6]
        total_fail_vtx = result[7]
        total_fail_kpt = result[8]
        total_fail_mu = result[9]
        total_mu_fill_ok = result[10]
        total_mu_fill_fail = result[11]
        total_mu_size_zero = result[12]
        total_mu_size_min = result[13]
        total_mu_size_max = result[14]
        total_mu_idx_max = result[15]
        total_mu_idx_neg = result[16]
        total_mu_idx_ge_size = result[17]
    else:
        with multiprocessing.Pool(processes=len(batch_files)) as pool:
            results = pool.starmap(process_file_batch,
                                   [(b, per_batch_events, muon_id, TREE_NAME) for b in batch_files])
        temp_files = [r[0] for r in results]
        total_processed = sum(r[1] for r in results)
        total_in_window = sum(r[2] for r in results)
        total_has_cand = sum(r[3] for r in results)
        total_pass_baseline = sum(r[4] for r in results)
        total_fail_mass = sum(r[5] for r in results)
        total_fail_pt = sum(r[6] for r in results)
        total_fail_vtx = sum(r[7] for r in results)
        total_fail_kpt = sum(r[8] for r in results)
        total_fail_mu = sum(r[9] for r in results)
        total_mu_fill_ok = sum(r[10] for r in results)
        total_mu_fill_fail = sum(r[11] for r in results)
        total_mu_size_zero = sum(r[12] for r in results)
        total_mu_size_min = min(r[13] for r in results if r[13] > 0) if results else 0
        total_mu_size_max = max(r[14] for r in results) if results else 0
        total_mu_idx_max = max(r[15] for r in results) if results else -1
        total_mu_idx_neg = sum(r[16] for r in results)
        total_mu_idx_ge_size = sum(r[17] for r in results)

    histos = create_histograms()
    for tf in temp_files:
        fin = TFile.Open(tf)
        if fin and not fin.IsZombie():
            merge_histograms(histos, fin)
        fin.Close()
        os.remove(tf)

    elapsed = time.time() - start_time
    n_to_process = total_processed if max_events < 0 else min(max_events, total_processed)
    print(f"[INFO] 处理事件总数: {total_processed}")
    print(f"[INFO] 有候选的事件数: {total_has_cand}")
    print(f"[INFO] 通过基线cuts的事件数: {total_pass_baseline}")
    print(f"[INFO] 窗口内事件数: {total_in_window}")
    print(f"[INFO] 失败统计: mass={total_fail_mass}, pt={total_fail_pt}, vtx={total_fail_vtx}, kpt={total_fail_kpt}, muID={total_fail_mu}")
    print(f"[INFO] Mu填充: ok={total_mu_fill_ok}, fail={total_mu_fill_fail}, neg_idx={total_mu_idx_neg}, ge_size={total_mu_idx_ge_size}")
    print(f"[INFO] Mu尺寸: size(min,max)=({total_mu_size_min},{total_mu_size_max}), idx_max={total_mu_idx_max}, size0={total_mu_size_zero}")
    rate = 0 if elapsed <= 0 else total_processed / elapsed
    print(f"[INFO] 耗时: {elapsed:.1f}s ({rate:.0f} evt/s)")

    if output_file is None:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        output_file = os.path.join(OUTPUT_DIR, "jjp_window_region.root")

    fout = TFile(output_file, "RECREATE")
    for h in histos.values():
        h.Write()
    fout.Close()
    print(f"[INFO] 输出保存到: {output_file}")

    # Plotting
    final_plot_dir = plot_dir
    if final_plot_dir is None:
        base_dir = os.path.dirname(output_file)
        final_plot_dir = os.path.join(base_dir if base_dir else '.', 'plots_JJP_window')
    print(f"[INFO] 保存图像到: {final_plot_dir}")
    save_plots(histos, final_plot_dir)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='JJP DPS window study (0.3<DeltaPhi<0.7, 0<DeltaY<0.4 for Jpsi2-Phi)')
    parser.add_argument('-n', '--max-events', type=int, default=-1, help='最大处理事件数 (-1=全部)')
    parser.add_argument('--muon-id', type=str, default='soft', choices=['loose', 'medium', 'tight', 'soft', 'none'], help='Muon ID要求')
    parser.add_argument('-i', '--input-dir', type=str, default=None, help='输入Ntuple目录')
    parser.add_argument('-o', '--output', type=str, default=None, help='输出ROOT文件路径')
    parser.add_argument('-p', '--plot-dir', type=str, default=None, help='输出图像目录')
    parser.add_argument('-j', '--jobs', type=int, default=1, help='并行进程数 (默认1)')
    args = parser.parse_args()

    global start_time
    start_time = time.time()

    mu_id = None if args.muon_id == 'none' else args.muon_id

    analyze_jjp_window(max_events=args.max_events,
                       muon_id=mu_id,
                       input_dir=args.input_dir,
                       output_file=args.output,
                       plot_dir=args.plot_dir,
                       n_workers=max(1, args.jobs))


if __name__ == '__main__':
    main()
