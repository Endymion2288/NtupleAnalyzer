#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
J/psi + J/psi + Phi (JJP) Ntuple角度关联分析

功能:
1. 从Ntuple加载数据并应用事件选择cuts
2. 计算粒子之间的角度关联 (Δy, Δφ)
3. 填充直方图并保存

使用方法:
    python analyze_ntuple_JJP.py -o output_jjp.root
    python analyze_ntuple_JJP.py -n 10000  # 处理前10000个事件
    python analyze_ntuple_JJP.py --muon-id soft  # 指定muon ID
"""

import ROOT
from ROOT import TFile, TChain, TH1F, TH2F, TLorentzVector
import os
import math
import argparse
import glob
import time

# =============================================================================
# 配置参数
# =============================================================================

# 数据路径
JJP_DATA_PATH = "/eos/user/x/xcheng/JpsiJpsiPhi_muon_pt_cut/merged_rootNtuple/"
JJP_DATA_FILES = [os.path.basename(f) for f in glob.glob(JJP_DATA_PATH + "*.root")]
TREE_NAME = "mkcands/X_data"

# 质量窗口
JPSI1_MASS_MIN, JPSI1_MASS_MAX = 2.9, 3.3
JPSI2_MASS_MIN, JPSI2_MASS_MAX = 2.9, 3.3
PHI_MASS_MIN, PHI_MASS_MAX = 0.99, 1.10

# Cuts 参数
JPSI_PT_MIN = 3.0
JPSI_VTXPROB_MIN = 0.05
PHI_PT_MIN = 2.0
PHI_VTXPROB_MIN = 0.05
PHI_K_PT_MIN = 2.0

# Muon ID 选择
JPSI_MUON_ID = 'soft'

# 输出目录
OUTPUT_DIR = "/eos/user/x/xcheng/CMSSW_14_0_18/src/NtupleAnalyzer/output/"


# =============================================================================
# 辅助函数
# =============================================================================

def setup_root():
    """配置ROOT环境"""
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)


def delta_phi(phi1, phi2):
    """计算Δφ，包装到[-π, π]"""
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2 * math.pi
    while dphi < -math.pi:
        dphi += 2 * math.pi
    return dphi


def rapidity_from_4vec(pt, eta, phi, mass):
    """从pT, η, φ, m计算快度y"""
    vec = TLorentzVector()
    vec.SetPtEtaPhiM(pt, eta, phi, mass)
    return vec.Rapidity()


def check_muon_id(chain, mu_idx, id_type):
    """检查muon是否通过ID选择"""
    if id_type is None:
        return True
    
    try:
        idx = int(mu_idx)
        if idx < 0:
            return False
        
        if id_type == 'loose':
            return chain.muIsPatLooseMuon.at(idx)
        elif id_type == 'medium':
            return chain.muIsPatMediumMuon.at(idx)
        elif id_type == 'tight':
            return chain.muIsPatTightMuon.at(idx)
        elif id_type == 'soft':
            return chain.muIsPatSoftMuon.at(idx)
        else:
            return True
    except:
        return False


def create_histograms():
    """创建所有直方图"""
    histograms = {}
    
    # 1D: Δy (快度差)
    histograms['h_dy_jpsi1_jpsi2'] = TH1F("h_dy_jpsi1_jpsi2",
        "#Delta y (J/#psi_{1} - J/#psi_{2});|#Delta y|;Events", 50, 0, 5)
    histograms['h_dy_jpsi1_phi'] = TH1F("h_dy_jpsi1_phi",
        "#Delta y (J/#psi_{1} - #phi);|#Delta y|;Events", 50, 0, 5)
    histograms['h_dy_jpsi2_phi'] = TH1F("h_dy_jpsi2_phi",
        "#Delta y (J/#psi_{2} - #phi);|#Delta y|;Events", 50, 0, 5)
    
    # 1D: Δφ (方位角差)
    histograms['h_dphi_jpsi1_jpsi2'] = TH1F("h_dphi_jpsi1_jpsi2",
        "#Delta#phi (J/#psi_{1} - J/#psi_{2});|#Delta#phi|;Events", 50, 0, math.pi)
    histograms['h_dphi_jpsi1_phi'] = TH1F("h_dphi_jpsi1_phi",
        "#Delta#phi (J/#psi_{1} - #phi);|#Delta#phi|;Events", 50, 0, math.pi)
    histograms['h_dphi_jpsi2_phi'] = TH1F("h_dphi_jpsi2_phi",
        "#Delta#phi (J/#psi_{2} - #phi);|#Delta#phi|;Events", 50, 0, math.pi)
    
    # 2D: Δy vs Δφ
    histograms['h2_dy_dphi_jpsi1_jpsi2'] = TH2F("h2_dy_dphi_jpsi1_jpsi2",
        "J/#psi_{1} - J/#psi_{2};|#Delta y|;|#Delta#phi|", 50, 0, 5, 50, 0, math.pi)
    histograms['h2_dy_dphi_jpsi1_phi'] = TH2F("h2_dy_dphi_jpsi1_phi",
        "J/#psi_{1} - #phi;|#Delta y|;|#Delta#phi|", 50, 0, 5, 50, 0, math.pi)
    histograms['h2_dy_dphi_jpsi2_phi'] = TH2F("h2_dy_dphi_jpsi2_phi",
        "J/#psi_{2} - #phi;|#Delta y|;|#Delta#phi|", 50, 0, 5, 50, 0, math.pi)
    
    # 运动学分布: pT
    histograms['h_jpsi1_pt'] = TH1F("h_jpsi1_pt",
        "J/#psi_{1} p_{T};p_{T} [GeV];Events", 100, 0, 50)
    histograms['h_jpsi2_pt'] = TH1F("h_jpsi2_pt",
        "J/#psi_{2} p_{T};p_{T} [GeV];Events", 100, 0, 50)
    histograms['h_phi_pt'] = TH1F("h_phi_pt",
        "#phi p_{T};p_{T} [GeV];Events", 100, 0, 50)
    
    # 运动学分布: η
    histograms['h_jpsi1_eta'] = TH1F("h_jpsi1_eta",
        "J/#psi_{1} #eta;#eta;Events", 60, -3, 3)
    histograms['h_jpsi2_eta'] = TH1F("h_jpsi2_eta",
        "J/#psi_{2} #eta;#eta;Events", 60, -3, 3)
    histograms['h_phi_eta'] = TH1F("h_phi_eta",
        "#phi #eta;#eta;Events", 60, -3, 3)
    
    # 运动学分布: y (快度)
    histograms['h_jpsi1_y'] = TH1F("h_jpsi1_y",
        "J/#psi_{1} y;y;Events", 60, -3, 3)
    histograms['h_jpsi2_y'] = TH1F("h_jpsi2_y",
        "J/#psi_{2} y;y;Events", 60, -3, 3)
    histograms['h_phi_y'] = TH1F("h_phi_y",
        "#phi y;y;Events", 60, -3, 3)
    
    # 运动学分布: φ (方位角)
    histograms['h_jpsi1_phi'] = TH1F("h_jpsi1_phi",
        "J/#psi_{1} #phi;#phi;Events", 60, -math.pi, math.pi)
    histograms['h_jpsi2_phi'] = TH1F("h_jpsi2_phi",
        "J/#psi_{2} #phi;#phi;Events", 60, -math.pi, math.pi)
    histograms['h_phi_phi'] = TH1F("h_phi_phi",
        "#phi #phi;#phi;Events", 60, -math.pi, math.pi)
    
    # 不变质量分布
    histograms['h_mass_jpsi1_jpsi2'] = TH1F("h_mass_jpsi1_jpsi2",
        "M(J/#psi_{1} + J/#psi_{2});M [GeV];Events", 100, 6, 30)
    histograms['h_mass_jpsi1_phi'] = TH1F("h_mass_jpsi1_phi",
        "M(J/#psi_{1} + #phi);M [GeV];Events", 100, 3, 20)
    histograms['h_mass_jpsi2_phi'] = TH1F("h_mass_jpsi2_phi",
        "M(J/#psi_{2} + #phi);M [GeV];Events", 100, 3, 20)
    histograms['h_mass_all'] = TH1F("h_mass_all",
        "M(J/#psi_{1} + J/#psi_{2} + #phi);M [GeV];Events", 100, 7, 40)
    
    return histograms


def fill_histograms(histograms, jpsi1_4vec, jpsi2_4vec, phi_4vec):
    """填充所有直方图"""
    # 计算快度
    y_jpsi1 = jpsi1_4vec.Rapidity()
    y_jpsi2 = jpsi2_4vec.Rapidity()
    y_phi = phi_4vec.Rapidity()
    
    # 计算Δy
    dy_jpsi1_jpsi2 = abs(y_jpsi1 - y_jpsi2)
    dy_jpsi1_phi = abs(y_jpsi1 - y_phi)
    dy_jpsi2_phi = abs(y_jpsi2 - y_phi)
    
    # 计算Δφ
    dphi_jpsi1_jpsi2 = abs(delta_phi(jpsi1_4vec.Phi(), jpsi2_4vec.Phi()))
    dphi_jpsi1_phi = abs(delta_phi(jpsi1_4vec.Phi(), phi_4vec.Phi()))
    dphi_jpsi2_phi = abs(delta_phi(jpsi2_4vec.Phi(), phi_4vec.Phi()))
    
    # 填充1D直方图
    histograms['h_dy_jpsi1_jpsi2'].Fill(dy_jpsi1_jpsi2)
    histograms['h_dy_jpsi1_phi'].Fill(dy_jpsi1_phi)
    histograms['h_dy_jpsi2_phi'].Fill(dy_jpsi2_phi)
    
    histograms['h_dphi_jpsi1_jpsi2'].Fill(dphi_jpsi1_jpsi2)
    histograms['h_dphi_jpsi1_phi'].Fill(dphi_jpsi1_phi)
    histograms['h_dphi_jpsi2_phi'].Fill(dphi_jpsi2_phi)
    
    # 填充2D直方图
    histograms['h2_dy_dphi_jpsi1_jpsi2'].Fill(dy_jpsi1_jpsi2, dphi_jpsi1_jpsi2)
    histograms['h2_dy_dphi_jpsi1_phi'].Fill(dy_jpsi1_phi, dphi_jpsi1_phi)
    histograms['h2_dy_dphi_jpsi2_phi'].Fill(dy_jpsi2_phi, dphi_jpsi2_phi)
    
    # 填充运动学直方图
    histograms['h_jpsi1_pt'].Fill(jpsi1_4vec.Pt())
    histograms['h_jpsi2_pt'].Fill(jpsi2_4vec.Pt())
    histograms['h_phi_pt'].Fill(phi_4vec.Pt())
    
    histograms['h_jpsi1_eta'].Fill(jpsi1_4vec.Eta())
    histograms['h_jpsi2_eta'].Fill(jpsi2_4vec.Eta())
    histograms['h_phi_eta'].Fill(phi_4vec.Eta())
    
    histograms['h_jpsi1_y'].Fill(y_jpsi1)
    histograms['h_jpsi2_y'].Fill(y_jpsi2)
    histograms['h_phi_y'].Fill(y_phi)
    
    histograms['h_jpsi1_phi'].Fill(jpsi1_4vec.Phi())
    histograms['h_jpsi2_phi'].Fill(jpsi2_4vec.Phi())
    histograms['h_phi_phi'].Fill(phi_4vec.Phi())
    
    # 填充不变质量直方图
    histograms['h_mass_jpsi1_jpsi2'].Fill((jpsi1_4vec + jpsi2_4vec).M())
    histograms['h_mass_jpsi1_phi'].Fill((jpsi1_4vec + phi_4vec).M())
    histograms['h_mass_jpsi2_phi'].Fill((jpsi2_4vec + phi_4vec).M())
    histograms['h_mass_all'].Fill((jpsi1_4vec + jpsi2_4vec + phi_4vec).M())


def analyze_jjp_ntuple(max_events=-1, muon_id='soft', output_file=None):
    """
    分析JJP Ntuple
    
    Args:
        max_events: 最大处理事件数 (-1表示全部)
        muon_id: muon ID类型 ('loose', 'medium', 'tight', 'soft', None)
        output_file: 输出ROOT文件路径
    """
    global JPSI_MUON_ID
    JPSI_MUON_ID = muon_id
    
    print("\n" + "="*60)
    print("J/psi + J/psi + Phi Ntuple 角度关联分析")
    print("="*60)
    
    start_time = time.time()
    
    # 创建TChain
    chain = TChain(TREE_NAME)
    n_files = 0
    for f in JJP_DATA_FILES:
        filepath = JJP_DATA_PATH + f
        if os.path.exists(filepath):
            chain.Add(filepath)
            n_files += 1
    
    print(f"[INFO] 加载 {n_files} 个文件")
    
    n_total = chain.GetEntries()
    n_to_process = n_total if max_events < 0 else min(max_events, n_total)
    print(f"[INFO] 总事件数: {n_total}, 处理: {n_to_process}")
    print(f"[INFO] Muon ID 要求: {muon_id}")
    
    # 创建直方图
    histograms = create_histograms()
    
    n_passed = 0
    
    for i_evt in range(n_to_process):
        if i_evt % 10000 == 0:
            elapsed = time.time() - start_time
            rate = i_evt / elapsed if elapsed > 0 else 0
            print(f"[INFO] 处理事件 {i_evt}/{n_to_process} ({rate:.0f} evt/s)")
        
        chain.GetEntry(i_evt)
        
        # 获取候选数
        try:
            n_cand = chain.Jpsi_1_mass.size()
        except:
            continue
        
        if n_cand == 0:
            continue
        
        # 选择最佳候选
        best_cand = None
        best_score = -1
        
        for i_cand in range(n_cand):
            try:
                # 获取质量
                jpsi1_mass = chain.Jpsi_1_mass.at(i_cand)
                jpsi2_mass = chain.Jpsi_2_mass.at(i_cand)
                phi_mass = chain.Phi_mass.at(i_cand)
                
                # 质量窗口cut
                if not (JPSI1_MASS_MIN < jpsi1_mass < JPSI1_MASS_MAX):
                    continue
                if not (JPSI2_MASS_MIN < jpsi2_mass < JPSI2_MASS_MAX):
                    continue
                if not (PHI_MASS_MIN < phi_mass < PHI_MASS_MAX):
                    continue
                
                # 运动学cuts
                jpsi1_pt = chain.Jpsi_1_pt.at(i_cand)
                jpsi2_pt = chain.Jpsi_2_pt.at(i_cand)
                phi_pt = chain.Phi_pt.at(i_cand)
                
                if jpsi1_pt < JPSI_PT_MIN:
                    continue
                if jpsi2_pt < JPSI_PT_MIN:
                    continue
                if phi_pt < PHI_PT_MIN:
                    continue
                
                # 顶点概率cuts
                jpsi1_vtxprob = chain.Jpsi_1_VtxProb.at(i_cand)
                jpsi2_vtxprob = chain.Jpsi_2_VtxProb.at(i_cand)
                phi_vtxprob = chain.Phi_VtxProb.at(i_cand)
                
                if jpsi1_vtxprob < JPSI_VTXPROB_MIN:
                    continue
                if jpsi2_vtxprob < JPSI_VTXPROB_MIN:
                    continue
                if phi_vtxprob < PHI_VTXPROB_MIN:
                    continue
                
                # K介子pT cut
                phi_k1_pt = chain.Phi_K_1_pt.at(i_cand)
                phi_k2_pt = chain.Phi_K_2_pt.at(i_cand)
                if phi_k1_pt < PHI_K_PT_MIN or phi_k2_pt < PHI_K_PT_MIN:
                    continue
                
                # Muon ID cuts
                if JPSI_MUON_ID:
                    jpsi1_mu1_idx = chain.Jpsi_1_mu_1_Idx.at(i_cand)
                    jpsi1_mu2_idx = chain.Jpsi_1_mu_2_Idx.at(i_cand)
                    if not check_muon_id(chain, jpsi1_mu1_idx, JPSI_MUON_ID):
                        continue
                    if not check_muon_id(chain, jpsi1_mu2_idx, JPSI_MUON_ID):
                        continue
                    
                    jpsi2_mu1_idx = chain.Jpsi_2_mu_1_Idx.at(i_cand)
                    jpsi2_mu2_idx = chain.Jpsi_2_mu_2_Idx.at(i_cand)
                    if not check_muon_id(chain, jpsi2_mu1_idx, JPSI_MUON_ID):
                        continue
                    if not check_muon_id(chain, jpsi2_mu2_idx, JPSI_MUON_ID):
                        continue
                
                # 获取运动学信息
                jpsi1_eta = chain.Jpsi_1_eta.at(i_cand)
                jpsi1_phi = chain.Jpsi_1_phi.at(i_cand)
                jpsi2_eta = chain.Jpsi_2_eta.at(i_cand)
                jpsi2_phi = chain.Jpsi_2_phi.at(i_cand)
                phi_eta = chain.Phi_eta.at(i_cand)
                phi_phi = chain.Phi_phi.at(i_cand)
                
                # 计算pT score
                score = math.sqrt(jpsi1_pt**2 + jpsi2_pt**2 + phi_pt**2)
                
                if score > best_score:
                    best_score = score
                    best_cand = {
                        'jpsi1_pt': jpsi1_pt, 'jpsi1_eta': jpsi1_eta,
                        'jpsi1_phi': jpsi1_phi, 'jpsi1_mass': jpsi1_mass,
                        'jpsi2_pt': jpsi2_pt, 'jpsi2_eta': jpsi2_eta,
                        'jpsi2_phi': jpsi2_phi, 'jpsi2_mass': jpsi2_mass,
                        'phi_pt': phi_pt, 'phi_eta': phi_eta,
                        'phi_phi': phi_phi, 'phi_mass': phi_mass,
                    }
                    
            except Exception as e:
                continue
        
        # 填充最佳候选
        if best_cand is not None:
            jpsi1_4vec = TLorentzVector()
            jpsi2_4vec = TLorentzVector()
            phi_4vec = TLorentzVector()
            
            jpsi1_4vec.SetPtEtaPhiM(best_cand['jpsi1_pt'], best_cand['jpsi1_eta'],
                                     best_cand['jpsi1_phi'], best_cand['jpsi1_mass'])
            jpsi2_4vec.SetPtEtaPhiM(best_cand['jpsi2_pt'], best_cand['jpsi2_eta'],
                                     best_cand['jpsi2_phi'], best_cand['jpsi2_mass'])
            phi_4vec.SetPtEtaPhiM(best_cand['phi_pt'], best_cand['phi_eta'],
                                   best_cand['phi_phi'], best_cand['phi_mass'])
            
            fill_histograms(histograms, jpsi1_4vec, jpsi2_4vec, phi_4vec)
            n_passed += 1
    
    elapsed = time.time() - start_time
    print(f"\n[INFO] 处理完成!")
    print(f"[INFO] 通过选择的事件数: {n_passed}/{n_to_process}")
    print(f"[INFO] 选择效率: {100*n_passed/n_to_process:.2f}%")
    print(f"[INFO] 耗时: {elapsed:.1f}s ({n_to_process/elapsed:.0f} evt/s)")
    
    # 保存直方图
    if output_file is None:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        output_file = os.path.join(OUTPUT_DIR, "jjp_ntuple_correlations.root")
    
    fout = TFile(output_file, "RECREATE")
    for h in histograms.values():
        h.Write()
    fout.Close()
    
    print(f"[INFO] 输出保存到: {output_file}")
    
    return histograms


def main():
    parser = argparse.ArgumentParser(description='JJP Ntuple角度关联分析')
    parser.add_argument('-n', '--max-events', type=int, default=-1,
                        help='最大处理事件数 (-1=全部)')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='输出ROOT文件路径')
    parser.add_argument('--muon-id', type=str, default='soft',
                        choices=['loose', 'medium', 'tight', 'soft', 'none'],
                        help='Muon ID要求 (默认: soft)')
    
    args = parser.parse_args()
    
    setup_root()
    
    muon_id = None if args.muon_id == 'none' else args.muon_id
    
    analyze_jjp_ntuple(
        max_events=args.max_events,
        muon_id=muon_id,
        output_file=args.output
    )


if __name__ == '__main__':
    main()
