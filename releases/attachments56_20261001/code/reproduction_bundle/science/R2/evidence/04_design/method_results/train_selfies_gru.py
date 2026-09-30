"""Single-CPU SELFIES GRU with optional first-round-neighbor fine-tuning.

Trained here from explicit public-library structures. No external weights,
affinity labels, reinforcement learning, or docking enter training.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import os
import random
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from portable_config import required_path

SCRATCH = required_path("TMC1_WORK_ROOT") / "04_design"

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import selfies as sf
import torch
from torch import nn
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Contrib.SA_Score import sascorer

import train_smiles_gru as common

HERE = Path(__file__).resolve().parent
rdBase.DisableLog("rdApp.error")
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
PARENT_ROWS = common.read_rows(common.PARENT_INPUT)
PARENT_FPS = [(r["parent_id"], FPGEN.GetFingerprint(Chem.MolFromSmiles(r["source_smiles"]))) for r in PARENT_ROWS]


def token_batches(seqs: list[tuple[str, ...]], stoi: dict, size: int, seed: int, shuffle: bool):
    idx = list(range(len(seqs)))
    if shuffle:
        random.Random(seed).shuffle(idx)
    for b in range(0, len(idx), size):
        block = [seqs[i] for i in idx[b:b+size]]
        codes = [[stoi["^"]] + [stoi[t] for t in tokens] + [stoi["$"]] for tokens in block]
        length = max(map(len, codes))
        x = torch.full((len(codes), length-1), stoi["_"], dtype=torch.long)
        y = torch.full_like(x, stoi["_"])
        for i, code in enumerate(codes):
            x[i, :len(code)-1] = torch.tensor(code[:-1])
            y[i, :len(code)-1] = torch.tensor(code[1:])
        yield x, y


def fit(model, seqs, val, stoi, epochs, start_epoch, lr, batch_size, seed, start_time):
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    loss_fn = nn.CrossEntropyLoss(ignore_index=stoi["_"], reduction="sum")
    history = []
    for epoch in range(start_epoch, start_epoch+epochs):
        model.train()
        total_loss, total_tokens = 0.0, 0
        for x, y in token_batches(seqs, stoi, batch_size, seed+epoch, True):
            opt.zero_grad()
            z, _ = model(x)
            loss = loss_fn(z.reshape(-1, len(stoi)), y.reshape(-1))
            tokens = int(y.ne(stoi["_"]).sum().item())
            (loss/tokens).backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            total_loss += loss.item()
            total_tokens += tokens
        model.eval()
        val_loss, val_tokens = 0.0, 0
        with torch.no_grad():
            for x, y in token_batches(val, stoi, batch_size, 0, False):
                z, _ = model(x)
                val_loss += loss_fn(z.reshape(-1, len(stoi)), y.reshape(-1)).item()
                val_tokens += int(y.ne(stoi["_"]).sum().item())
        row = {"epoch": epoch, "train_token_loss": round(total_loss/total_tokens, 5),
               "scaffold_heldout_token_loss": round(val_loss/val_tokens, 5),
               "wall_seconds_cumulative": round(time.monotonic()-start_time, 1)}
        history.append(row)
        print(json.dumps(row), flush=True)
    return history


def properties(mol):
    return {
        "mw": round(Descriptors.MolWt(mol), 3), "clogp": round(Descriptors.MolLogP(mol), 3),
        "tpsa_A2": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "hbd": rdMolDescriptors.CalcNumHBD(mol), "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "sa_score": round(sascorer.calculateScore(mol), 3), "qed": round(QED.qed(mol), 4),
        "pains": common.alerts(mol, common.PAINS), "brenk": common.alerts(mol, common.BRENK),
    }


def assess(model, label, itos, stoi, args, run_dir, input_keys, train_keys,
           legacy_screen_keys, enumerated_keys,
           library_hash, sample_seed):
    raw_samples = common.sample(model, itos, stoi, args.samples, 256, args.max_tokens,
                                args.temperature, sample_seed)
    raw_path = run_dir / f"{label}_raw_generated.csv.gz"
    seen = set()
    counts = Counter()
    sim_yield = Counter()
    eligible = []
    with gzip.open(raw_path, "wt", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_index", "sampled_selfies", "decoded_smiles", "canonical_smiles", "inchikey", "status", "reasons", "nearest_first_round_ecfp4_tanimoto"])
        writer.writeheader()
        for i, sampled in enumerate(raw_samples, 1):
            row = {"sample_index": i, "sampled_selfies": sampled, "decoded_smiles": "", "canonical_smiles": "", "inchikey": "", "status": "", "reasons": "", "nearest_first_round_ecfp4_tanimoto": ""}
            try:
                decoded = sf.decoder(sampled)
                mol = Chem.MolFromSmiles(decoded) if decoded else None
            except Exception:
                decoded, mol = "", None
            row["decoded_smiles"] = decoded or ""
            if mol is None:
                row["status"] = "invalid_decode_or_rdkit"
                counts[row["status"]] += 1
                writer.writerow(row)
                continue
            counts["valid_samples"] += 1
            smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            key = Chem.MolToInchiKey(mol)
            row.update({"canonical_smiles": smi, "inchikey": key})
            if key in seen:
                row["status"] = "duplicate_generated"
            else:
                seen.add(key)
                counts["unique_valid"] += 1
                if key in train_keys:
                    row["status"] = "seen_in_training"
                elif key in input_keys:
                    row["status"] = "seen_in_input_library"
                elif key in legacy_screen_keys:
                    row["status"] = "seen_in_first_round_screen"
                elif key in enumerated_keys:
                    row["status"] = "seen_in_rule_enumeration"
                else:
                    prop = properties(mol)
                    reasons = []
                    if not 150 <= prop["mw"] <= 550: reasons.append("mw")
                    if not -1 <= prop["clogp"] <= 5: reasons.append("clogp")
                    if not 30 <= prop["tpsa_A2"] <= 140: reasons.append("tpsa")
                    if prop["hbd"] > 5: reasons.append("hbd")
                    if prop["hba"] > 10: reasons.append("hba")
                    if prop["rotatable_bonds"] > 10: reasons.append("rotatable_bonds")
                    if prop["sa_score"] > 5.5: reasons.append("sa_score")
                    if Chem.GetFormalCharge(mol) != 0: reasons.append("charge")
                    if prop["pains"]: reasons.append("PAINS")
                    if prop["brenk"]: reasons.append("Brenk")
                    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()): reasons.append("radical")
                    if any(str(s.specified) == "Unspecified" for s in Chem.FindPotentialStereo(mol)):
                        reasons.append("unassigned_stereochemistry")
                    if "." in smi: reasons.append("disconnected")
                    if any(a.GetAtomicNum() not in (1, 6, 7, 8, 9, 15, 16, 17, 35) for a in mol.GetAtoms()): reasons.append("element")
                    fp = FPGEN.GetFingerprint(mol)
                    closest, sim = max(((pid, DataStructs.TanimotoSimilarity(fp, pfp)) for pid, pfp in PARENT_FPS), key=lambda x: x[1])
                    row["nearest_first_round_ecfp4_tanimoto"] = round(sim, 4)
                    if not reasons:
                        counts["property_and_alert_pass"] += 1
                        for t in (0.15, 0.18, 0.2, 0.22, 0.25, 0.3):
                            if sim >= t: sim_yield[str(t)] += 1
                    if not 0.15 <= sim <= 0.8: reasons.append("first_round_similarity")
                    if reasons:
                        row["status"] = "constraint_excluded"
                        row["reasons"] = "|".join(reasons)
                    else:
                        row["status"] = "eligible"
                        eligible.append({
                            "compound_id": "GEN_"+key[:14], "parent_id": "",
                            "nearest_first_round_candidate_id": closest,
                            "source_smiles": decoded, "standardized_isomeric_smiles": smi,
                            "inchikey": key, "method": "locally_trained_SELFIES_GRU_sampling",
                            "model_version": f"selfies_gru_{label}_v1", "training_library_hash": library_hash,
                            "novelty_to_training_set": "exact_InChIKey_not_found",
                            "novelty_to_input_library": "exact_InChIKey_not_found",
                            "scaffold": common.scaffold(mol),
                            "nearest_first_round_ecfp4_tanimoto": round(sim, 4),
                            "mw": prop["mw"], "clogp": prop["clogp"], "tpsa_A2": prop["tpsa_A2"],
                            "hbd": prop["hbd"], "hba": prop["hba"],
                            "rotatable_bonds": prop["rotatable_bonds"], "sa_score": prop["sa_score"],
                            "qed": prop["qed"], "pains_alerts": "", "brenk_alerts": "",
                            "sample_index": i, "sampled_selfies": sampled,
                        })
            counts[row["status"]] += 1
            writer.writerow(row)
    eligible.sort(key=lambda r: (-r["nearest_first_round_ecfp4_tanimoto"], -r["qed"], r["sa_score"], r["compound_id"]))
    fixed, fixed_fps = [], []
    scaffold_counts = Counter()
    for r in eligible:
        if len(fixed) >= 150: break
        if scaffold_counts[r["scaffold"]] >= 2: continue
        fp = FPGEN.GetFingerprint(Chem.MolFromSmiles(r["standardized_isomeric_smiles"]))
        if fixed_fps and max(DataStructs.BulkTanimotoSimilarity(fp, fixed_fps)) > 0.7: continue
        fixed.append(r)
        fixed_fps.append(fp)
        scaffold_counts[r["scaffold"]] += 1
    for i, r in enumerate(fixed, 1): r["selection_rank"] = i
    keys = ["compound_id", "parent_id", "nearest_first_round_candidate_id", "source_smiles",
            "standardized_isomeric_smiles", "inchikey", "method", "model_version",
            "training_library_hash", "novelty_to_training_set", "novelty_to_input_library",
            "scaffold", "nearest_first_round_ecfp4_tanimoto", "mw", "clogp", "tpsa_A2",
            "hbd", "hba", "rotatable_bonds", "sa_score", "qed", "pains_alerts",
            "brenk_alerts", "sample_index", "sampled_selfies", "selection_rank"]
    common.write_csv(HERE/f"{args.run_id}_{label}_all_eligible.csv", eligible, keys)
    common.write_csv(HERE/f"{args.run_id}_{label}_fixed_list.csv", fixed, keys)
    return {"raw_path": str(raw_path), "raw_sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
            "sample_counts": dict(counts), "similarity_yield_after_properties": dict(sim_yield),
            "eligible": len(eligible), "fixed_diverse": len(fixed)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--library", required=True)
    p.add_argument("--smiles-column", required=True)
    p.add_argument("--run-id", required=True)
    p.add_argument("--base-epochs", type=int, default=10)
    p.add_argument("--finetune-epochs", type=int, default=3)
    p.add_argument("--samples", type=int, default=10000)
    p.add_argument("--max-tokens", type=int, default=140)
    p.add_argument("--temperature", type=float, default=0.85)
    p.add_argument("--batch-size", type=int, default=64)
    args = p.parse_args()
    torch.set_num_threads(1)
    random.seed(20260929)
    torch.manual_seed(20260929)
    start = time.monotonic()
    run_dir = SCRATCH/args.run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    lib = Path(args.library)
    lib_hash = hashlib.sha256(lib.read_bytes()).hexdigest()
    raw = common.read_rows(lib)
    rejections = Counter()
    records = []
    seen = set()
    input_keys = set()
    for r in raw:
        original = r.get(args.smiles_column, "")
        mol = Chem.MolFromSmiles(original) if original else None
        if mol is None:
            rejections["input_invalid"] += 1
            continue
        key = Chem.MolToInchiKey(mol)
        input_keys.add(key)
        if key in seen:
            rejections["duplicate_input"] += 1
            continue
        seen.add(key)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        if "." in smiles:
            rejections["disconnected"] += 1
            continue
        try:
            tokens = tuple(sf.split_selfies(sf.encoder(smiles)))
        except Exception:
            rejections["selfies_encoding_failed"] += 1
            continue
        if not 15 <= len(tokens) <= args.max_tokens-2:
            rejections["token_length"] += 1
            continue
        if any(a.GetAtomicNum() not in (1, 6, 7, 8, 9, 15, 16, 17, 35) for a in mol.GetAtoms()):
            rejections["element"] += 1
            continue
        fp = FPGEN.GetFingerprint(mol)
        sim = max(DataStructs.TanimotoSimilarity(fp, pfp) for _, pfp in PARENT_FPS)
        records.append({"inchikey": key, "smiles": smiles, "tokens": tokens,
                        "scaffold": common.scaffold(mol), "nearest_similarity": sim})
    assert len(records) >= 300
    legacy_path = required_path("TMC1_R1_ROOT") / "screen_v2/inputs/library/eligible_pool.csv"
    legacy_screen_keys = set()
    for legacy in common.read_rows(legacy_path):
        legacy_mol = Chem.MolFromSmiles(legacy["source_smiles"])
        if legacy_mol is not None:
            legacy_screen_keys.add(Chem.MolToInchiKey(legacy_mol))
    enumeration_path = HERE/"design_constraints_all.csv"
    if hashlib.sha256(enumeration_path.read_bytes()).hexdigest() != "2a46fc4c47441b6773cc95c72311d4ff1b7e6ea1399b79657d83bdf853629f9b":
        raise ValueError("Frozen rule-enumeration SHA256 mismatch")
    enumerated_keys = {r["inchikey"] for r in common.read_rows(enumeration_path)}
    train_records = [r for r in records if int(hashlib.sha256(r["scaffold"].encode()).hexdigest()[:8], 16)%10 != 0]
    val_records = [r for r in records if int(hashlib.sha256(r["scaffold"].encode()).hexdigest()[:8], 16)%10 == 0]
    assert train_records and val_records
    focused = [r for r in train_records if r["nearest_similarity"] >= 0.25]
    assert len(focused) >= 20
    train = [r["tokens"] for r in train_records]
    val = [r["tokens"] for r in val_records]
    focus = [r["tokens"] for r in focused]
    train_keys = {r["inchikey"] for r in train_records}
    alphabet = sorted({token for r in records for token in r["tokens"]})
    itos = ["_", "^", "$"] + alphabet
    stoi = {token: idx for idx, token in enumerate(itos)}
    (run_dir/"vocab.json").write_text(json.dumps(itos), encoding="utf-8")
    model = common.SmilesGRU(len(itos), 48, 96)
    base_history = fit(model, train, val, stoi, args.base_epochs, 1, 0.002, args.batch_size, 20260929, start)
    base_weights = run_dir/"base_model_state.pt"
    torch.save(model.state_dict(), base_weights)
    base_result = assess(model, "base", itos, stoi, args, run_dir, input_keys, train_keys,
                         legacy_screen_keys, enumerated_keys, lib_hash, 20260930)
    if args.finetune_epochs > 0:
        focused_model = copy.deepcopy(model)
        focus_history = fit(focused_model, focus, val, stoi, args.finetune_epochs, args.base_epochs+1,
                            0.0005, min(args.batch_size, 32), 20260929, start)
        focus_weights = run_dir/"focused_model_state.pt"
        torch.save(focused_model.state_dict(), focus_weights)
        focus_result = assess(focused_model, "focused", itos, stoi, args, run_dir, input_keys, train_keys,
                              legacy_screen_keys, enumerated_keys, lib_hash, 20260931)
    else:
        focus_history, focus_weights, focus_result = [], None, None
    manifest = {
        "run_id": args.run_id, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "method": "locally trained SELFIES single-layer GRU with structural-neighbor fine-tuning",
        "not_pretrained": True, "not_TMC1_affinity_conditioned": True,
        "torch_version": torch.__version__, "selfies_package_metadata_version": "2.2.0",
        "selfies_runtime_version_reported": sf.__version__, "rdkit_version": rdBase.rdkitVersion,
        "library_path": str(lib), "library_sha256": lib_hash, "raw_input_rows": len(raw),
        "first_round_screen_path": str(legacy_path),
        "first_round_screen_sha256": hashlib.sha256(legacy_path.read_bytes()).hexdigest(),
        "rule_enumeration_path": str(enumeration_path),
        "rule_enumeration_sha256": hashlib.sha256(enumeration_path.read_bytes()).hexdigest(),
        "unique_structures_in_model_scope": len(records), "smiles_augmented_examples": 0,
        "base_training_structures": len(train_records), "scaffold_heldout_structures": len(val_records),
        "focused_finetune_structures": len(focused), "focused_similarity_threshold": 0.25,
        "focused_examples_from_training_only": True,
        "input_rejections": dict(rejections), "vocab_size": len(itos),
        "model_parameters": sum(p.numel() for p in model.parameters()),
        "architecture": {"embedding": 48, "hidden": 96, "layers": 1},
        "base_optimizer": "Adam lr=0.002", "focus_optimizer": "Adam lr=0.0005",
        "base_history": base_history, "focused_history": focus_history,
        "seed_training": 20260929, "seed_sample_base": 20260930,
        "seed_sample_focused": 20260931, "temperature": args.temperature,
        "sample_max_tokens": args.max_tokens, "sample_requested_each_model": args.samples,
        "base_weights_path": str(base_weights), "base_weights_sha256": hashlib.sha256(base_weights.read_bytes()).hexdigest(),
        "focused_weights_path": str(focus_weights) if focus_weights else None,
        "focused_weights_sha256": hashlib.sha256(focus_weights.read_bytes()).hexdigest() if focus_weights else None,
        "base_result": base_result, "focused_result": focus_result,
        "wall_seconds": round(time.monotonic()-start, 1),
        "selection_note": "post-sampling physicochemical/PAINS/Brenk filters; nearest first-round similarity 0.15-0.80 is exploratory, not a parent-child synthetic lineage; <=2 per scaffold, pairwise ECFP4 <=0.7; no docking scores used",
    }
    (HERE/f"{args.run_id}_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"run_id": args.run_id, "base_result": base_result, "focused_result": focus_result,
                      "wall_seconds": manifest["wall_seconds"]}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
