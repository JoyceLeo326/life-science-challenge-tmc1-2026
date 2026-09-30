"""Train and sample a small SMILES GRU language model on an explicit library.

This is a locally trained generative neural model, not a pretrained TMC1 model.
There is no affinity conditioning or biological activity claim. Intermediate
weights and all raw samples stay in the explicitly configured work directory.
"""

from __future__ import annotations

import argparse
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
import torch
from torch import nn
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Contrib.SA_Score import sascorer

HERE = Path(__file__).resolve().parent
rdBase.DisableLog("rdApp.error")
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
PARENT_INPUT = required_path("TMC1_R1_ROOT") / "confirmation_v1/inputs/library/selected_candidates.csv"
TRAINED_MODEL_NAME = "single_layer_char_gru_v1"


class SmilesGRU(nn.Module):
    def __init__(self, vocab: int, embed: int = 48, hidden: int = 96):
        super().__init__()
        self.embedding = nn.Embedding(vocab, embed)
        self.gru = nn.GRU(embed, hidden, batch_first=True)
        self.output = nn.Linear(hidden, vocab)

    def forward(self, x, state=None):
        y, state = self.gru(self.embedding(x), state)
        return self.output(y), state


def read_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], keys: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def scaffold(mol) -> str:
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False) or "ACYCLIC"


def catalog(kind):
    p = FilterCatalogParams()
    p.AddCatalog(kind)
    return FilterCatalog(p)


PAINS = catalog(FilterCatalogParams.FilterCatalogs.PAINS)
BRENK = catalog(FilterCatalogParams.FilterCatalogs.BRENK)


def alerts(mol, cat) -> list[str]:
    return sorted({hit.GetDescription() for hit in cat.GetMatches(mol)})


def example_batches(items: list[str], stoi: dict[str, int], batch_size: int, seed: int, shuffle: bool):
    indices = list(range(len(items)))
    if shuffle:
        random.Random(seed).shuffle(indices)
    for start in range(0, len(indices), batch_size):
        samples = [items[i] for i in indices[start:start+batch_size]]
        tokens = [[stoi["^"]] + [stoi[c] for c in s] + [stoi["$"]] for s in samples]
        length = max(len(t) for t in tokens)
        x = torch.full((len(tokens), length-1), stoi["_"], dtype=torch.long)
        y = torch.full_like(x, stoi["_"])
        for i, t in enumerate(tokens):
            x[i, :len(t)-1] = torch.tensor(t[:-1])
            y[i, :len(t)-1] = torch.tensor(t[1:])
        yield x, y


@torch.no_grad()
def sample(model: SmilesGRU, itos: list[str], stoi: dict[str, int], n: int, batch: int, max_len: int,
           temperature: float, seed: int) -> list[str]:
    torch.manual_seed(seed)
    model.eval()
    output = []
    for begin in range(0, n, batch):
        count = min(batch, n-begin)
        prev = torch.full((count, 1), stoi["^"], dtype=torch.long)
        state = None
        strings = [""]*count
        ended = [False]*count
        for _ in range(max_len):
            logits, state = model(prev, state)
            scores = logits[:, -1, :]/temperature
            scores[:, stoi["^"]] = -1e9
            scores[:, stoi["_"]] = -1e9
            next_ids = torch.multinomial(torch.softmax(scores, dim=-1), 1).squeeze(1)
            for i, idx in enumerate(next_ids.tolist()):
                if ended[i]:
                    continue
                char = itos[idx]
                if char == "$":
                    ended[i] = True
                else:
                    strings[i] += char
            prev = next_ids.unsqueeze(1)
            if all(ended):
                break
        output.extend(strings)
    return output


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", required=True)
    ap.add_argument("--smiles-column", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--epochs", type=int, default=3)
    ap.add_argument("--samples", type=int, default=6000)
    ap.add_argument("--temperature", type=float, default=0.85)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--max-len", type=int, default=110)
    args = ap.parse_args()
    assert 0.4 <= args.temperature <= 1.5
    assert 1 <= args.epochs <= 20
    torch.set_num_threads(1)
    random.seed(20260929)
    torch.manual_seed(20260929)
    start = time.monotonic()
    run_dir = SCRATCH / args.run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    lib_path = Path(args.library)
    input_sha256 = hashlib.sha256(lib_path.read_bytes()).hexdigest()
    raw_rows = read_rows(lib_path)
    training_records = []
    all_input_keys = set()
    rejected = Counter()
    seen = set()
    for r in raw_rows:
        raw = r.get(args.smiles_column, "")
        mol = Chem.MolFromSmiles(raw) if raw else None
        if mol is None:
            rejected["invalid_input"] += 1
            continue
        all_input_keys.add(Chem.MolToInchiKey(mol))
        s = Chem.MolToSmiles(mol, isomericSmiles=True)
        if s in seen:
            rejected["duplicate_input"] += 1
            continue
        if not 15 <= len(s) <= args.max_len - 2:
            rejected["length_out_of_scope"] += 1
            continue
        if "." in s or any(a.GetAtomicNum() not in (1, 6, 7, 8, 9, 15, 16, 17, 35) for a in mol.GetAtoms()):
            rejected["salt_or_element_out_of_scope"] += 1
            continue
        seen.add(s)
        training_records.append({"smiles": s, "scaffold": scaffold(mol), "inchikey": Chem.MolToInchiKey(mol)})
    assert len(training_records) >= 300
    # Scaffold-group validation removes exact-scaffold overlap from the model
    # validation set. This is a language-model diagnostic, not affinity testing.
    train = [r["smiles"] for r in training_records
             if int(hashlib.sha256(r["scaffold"].encode()).hexdigest()[:8], 16) % 10 != 0]
    val = [r["smiles"] for r in training_records
           if int(hashlib.sha256(r["scaffold"].encode()).hexdigest()[:8], 16) % 10 == 0]
    if not val:
        val = train[-max(1, len(train)//10):]
        train = train[:-len(val)]
    train_set = set(train)
    train_keys = {r["inchikey"] for r in training_records if r["smiles"] in train_set}
    lib_keys = all_input_keys
    alphabet = sorted(set("".join(r["smiles"] for r in training_records)))
    itos = ["_", "^", "$"] + alphabet
    stoi = {c: i for i, c in enumerate(itos)}
    model = SmilesGRU(len(itos))
    optimizer = torch.optim.Adam(model.parameters(), lr=0.002)
    loss_fn = nn.CrossEntropyLoss(ignore_index=stoi["_"], reduction="sum")
    history = []
    for epoch in range(1, args.epochs+1):
        model.train()
        totals = [0.0, 0]
        for x, y in example_batches(train, stoi, args.batch_size, 20260929+epoch, True):
            optimizer.zero_grad()
            logits, _ = model(x)
            loss = loss_fn(logits.reshape(-1, len(itos)), y.reshape(-1))
            (loss/y.ne(stoi["_"]).sum()).backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            totals[0] += loss.item()
            totals[1] += int(y.ne(stoi["_"]).sum().item())
        model.eval()
        vtotal = [0.0, 0]
        with torch.no_grad():
            for x, y in example_batches(val, stoi, args.batch_size, 0, False):
                logits, _ = model(x)
                vtotal[0] += loss_fn(logits.reshape(-1, len(itos)), y.reshape(-1)).item()
                vtotal[1] += int(y.ne(stoi["_"]).sum().item())
        row = {"epoch": epoch, "train_token_loss": round(totals[0]/totals[1], 5),
               "val_token_loss": round(vtotal[0]/vtotal[1], 5),
               "wall_seconds_cumulative": round(time.monotonic()-start, 1)}
        history.append(row)
        print(json.dumps(row), flush=True)
    weights = run_dir / "model_state.pt"
    torch.save(model.state_dict(), weights)
    (run_dir / "vocab.json").write_text(json.dumps(itos), encoding="utf-8")
    generated = sample(model, itos, stoi, args.samples, 256, args.max_len, args.temperature, 20260930)
    parent_rows = read_rows(PARENT_INPUT)
    parents = [(r["parent_id"], Chem.MolFromSmiles(r["source_smiles"])) for r in parent_rows]
    parent_fps = [(pid, FPGEN.GetFingerprint(m)) for pid, m in parents]
    fixed_ids = {r["inchikey"] for r in read_rows(HERE / "design_constraints_all.csv")}
    raw_output = run_dir / "raw_generated.csv.gz"
    candidates = []
    counts = Counter()
    valid_keys = set()
    eligible_keys = set()
    with gzip.open(raw_output, "wt", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_index", "sampled_smiles", "status", "canonical_smiles", "inchikey", "reason"])
        writer.writeheader()
        for idx, s in enumerate(generated, 1):
            mol = Chem.MolFromSmiles(s) if s else None
            row = {"sample_index": idx, "sampled_smiles": s, "status": "", "canonical_smiles": "", "inchikey": "", "reason": ""}
            if mol is None:
                row["status"] = "invalid_smiles"
                counts[row["status"]] += 1
                writer.writerow(row)
                continue
            canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
            key = Chem.MolToInchiKey(mol)
            row.update({"canonical_smiles": canonical, "inchikey": key})
            counts["valid_samples"] += 1
            if key in valid_keys:
                row["status"] = "duplicate_generated"
            else:
                valid_keys.add(key)
                counts["unique_valid"] += 1
                if key in train_keys:
                    row["status"] = "seen_in_training"
                elif key in lib_keys:
                    row["status"] = "seen_in_library_validation"
                elif key in fixed_ids:
                    row["status"] = "matches_enumerated_analogue"
                else:
                    reasons = []
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    tpsa = rdMolDescriptors.CalcTPSA(mol)
                    hbd = rdMolDescriptors.CalcNumHBD(mol)
                    hba = rdMolDescriptors.CalcNumHBA(mol)
                    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
                    sa = sascorer.calculateScore(mol)
                    pains = alerts(mol, PAINS)
                    brenk = alerts(mol, BRENK)
                    if not 150 <= mw <= 550: reasons.append("mw")
                    if not -1 <= logp <= 5: reasons.append("clogp")
                    if not 30 <= tpsa <= 140: reasons.append("tpsa")
                    if hbd > 5: reasons.append("hbd")
                    if hba > 10: reasons.append("hba")
                    if rb > 10: reasons.append("rotatable_bonds")
                    if sa > 5.5: reasons.append("sa_score")
                    if Chem.GetFormalCharge(mol) != 0: reasons.append("charge")
                    if pains: reasons.append("PAINS")
                    if brenk: reasons.append("Brenk")
                    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()): reasons.append("radical")
                    if any(str(s.specified) == "Unspecified" for s in Chem.FindPotentialStereo(mol)):
                        reasons.append("unassigned_stereochemistry")
                    if any(a.GetAtomicNum() not in (1, 6, 7, 8, 9, 15, 16, 17, 35) for a in mol.GetAtoms()): reasons.append("element")
                    if "." in canonical: reasons.append("disconnected")
                    p_id, sim = max(((pid, DataStructs.TanimotoSimilarity(FPGEN.GetFingerprint(mol), fp)) for pid, fp in parent_fps), key=lambda pair: pair[1])
                    # Yield-calibrated on the 5067-molecule pre-docking pilot:
                    # 0.15 gives exploratory structures near the first-round
                    # panel, without pretending they share a synthetic parent.
                    if not 0.15 <= sim <= 0.8: reasons.append("first_round_similarity")
                    if reasons:
                        row["status"] = "constraint_excluded"
                        row["reason"] = "|".join(reasons)
                    else:
                        row["status"] = "eligible"
                        eligible_keys.add(key)
                        candidates.append({
                            "compound_id": "GEN_"+key[:14], "parent_id": "",
                            "nearest_first_round_candidate_id": p_id,
                            "source_smiles": s, "standardized_isomeric_smiles": canonical,
                            "inchikey": key, "method": "locally_trained_char_GRU_sampling",
                            "model_version": TRAINED_MODEL_NAME, "training_library_hash": input_sha256,
                            "novelty_to_training_set": "exact_InChIKey_not_found",
                            "novelty_to_input_library": "exact_InChIKey_not_found",
                            "scaffold": scaffold(mol), "nearest_first_round_ecfp4_tanimoto": round(sim, 4),
                            "mw": round(mw, 3), "clogp": round(logp, 3), "tpsa_A2": round(tpsa, 2),
                            "hbd": hbd, "hba": hba, "rotatable_bonds": rb,
                            "sa_score": round(sa, 3), "qed": round(QED.qed(mol), 4),
                            "pains_alerts": "", "brenk_alerts": "",
                            "sample_index": idx, "model_weight_path": str(weights),
                        })
            counts[row["status"]] += 1
            writer.writerow(row)
    # Diversity-greedy fixed list, still pre-docking and independent of scores.
    candidates.sort(key=lambda r: (-r["nearest_first_round_ecfp4_tanimoto"], -r["qed"], r["sa_score"], r["compound_id"]))
    selected = []
    selected_fps = []
    scaffold_counts = Counter()
    for r in candidates:
        if len(selected) >= 150:
            break
        if scaffold_counts[r["scaffold"]] >= 2:
            continue
        fp = FPGEN.GetFingerprint(Chem.MolFromSmiles(r["standardized_isomeric_smiles"]))
        if selected_fps and max(DataStructs.BulkTanimotoSimilarity(fp, selected_fps)) > 0.7:
            continue
        selected.append(r)
        selected_fps.append(fp)
        scaffold_counts[r["scaffold"]] += 1
    for i, r in enumerate(selected, 1):
        r["selection_rank"] = i
    keys = ["compound_id", "parent_id", "nearest_first_round_candidate_id",
            "source_smiles", "standardized_isomeric_smiles",
            "inchikey", "method", "model_version", "training_library_hash",
            "novelty_to_training_set", "novelty_to_input_library", "scaffold",
            "nearest_first_round_ecfp4_tanimoto", "mw", "clogp", "tpsa_A2", "hbd", "hba",
            "rotatable_bonds", "sa_score", "qed", "pains_alerts", "brenk_alerts",
            "sample_index", "model_weight_path", "selection_rank"]
    write_csv(HERE / f"{args.run_id}_all_eligible.csv", candidates, keys)
    write_csv(HERE / f"{args.run_id}_fixed_list.csv", selected, keys)
    manifest = {
        "run_id": args.run_id, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "model": TRAINED_MODEL_NAME, "implementation": "local PyTorch GRU; inspired by published SMILES RNN design but not REINVENT implementation or weights",
        "torch_version": torch.__version__, "rdkit_version": rdBase.rdkitVersion,
        "library_path": str(lib_path), "library_sha256": input_sha256,
        "raw_input_rows": len(raw_rows), "dedup_valid_model_scope": len(training_records),
        "training_molecules": len(train), "heldout_scaffold_molecules": len(val),
        "input_rejections": dict(rejected), "vocab_size": len(itos),
        "model_parameters": sum(p.numel() for p in model.parameters()),
        "architecture": {"embedding": 48, "hidden": 96, "layers": 1},
        "optimizer": "Adam lr=0.002", "epoch_history": history,
        "train_seed": 20260929, "sample_seed": 20260930,
        "sample_temperature": args.temperature, "sample_max_chars": args.max_len,
        "sample_requested": args.samples, "sample_counts": dict(counts),
        "unique_eligible": len(eligible_keys), "fixed_diverse_list": len(selected),
        "weights_path": str(weights), "weights_sha256": hashlib.sha256(weights.read_bytes()).hexdigest(),
        "raw_generated_log": str(raw_output), "raw_generated_log_sha256": hashlib.sha256(raw_output.read_bytes()).hexdigest(),
        "wall_seconds": round(time.monotonic()-start, 1),
        "evidence_limit": "No TMC1 affinity/function conditioning, docking, synthesis, PubChem novelty, or patent search performed by this model run",
    }
    (HERE / f"{args.run_id}_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in ("run_id", "training_molecules", "heldout_scaffold_molecules", "model_parameters", "sample_counts", "unique_eligible", "fixed_diverse_list", "wall_seconds")}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
