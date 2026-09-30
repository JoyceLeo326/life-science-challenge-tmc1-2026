"""Reproduce frozen structural-alert diagnostics, retaining every domain row.

Run from any directory with RDKit 2025.09.6 and matplotlib installed.
All release paths are relative to this script. No filtering or reranking occurs.
"""
import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"
import csv
import ctypes
import hashlib
import json
import platform
import random
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "diagnostics"
OUT.mkdir(exist_ok=True)
CST = timezone(timedelta(hours=8))
CANDIDATES = {"LIB_YWDIXVZCRMMQQY", "LIB_MYAFDHCCQBIZLC", "LIB_ZNDXPRCORYDEAP", "LIB_FCRHJSNRFAAYQD", "LIB_DNDVGFIEWOHGAH", "LIB_ZVHBDYFMHCQADD"}
EXPECTED_INPUT = "79b775967bfc87ded9a07aae03549c9233691ec4fc3375a2c2662b64ee4b2841"
FAMILIES = ("pains", "pains_A", "pains_B", "pains_C", "brenk", "reactive")
REPRESENTATIONS = {"primary": "rule_smiles", "sensitivity": "standardized_isomeric_smiles"}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def write_csv(path, rows, fields=None):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

def iso():
    return datetime.now(CST).isoformat()

def set_one_cpu():
    if os.name == "nt":
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.GetCurrentProcess.restype = ctypes.c_void_p
        kernel.GetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
        kernel.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        handle = kernel.GetCurrentProcess()
        before, system = ctypes.c_size_t(), ctypes.c_size_t()
        assert kernel.GetProcessAffinityMask(handle, ctypes.byref(before), ctypes.byref(system))
        selected = before.value & -before.value
        assert selected and kernel.SetProcessAffinityMask(handle, selected)
        after = ctypes.c_size_t()
        assert kernel.GetProcessAffinityMask(handle, ctypes.byref(after), ctypes.byref(system))
        assert after.value == selected and selected.bit_count() == 1
        return {"method": "Windows SetProcessAffinityMask", "before_mask_hex": hex(before.value), "actual_mask_hex": hex(after.value), "logical_cpus_allowed": 1, "worker_processes": 1}
    available = os.sched_getaffinity(0)
    os.sched_setaffinity(0, {min(available)})
    assert len(os.sched_getaffinity(0)) == 1
    return {"method": "sched_setaffinity", "logical_cpus_allowed": 1, "worker_processes": 1}

def official_entries(text, family):
    block = re.search(r"const FilterData_t " + family + r"\[\]\s*=\s*\{(.*?)\};", text, re.S).group(1)
    string = r'"(?:\\.|[^"\\])*"'
    expression = r'(?:' + string + r'\s*)+'
    pattern = r'\{\s*(' + string + r')\s*,\s*(' + expression + r')\s*,\s*(\d+)\s*,'
    entries = []
    for match in re.finditer(pattern, block, re.S):
        name = json.loads(match.group(1))
        smarts = "".join(json.loads(s) for s in re.findall(string, match.group(2)))
        entries.append((name, smarts, int(match.group(3))))
    return entries

def main():
    started = iso()
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    affinity = set_one_cpu()
    from rdkit import Chem, rdBase, RDLogger
    from rdkit.Chem import FilterCatalog
    RDLogger.DisableLog("rdApp.warning")
    rules_path = ROOT / "rules/RULES_FROZEN.json"
    catalog_path = ROOT / "rules/rule_catalog.csv"
    rules = json.loads(rules_path.read_text(encoding="utf-8-sig"))
    domain = json.loads((ROOT / "frozen_domain/DOMAIN_FREEZE.json").read_text(encoding="utf-8-sig"))
    assert datetime.fromisoformat(started) > datetime.fromisoformat(rules["frozen_at_cst"])
    assert rdBase.rdkitVersion == rules["rdkit_version"] == "2025.09.6"
    assert sha(catalog_path) == rules["catalog_csv_sha256"]
    assert sha(ROOT / "frozen_domain/all_predictions.csv") == EXPECTED_INPUT
    assert sha(ROOT / "frozen_domain/eligible_9098.csv") == domain["eligible_csv_sha256"]
    for source in rules["source_files"]:
        assert sha(ROOT / "rules" / source["file"]) == source["sha256"], source["file"]
    frozen_rules = read_csv(catalog_path)
    assert len(frozen_rules) == 585
    by_description = {r["description"]: r for r in frozen_rules}
    assert len(by_description) == 585
    reactive = {r["description"] for r in rules["reactive_subset"]}
    assert len(reactive) == 10
    for r in rules["reactive_subset"]:
        assert by_description[r["description"]]["family"] == "BRENK"
        assert by_description[r["description"]]["smarts"] == r["smarts"]
    cpp = (ROOT / "rules/official_sources/Filters.cpp").read_text(encoding="utf-8")
    source_counts = {}
    for family in ("PAINS_A", "PAINS_B", "PAINS_C", "BRENK"):
        source_text = cpp if family == "BRENK" else (ROOT / "rules/official_sources" / (family.lower() + ".in")).read_text(encoding="utf-8")
        entries = official_entries(source_text, family)
        expected = [(r["description"], r["smarts"], 0) for r in frozen_rules if r["family"] == family]
        assert entries == expected, family
        source_counts[family] = len(entries)
    cats = {}
    built_entries = {}
    query_by_name = {}
    checks = []
    for family, enum in (("PAINS", FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS), ("BRENK", FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)):
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(enum)
        cat = FilterCatalog.FilterCatalog(params)
        assert cat.GetNumEntries() == (480 if family == "PAINS" else 105)
        cats[family.lower()] = cat
        for i in range(cat.GetNumEntries()):
            actual = cat.GetEntryWithIdx(i)
            name = actual.GetDescription()
            r = by_description[name]
            assert r["family"].startswith(family)
            assert actual.GetProp("FilterSet") == ("Brenk" if r["family"] == "BRENK" else r["family"])
            query = Chem.MolFromSmarts(r["smarts"], mergeHs=True)
            assert query is not None
            matcher = FilterCatalog.SmartsMatcher(name, query, int(r["min_count"]), int(r["max_count"]))
            reconstructed = FilterCatalog.FilterCatalogEntry(name, matcher)
            for prop in actual.GetPropList():
                reconstructed.SetProp(prop, actual.GetProp(prop))
            # Byte-identical full serialized entry verifies built-in query, description,
            # match bounds, properties and source-defined hydrogen-merging semantics.
            assert actual.Serialize() == reconstructed.Serialize(), name
            query_by_name[name] = query
            built_entries[name] = actual
            checks.append({"family": r["family"], "rule_id": name, "description_exact": True, "smarts_entry_serialization_exact": True, "min_count": matcher.GetMinCount(), "max_count": matcher.GetMaxCount(), "official_source_exact": True})
    write_csv(OUT / "catalog_runtime_checks.csv", checks)
    rows = read_csv(ROOT / "frozen_domain/eligible_9098.csv")
    all_rows = read_csv(ROOT / "frozen_domain/all_predictions.csv")
    expected_ids = [r["compound_id"] for r in all_rows if r["docking_eligibility"] == "eligible"]
    assert len(all_rows) == 16698 and len(rows) == len(expected_ids) == 9098
    assert [r["compound_id"] for r in rows] == expected_ids
    assert len({r["compound_id"] for r in rows}) == len({r["parent_key"] for r in rows}) == 9098
    assert CANDIDATES <= set(expected_ids)
    sources = defaultdict(set)
    for r in read_csv(ROOT / "library/all_source_records.csv"):
        sources[r["compound_id"]].add(r["source_id"])
    assert all(r["source_id"] in sources[r["compound_id"]] for r in rows)
    cache = {}
    def evaluate(smiles):
        if smiles in cache:
            return cache[smiles]
        hits = {f: [] for f in FAMILIES}
        try:
            if not smiles:
                raise ValueError("missing_smiles")
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                raise ValueError("SMILES_parse_or_sanitize_failed")
            for family in ("pains", "brenk"):
                hits[family] = sorted(e.GetDescription() for e in cats[family].GetMatches(molecule))
            for sf in ("A", "B", "C"):
                hits["pains_" + sf] = [n for n in hits["pains"] if by_description[n]["family"] == "PAINS_" + sf]
            hits["reactive"] = sorted(set(hits["brenk"]) & reactive)
            result = {"status": "OK", "error": "", "hits": hits}
        except Exception as exc:
            result = {"status": "FAILED", "error": type(exc).__name__ + ": " + str(exc), "hits": hits}
        cache[smiles] = result
        return result
    diagnostics = []
    for index, row in enumerate(rows, 1):
        out = {k: row[k] for k in ("compound_id", "source_id", "parent_key", "inchi_key", "rule_inchikey", "rule_smiles", "standardized_isomeric_smiles")}
        out["source_ids_json"] = json.dumps(sorted(sources[row["compound_id"]]), separators=(",", ":"))
        out["is_six_candidate"] = int(row["compound_id"] in CANDIDATES)
        results = {}
        for rep, field in REPRESENTATIONS.items():
            result = evaluate(row[field])
            results[rep] = result
            out[rep + "_status"] = result["status"]
            out[rep + "_error"] = result["error"]
            for family in FAMILIES:
                names = result["hits"][family]
                out[f"{rep}_{family}_rule_count"] = len(names)
                out[f"{rep}_{family}_rules_json"] = json.dumps(names, separators=(",", ":"))
        ok_count = sum(r["status"] == "OK" for r in results.values())
        out["overall_status"] = "OK" if ok_count == 2 else ("PARTIAL_FAILURE" if ok_count == 1 else "FAILED")
        out["union_complete"] = int(ok_count == 2)
        for family in FAMILIES:
            names = sorted(set(results["primary"]["hits"][family]) | set(results["sensitivity"]["hits"][family]))
            out[f"union_{family}_rule_count"] = len(names)
            out[f"union_{family}_rules_json"] = json.dumps(names, separators=(",", ":"))
        out["representation_alert_set_differs"] = int(any(results["primary"]["hits"][f] != results["sensitivity"]["hits"][f] for f in FAMILIES))
        diagnostics.append(out)
        if index % 500 == 0:
            print(f"Evaluated {index}/9098; unique representations={len(cache)}; wall={time.perf_counter()-wall_start:.1f}s", flush=True)
    diag_path = OUT / "chemical_alerts_9098.csv"
    write_csv(diag_path, diagnostics)
    # Independent re-read and recount use CSV output only, retaining failed attempts.
    reread = read_csv(diag_path)
    assert [r["compound_id"] for r in reread] == expected_ids
    assert len(reread) == len({r["compound_id"] for r in reread}) == len({r["parent_key"] for r in reread}) == 9098
    stats = []
    frequencies = Counter()
    overlap = Counter()
    for rep in ("primary", "sensitivity", "union"):
        evaluable = sum(r[rep + "_status"] == "OK" for r in reread) if rep != "union" else sum(int(r["union_complete"]) for r in reread)
        for family in FAMILIES:
            count = 0
            for row in reread:
                names = json.loads(row[f"{rep}_{family}_rules_json"])
                assert len(names) == len(set(names)) == int(row[f"{rep}_{family}_rule_count"])
                count += bool(names)
                frequencies.update((rep, family, name) for name in names)
                if rep == "union":
                    assert set(names) == set(json.loads(row[f"primary_{family}_rules_json"])) | set(json.loads(row[f"sensitivity_{family}_rules_json"]))
            stats.append({"representation": rep, "family": family, "hit_molecules": count, "denominator_all_attempted": 9098, "percent_all_attempted": round(count / 9098 * 100, 6), "evaluable_molecules": evaluable, "failed_or_incomplete_molecules": 9098-evaluable})
    for row in reread:
        p, b, r = [bool(int(row[f"union_{f}_rule_count"])) for f in ("pains", "brenk", "reactive")]
        assert not r or b
        key = (int(p), int(b), int(r), int(row["union_complete"]))
        overlap[key] += 1
    overlap_rows = [{"union_pains_hit": p, "union_brenk_hit": b, "union_reactive_hit": r, "union_complete": complete, "molecules": count, "denominator_all_attempted": 9098, "percent_all_attempted": round(count/9098*100, 6)} for (p,b,r,complete), count in sorted(overlap.items())]
    assert sum(x["molecules"] for x in overlap_rows) == 9098
    frequency_rows = []
    for rep in ("primary", "sensitivity", "union"):
        for family in FAMILIES:
            names = [r["description"] for r in frozen_rules if (family == "pains" and r["family"].startswith("PAINS")) or (family.startswith("pains_") and r["family"] == "PAINS_" + family[-1]) or (family == "brenk" and r["family"] == "BRENK") or (family == "reactive" and r["description"] in reactive)]
            for name in names:
                frequency_rows.append({"representation": rep, "family": family, "rule_id": name, "hit_molecules": frequencies[(rep,family,name)], "denominator_all_attempted": 9098, "source_url": by_description[name]["source_url"]})
    write_csv(OUT / "alert_summary.csv", stats)
    write_csv(OUT / "rule_frequencies.csv", frequency_rows)
    write_csv(OUT / "overlap_patterns.csv", overlap_rows)
    candidate_rows = [r for r in reread if int(r["is_six_candidate"])]
    assert len(candidate_rows) == 6 and {r["compound_id"] for r in candidate_rows} == CANDIDATES
    write_csv(OUT / "six_candidate_alerts.csv", candidate_rows)
    # Deterministic random molecules and all six candidates: direct official SMARTS
    # matching with mergeHs=True, using every frozen rule, independently of GetMatches.
    selected_ids = set(random.Random(20260930).sample(expected_ids, 12)) | CANDIDATES
    direct_checks = []
    for row in reread:
        if row["compound_id"] not in selected_ids:
            continue
        for rep, field in REPRESENTATIONS.items():
            if row[rep + "_status"] != "OK":
                direct_checks.append({"compound_id": row["compound_id"], "representation": rep, "status": "FAILED_RETAINED", "rules_compared": 0, "mismatches": 0})
                continue
            molecule = Chem.MolFromSmiles(row[field])
            direct = {n for n,q in query_by_name.items() if molecule.HasSubstructMatch(q)}
            catalog_hits = set(json.loads(row[f"{rep}_pains_rules_json"])) | set(json.loads(row[f"{rep}_brenk_rules_json"]))
            mismatch = direct ^ catalog_hits
            assert not mismatch, (row["compound_id"], rep, sorted(mismatch))
            direct_checks.append({"compound_id": row["compound_id"], "representation": rep, "status": "PASS", "rules_compared": 585, "mismatches": 0})
    write_csv(OUT / "direct_smarts_spotchecks.csv", direct_checks)
    failed = Counter(r["overall_status"] for r in reread)
    summary = {"status": "DIAGNOSTIC_COMPUTE_COMPLETE", "domain_rows": 9098, "unique_compound_ids": 9098, "unique_parent_keys": 9098, "representation_status_counts": {rep: dict(Counter(r[rep+"_status"] for r in reread)) for rep in REPRESENTATIONS}, "overall_status_counts": dict(failed), "representation_alert_set_differs_molecules": sum(int(r["representation_alert_set_differs"]) for r in reread), "summary": stats, "overlap_patterns": overlap_rows, "six_candidates": [{"compound_id": r["compound_id"], "source_id": r["source_id"], "status": r["overall_status"], "union_pains_rules": json.loads(r["union_pains_rules_json"]), "union_brenk_rules": json.loads(r["union_brenk_rules_json"]), "union_reactive_rules": json.loads(r["union_reactive_rules_json"])} for r in candidate_rows], "diagnostic_only": True, "changes_to_library_or_original_rank": False, "primary_representation": "rule_smiles", "sensitivity_representation": "standardized_isomeric_smiles", "failed_attempts_kept_in_all_denominators": True, "union_is_not_sum_of_family_counts": True}
    dump(OUT / "summary.json", summary)
    print("COMPUTE_COMPLETE " + json.dumps({"rows":9098, "overall_status_counts":dict(failed), "counts":[s for s in stats if s["family"] in ("pains","brenk","reactive")], "six_candidates":summary["six_candidates"]}, ensure_ascii=False), flush=True)
    # One publication-style two-panel figure. Bars are separate overlapping alerts;
    # right panel gives mutually exclusive union patterns, avoiding double counting.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font_candidates = [f.name for f in font_manager.fontManager.ttflist if f.name in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC")]
    plt.rcParams.update({"font.family": font_candidates[0] if font_candidates else "DejaVu Sans", "font.size": 10, "axes.unicode_minus": False, "pdf.fonttype": 42})
    figure_source = []
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.5), gridspec_kw={"width_ratios":[1,1.35]})
    colors = ["#4477AA", "#EE6677", "#228833"]
    labels = ["PAINS", "BRENK", "反应性提示子集"]
    for j, rep in enumerate(("primary", "sensitivity", "union")):
        for i,family in enumerate(("pains","brenk","reactive")):
            s = next(x for x in stats if x["representation"] == rep and x["family"] == family)
            x = i+(j-1)*0.25
            axes[0].bar(x, s["percent_all_attempted"], width=0.23, color=colors[j], label=["主表示","源表示敏感性","两表示并集"][j] if i==0 else None)
            axes[0].text(x, s["percent_all_attempted"]+0.7, f'{s["hit_molecules"]}/9098', rotation=90, ha="center", va="bottom", fontsize=8)
            figure_source.append({"panel":"A", "representation":rep,"category":family,"molecules":s["hit_molecules"],"denominator":9098,"percent":s["percent_all_attempted"]})
    axes[0].set_xticks(range(3), labels)
    axes[0].set_ylabel("命中占全部9098个尝试分子 (%)")
    axes[0].set_title("A 结构提示：各组有重叠")
    axes[0].set_ylim(0, min(105,max(x["percent_all_attempted"] for x in stats if x["family"] in ("pains","brenk","reactive"))+17))
    axes[0].legend(loc="upper right", frameon=False, fontsize=9)
    pattern_labels, pattern_counts = [], []
    for r in overlap_rows:
        if not r["union_complete"]:
            label = "表示不完整（保留分母）"
        elif not r["union_pains_hit"] and not r["union_brenk_hit"]:
            label = "无上述结构提示"
        else:
            parts = []
            if r["union_pains_hit"]: parts.append("PAINS")
            if r["union_brenk_hit"]: parts.append("BRENK")
            if r["union_reactive_hit"]: parts.append("含反应性子集")
            label = " + ".join(parts)
        pattern_labels.append(label)
        pattern_counts.append(r["molecules"])
        figure_source.append({"panel":"B","representation":"union","category":label,"molecules":r["molecules"],"denominator":9098,"percent":r["percent_all_attempted"]})
    axes[1].barh(range(len(pattern_counts)), [n/9098*100 for n in pattern_counts], color="#4477AA")
    axes[1].set_yticks(range(len(pattern_counts)), pattern_labels)
    axes[1].invert_yaxis()
    for i,n in enumerate(pattern_counts):
        axes[1].text(n/9098*100+0.7,i,f"{n}/9098 ({n/9098*100:.1f}%)",va="center",fontsize=9)
    axes[1].set_xlim(0, min(105,max(pattern_counts)/9098*100+22))
    axes[1].set_xlabel("占全部9098个尝试分子 (%)")
    axes[1].set_title("B 并集的互斥重叠模式（合计9098）")
    for ax in axes:
        ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("冻结9098分子库的结构提示诊断 · RDKit 2025.09.6",fontsize=14)
    fig.text(0.02,0.035,"反应性提示是冻结BRENK的10条规则子集；命中不等于实验反应性、毒性、干扰或TMC1活性。",fontsize=9)
    fig.tight_layout(rect=(0,0.08,1,0.94),w_pad=2.5)
    fig.savefig(OUT / "chemical_alerts_9098.png",dpi=300,bbox_inches="tight")
    fig.savefig(OUT / "chemical_alerts_9098.pdf",bbox_inches="tight")
    plt.close(fig)
    write_csv(OUT / "figure_source.csv",figure_source)
    # Recompute again from the saved CSV and compare every summary record.
    for s in read_csv(OUT / "alert_summary.csv"):
        count = sum(int(r[f'{s["representation"]}_{s["family"]}_rule_count']) > 0 for r in reread)
        assert count == int(s["hit_molecules"]) and int(s["denominator_all_attempted"]) == 9098
    assert sum(int(r["molecules"]) for r in read_csv(OUT / "overlap_patterns.csv")) == 9098
    assert all(c["status"] == "PASS" for c in direct_checks)
    ended = iso()
    receipt = {"started_at_cst":started,"finished_at_cst":ended,"wall_seconds":round(time.perf_counter()-wall_start,6),"process_cpu_seconds":round(time.process_time()-cpu_start,6),"cpu_affinity":affinity,"thread_environment":{k:os.environ[k] for k in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS")},"software":{"python":platform.python_version(),"rdkit":rdBase.rdkitVersion,"matplotlib":matplotlib.__version__,"operating_system":platform.system(),"machine":platform.machine()},"logical_representation_attempts":18196,"unique_smiles_catalog_evaluations":len(cache),"identical_smiles_reuse":18196-len(cache),"frozen_rules_sha256":sha(rules_path),"rule_catalog_sha256":sha(catalog_path),"original_all_predictions_sha256":EXPECTED_INPUT,"eligible_9098_csv_sha256":sha(ROOT/"frozen_domain/eligible_9098.csv"),"domain_freeze_sha256":sha(ROOT/"frozen_domain/DOMAIN_FREEZE.json"),"script_sha256":sha(Path(__file__)),"official_source_entry_counts":source_counts,"all_585_runtime_entries_serialize_identically_to_frozen_official_smarts":True,"direct_smarts_spotcheck":{"seed":20260930,"random_molecules":12,"plus_six_candidates":True,"molecules":len(selected_ids),"representations":len(direct_checks),"rule_comparisons":sum(r["rules_compared"] for r in direct_checks),"mismatches":0},"independent_csv_recounts_pass":True,"all_9098_ids_order_and_parent_uniqueness_pass":True,"denominator_contains_every_attempt_including_failures":True}
    dump(OUT / "execution_receipt.json",receipt)
    def hit(rep,fam): return next(s["hit_molecules"] for s in stats if s["representation"]==rep and s["family"]==fam)
    candidate_lines = "\n".join(f'| {r["compound_id"]} | {len(r["union_pains_rules"])} | {len(r["union_brenk_rules"])} | {len(r["union_reactive_rules"])} |' for r in summary["six_candidates"])
    readme = f'''# 冻结9098分子结构提示诊断

已实际执行RDKit 2025.09.6官方FilterCatalog：PAINS 480条（A=16、B=55、C=409），BRENK 105条。反应性提示为冻结JSON指定的10条BRENK子集。全部585条运行条目与冻结官方SMARTS重建后的完整序列化逐一一致。

输入来自冻结R2的16698行预测表，严格取 `docking_eligibility == eligible` 的9098行。compound_id与parent_key均唯一，顺序与原输入逐一一致。主表示为原冻结的 `rule_smiles`；敏感性表示为 `standardized_isomeric_smiles`；同时报告两表示命中并集，不改变库或原排名。

| 表示 | PAINS命中/9098 | BRENK命中/9098 | 反应性子集命中/9098 |
|---|---:|---:|---:|
| 主表示 | {hit("primary","pains")} | {hit("primary","brenk")} | {hit("primary","reactive")} |
| 源表示敏感性 | {hit("sensitivity","pains")} | {hit("sensitivity","brenk")} | {hit("sensitivity","reactive")} |
| 两表示并集 | {hit("union","pains")} | {hit("union","brenk")} | {hit("union","reactive")} |

全部尝试行保留；主表示失败{9098-next(s["evaluable_molecules"] for s in stats if s["representation"]=="primary")}，源表示失败{9098-next(s["evaluable_molecules"] for s in stats if s["representation"]=="sensitivity")}。两表示规则集合不同的分子数为{summary["representation_alert_set_differs_molecules"]}。分母始终为9098，另列可解析分母；失败不能解释为无提示。

| 六候选 | PAINS并集规则数 | BRENK并集规则数 | 反应性并集规则数 |
|---|---:|---:|---:|
{candidate_lines}

具体规则见 `chemical_alerts_9098.csv` 与 `six_candidate_alerts.csv`。每条规则频次（含零命中）见 `rule_frequencies.csv`；所有比例与失败分母见 `alert_summary.csv`；互斥并集模式见 `overlap_patterns.csv`。图 `chemical_alerts_9098.png` / `.pdf` 使用独立条形展示重叠提示，并另列合计9098的互斥模式；可复算图源见 `figure_source.csv`。

自检从产出CSV独立重读重算汇总、规则并集与全部身份。固定种子随机12分子加六候选的两个表示用全部585条官方SMARTS（mergeHs=True）直接核对目录，共{receipt["direct_smarts_spotcheck"]["rule_comparisons"]}次比较，零差异。条目检查见 `catalog_runtime_checks.csv`，抽查见 `direct_smarts_spotchecks.csv`。输入、规则、脚本SHA256及真实软件、时间和单CPU亲和性回执见 `execution_receipt.json`。

结构命中是复核提示，不能当作实验反应性、毒性、测定干扰、TMC1生物学活性或疗效证据。BRENK子集包含广义三元杂环，不能仅解释为环氧化物。各提示集合有重叠，禁止直接相加得到分子总数。规则来源、许可与原始引用完整保留于 `../rules/README.md`、`../rules/RULES_FROZEN.json` 和 `../rules/official_sources/`。

复现：使用RDKit 2025.09.6及matplotlib运行 `../scripts/run_chemical_alerts.py`。程序主动限制为1个逻辑CPU与1个工作进程，并按冻结SHA校验输入和规则；不会写回R2源文件。
'''
    (OUT/"README.md").write_text(readme,encoding="utf-8")
    hashes = {p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name != "DIAGNOSTIC_READY.json"}
    ready = {"status":"DIAGNOSTIC_READY","ready_at_cst":iso(),"rows":9098,"all_checks_pass":True,"failed_primary":9098-next(s["evaluable_molecules"] for s in stats if s["representation"]=="primary"),"failed_sensitivity":9098-next(s["evaluable_molecules"] for s in stats if s["representation"]=="sensitivity"),"union_alert_counts":{f:hit("union",f) for f in ("pains","brenk","reactive")},"public_outputs_sha256":hashes,"no_changes_to_original_library_or_ranking":True,"scope":"Frozen structural-alert diagnostic only; not measured biology, reactivity or toxicity."}
    dump(OUT/"DIAGNOSTIC_READY.json",ready)
    print("DIAGNOSTIC_READY " + json.dumps(ready,ensure_ascii=False),flush=True)

if __name__ == "__main__":
    main()
