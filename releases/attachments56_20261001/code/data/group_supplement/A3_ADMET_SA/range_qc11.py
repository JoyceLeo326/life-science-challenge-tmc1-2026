import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
with (HERE / "ADMET_AI_v2_endpoint_metadata.csv").open(encoding="utf-8-sig", newline="") as f:
    info = {r["id"]: r for r in csv.DictReader(f)}
with (HERE / "ADMET_AI_v2_11_full_predictions.csv").open(encoding="utf-8-sig", newline="") as f:
    predictions = list(csv.DictReader(f))

audit = []
for row in predictions:
    for endpoint, meta in info.items():
        if endpoint not in row:
            continue
        value = float(row[endpoint])
        minimum, maximum = float(meta["minimum"]), float(meta["maximum"])
        audit.append({
            "id": row["id"], "endpoint": endpoint, "value": row[endpoint],
            "name": meta["name"], "task_type": meta["task_type"], "units": meta["units"],
            "metadata_minimum": meta["minimum"], "metadata_maximum": meta["maximum"],
            "outside_metadata_range": (not math.isfinite(value)) or value < minimum or value > maximum,
        })

fields = list(audit[0])
for name, rows in [
    ("ADMET_AI_v2_range_qc_all.csv", audit),
    ("ADMET_AI_v2_out_of_range.csv", [r for r in audit if r["outside_metadata_range"]]),
]:
    with (HERE / name).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
print(f"{len(audit)} endpoint values; {sum(x['outside_metadata_range'] for x in audit)} outside metadata range")

