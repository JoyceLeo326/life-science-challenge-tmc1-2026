# 历史审计与现行证据的范围

本目录保留来源历史，不能将不同快照的数目混用。

| 记录 | 对应范围 | 本轮使用方式 |
|---|---|---|
| summary.json、quality_profile.json、snapshot_manifest.json | 23475来源记录及17645理化合格母体的原快照 | 来源、去重、过滤和历史SHA追踪 |
| final_library_independent_audit.json | 17645行CSV与早期SDF抽查 | CSV身份和母体/前向池审计；其中17645条旧SDF数不是本次exact SDF条数 |
| repaired_sdf_independent_audit.json | SHA46d107…的17632条exact SDF与13条隔离 | 本次正式SDF的全量图/完整立体键既有独立审计 |
| library_split_audit.json | 早期5067行试验库与90行试选、81个独立分子 | 历史问题记录，不能代表17645/16698/9098正式分母或最终三臂结果 |
| historical_overlap_summary.json | 971个历史身份、935个历史QC标签 | 按来源ID或非立体母体键排除历史交集 |
| ../frozen_domain/DOMAIN_FREEZE.json | 16698前向行中固定9098可对接行 | 当前结构警示诊断的唯一分母 |

SOURCE_PUBLIC_HASH_MAP.csv保留original_sha256与public_sha256。仅对绝对或旧逻辑路径、目录导航和说明文字作公开副本调整；科学CSV内容、原始API页、正式SDF及历史审计科学数值不改。原审计中旧文件SHA仍表示其原快照，现公开文件字节以本包MANIFEST.json为准。
