# 第三轮阶段二 P0 库证据

本目录补齐阶段二要求的可追溯证据，不改动第二轮筛选库和模型结果。

## 分母链

| 层级 | 数量 | 说明 |
|---|---:|---|
| ChEMBL来源记录 | 23,475 | 第一轮只读3,475条加本轮公开API新增20,000条 |
| 去重来源ID | 23,403 | 按来源ID去重 |
| 标准化独立母体 | 21,356 | 去盐、标准化及母体归并 |
| 理化过滤后 | 17,645 | 按预设性质范围 |
| 前瞻未标注池 | 16,698 | 排除历史身份后的未标注池 |
| 可进入实际对接 | 9,098 | 事前冻结化学与立体门禁后 |

## 证据成员

- `all_source_records.csv` 和 `all_unique_parents.csv` 保留来源和母体层
- `filtered_for_screening.csv` 与 `forward_unlabeled_pool.csv` 保留17,645和16,698层
- `filtered_for_screening_stereo_safe_2d.sdf` 是安全二维结构文件
- 四个 `*_exclusions.csv` 是立体信息隔离清单，明确哪些记录未进入安全SDF表示
- `final_library_independent_audit.json` 和 `library_split_audit.json` 是独立逐条核查结果

## 边界

这些文件证明分母和输入链可以重建，不等于分子具有TMC1活性，也不等于已经完成湿实验。
