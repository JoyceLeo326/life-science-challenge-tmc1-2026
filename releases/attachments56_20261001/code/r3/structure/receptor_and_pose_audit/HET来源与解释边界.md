# 7USX 非蛋白条目的来源与解释边界

## 化学名称与分类证据

7USX 官方 PDB 的 HETNAM 和 RCSB Chemical Component Dictionary（CCD）名称相符。CCD 定义模型条目的化学身份，不能单凭 CCD 名称判断实验样本中每一个位置的分子来源或是否属于生理膜成分。

| CCD代码 | 7USX/CCD化学名称 | 可确定的化学类别 | 本结构来源性质的证据界限 |
|---|---|---|---|
| D12 | dodecane，十二烷 | 烃 | 实验模型的烃占位条目；逐项来源未确认，不能直接认定为生理脂质或特定去污剂片段。[CCD](https://www.rcsb.org/ligand/D12) |
| R16 | hexadecane，十六烷 | 烃 | 实验模型的烃占位条目；逐项来源未确认，不能直接认定为生理脂质或特定去污剂片段。[CCD](https://www.rcsb.org/ligand/R16) |
| 3PE | 1,2-distearoyl-sn-glycerophosphoethanolamine | 磷脂酰乙醇胺 | 模型具有磷脂化学身份；具体来源和脂链组成不能由 CCD 名称作为实验化学鉴定。[CCD](https://www.rcsb.org/ligand/3PE) |
| PEE | 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine | 磷脂酰乙醇胺 | 模型具有磷脂化学身份；具体来源和脂链组成不能由 CCD 名称作为实验化学鉴定。[CCD](https://www.rcsb.org/ligand/PEE) |
| PLM | palmitic acid，棕榈酸 | 饱和脂肪酸/模型中的棕榈酰连接 | 7USX LINK 明确连接 TMIE D/F CYS43/44；原论文讨论 TMIE C44 棕榈酰化。[CCD](https://www.rcsb.org/ligand/PLM) |
| CA | calcium ion，钙离子 | 无机离子 | 原论文纯化方法加入 3 mM CaCl2；不能由此判定每个模型钙离子的天然来源。[CCD](https://www.rcsb.org/ligand/CA) |
| CLR | cholesterol，胆固醇 | 甾醇 | 官方模型名称；逐项天然来源没有在本审计中额外实证。[CCD](https://www.rcsb.org/ligand/CLR) |
| NAG | 2-acetamido-2-deoxy-beta-D-glucopyranose | 糖 | 7USX LINK 明确连接 TMC1 A/B ASN209，为模型的糖基化环境。[CCD](https://www.rcsb.org/ligand/NAG) |

## 原始论文的实验背景

2022 年原始研究的 Methods“Isolation of the native TMC-1 complex”记录用 2%（w/v）GDN 提取，并在包含 0.02%（w/v）GDN 和 3 mM CaCl2 的缓冲液中进行 SEC。正文“Lipid-mediated interactions of TMIE with TMC-1”将 TMIE/TMC1 腔体的占据物描述为去污剂或脂质。该证据支持样本包含去污剂环境，不能把每个 D12/R16 坐标逐项鉴定为特定去污剂或生理脂链。[原始论文](https://www.nature.com/articles/s41586-022-05314-8)

本审计实际读取本次保存的原始论文全文 XML，SHA256 为 `f62fcc6044aa8a39c38de48074fb454bfb1952c3cc8256016cd04d885fd5727a`。论文未在所核对内容中明确给出 D12/R16 各位置的具体化学来源；没有推测它们的来源。

## 共价及配位环境

官方 PDB 的 LINK 记录包含：D/F 链 CYS43、CYS44 的 SG 与 PLM C1 的四个连接；A/B 链 ASN209 ND2 与 NAG C1 的两个连接；另有 Ca 与酸性残基氧的配位记录。本蛋白单独计算条件省略这些 HET 条目，因此也省略了原模型中的脂化、糖基化及离子环境。它保持全部保留蛋白的实验重原子坐标，但不能宣称保留完整实验化学环境。[官方 7USX PDB](https://files.rcsb.org/download/7USX.pdb)

## 姿态重叠能够说明什么

全部 10 个线虫第一姿态都与原始模型的非蛋白重原子配置出现小于 1.5 Å 的几何重叠。该结果表示“与冻结实验模型的当前非蛋白原子配置不兼容”。它没有证明生理脂质兼容性，也没有证明实际实验中不能结合；分子置换、移动或其它环境重排是否发生尚未实证。

原始姿态、第一模式规则、框、门禁及所有分数保留不变。该解释修正不生成新计算条件，不进行新的 MD，不把环境标志转换成药效阈值或直接 TMC1 阴性结论。
