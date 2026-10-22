# Microglia atlas 主图代码

本仓库保存主图相关的分析与绘图代码快照，整理于 **2026-10-22**。Figure5、Figure6 使用本次确认的最新 PCSF 方案。

| 索引 | 代码范围 | 版本依据 |
|---|---|---|
| [Figure1](figures/Figure1.md) | 疾病覆盖、冻结 UMAP、样本级 pseudobulk DEG 与 GO | Disease17 主分析与图件收尾；全样本覆盖核查分支单独标注 |
| [Figure2](figures/Figure2.md) | NMF、参考签名、疾病 usage、平台内差值与 MP04 | 冻结 consensus **12 MP** |
| [Figure3](figures/Figure3.md) | 30 簇身份、M1/M2-like 连续评分、疾病组成 | 冻结 res0.6 及其现行图源 |
| [Figure4](figures/Figure4.md) | clean24 FUN8、CollecTRI ULM、转录相似性、pan14 轨迹与 PTG | 冻结轨迹与逐簇 TF Top3 |
| [Figure5](figures/Figure5.md) | PTG14 模型、PPI/PCSF、seed42 网络、权重、GO 与模型内扰动 | **b=10、w=3、g=3，seed42：82 节点、81 骨架边；76 个模型证据节点、6 个连接节点** |
| [Figure6](figures/Figure6.md) | Disease vs Ctrl 全基因模型、SHAP、逐病 PCSF、药物网络接近度与双支持 Top3 | **30–75 节点方案**；15 病有网络，HD 无合格网络 |
| [共享上游](figures/Shared.md) | 图谱构建来源、冻结颜色、签名契约和第三方许可 | 保留原始脚本与配置 |

Figure1–4 的分组沿用已有溯源索引，便于按分析内容查找；源文件中的历史 Fig2/Fig3 名称保留，不自动重命名。最终拼版 panel 字母由 Illustrator 排版决定。

## 文件组织与来源

- `source/project/` 对应原 `D:/data/microglia_new/` 下的代码。
- `source/outputs/` 对应原 `E:/micro_output/` 下随分析保存的代码和白名单配置。
- [code_manifest.csv](code_manifest.csv) 记录每个来源文件的归属、用途、仓库路径、原始路径、字节数及 SHA256。
- [dependency_inventory.csv](dependency_inventory.csv) 列出源码中静态识别的 R 包和 Python import 名称，包含标准库及本地模块，**不是安装清单或锁定环境**。

脚本按原字节复制，保留原始路径、参数、注释、来源哈希及版本。`.gitattributes` 禁用自动换行转换，使仓库中的脚本哈希可与本地来源核对。当前图件登记所引用的准确生产脚本也已收录；例如 Figure5B 的实际生成脚本位于原归档目录，这是当前登记确认的来源，不是按目录日期推定的替代版本。

## 使用条件

本快照用于代码归档和方法追溯。脚本保留原 Windows 路径，需在原项目环境及冻结输入可用时使用；直接克隆到其他机器不能立即复现全部分析。仓库不含单细胞表达、细胞/供体级元数据、大 RDS、模型权重、完整 SHAP、导师 PPI、药物数据库和成图。

本地 R 使用 `D:/R/R-4.4.3/bin/Rscript.exe`。Python 历史步骤使用各自运行环境；PCSF 记录为 OmicsIntegrator2 2.3.10 / pcst_fast 1.0.7。不要用当前系统 Python 替代冻结的训练环境。部分脚本含写入、清理中间文件或训练步骤，应先核对路径和输入；本次整理没有运行这些分析入口。

Figure5 的旧 `run_pan14_pcsf.py` 作为 PPI 映射和输入生产来源保留，其默认 b5/w3/g5 **不是当前主参数**；当前 b10/w3/g3 以参数评估配置和 seed42 来源准备脚本为准。旧多数28候选与 Figure6 OVR Top10 没有作为当前主分析入口。

Figure6 的 30–75 主运行复用了较早参数扫描产物，因此同时保留了显式依赖的扫描代码；历史扫描参数不替代当前 `selected_parameters.csv` 及其锁文件。当前药物图使用最终网络全部入选节点，包括连接节点；PPI 为 number≥2 的无向无权最大连通分量，逆转支持沿用冻结历史结果。

本次只检查代码归档、语法、当前版本映射和复制完整性，没有重算科学结果。部分原 R 运行记录为产物核验完成但进程异常退出，本次代码归档不将其改写为正常运行；第三方旧 docstring 的 Python 转义警告也保留原文。

## 第三方代码许可

两套 CauFinder 适配器按冻结版本分别保留，来源为 [CauFinder](https://github.com/ChengmingZhang-CAS/CauFinder-main)，沿用随原适配器保存的 [MIT License](source/outputs/分析审阅中心/待审阅/PTG02_CD_分支疾病与调控候选/data/pan14_disease/caufinder_adapter/LICENSE)。本次未给整个研究项目另行指定开源许可证。
