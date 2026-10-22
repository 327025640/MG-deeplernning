# fig5 代码

按图号组织；独立样本核查和输入核验脚本不收录。

| 文件 | 用途 |
|---|---|
| [01_export_ptg14.R](01_export_ptg14.R) | 冻结PTG14模型上游 |
| [03_train.py](03_train.py) | 冻结PTG14模型上游 |
| [CauFinder/LICENSE](CauFinder/LICENSE) | 共享上游/配置；CauFinder原始MIT许可；适用于本仓库CauFinder适配器 |
| [CauFinder/__init__.py](CauFinder/__init__.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/benchmark.py](CauFinder/benchmark.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/benchmark_crn.py](CauFinder/benchmark_crn.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/caufinder_main.py](CauFinder/caufinder_main.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/causal_effect.py](CauFinder/causal_effect.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/dataloader.py](CauFinder/dataloader.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/dataloader_crn.py](CauFinder/dataloader_crn.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/driver_regulators.py](CauFinder/driver_regulators.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/layers.py](CauFinder/layers.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/models.py](CauFinder/models.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/utils.py](CauFinder/utils.py) | 本项目冻结CauFinder适配器；沿用原MIT许可 |
| [assess_pan14_seed42_representativeness.py](assess_pan14_seed42_representativeness.py) | seed42及跨seed支持核对 |
| [config/assessment_config.json](config/assessment_config.json) | 当前seed42参数和跨seed设置 |
| [config/current_artifact_registry.csv](config/current_artifact_registry.csv) | 当前图件与脚本溯源登记；不是本次重跑记录 |
| [config/run_config.json](config/run_config.json) | 原始PCSF网格配置 |
| [config/script_sha256_manifest.csv](config/script_sha256_manifest.csv) | 原始脚本哈希记录 |
| [fig5e_seed42_expression_perturbation.py](fig5e_seed42_expression_perturbation.py) | 当前seed42主图入口 |
| [pan14_pcsf_parameter_adaptation.py](pan14_pcsf_parameter_adaptation.py) | b10/w3/g3来源参数网格 |
| [prepare_fig5_seed42_sources.R](prepare_fig5_seed42_sources.R) | 当前seed42主图入口 |
| [ptg02_p14d_07_network_fulltrain_stability_20260914.py](ptg02_p14d_07_network_fulltrain_stability_20260914.py) | 训练脚本显式引用的网络辅助函数 |
| [render_fig5_seed42_BCD.R](render_fig5_seed42_BCD.R) | 当前seed42主图入口 |
| [render_fig5_seed42_BCD_pre_D_name_fix.R](render_fig5_seed42_BCD_pre_D_name_fix.R) | 当前图件登记的准确生产脚本 |
| [render_fig5_seed42_C_final.R](render_fig5_seed42_C_final.R) | 当前seed42主图入口 |
| [render_fig5_seed42_D_only.R](render_fig5_seed42_D_only.R) | 当前seed42主图入口 |
| [run_pan14_pcsf.py](run_pan14_pcsf.py) | PPI映射和原始输入生产；其旧默认参数不是当前主参数 |
