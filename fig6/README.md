# fig6 代码

按图号组织；独立样本核查和输入核验脚本不收录。

| 文件 | 用途 |
|---|---|
| [01_export_full_gene_inputs.R](01_export_full_gene_inputs.R) | Disease-vs-Ctrl全基因模型和SHAP上游 |
| [02_train_disease_vs_ctrl.py](02_train_disease_vs_ctrl.py) | Disease-vs-Ctrl全基因模型和SHAP上游 |
| [CauFinder/__init__.py](CauFinder/__init__.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/benchmark.py](CauFinder/benchmark.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/benchmark_crn.py](CauFinder/benchmark_crn.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/caufinder_main.py](CauFinder/caufinder_main.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/causal_effect.py](CauFinder/causal_effect.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/dataloader.py](CauFinder/dataloader.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/dataloader_crn.py](CauFinder/dataloader_crn.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/driver_regulators.py](CauFinder/driver_regulators.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/layers.py](CauFinder/layers.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/models.py](CauFinder/models.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [CauFinder/utils.py](CauFinder/utils.py) | Disease-vs-Ctrl冻结CauFinder适配器；沿用原MIT许可 |
| [Fig6_pcsf_drug_accessibility_30_75_all_nodes_ge2.R](Fig6_pcsf_drug_accessibility_30_75_all_nodes_ge2.R) | 全节点药物网络接近度及双支持Top3 |
| [assemble_fig6_disease_compact_le50_base_grid.py](assemble_fig6_disease_compact_le50_base_grid.py) | 当前运行显式引用或复用的参数扫描来源 |
| [config/all_disease_training_configurations.json](config/all_disease_training_configurations.json) | 冻结逐病模型训练配置 |
| [config/current_figure_registry.csv](config/current_figure_registry.csv) | 当前图件与脚本溯源登记；不是本次重跑记录 |
| [config/run_config.json](config/run_config.json) | 30–75节点运行配置与逐病参数锁 |
| [config/selected_parameters.csv](config/selected_parameters.csv) | 冻结逐病PCSF参数 |
| [config/selected_parameters.lock.json](config/selected_parameters.lock.json) | 30–75节点运行配置与逐病参数锁 |
| [dependencies/fig6_30_75_final_closeout.R](dependencies/fig6_30_75_final_closeout.R) | 由归档脚本的现存绝对代码路径引用补齐 |
| [dependencies/render_Fig6_Top10_PPI_compact.R](dependencies/render_Fig6_Top10_PPI_compact.R) | 由归档脚本的现存绝对代码路径引用补齐 |
| [dependencies/render_fig6_pcsf_drug_dual_support_top3_12d.R](dependencies/render_fig6_pcsf_drug_dual_support_top3_12d.R) | 由归档脚本的现存绝对代码路径引用补齐 |
| [dependencies/render_fig6_pcsf_networks_30_75_no_modules.R](dependencies/render_fig6_pcsf_networks_30_75_no_modules.R) | 由归档脚本的现存绝对代码路径引用补齐 |
| [dependencies/render_fig6_pcsf_networks_layout_v3.R](dependencies/render_fig6_pcsf_networks_layout_v3.R) | 由归档脚本的现存绝对代码路径引用补齐 |
| [export_orgdb_symbol_entrez.R](export_orgdb_symbol_entrez.R) | SHAP到PPI映射输入 |
| [fig6_30_75_final_closeout.R](fig6_30_75_final_closeout.R) | 当前图件登记的准确生产脚本 |
| [fig6_disease_compact_30_75_parameter_selection.py](fig6_disease_compact_30_75_parameter_selection.py) | 当前运行显式引用或复用的参数扫描来源 |
| [fig6_disease_compact_30_75_parameter_selection_1291b34da502.py](fig6_disease_compact_30_75_parameter_selection_1291b34da502.py) | 30–75节点正式运行来源 |
| [fig6_disease_compact_le50_b2g4_base_addendum.py](fig6_disease_compact_le50_b2g4_base_addendum.py) | 当前运行显式引用或复用的参数扫描来源 |
| [fig6_disease_compact_le50_lowb_base_grid.py](fig6_disease_compact_le50_lowb_base_grid.py) | 当前运行显式引用或复用的参数扫描来源 |
| [fig6_disease_compact_le50_parameter_selection.py](fig6_disease_compact_le50_parameter_selection.py) | 当前运行显式引用或复用的参数扫描来源 |
| [fig6_disease_pcsf_size_target_scan.py](fig6_disease_pcsf_size_target_scan.py) | 当前运行显式引用或复用的参数扫描来源 |
| [fig6_pcsf_stage_summary_compat.py](fig6_pcsf_stage_summary_compat.py) | 30–75节点正式运行来源 |
| [pcsf_runtime_patched.py](pcsf_runtime_patched.py) | 30–75节点正式运行来源 |
| [prepare_pcsf_inputs.py](prepare_pcsf_inputs.py) | SHAP到PPI映射输入 |
| [render_fig6_pcsf_drug_dual_support_top3_30_75.R](render_fig6_pcsf_drug_dual_support_top3_30_75.R) | 全节点药物网络接近度及双支持Top3 |
| [render_fig6_pcsf_networks_30_75_no_modules.R](render_fig6_pcsf_networks_30_75_no_modules.R) | 30–75节点正式运行来源 |
| [render_sc_30_75_only.R](render_sc_30_75_only.R) | 当前图件登记的准确生产脚本 |
