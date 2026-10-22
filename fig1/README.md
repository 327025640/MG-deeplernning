# fig1 代码

按图号组织；独立样本核查和输入核验脚本不收录。

| 文件 | 用途 |
|---|---|
| [Figure1A_all_samples_coverage.R](Figure1A_all_samples_coverage.R) | 最终样式入口 |
| [Figure1B_umap_layers.R](Figure1B_umap_layers.R) | 最终样式入口 |
| [Figure1F_panel_combined.R](Figure1F_panel_combined.R) | 最终样式入口 |
| [Figure1_Disease_V3_Disease17_all_in_one.R](Figure1_Disease_V3_Disease17_all_in_one.R) | DEG、GO及主图生产脚本 |
| [config/cluster_color_unified_res0.6.csv](config/cluster_color_unified_res0.6.csv) | 共享上游/配置；现行30簇res0.6颜色；绘图脚本实际读取的映射 |
| [config/dataset_palette_contract.csv](config/dataset_palette_contract.csv) | 共享上游/配置；冻结疾病与数据集颜色契约 |
| [config/disease_palette_contract.csv](config/disease_palette_contract.csv) | 共享上游/配置；冻结疾病与数据集颜色契约 |
| [upstream/01_counts_only_merge_v3.R](upstream/01_counts_only_merge_v3.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
| [upstream/02_normalize_hvg_pca_v3.R](upstream/02_normalize_hvg_pca_v3.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
| [upstream/02b_build_metadata_injection_v3.R](upstream/02b_build_metadata_injection_v3.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
| [upstream/03d_harmony_dims16_preview_v3.R](upstream/03d_harmony_dims16_preview_v3.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
| [upstream/07_update_atlas_disease_metadata_v3.R](upstream/07_update_atlas_disease_metadata_v3.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
| [upstream/08_multires_scan_v3.R](upstream/08_multires_scan_v3.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
| [upstream/19_create_merged_obj_res06_protein_coding.R](upstream/19_create_merged_obj_res06_protein_coding.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
| [upstream/19c_add_disease_analyse_alias.R](upstream/19c_add_disease_analyse_alias.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
| [upstream/20_findallmarkers_res06_from_merged_obj.R](upstream/20_findallmarkers_res06_from_merged_obj.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
| [upstream/20b_resume_findallmarkers_res06_stage1_publish.R](upstream/20b_resume_findallmarkers_res06_stage1_publish.R) | 共享上游/配置；上游来源；不作为本次重算指令 |
