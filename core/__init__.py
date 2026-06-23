#!/usr/bin/env python3
"""
core/ — Fuxi (伏羲) 公共基础设施

统一 scRNA-seq 和 scATAC-seq 管线的共享代码:
  - utils.py:      safe_write, safe_plot, setup_logger, resolve_config, validate_adata, monitor_performance
  - ai_caller.py:  统一 LLM 调用 (重试、思考模式、磁盘缓存、模型自发现)
  - ai_prompts.py: 统一 prompt 模板 (RNA + ATAC + modality-specific)
  - config.py:     统一配置 (BaseConfig + modality-specific fields)
  - run_pipeline.py: 统一 CLI (--modality rna|atac|spatial)
  - dataset_schema.py:   dataset.yaml Python 数据模型
  - dataset_detector.py: 自动检测数据组学类型
"""
