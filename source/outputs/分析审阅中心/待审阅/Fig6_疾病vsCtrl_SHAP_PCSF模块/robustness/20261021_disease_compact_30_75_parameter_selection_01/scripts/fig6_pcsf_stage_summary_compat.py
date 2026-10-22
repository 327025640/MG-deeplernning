#!/usr/bin/env python3
"""Compatibility entry point for the frozen 30-75 PCSF driver.

The frozen driver chooses repeat-record paths from runtime_module.CURRENT_PHASE.
During stage-summary reconstruction it did not set that value to the stage
being read. This wrapper supplies that stage locally without changing either
the frozen driver/runtime files or their SHA-locked run configuration.
"""
from __future__ import print_function

import importlib.util
import os
import sys


DRIVER_PATH = r"D:\data\microglia_new\scripts\fig6_disease_compact_30_75_parameter_selection.py"


def main():
    if not os.path.isfile(DRIVER_PATH):
        raise FileNotFoundError(DRIVER_PATH)
    spec = importlib.util.spec_from_file_location("fig6_compact_30_75_driver", DRIVER_PATH)
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)

    original_load_runtime = driver.load_runtime
    original_summary = driver.save_summary_for_stage

    def load_runtime_with_stage_summary(path):
        runtime_module = original_load_runtime(path)

        def save_summary_with_explicit_stage(mod, stage, *args, **kwargs):
            mod.CURRENT_PHASE = str(stage)
            return original_summary(mod, stage, *args, **kwargs)

        driver.save_summary_for_stage = save_summary_with_explicit_stage
        return runtime_module

    driver.load_runtime = load_runtime_with_stage_summary
    sys.argv = [DRIVER_PATH] + sys.argv[1:]
    driver.main()


if __name__ == "__main__":
    main()
