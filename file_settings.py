"""File settings for running multiple analyses across datasets"""
from enum import Enum


class FileSettings(object):

    DATA_FOLDER = "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data"
    IGNORED = {"esc_unusable"}
    TABLES_FOLDER = (
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables"
    )
    FIGURES_FOLDER = (
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/figures"
    )
    SELECTED_CONDITION = ("100%", " 1 ms")

    THRESHOLD_LIST = ["nothresh", 1, 2, 4]

