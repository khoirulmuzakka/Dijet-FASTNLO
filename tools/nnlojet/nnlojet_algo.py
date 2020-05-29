#!/usr/bin/env python3

import numpy as np


def is_outlier_MAD(points, thresh=3.5):

    median = np.nanmedian(points)
    abs_diff = np.abs(points - median)

    mad = np.nanmedian(abs_diff)
    if mad == 0.:
        modified_z_score = 0. * abs_diff
    else:
        modified_z_score = 0.6745 * abs_diff / mad

    return modified_z_score > thresh


def is_outlier_doubleMAD(points, thresh=3.5):

    median = np.nanmedian(points)
    abs_diff = np.abs(points - median)

    left_mad = np.nanmedian(abs_diff[points <= median])
    right_mad = np.nanmedian(abs_diff[points >= median])

    if left_mad == 0. or right_mad == 0.:
        modified_z_score = 0. * abs_diff

    else: 
        y_mad = left_mad * np.ones(len(points))
        y_mad[points > median] = right_mad
        modified_z_score = 0.6745 * abs_diff / y_mad
        modified_z_score[points == median] = 0.

    return modified_z_score > thresh


def is_outlier_IQR(points, thresh=1.5):

    median = np.nanmedian(points)
    lower_quartile = np.nanpercentile(points, 25.)
    upper_quartile = np.nanpercentile(points, 75.)
    
    # the inter-quartile-range
    iqr = upper_quartile - lower_quartile

    # fences
    lower_fence = lower_quartile - thresh * iqr
    upper_fence = upper_quartile + thresh * iqr

    # return (points < lower_fence) | (points > upper_fence)
    return np.array([ (pt < lower_fence) | (pt > upper_fence) if np.isfinite(pt) else False for pt in points ])
