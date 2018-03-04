from collections import namedtuple

FontSettings = namedtuple('FontSettings',
    [
    'legend_size',
    'label_size',
    'title_size',
    'tick_size',
    ])

legend = 16
label   = 20
title  = 25
tick = legend

PlotSettings = FontSettings(
    legend_size = legend,
    label_size  = label,
    title_size  = title,
    tick_size   = tick,
    )
