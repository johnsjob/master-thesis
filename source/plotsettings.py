from collections import namedtuple

FontSettings = namedtuple('FontSettings',
    [
    'legend_size',
    'label_size',
    'title_size',
    'tick_size',
    ])

title   = 25
label   = 20
legend  = 25
tick    = 16

PlotSettings = FontSettings(
    legend_size = legend,
    label_size  = label,
    title_size  = title,
    tick_size   = tick,
    )
