# -*- coding: utf-8 -*-
def classFactory(iface):
    from .plugin import OnshoreSpillPlugin
    return OnshoreSpillPlugin(iface)
