"""Compatibility helpers so one code base runs on QGIS 3.x (PyQt5) and QGIS 4 (PyQt6)."""
from qgis.PyQt.QtCore import QVariant
from qgis.core import QgsProcessing, QgsWkbTypes

try:
    from qgis.PyQt.QtCore import QMetaType
except ImportError:  # very old builds
    QMetaType = None

try:
    from qgis.core import Qgis
except ImportError:
    Qgis = None


def _field_type(name):
    """Return the field type constant accepted by QgsField on this QGIS build."""
    try:
        return getattr(QVariant, name)  # PyQt5 builds
    except AttributeError:
        return getattr(QMetaType.Type, name)  # PyQt6 builds


class VT:
    """Field type constants for QgsField."""
    Double = _field_type('Double')
    Int = _field_type('Int')
    String = _field_type('String')


def _geometry(new_name, old_name):
    enum = getattr(Qgis, 'GeometryType', None) if Qgis is not None else None
    if enum is not None and hasattr(enum, new_name):
        return getattr(enum, new_name)
    return getattr(QgsWkbTypes, old_name)


GEOM_POINT = _geometry('Point', 'PointGeometry')
GEOM_LINE = _geometry('Line', 'LineGeometry')
GEOM_POLYGON = _geometry('Polygon', 'PolygonGeometry')


def _source(new_name, old_name):
    enum = getattr(Qgis, 'ProcessingSourceType', None) if Qgis is not None else None
    if enum is not None and hasattr(enum, new_name):
        return getattr(enum, new_name)
    return getattr(QgsProcessing, old_name)


SRC_POINT = _source('VectorPoint', 'TypeVectorPoint')
SRC_LINE = _source('VectorLine', 'TypeVectorLine')
SRC_POLYGON = _source('VectorPolygon', 'TypeVectorPolygon')
SRC_ANY = _source('VectorAnyGeometry', 'TypeVectorAnyGeometry')
