from .about import about
from .management import init, remove
from .library import lib_info, lib_edit, lib_generate
from .diversityPoints import pnt_list, pnt_show, pnt_add, pnt_del, pnt_edit, pnt_clear
from .diversitySuperPoints import ssp_add, ssp_pnt_list, ssp_pnt_add, ssp_pnt_del, ssp_pnt_reindex
from .diversityElements import elm_list, elm_show, elm_add, elm_del, elm_import, elm_copy

commands = [
    about,
    init, remove,
    lib_info, lib_edit, lib_generate,
    pnt_list, pnt_show, pnt_add, pnt_del, pnt_edit, pnt_clear,
    elm_list, elm_show, elm_add, elm_del, elm_import, elm_copy,
    ssp_add, ssp_pnt_list, ssp_pnt_add, ssp_pnt_del, ssp_pnt_reindex,
]
