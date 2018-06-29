from .about import about
from .management import init, remove
from .library import lib_info, lib_edit, lib_generate
from .categories import cat_list, cat_show, cat_add, cat_del, cat_edit, cat_clear
from .elements import elm_list, elm_show, elm_add, elm_del, elm_import

commands = [
    about,
    init, remove,
    cat_list, cat_show, cat_add, cat_del, cat_edit, cat_clear,
    elm_list, elm_show, elm_add, elm_del, elm_import,
    lib_info, lib_edit, lib_generate,
]
