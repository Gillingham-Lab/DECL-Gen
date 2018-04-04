from .about import about
from .management import init, remove
from .categories import cat_add, cat_clear, cat_del, cat_edit, cat_list, cat_show
from .library import lib_info

commands = [
    about,
    init, remove,
    cat_add, cat_clear, cat_del, cat_edit, cat_list, cat_show,
    lib_info,
]
